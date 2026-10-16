/* Machine-precision ProductLog (Lambert W0): the shared double kernels of
 * src/special_functions/sf_machine_plog.c behind the scalar builtin, the
 * NDArray element-wise path and Compile[].
 *
 * Covers: spot values (real and complex) against mpmath 1.3 lambertw at 30
 * digits (computed for this file, not recalled); the branch point -1/e and the
 * double nearest it (which lies BELOW -1/e, so its W0 is complex); the real
 * axis below -1/e (complex result); the regime boundaries of the kernel (the
 * branch-point u-iteration, the Pade seed, the asymptotic seed); declines (a
 * branch k != 0 goes to the MPFR core and must still answer); exact special
 * values unchanged; and the bit-match
 *     scalar == unpacked List == packed List == NDArray == Compile
 * for real and complex arguments.
 */
#include "eval.h"
#include "parse.h"
#include "expr.h"
#include "symtab.h"
#include "core.h"
#include "arithmetic.h"
#include "test_utils.h"
#include "special_functions/sf_machine.h"
#include "special_functions/sf_machine_plog.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ---- helpers ---------------------------------------------------------- */

/* Read a machine number (Real, or Complex[Real, Real]) out of `e`. */
static bool leaf_value(const Expr* e, double* re, double* im) {
    if (e->type == EXPR_REAL) { *re = e->data.real; *im = 0.0; return true; }
    Expr *a, *b;
    if (is_complex((Expr*)e, &a, &b) && a->type == EXPR_REAL && b->type == EXPR_REAL) {
        *re = a->data.real; *im = b->data.real; return true;
    }
    return false;
}

static Expr* eval_str(const char* s) {
    Expr* e = parse_expression(s);
    ASSERT_MSG(e != NULL, "parse failed: %s", s);
    Expr* r = evaluate(e);
    expr_free(e);
    return r;
}

/* Evaluate `s` to a machine scalar. */
static void eval_num(const char* s, double* re, double* im) {
    Expr* r = eval_str(s);
    ASSERT_MSG(leaf_value(r, re, im), "%s: not a machine number", s);
    expr_free(r);
}

/* Evaluate `s` to a List of machine numbers. */
static void eval_list(const char* s, int n, double* re, double* im) {
    Expr* r = eval_str(s);
    ASSERT_MSG(r->type == EXPR_FUNCTION && (int)r->data.function.arg_count == n,
               "%s: expected a list of %d", s, n);
    for (int i = 0; i < n; i++)
        ASSERT_MSG(leaf_value(r->data.function.args[i], &re[i], &im[i]),
                   "%s: element %d not a machine number", s, i + 1);
    expr_free(r);
}

static bool same_bits(double a, double b) { return memcmp(&a, &b, sizeof a) == 0; }

static void assert_rel(const char* in, double want_re, double want_im, double tol) {
    double r, i;
    eval_num(in, &r, &i);
    double err = hypot(r - want_re, i - want_im), mag = hypot(want_re, want_im);
    ASSERT_MSG(err <= tol * mag, "%s: got %.17g + %.17g I, want %.17g + %.17g I",
               in, r, i, want_re, want_im);
}

/* ---- spot values (mpmath lambertw, 30 digits) ------------------------- */

void test_plog_real_spot() {
    assert_rel("ProductLog[1.`]",    0.567143290409783872999968662210, 0, 4e-16);
    assert_rel("ProductLog[10.`]",   1.74552800274069938307430126488, 0, 4e-16);
    assert_rel("ProductLog[-0.3`]", -0.489402227180214933565021502577, 0, 4e-16);
    assert_rel("ProductLog[0.001`]", 0.00099900149733853091073296967597, 0, 4e-16);
    assert_rel("ProductLog[1.`*^300]", 684.247208629760849292015760652, 0, 4e-16);
    /* Inside the branch-point u-iteration (e (x + 1/e) < 1/2). */
    assert_rel("ProductLog[-0.36`]", -0.806084315970817624450041805922, 0, 8e-16);
    assert_rel("ProductLog[-0.3678`]", -0.979360714957830517130292458791, 0, 8e-16);
}

void test_plog_complex_spot() {
    assert_rel("ProductLog[Complex[1.`, 1.`]]",
               0.656966069230436405873935189615, 0.325450339413415029989281951389, 4e-16);
    assert_rel("ProductLog[Complex[-2.`, 0.5`]]",
               0.290053967513204697800529012275, 1.51498512499366943360199373533, 4e-16);
    assert_rel("ProductLog[Complex[100.`, -50.`]]",
               3.46784445600891799740913332159, -0.360161243768948779000070295568, 4e-16);
    /* Real axis below -1/e: a complex result, the cut's upper side. */
    assert_rel("ProductLog[-0.5`]",
               -0.794023632344689367963015321901, 0.770111750510379109681313077405, 4e-16);
    assert_rel("ProductLog[-1.`*^6]",
               11.3544677771688590841008960817, 2.8921791489948752450172209892, 4e-16);
}

/* ---- the branch point ------------------------------------------------- */

void test_plog_branch_point() {
    /* The double nearest -1/e is -1/e - 1.24e-17, BELOW the branch point: its
     * W0 is complex, -1 + i sqrt(2 e * 1.24e-17) (to first order). */
    double r, i;
    eval_num("ProductLog[-0.36787944117144233`]", &r, &i);
    ASSERT_MSG(fabs(r + 1.0) < 1e-15 && fabs(i - 8.22007971483661770e-09) < 1e-22,
               "W0 at the double nearest -1/e: %.17g + %.17g I", r, i);
    /* The real kernel declines it (the value is complex) ... */
    double w;
    ASSERT_MSG(!sf_machine_productlog(-0.36787944117144233, &w), "real kernel must decline");
    /* ... and answers just above it with full relative accuracy in w + 1. */
    double x = nextafter(-0.36787944117144233, 0.0);   /* -1/e + 4.31e-17 */
    ASSERT(sf_machine_productlog(x, &w));
    /* mpmath: W0(x) + 1 = 1.53042541284951713e-8 */
    ASSERT_MSG(fabs((w + 1.0) - 1.53042541284951713e-08) < 2.3e-16, "w + 1 = %.17g", w + 1.0);
}

/* ---- exact values / declines are unchanged ----------------------------- */

void test_plog_exact_and_declines() {
    assert_eval_eq("ProductLog[0]", "0", 0);
    assert_eval_eq("ProductLog[E]", "1", 0);
    assert_eval_eq("ProductLog[-1/E]", "-1", 0);
    assert_eval_eq("ProductLog[-Pi/2]", "(1/2*I) Pi", 0);
    assert_eval_eq("ProductLog[x]", "ProductLog[x]", 0);
    assert_eval_eq("ProductLog[0.`]", "0.0", 0);
    /* Branch k != 0 is not the machine kernel's; the MPFR core still answers. */
    assert_rel("ProductLog[-1, -0.2`]", -2.54264135777352633279817223827, 0, 4e-16);
    /* Arbitrary precision is untouched. */
    assert_eval_startswith("N[ProductLog[2], 30]", "0.852605502013725491346472414");
}

/* ---- five surfaces agree bit-for-bit ----------------------------------- */

#define NPT 8
void test_plog_surfaces_real() {
    static const char* xs = "{-0.367`, -0.3`, -0.1`, 0.25`, 1.`, 2.718`, 37.5`, 1.`*^8}";
    double s_re[NPT], s_im[NPT], a_re[NPT], a_im[NPT];
    char buf[512];
    /* scalar, one call each */
    Expr* l = eval_str(xs);
    for (int k = 0; k < NPT; k++) {
        snprintf(buf, sizeof buf, "ProductLog[%s[[%d]]]", xs, k + 1);
        eval_num(buf, &s_re[k], &s_im[k]);
        double w; ASSERT(sf_machine_productlog(l->data.function.args[k]->data.real, &w));
        ASSERT_MSG(same_bits(w, s_re[k]), "scalar != kernel at %d", k);
    }
    expr_free(l);
    const char* forms[] = {
        "ProductLog[FromPackedArray[%s]]",
        "FromPackedArray[ProductLog[ToPackedArray[%s]]]",
        "Normal[ProductLog[NDArray[%s]]]",
        "Compile[{{x, _Real}}, ProductLog[x]] /@ %s",
        "FromPackedArray[Compile[{{v, _Real, 1}}, ProductLog[v]][%s]]",
    };
    for (size_t f = 0; f < sizeof forms / sizeof forms[0]; f++) {
        snprintf(buf, sizeof buf, forms[f], xs);
        eval_list(buf, NPT, a_re, a_im);
        for (int k = 0; k < NPT; k++)
            ASSERT_MSG(same_bits(a_re[k], s_re[k]) && a_im[k] == 0.0,
                       "%s: element %d %.17g != scalar %.17g", buf, k + 1, a_re[k], s_re[k]);
    }
}

void test_plog_surfaces_complex() {
    static const char* zs = "{Complex[-0.36`, 0.01`], Complex[0.1`, -0.2`], "
                            "Complex[1.`, 1.`], Complex[-2.`, 0.5`], Complex[30.`, -40.`], "
                            "Complex[-0.5`, 0.`], Complex[3.`, 0.`], Complex[-1.`*^5, 1.`]}";
    double s_re[NPT], s_im[NPT], a_re[NPT], a_im[NPT];
    char buf[640];
    for (int k = 0; k < NPT; k++) {
        snprintf(buf, sizeof buf, "ProductLog[%s[[%d]]]", zs, k + 1);
        eval_num(buf, &s_re[k], &s_im[k]);
    }
    const char* forms[] = {
        "ProductLog[%s]",
        "Normal[ProductLog[NDArray[%s]]]",
        "Compile[{{z, _Complex}}, ProductLog[z]] /@ %s",
        "FromPackedArray[Compile[{{v, _Complex, 1}}, ProductLog[v]][%s]]",
    };
    for (size_t f = 0; f < sizeof forms / sizeof forms[0]; f++) {
        snprintf(buf, sizeof buf, forms[f], zs);
        eval_list(buf, NPT, a_re, a_im);
        for (int k = 0; k < NPT; k++)
            ASSERT_MSG(same_bits(a_re[k], s_re[k]) && same_bits(a_im[k], s_im[k]),
                       "%s: element %d differs from scalar", buf, k + 1);
    }
}

int main() {
    symtab_init();
    core_init();
    TEST(test_plog_real_spot);
    TEST(test_plog_complex_spot);
    TEST(test_plog_branch_point);
    TEST(test_plog_exact_and_declines);
    TEST(test_plog_surfaces_real);
    TEST(test_plog_surfaces_complex);
    printf("All productlog_machine tests passed!\n");
    return 0;
}
