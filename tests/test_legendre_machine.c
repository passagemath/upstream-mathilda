/* Machine-precision LegendreP[n, x] / LegendreQ[n, x] (integer degree): the
 * shared double helpers of src/special_functions/sf_machine_plog.c behind the
 * scalar builtin, the NDArray binary kernels and Compile[].
 *
 * Covers: spot values against mpmath 1.3 at 30 digits (computed for this
 * file); LegendreQ type 1 off the cut and on the real axis beyond +-1 (the
 * closed form's branch: Log[1 - x] = log(x - 1) + i Pi for x > 1); the small-|z|
 * Q_0 (no cancellation); negative degree P_{-n-1} = P_n; declines (Q at x = 1,
 * non-integer degree, exact input stays exact); and the bit-match
 *     scalar == unpacked List == packed List == NDArray == Compile.
 */
#include "eval.h"
#include "parse.h"
#include "expr.h"
#include "symtab.h"
#include "core.h"
#include "arithmetic.h"
#include "test_utils.h"
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

/* ---- spot values ---------------------------------------------------- */

void test_leg_p_spot() {
    assert_rel("LegendreP[3, 0.3`]", -0.382499999999999990840660046842, 0, 4e-16);
    assert_rel("LegendreP[5, -0.6`]", 0.152639999999999945110573662532, 0, 8e-16);
    assert_rel("LegendreP[3, 2.5`]", 35.3125, 0, 4e-16);
    assert_rel("LegendreP[-4, 2.5`]", 35.3125, 0, 4e-16);         /* P_{-4} = P_3 */
    assert_rel("LegendreP[3, Complex[0.5`, 0.25`]]", -0.671875, 0.0546875, 4e-16);
    assert_rel("LegendreP[10, Complex[-1.2`, 0.7`]]",
               3088.19030931748874062720487263, -745.975390905215263646180161893, 1e-15);
    assert_rel("LegendreP[3., 0.3`]", -0.382499999999999990840660046842, 0, 4e-16);
}

void test_leg_q_spot() {
    assert_rel("LegendreQ[2, 0.5`]", -0.818663268041756855712202827308, 0, 4e-16);
    assert_rel("LegendreQ[3, -0.3`]", 0.323275418058976459652773331422, 0, 4e-16);
    assert_rel("LegendreQ[2, 1.5`]",
               0.0635669991240192884885915415129, -4.51603943953532778029004986346, 4e-16);
    assert_rel("LegendreQ[2, -1.5`]",
               -0.0635669991240192884885915415129, 4.51603943953532778029004986346, 4e-16);
    assert_rel("LegendreQ[2, Complex[0.5`, 0.25`]]",
               -0.977355241214384839032910803681, -0.256135843806988847857026150327, 4e-16);
    assert_rel("LegendreQ[4, Complex[2.`, -3.`]]",
               895.353917929091492522336129846, 787.754319679240479570819765857, 1e-15);
    /* Q_0 near the origin: Kahan's formula keeps full relative accuracy. */
    assert_rel("LegendreQ[0, Complex[1.`*^-8, 1.`*^-9]]",
               1.00000000000000000053e-8, 1.00000000000000001619e-9, 4e-16);
}

void test_leg_exact_and_declines() {
    assert_eval_eq("LegendreP[3, x]", "-3/2 x + 5/2 x^3", 0);
    assert_eval_eq("LegendreP[3, 1/2]", "-7/16", 0);
    assert_eval_eq("LegendreP[3., 1]", "1", 0);
    assert_eval_eq("LegendreQ[2, 1.`]", "Indeterminate", 0);
    /* Non-integer degree: not the machine kernel's; the 2F1 path answers. */
    double r, i;
    eval_num("LegendreP[2.5`, 0.3`]", &r, &i);
    ASSERT_MSG(fabs(r - (-0.463661013457223193)) < 1e-14, "P_2.5(0.3) = %.17g", r);
    double rr, ri;
    ASSERT(!sf_machine_legendre_q_c(2.0, 1.0, 0.0, &rr, &ri));
    ASSERT(!sf_machine_legendre_p_c(2.5, 0.3, 0.0, &rr, &ri));
}

/* ---- five surfaces agree bit-for-bit ----------------------------------- */

#define NPT 7
static void surfaces(const char* head, const char* deg, const char* xs, bool cplx) {
    double s_re[NPT], s_im[NPT], a_re[NPT], a_im[NPT];
    char buf[700];
    for (int k = 0; k < NPT; k++) {
        snprintf(buf, sizeof buf, "%s[%s, %s[[%d]]]", head, deg, xs, k + 1);
        eval_num(buf, &s_re[k], &s_im[k]);
    }
    const char* rforms[] = {
        "%s[%s, FromPackedArray[%s]]",
        "FromPackedArray[%s[%s, ToPackedArray[%s]]]",
        "Normal[%s[%s, NDArray[%s]]]",
        "Compile[{{x, _Real}}, %s[%s, x]] /@ %s",
        "FromPackedArray[Compile[{{v, _Real, 1}}, %s[%s, v]][%s]]",
    };
    const char* cforms[] = {
        "%s[%s, %s]",
        "Normal[%s[%s, NDArray[%s]]]",
        "Compile[{{z, _Complex}}, %s[%s, z]] /@ %s",
        "FromPackedArray[Compile[{{v, _Complex, 1}}, %s[%s, v]][%s]]",
    };
    const char** forms = cplx ? cforms : rforms;
    size_t nf = cplx ? 4 : 5;
    for (size_t f = 0; f < nf; f++) {
        snprintf(buf, sizeof buf, forms[f], head, deg, xs);
        eval_list(buf, NPT, a_re, a_im);
        for (int k = 0; k < NPT; k++)
            ASSERT_MSG(same_bits(a_re[k], s_re[k]) && same_bits(a_im[k], s_im[k]),
                       "%s: element %d (%.17g, %.17g) != scalar (%.17g, %.17g)",
                       buf, k + 1, a_re[k], a_im[k], s_re[k], s_im[k]);
    }
}

void test_leg_surfaces() {
    surfaces("LegendreP", "3", "{-0.95`, -0.7745`, -0.2`, 0.`, 0.4`, 0.91`, 7.5`}", false);
    surfaces("LegendreQ", "2", "{-0.9`, -0.57`, -0.1`, 0.`, 0.3`, 0.66`, 0.999`}", false);
    surfaces("LegendreP", "3", "{Complex[0.5`, 0.25`], Complex[-0.9`, 0.9`], Complex[3.`, -2.`], "
             "Complex[0.`, 1.`], Complex[-0.7`, -0.1`], Complex[0.1`, 0.`], Complex[8.`, 8.`]}", true);
    surfaces("LegendreQ", "2", "{Complex[0.5`, 0.25`], Complex[-0.8`, 0.8`], Complex[3.`, -2.`], "
             "Complex[0.`, 1.`], Complex[-0.7`, -0.1`], Complex[1.5`, 0.5`], Complex[-4.`, 0.2`]}", true);
}

int main() {
    symtab_init();
    core_init();
    TEST(test_leg_p_spot);
    TEST(test_leg_q_spot);
    TEST(test_leg_exact_and_declines);
    TEST(test_leg_surfaces);
    printf("All legendre_machine tests passed!\n");
    return 0;
}
