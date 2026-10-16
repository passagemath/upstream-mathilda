/* Tests for the machine kernels of the COMPLETE elliptic integrals
 * EllipticK[m] and EllipticE[m] (src/special_functions/sf_machine_ellk.c).
 *
 * Covers: real spot values across the regimes (0 <= m < 1, m near 1, m < 0
 * through the imaginary-modulus transformation, huge negative m); machine
 * complex values (complex AGM) including the cut m > 1, where a machine REAL
 * argument must give a machine Complex; the declines (the pole K(1), m = 1
 * complex) that fall through to the exact/Arb path; machine-in => machine-out
 * (no 53-bit MPFR leaks from the Arb bridge); exact input staying exact; and
 * the five-surface bit-match: scalar == unpacked List == packed List ==
 * NDArray == Compile, for real and complex arguments.
 *
 * Reference values: mpmath 1.3 (mp.ellipk / mp.ellipe) at 25-50 digits, for the
 * exact double argument, computed for this file (not recalled).
 */

#include "eval.h"
#include "parse.h"
#include "expr.h"
#include "symtab.h"
#include "attr.h"
#include "core.h"
#include "print.h"
#include "test_utils.h"
#include "arithmetic.h"
#include "special_functions/sf_machine_ellk.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static Expr* ev(const char* input) {
    Expr* e = parse_expression(input);
    ASSERT(e != NULL);
    Expr* r = evaluate(e);
    expr_free(e);
    return r;
}

/* Relative closeness within `ulps` units of 2^-52. */
static void assert_real(const char* input, double expected, double ulps) {
    Expr* r = ev(input);
    ASSERT_MSG(r->type == EXPR_REAL, "%s: expected a machine Real", input);
    double v = r->data.real;
    expr_free(r);
    ASSERT_MSG(fabs(v - expected) <= ulps * 2.220446049250313e-16 * fabs(expected),
               "%s: expected %.17g, got %.17g", input, expected, v);
}

static void assert_cplx(const char* input, double er, double ei, double ulps) {
    Expr* r = ev(input);
    Expr *re, *im;
    ASSERT_MSG(is_complex(r, &re, &im), "%s: expected Complex", input);
    ASSERT_MSG(re->type == EXPR_REAL && im->type == EXPR_REAL,
               "%s: expected MACHINE Complex parts", input);
    double dr = re->data.real - er, di = im->data.real - ei;
    double err = sqrt(dr * dr + di * di), mag = sqrt(er * er + ei * ei);
    expr_free(r);
    ASSERT_MSG(err <= ulps * 2.220446049250313e-16 * mag,
               "%s: expected %.17g + %.17g I", input, er, ei);
}

static void assert_true(const char* input) {
    Expr* r = ev(input);
    bool ok = r->type == EXPR_SYMBOL && strcmp(r->data.symbol.name, "True") == 0;
    expr_free(r);
    ASSERT_MSG(ok, "%s: expected True", input);
}

/* ---- real kernel -------------------------------------------------------- */

void test_ellk_real_values() {
    assert_real("EllipticK[0.5]", 1.8540746773013719184, 2);
    assert_real("EllipticE[0.5]", 1.3506438810476755025, 2);
    assert_real("EllipticK[0.9]", 2.5780921133481732927, 2);
    assert_real("EllipticE[0.9]", 1.1047747327040733079, 2);
    assert_real("EllipticK[0.999999]", 8.2940514636010622019, 2);
    assert_real("EllipticE[0.999999]", 1.000003897026172166, 2);
    assert_real("EllipticK[0.9999999999999999`]", 19.754694645958441839, 2);
    assert_real("EllipticE[0.9999999999999999`]", 1.0000000000000010689, 2);
    assert_real("EllipticK[1.*^-300]", 1.5707963267948966192, 1);
    assert_real("EllipticK[0.]", 1.5707963267948966192, 1);
    assert_real("EllipticE[0.]", 1.5707963267948966192, 1);
    assert_real("EllipticE[1.]", 1.0, 0);
    /* m < 0: imaginary-modulus transformation */
    assert_real("EllipticK[-0.5]", 1.4157372084259561989, 2);
    assert_real("EllipticE[-0.5]", 1.751771275694817862, 2);
    assert_real("EllipticK[-3.]", 1.0782578237498216177, 2);
    assert_real("EllipticE[-3.]", 2.4221120551369190496, 2);
    assert_real("EllipticK[-1.*^6]", 0.0082940478165906199329, 3);
    assert_real("EllipticE[-1.*^6]", 1000.0043970243485481, 3);
    assert_real("EllipticK[-1.*^300]", 3.4677405831022673414e-148, 3);
    assert_real("EllipticE[-1.*^300]", 1.0000000000000000263e+150, 3);

    /* the C entry points' contract */
    double v;
    ASSERT(!sf_machine_ellk_k(1.0, &v));
    ASSERT(!sf_machine_ellk_k(2.0, &v));
    ASSERT(!sf_machine_ellk_e(1.5, &v));
    ASSERT(!sf_machine_ellk_k(NAN, &v));
    ASSERT(!sf_machine_ellk_k(-INFINITY, &v));
    ASSERT(sf_machine_ellk_e(1.0, &v) && v == 1.0);
}

/* ---- complex kernel ----------------------------------------------------- */

void test_ellk_complex_values() {
    assert_cplx("EllipticK[Complex[0.3, 0.2]]", 1.6961165546095899446, 0.1134950278395376944, 4);
    assert_cplx("EllipticE[Complex[0.3, 0.2]]", 1.4498607813767230317, -0.088990096584979846297, 4);
    assert_cplx("EllipticK[Complex[-0.5, 0.7]]", 1.3741460694006835397, 0.16015821224245015307, 4);
    assert_cplx("EllipticE[Complex[-0.5, 0.7]]", 1.772346583376868449, -0.23121783571928993103, 4);
    assert_cplx("EllipticK[Complex[2., 0.1]]", 1.3260612688769483745, 1.2629470674171812187, 4);
    assert_cplx("EllipticE[Complex[2., 0.1]]", 0.64598096611236153401, -0.61764424511657842121, 4);
    assert_cplx("EllipticK[Complex[5., -3.]]", 0.80441615395715038693, -0.75798242705210756555, 4);
    assert_cplx("EllipticE[Complex[5., -3.]]", 1.0664019827122897894, 1.846546074710714728, 4);
    /* On the cut, from below (Mathematica / mpmath): K(2) = 1.311 - 1.311 I.
     * A machine REAL m > 1 gives a machine Complex. */
    assert_cplx("EllipticK[2.]", 1.3110287771460599052, -1.3110287771460599052, 4);
    assert_cplx("EllipticE[2.]", 0.59907011736779610372, 0.59907011736779610372, 4);
    assert_cplx("EllipticK[Complex[2., 0.]]", 1.3110287771460599052, -1.3110287771460599052, 4);
    /* machine in => machine out (was a 53-bit MPFR from the Arb bridge) */
    assert_eval_eq("Precision[EllipticK[Complex[0.6, 0.5]]]", "MachinePrecision", 0);
    assert_eval_eq("Precision[EllipticE[Complex[0.6, 0.5]]]", "MachinePrecision", 0);

    double r, i;
    ASSERT(!sf_machine_ellk_k_c(1.0, 0.0, &r, &i));        /* the pole */
    ASSERT(sf_machine_ellk_e_c(1.0, 0.0, &r, &i) && r == 1.0 && i == 0.0);
    ASSERT(!sf_machine_ellk_k_c(NAN, 0.0, &r, &i));
    ASSERT(!sf_machine_ellk_k_c(1e305, 1.0, &r, &i));
}

/* ---- declines, exact input, special values ------------------------------ */

void test_ellk_special() {
    assert_eval_eq("EllipticK[1.]", "ComplexInfinity", 0);
    assert_eval_eq("EllipticK[1]", "ComplexInfinity", 0);
    assert_eval_eq("EllipticK[0]", "Times[Rational[1, 2], Pi]", 1);
    assert_eval_eq("EllipticE[1]", "1", 0);
    assert_eval_eq("EllipticK[1/3]", "EllipticK[1/3]", 0);
    assert_eval_eq("EllipticE[m]", "EllipticE[m]", 0);
    /* arbitrary precision is untouched by the machine path */
    assert_true("Precision[EllipticK[0.5`30]] > 29");
    assert_true("Precision[EllipticE[SetPrecision[2, 30]]] > 29");
}

/* ---- five surfaces: scalar == List == packed == NDArray == Compile ------- */

void test_ellk_surfaces() {
    const char* defs =
        "xs$ = {0.3, -2.5, 0.95, 0.999999, 0., -1.*^6, 0.62}; "
        "zs$ = {Complex[0.3, 0.2], Complex[-0.5, 0.7], Complex[2., 0.1], "
        "Complex[-7., -3.], Complex[0.9, 0.]};";
    Expr* r = ev(defs); expr_free(r);
    const char* heads[2] = { "EllipticK", "EllipticE" };
    char buf[1024];
    for (int h = 0; h < 2; h++) {
        const char* f = heads[h];
        snprintf(buf, sizeof buf, "s$ = Map[%s, xs$]; "
                 "cr$ = Compile[{{x, _Real}}, %s[x]]; cv$ = Compile[{{v, _Real, 1}}, %s[v]];"
                 "cz$ = Compile[{{z, _Complex}}, %s[z]]; czv$ = Compile[{{v, _Complex, 1}}, %s[v]];",
                 f, f, f, f, f);
        r = ev(buf); expr_free(r);
        snprintf(buf, sizeof buf, "s$ === %s[FromPackedArray[xs$]]", f);         assert_true(buf);
        snprintf(buf, sizeof buf, "s$ === %s[ToPackedArray[xs$]]", f);           assert_true(buf);
        snprintf(buf, sizeof buf, "PackedArrayQ[%s[ToPackedArray[xs$]]]", f);   assert_true(buf);
        snprintf(buf, sizeof buf, "s$ === Normal[%s[NDArray[xs$]]]", f);        assert_true(buf);
        snprintf(buf, sizeof buf, "Head[%s[NDArray[xs$]]] === NDArray", f);      assert_true(buf);
        assert_true("s$ === Map[cr$, xs$]");
        assert_true("s$ === Normal[cv$[xs$]]");
        snprintf(buf, sizeof buf, "c$ = Map[%s, zs$];", f); r = ev(buf); expr_free(r);
        snprintf(buf, sizeof buf, "c$ === %s[zs$]", f);                          assert_true(buf);
        snprintf(buf, sizeof buf, "c$ === Normal[%s[NDArray[zs$]]]", f);         assert_true(buf);
        snprintf(buf, sizeof buf, "Head[%s[NDArray[zs$]]] === NDArray", f);      assert_true(buf);
        assert_true("c$ === Map[cz$, zs$]");
        assert_true("c$ === Normal[czv$[zs$]]");
        snprintf(buf, sizeof buf, "CompileDiagnostics[{{x, _Real}}, %s[x]][[1, 2]]", f);
        assert_true(buf);
        snprintf(buf, sizeof buf, "CompileDiagnostics[{{v, _Real, 1}}, %s[v]][[1, 2]]", f);
        assert_true(buf);
        snprintf(buf, sizeof buf, "CompileDiagnostics[{{z, _Complex}}, %s[z]][[1, 2]]", f);
        assert_true(buf);
    }
    /* a real buffer with an element on the cut (m > 1) falls back, and still
     * agrees with the scalar (complex) answer element by element */
    assert_true("Map[EllipticK, {0.5, 2.}] === Normal[EllipticK[NDArray[{0.5, 2.}]]]");
}

int main() {
    symtab_init();
    core_init();

    TEST(test_ellk_real_values);
    TEST(test_ellk_complex_values);
    TEST(test_ellk_special);
    TEST(test_ellk_surfaces);

    printf("All elliptic machine-kernel tests passed.\n");
    return 0;
}
