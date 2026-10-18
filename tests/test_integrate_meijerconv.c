/* test_integrate_meijerconv.c
 *
 * The MeijerConvolution definite-integration stage
 * (src/internal/MeijerConvolution.m, wired in src/calculus/integrate.c):
 * C x^alpha K1(x) K2(x) on (0, Infinity) by the Meijer G convolution theorem,
 * kernels BesselJ / BesselY / BesselK, E^(-p x), E^(-p x^2), Cos, Sin,
 * (a^2 + x^2)^-r at arbitrary scales.  Reference values are Mathematica
 * NIntegrate at 30 digits (closed forms where Mathematica gives them).
 */

#include "core.h"
#include "test_utils.h"
#include "expr.h"
#include "eval.h"
#include "parse.h"
#include "symtab.h"

#include <stdio.h>
#include <string.h>

static void check_eq(const char* input, const char* expected) {
    Expr* p = parse_expression(input);
    ASSERT(p != NULL);
    Expr* r = evaluate(p);
    char* s = expr_to_string(r);
    if (strcmp(s, expected) != 0) {
        fprintf(stderr, "FAIL: %s\n  expected: %s\n  actual:   %s\n",
                input, expected, s);
    }
    ASSERT_STR_EQ(s, expected);
    free(s);
    expr_free(p);
    expr_free(r);
}

static void check_value(const char* integral, const char* ref) {
    char buf[1024];
    snprintf(buf, sizeof buf, "Chop[N[(%s) - (%s)], 10^-9]", integral, ref);
    check_eq(buf, "0");
}

static void test_two_kernels(void) {
    check_value("Integrate[x^(3/5) E^(-13/10 x) BesselJ[7/10, 9/10 x], {x, 0, Infinity}]",
                "0.2592479223723994");
    check_value("Integrate[x^2 BesselK[0, x] BesselK[1, 2 x], {x, 0, Infinity}]",
                "0.1413987469155312");
    /* different frequencies, amplitude ~1/x: conditionally convergent; = 0 */
    check_value("Integrate[BesselJ[1, x] BesselY[1, 3 x], {x, 0, Infinity}]", "0");
    /* x J1(x) Y1(3x) has a non-decaying amplitude: not convergent, must decline */
    check_eq("Head[Integrate[x BesselJ[1, x] BesselY[1, 3 x], {x, 0, Infinity}, "
             "Method -> \"MeijerConvolution\"]]", "Integrate");
    check_value("Integrate[x BesselJ[0, x] E^(-x^2/2), {x, 0, Infinity}]", "E^(-1/2)");
    check_value("Integrate[BesselJ[0, 2 x] Cos[x], {x, 0, Infinity}]", "1/Sqrt[3]");
    check_value("Integrate[x BesselJ[0, x]/(1 + x^2), {x, 0, Infinity}]", "BesselK[0, 1]");
}

static void test_symbolic(void) {
    /* Laplace transform of J_nu: a^nu/(Sqrt[a^2+p^2] (p + Sqrt[a^2+p^2])^nu) */
    check_eq("Chop[N[(Integrate[BesselJ[nu, a x] E^(-p x), {x, 0, Infinity}, "
             "Assumptions -> a > 0 && p > 0 && nu > 0] - a^nu/(Sqrt[a^2 + p^2] "
             "(p + Sqrt[a^2 + p^2])^nu)) /. {a -> 2/3, p -> 5/4, nu -> 3/7}], 10^-9]", "0");
}

static void test_declines(void) {
    /* divergent at 0: x^-2 J0 */
    check_eq("Head[Integrate[BesselJ[0, x] BesselJ[0, 2 x]/x^2, {x, 0, Infinity}, "
             "Method -> \"MeijerConvolution\"]]", "Integrate");
    /* equal-frequency product without enough decay */
    check_eq("Head[Integrate[BesselJ[0, x]^2, {x, 0, Infinity}, "
             "Method -> \"MeijerConvolution\"]]", "Integrate");
    /* (x^2 - b^2)^-1 has real poles: the (a + x^2)^-r kernel needs Re[a] > 0,
       and the refuted side condition must decline (it returned a wrong
       finite value for these principal-value integrals before) */
    check_eq("Head[Integrate[Sin[a x]/(x^2 - b^2), {x, 0, Infinity}, "
             "Assumptions -> a > 0 && b > 0, Method -> \"MeijerConvolution\"]]", "Integrate");
    check_eq("Head[Integrate[x Cos[a x]/(x^2 - b^2), {x, 0, Infinity}, "
             "Assumptions -> a > 0 && b > 0, Method -> \"MeijerConvolution\"]]", "Integrate");
    /* ... while Re[a] > 0 still closes */
    check_eq("Chop[N[(Integrate[Cos[a x]/(x^2 + b^2), {x, 0, Infinity}, "
             "Assumptions -> a > 0 && b > 0, Method -> \"MeijerConvolution\"] - "
             "Pi E^(-a b)/(2 b)) /. {a -> 2, b -> 3}, 20], 10^-12]", "0");
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_two_kernels);
    TEST(test_symbolic);
    TEST(test_declines);

    printf("All Integrate MeijerConvolution tests passed!\n");
    return 0;
}
