/* test_integrate_bessel_hyper.c
 *
 * The BesselHypergeometric indefinite-integration stage
 * (src/internal/BesselHyper.m, wired in src/calculus/integrate.c):
 *   - the series route: one Bessel factor, or J_mu J_nu / I_mu I_nu (the
 *     Chu-Vandermonde 2F3), integrated term by term into a HypergeometricPFQ,
 *     including the Log[x] split when the integrated exponent hits zero;
 *   - Y / K of non-integer order via Y = (J Cos - J_-nu)/Sin, K likewise;
 *   - integer and symbolic Y / K orders via Mellin-Barnes MeijerG forms;
 *   - integer order differences reduced to diagonal products by the
 *     recurrences and the Wronskian, and Int D_n/x by telescoping;
 *   - declines: different arguments, mixed ordinary/modified families, three
 *     factors;
 *   - the cascade: Holonomic's Bessel-product forms still win.
 *
 * Antiderivatives are checked by differentiation at three numeric points;
 * their surface form is not pinned except where noted.
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

/* Integrate[f, x] closes (pinned to the stage) and differentiates back to f
 * at x = 37/100, 13/10, 27/10 with nu = 7/4, mu = 5/6, n = 3. */
static void check_antideriv(const char* f) {
    char buf[1024];
    snprintf(buf, sizeof buf,
             "Module[{F = Integrate[%s, x, Method -> \"BesselHypergeometric\"]}, "
             "If[!FreeQ[F, Integrate], \"open\", "
             "Chop[Max[Table[Abs[N[(D[F, x] - (%s)) /. {x -> x0, nu -> 7/4, mu -> 5/6, n -> 3}, 30]], "
             "{x0, {37/100, 13/10, 27/10}}]], 10^-20]]]",
             f, f);
    check_eq(buf, "0");
}

static void test_series_route(void) {
    check_antideriv("BesselJ[1/3, x]^2");
    check_antideriv("BesselJ[0, x]^2");
    check_antideriv("BesselJ[nu, x]^2");
    check_antideriv("BesselJ[1/3, x] BesselJ[2/3, x]");
    check_antideriv("BesselI[1/3, x] BesselI[5/4, x]");
    check_antideriv("BesselJ[1/3, x]");
    check_antideriv("Sqrt[x] BesselJ[0, x]");
    check_antideriv("BesselJ[1/3, 2 x]^2");
    check_antideriv("BesselJ[nu, x^2]");
    check_antideriv("x^(1/3) BesselJ[1/3, 3 x^(2/3)]");
    check_eq("Integrate[BesselJ[0, x]^2, x]",
             "x HypergeometricPFQ[{1/2, 1/2}, {1, 1, 3/2}, -x^2]");
}

static void test_log_case(void) {
    check_antideriv("BesselJ[0, x]/x");
    check_antideriv("BesselJ[1, x]/x^2");
    check_antideriv("BesselI[2, x]/x^5");
    check_antideriv("BesselJ[1/3, x] BesselJ[-1/3, x]/x");
}

static void test_noninteger_y_k(void) {
    check_antideriv("BesselK[1/3, x]");
    check_antideriv("BesselK[1/3, x]^2");
    check_antideriv("BesselJ[1/3, x] BesselY[1/3, x]");
    check_antideriv("BesselJ[1/3, x] BesselY[0, x]");
    check_antideriv("x^2 BesselI[1/3, x] BesselK[0, x]");
}

static void test_meijer_route(void) {
    check_antideriv("BesselK[0, x]");
    check_antideriv("BesselY[0, x]");
    check_antideriv("BesselY[nu, x]");
    check_antideriv("BesselK[0, x]^2");
    check_antideriv("BesselJ[0, x] BesselY[0, x]");
    check_antideriv("BesselI[0, x] BesselK[0, x]");
    check_antideriv("BesselY[1, x]^2");
    check_antideriv("BesselY[nu, x]^2");
    check_antideriv("BesselK[nu, x] BesselK[mu, x]");
}

static void test_integer_difference(void) {
    check_antideriv("BesselJ[0, x] BesselY[1, x]");
    check_antideriv("BesselJ[1, x] BesselY[2, x]");
    check_antideriv("BesselI[0, x] BesselK[1, x]");
    check_antideriv("BesselY[0, x] BesselY[2, x]");
    check_antideriv("BesselY[1, x] BesselY[2, x]");
    check_antideriv("BesselJ[3, x] BesselY[3, x]/x");
    check_antideriv("BesselK[2, x]^2/x");
    check_antideriv("BesselI[1, x] BesselK[1, x]/x");
}

static void test_declines(void) {
    /* different arguments: an Appell F4 double series */
    check_eq("Head[Integrate[BesselJ[0, x] BesselJ[0, 2 x], x, Method -> \"BesselHypergeometric\"]]", "Integrate");
    /* mixed ordinary / modified families */
    check_eq("Head[Integrate[BesselJ[0, x] BesselI[1, x], x, Method -> \"BesselHypergeometric\"]]", "Integrate");
    /* three factors */
    check_eq("Head[Integrate[BesselJ[0, x]^3, x, Method -> \"BesselHypergeometric\"]]", "Integrate");
    /* K_0^2/x: the log-degenerate Mellin-Barnes case with nothing to telescope */
    check_eq("Head[Integrate[BesselK[0, x]^2/x, x, Method -> \"BesselHypergeometric\"]]", "Integrate");
}

static void test_cascade(void) {
    /* Holonomic runs first: closed Bessel-product forms are kept */
    check_eq("FreeQ[Integrate[x^2 BesselK[1/3, x]^3, x], HypergeometricPFQ | MeijerG]", "True");
    check_eq("FreeQ[Integrate[x BesselJ[1/3, x]^2, x], HypergeometricPFQ | MeijerG]", "True");
    check_eq("FreeQ[Integrate[BesselJ[1/3, x]^2, x], Integrate]", "True");
    /* the definite integral through the antiderivative */
    check_eq("Chop[N[Integrate[BesselJ[1/3, x]^2, {x, 0, 1}] - NIntegrate[BesselJ[1/3, x]^2, {x, 0, 1}]], 10^-10]", "0");
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_series_route);
    TEST(test_log_case);
    TEST(test_noninteger_y_k);
    TEST(test_meijer_route);
    TEST(test_integer_difference);
    TEST(test_declines);
    TEST(test_cascade);

    printf("All Integrate BesselHypergeometric tests passed!\n");
    return 0;
}
