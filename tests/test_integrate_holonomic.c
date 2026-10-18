/* test_integrate_holonomic.c
 *
 * The Holonomic indefinite-integration stage (src/internal/Holonomic.m, wired
 * in src/calculus/integrate.c):
 *   - the module ansatz: x^alpha L(x) times monomials in Bessel J/Y/I/K,
 *     Airy and Sin/Cos kernels, solved as one linear system;
 *   - the Wronskian route: W(x) g(y1, y2) with g homogeneous of degree -2,
 *     by t = y2/y1;
 *   - declines: x BesselJ[1, x] (needs Struve H, outside the module).
 *
 * Antiderivatives are checked by differentiation at a numeric point; their
 * surface form is not pinned.
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
 * at x = 37/100 with nu = 7/4, mu = 5/6, n = 3. */
static void check_antideriv(const char* f) {
    char buf[1024];
    snprintf(buf, sizeof buf,
             "Module[{F = Integrate[%s, x, Method -> \"Holonomic\"]}, "
             "If[!FreeQ[F, Integrate], \"open\", "
             "Chop[N[(D[F, x] - (%s)) /. {x -> 37/100, nu -> 7/4, mu -> 5/6, n -> 3}, 20], 10^-12]]]",
             f, f);
    check_eq(buf, "0");
}

static void test_bessel_products(void) {
    check_antideriv("x BesselJ[0, x]^2");
    check_antideriv("BesselJ[1, x]^4/x");               /* Mathematica: open */
    check_antideriv("x^2 BesselJ[1/3, x]^3");           /* Mathematica: open */
    check_antideriv("BesselJ[3, x]^4/x^3");             /* Mathematica: open */
    check_antideriv("BesselJ[nu, x]^2/x^2");
    check_antideriv("BesselJ[nu, x]^2/x^4");            /* needs the D - alpha window */
    check_antideriv("BesselJ[mu, x] BesselY[nu, x]/x^3");
    check_antideriv("x^(1 + mu + nu) BesselJ[mu, x] BesselY[nu, x]");
    /* orders nu0 + k share one family; mu in the order must not be captured */
    check_antideriv("x^(-1 - mu - nu) BesselJ[1 + mu, x] BesselY[1 + nu, x]");
    check_antideriv("x BesselI[0, x] BesselK[0, x]");
    check_antideriv("x^3 BesselI[0, 2 x]");
}

static void test_mixed_kernels(void) {
    check_antideriv("BesselJ[nu, x] Sin[x]/x^(3/2)");
    check_antideriv("BesselJ[nu, x] Cos[x]/x^(3/2)");
    check_antideriv("x AiryAi[x]");
    check_antideriv("x AiryAi[x]^2");
}

static void test_wronskian(void) {
    check_antideriv("1/AiryAi[x]^2");
    check_antideriv("1/(AiryAi[x] AiryBi[x])");
    check_antideriv("1/(AiryAi[x] + AiryBi[x])^2");
    check_antideriv("AiryAi[x]^(-2 - n) AiryBi[x]^n");
    check_antideriv("AiryAi[x]^(n - 1) AiryBi[x]^(n - 1)/(AiryAi[x]^n + AiryBi[x]^n)^2");
    check_antideriv("1/(x BesselJ[nu, x]^2)");
    check_antideriv("1/(x BesselI[0, x] BesselK[0, x])");
    check_eq("Integrate[1/AiryAi[x]^2, x]", "(Pi AiryBi[x])/AiryAi[x]");
}

static void test_declines(void) {
    /* x J_1 needs Struve H: no antiderivative in the module */
    check_eq("Head[Integrate[x BesselJ[1, x], x, Method -> \"Holonomic\"]]", "Integrate");
    /* a pole away from 0 is outside the Laurent ansatz */
    check_eq("Head[Integrate[BesselJ[0, x]/(1 + x^2), x, Method -> \"Holonomic\"]]", "Integrate");
    /* no Bessel / Airy kernel: the stage's gate */
    check_eq("Head[Integrate[x Sin[x]^2 Log[x], x, Method -> \"Holonomic\"]]", "Integrate");
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_bessel_products);
    TEST(test_mixed_kernels);
    TEST(test_wronskian);
    TEST(test_declines);

    printf("All Integrate Holonomic tests passed!\n");
    return 0;
}
