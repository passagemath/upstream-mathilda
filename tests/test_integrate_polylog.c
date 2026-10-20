/* test_integrate_polylog.c
 *
 * The PolyLogLadder Integrate stage (src/internal/PolyLogLadder.m, wired in
 * src/calculus/integrate.c) and the pieces it relies on:
 *   - log-kernel ladder: R(x) Log[a x + b]^k with poles anywhere;
 *   - trig route: polynomial x rational trig/hyperbolic -> exponential ladder;
 *   - E^(I q Pi) for integer / half-integer q (power.c), PolyLog[2, +-I]
 *     (polylog.c), the Newton-Leibniz real-form finishing pass and its radical
 *     endpoint pole test (integrate_newton_leibniz.c).
 *
 * Antiderivatives are checked by differentiation at a numeric point (their
 * surface form is not pinned); definite values against independently known
 * closed forms.
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

/* D[Integrate[f, x], x] == f at x = 37/100 (and the stage closed it). */
static void check_antideriv(const char* f) {
    char buf[1024];
    snprintf(buf, sizeof buf,
             "Module[{F = Integrate[%s, x]}, If[!FreeQ[F, Integrate], \"open\", "
             "Chop[N[(D[F, x] - (%s)) /. x -> 37/100, 20], 10^-12]]]", f, f);
    check_eq(buf, "0");
}

static void test_log_ladder(void) {
    check_antideriv("Log[x]/(1 + x^2)");
    check_antideriv("Log[1 + x]/(1 + x^2)");
    check_antideriv("Log[x]/(x - I)");
    check_antideriv("Log[x]^2/(1 + x)");
    check_antideriv("Log[x]^3/(1 - x)");
    check_antideriv("x Log[x]/(1 + x)^2");
    check_antideriv("Log[x]/(1 + x)^3");
    check_antideriv("(1 + x^2) Log[2 x + 1]^2/(x^2 - 3)");
    check_antideriv("Log[x]/(x^3 + 1) + Log[1 - x]/x");
    /* symbolic parameter */
    check_eq("Module[{F = Integrate[Log[1 + t x]/(1 + x^2), x]}, "
             "Chop[N[(D[F, x] - Log[1 + t x]/(1 + x^2)) /. {x -> 37/100, t -> 1/3}, 20], 10^-12]]", "0");
}

static void test_trig_route(void) {
    check_antideriv("x^2 Csc[x]");
    check_antideriv("x Cot[x]");
    check_antideriv("x Tan[x]");
    check_antideriv("x^2 Sec[x]^2");
    check_antideriv("x^3 Csc[2 x]");
    check_antideriv("x Sin[x]/(Cos[x] + Sin[x])");
    check_antideriv("x Csch[x]");
    check_antideriv("x/(1 + Cosh[x])");
}

static void test_log_trig(void) {
    check_antideriv("Log[Sin[x]]");
    check_antideriv("x Log[Sin[x]]");
    check_antideriv("Log[Cot[x]]");
    check_antideriv("x^2 Log[2 Cos[x]]");
    check_antideriv("Log[Cosh[x]]");
    check_antideriv("Log[Tan[x]]");
    check_eq("Integrate[x Log[Sin[x]], {x, 0, Pi}]", "-1/2 Log[2] Pi^2");
    check_eq("Chop[N[Integrate[Log[Sin[x]], {x, 0, Pi/4}] + (2 Catalan + Pi Log[2])/4]]", "0");
    check_eq("Chop[N[Integrate[Log[Cot[x/2]^2], {x, 0, Pi/2}] - 4 Catalan]]", "0");
}

static void test_definite(void) {
    check_eq("Integrate[Log[x]/(1 + x^2), {x, 0, 1}]", "-Catalan");
    check_eq("Integrate[Log[x]^2/(1 + x), {x, 0, 1}]", "3/2 Zeta[3]");
    check_eq("Integrate[x Csc[x], {x, 0, Pi/2}]", "2 Catalan");
    check_eq("Chop[N[Integrate[Log[1 + x]/(1 + x^2), {x, 0, 1}] - Pi Log[2]/8]]", "0");
    /* radical endpoint: no false divergence */
    check_eq("Chop[N[Integrate[Log[x]/(1 + x^2), {x, 0, 2 - Sqrt[3]}] + 0.6106437294514793]]", "0");
    /* Values at roots of unity folded by Integrate's final pass: conjugate
     * pairs of PolyLog / Log[1 + w] cancel their imaginary parts exactly.
     * The pre-fold Newton-Leibniz values are checked to 30 digits too. */
    check_eq("Integrate[ArcTan[x]/x, {x, 0, 1}]", "Catalan");
    check_eq("Integrate[x Log[2 Cos[x]], {x, 0, Pi/2}]", "-7/16 Zeta[3]");
    check_eq("Abs[N[PolyLog[3, I] + PolyLog[3, -I] - Zeta[3]/4 + 7 Zeta[3]/16, 40]] < 10^-30", "True");
    check_eq("Abs[N[Integrate[ArcTan[x]/x, {x, 0, 1}, Method -> \"NewtonLeibniz\"] - Catalan, 40]]"
             " < 10^-30", "True");
    /* a lone complex Log value is not "simpler" expanded: left as it was */
    check_eq("Integrate[1/(x + I), {x, 0, 1}] === Log[1 + I] - I Pi/2", "True");
}

static void test_support(void) {
    check_eq("{E^(I Pi), E^(-I Pi/2), E^(2 I Pi), E^(3 I Pi/2)}", "{-1, -I, 1, -I}");
    check_eq("E^(I Pi/4)", "E^((1/4*I) Pi)");
    check_eq("PolyLog[2, I]", "I Catalan - 1/48 Pi^2");
    check_eq("PolyLog[2, -I]", "-I Catalan - 1/48 Pi^2");
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_log_ladder);
    TEST(test_trig_route);
    TEST(test_log_trig);
    TEST(test_definite);
    TEST(test_support);

    printf("All Integrate PolyLogLadder tests passed!\n");
    return 0;
}
