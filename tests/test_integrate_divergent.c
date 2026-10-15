/* test_integrate_divergent.c
 *
 * An ordinary (non-principal-value) definite integral whose integrand has a
 * NON-INTEGRABLE pole inside the closed interval DIVERGES and must be left
 * UNEVALUATED -- never reduced to a spurious finite value.  The regression that
 * motivated these tests was Integrate[1/x, {x, -1, 1}] returning 0: several
 * "trick" stages (King reflection, origin symmetry, differentiation under the
 * integral sign, and the Newton-Leibniz telescoping) applied an identity valid
 * only for a convergent integrand without first checking integrability at the
 * interval's poles, so the odd pole's two infinite one-sided pieces cancelled to
 * the Cauchy principal value (0) even though PrincipalValue -> True was NOT given.
 *
 * The guards live in:
 *   - integrate_king.c          (king_interior_pole)
 *   - integrate_symmetry.c      (pole_at_center)
 *   - integrate_newton_leibniz.c(nl_nonintegrable_at, interior + endpoint)
 *   - integrate_diffunderint.c  (finite_interior_pole)
 * with the integrability test Limit[(x - p) f, x -> p] != 0 (a pole of order >= 1).
 *
 * What must STILL work, and is tested here:
 *   - convergent integrals (including odd -> 0 and even -> 2*half by symmetry),
 *   - integrable endpoint singularities (Log[x], 1/Sqrt[x]) and removable ones
 *     (Sin[x]/x),
 *   - the genuine Cauchy principal value under PrincipalValue -> True.
 */

#include "core.h"
#include "test_utils.h"
#include "expr.h"
#include "eval.h"
#include "parse.h"
#include "symtab.h"

#include <stdio.h>
#include <string.h>

/* True if e contains an Integrate head anywhere -- i.e. the integral was left
 * unevaluated rather than reduced to a value. */
static bool has_integrate_head(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if (e->data.function.head && e->data.function.head->type == EXPR_SYMBOL &&
        strcmp(e->data.function.head->data.symbol.name, "Integrate") == 0)
        return true;
    if (has_integrate_head(e->data.function.head)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (has_integrate_head(e->data.function.args[i])) return true;
    return false;
}

/* The integral must stay UNEVALUATED (divergent, no PrincipalValue): the result
 * still carries an Integrate head and is free of a spurious finite value. */
static void check_unevaluated(const char* input) {
    Expr* p = parse_expression(input);
    ASSERT(p != NULL);
    Expr* r = evaluate(p);
    bool uneval = has_integrate_head(r);
    if (!uneval) {
        char* s = expr_to_string(r);
        fprintf(stderr, "FAIL: %s\n  expected: unevaluated (divergent)\n  actual:   %s\n",
                input, s);
        free(s);
    }
    ASSERT(uneval);
    expr_free(p);
    expr_free(r);
}

/* The integral must reduce to exactly `expected` (canonical print form). */
static void check_eq(const char* input, const char* expected) {
    Expr* p = parse_expression(input);
    ASSERT(p != NULL);
    Expr* r = evaluate(p);
    char* s = expr_to_string(r);
    if (strcmp(s, expected) != 0)
        fprintf(stderr, "FAIL: %s\n  expected: %s\n  actual:   %s\n", input, expected, s);
    ASSERT_STR_EQ(s, expected);
    free(s);
    expr_free(p);
    expr_free(r);
}

/* -------------------------------------------------------------------------
 * Non-integrable pole strictly inside the interval -> divergent -> unevaluated.
 * ---------------------------------------------------------------------- */
static void test_divergent_interior_pole(void) {
    check_unevaluated("Integrate[1/x, {x, -1, 1}]");          /* the motivating bug */
    check_unevaluated("Integrate[1/x^3, {x, -1, 1}]");        /* odd order-3 pole */
    check_unevaluated("Integrate[1/x^2, {x, -1, 1}]");        /* even order-2 pole */
    check_unevaluated("Integrate[1/x, {x, -2, 3}]");          /* non-symmetric interval */
    check_unevaluated("Integrate[1/(x - 1/2), {x, 0, 1}]");   /* off-centre simple pole */
    check_unevaluated("Integrate[1/(x^2 - 1), {x, 0, 2}]");   /* pole at an interior x = 1 */
}

/* -------------------------------------------------------------------------
 * Non-integrable pole AT a finite endpoint -> divergent -> unevaluated.
 * ---------------------------------------------------------------------- */
static void test_divergent_boundary_pole(void) {
    check_unevaluated("Integrate[1/x, {x, 0, 1}]");           /* boundary pole at 0 */
    check_unevaluated("Integrate[Tan[x], {x, 0, Pi}]");       /* Tan pole at Pi/2 */
}

/* -------------------------------------------------------------------------
 * Convergent integrals are unaffected -- including odd -> 0 and even reductions.
 * ---------------------------------------------------------------------- */
static void test_convergent_unaffected(void) {
    check_eq("Integrate[x^2, {x, -1, 1}]", "2/3");            /* even, no pole */
    check_eq("Integrate[x^3, {x, -2, 2}]", "0");              /* odd, no pole -> 0 */
    check_eq("Integrate[1/(x^2 + 1), {x, -1, 1}]", "1/2 Pi"); /* even, no REAL pole */
    check_eq("Integrate[x/(x^2 + 1), {x, -1, 1}]", "0");      /* odd, no REAL pole -> 0 */
}

/* -------------------------------------------------------------------------
 * Integrable / removable singularities still evaluate (the guard uses the
 * integrability test Limit[(x-p) f] = 0, not mere denominator vanishing).
 * ---------------------------------------------------------------------- */
static void test_integrable_singularity(void) {
    check_eq("Integrate[Log[x], {x, 0, 1}]", "-1");           /* Log singular, integrable */
    check_eq("Integrate[1/Sqrt[x], {x, 0, 1}]", "2");         /* x^(-1/2), integrable */
    /* Sin[x]/x has a removable singularity at 0; value 2 SinIntegral[1]. */
    check_eq("Chop[N[Integrate[Sin[x]/x, {x, -1, 1}] - 2 SinIntegral[1]]]", "0");
}

/* -------------------------------------------------------------------------
 * The genuine Cauchy principal value is still returned under PrincipalValue.
 * ---------------------------------------------------------------------- */
static void test_principal_value_still_works(void) {
    check_eq("Integrate[1/x, {x, -1, 1}, PrincipalValue -> True]", "0");
    check_eq("Integrate[1/(x - 1), {x, 0, 3}, PrincipalValue -> True]", "Log[2]");
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_divergent_interior_pole);
    TEST(test_divergent_boundary_pole);
    TEST(test_convergent_unaffected);
    TEST(test_integrable_singularity);
    TEST(test_principal_value_still_works);

    printf("All integrate_divergent tests passed.\n");
    return 0;
}
