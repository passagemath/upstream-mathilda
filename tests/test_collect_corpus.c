/*
 * test_collect_corpus.c -- Regression corpus for Collect.
 *
 * Captures the user-supplied corpus from the 2026-05-10 audit
 * (14 cases the user marked as "passing" against MMA + Bug 1 fix +
 * Bug 2 current behaviour). Each case pins Mathilda's canonical
 * output string, so any future change in Collect / Power / Times
 * canonicalisation that perturbs these forms surfaces here as a
 * test failure.
 *
 * Bug 1 -- Collect[a x^(2c) + b x^(2c), x^c]
 *
 *   Pre-fix: returned the input unchanged because decompose_to_bp +
 *   the single-base path could not match Power[x, 2c] (an atomic
 *   bp entry) against the target Power[x, c]. The fix adds a
 *   power-of-power fallback in collect_internal: when the target
 *   is Power[B, e_t] and a term factor decomposes to an atomic
 *   Power[B, e_term], compute k = e_term / e_t (via Cancel) and use
 *   it as the multiplier when k is a positive integer.
 *
 * Bug 2 -- Collect[D[f[Sqrt[x^2+1]], {x,3}], Derivative[_][f][_]]
 *
 *   The target is a pattern. Pre-v0.516 Collect compared keys with
 *   expr_eq only, so any pattern target (Power[_, 1/2], f[_], x^_)
 *   returned the input unchanged. builtin_collect now expands each
 *   pattern variable into the distinct matching subexpressions
 *   (levels 1..Infinity, canonical order), each a collection variable
 *   of its own.
 */

#include "test_utils.h"
#include "eval.h"
#include "parse.h"
#include "expr.h"
#include "symtab.h"
#include "core.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* --- Bug 1: Symbolic-exponent ratio in target ----------------------- */
static void test_bug1_symbolic_exp_ratio(void) {
    assert_eval_eq("Collect[a*x^(2*c) + b*x^(2*c), x^c]",
                   "(a + b) x^(2 c)", 0);
    /* Generalisation -- ratio 3 should also work. */
    assert_eval_eq("Collect[a*x^(3*c) + b*x^(3*c), x^c]",
                   "(a + b) x^(3 c)", 0);
    /* Mixed ratios (k=2 and k=1) should produce two groups. */
    assert_eval_eq("Collect[a*x^c + b*x^(2*c), x^c]",
                   "a x^c + b x^(2 c)", 0);
    /* Non-integer ratio (3/2) -> stays unevaluated against target x^c. */
    assert_eval_eq("Collect[a*x^(3 c/2), x^c]",
                   "a x^(3/2 c)", 0);
}

/* --- Bug 2: Pattern-target Collect --------------------------------- */
static void test_bug2_pattern_target(void) {
    assert_eval_eq(
        "e = D[f[Sqrt[x^2 + 1]], {x, 3}]; Collect[e, Derivative[_][f][_]]",
        "(x^3 Derivative[3][f][Sqrt[1 + x^2]])/(1 + x^2)^(3/2)"
        " + (-(3 x^3)/(1 + x^2)^2 + (3 x)/(1 + x^2)) Derivative[2][f][Sqrt[1 + x^2]]"
        " + (-(3 x)/(1 + x^2)^(3/2) + (3 x^3)/(1 + x^2)^(5/2)) Derivative[1][f][Sqrt[1 + x^2]]",
        0);
    assert_eval_eq(
        "Collect[1 + 3 Sqrt[x^2 + 1] + x Sqrt[x^2 + 1], Power[_, 1/2]]",
        "1 + (3 + x) Sqrt[1 + x^2]", 0);
    /* The Integrate output that surfaced the bug, directly and via //@. */
    assert_eval_eq(
        "u = 10001 + 3124 x^2 - 1408 x^3 + 54 x^4 - 128 x^5 + 20 x^6 + x^8"
        " + 781 Sqrt[-71 - 96 x + 10 x^2 + x^4] + x^6 Sqrt[-71 - 96 x + 10 x^2 + x^4]"
        " - 528 x Sqrt[-71 - 96 x + 10 x^2 + x^4] - 80 x^3 Sqrt[-71 - 96 x + 10 x^2 + x^4]"
        " + 15 x^4 Sqrt[-71 - 96 x + 10 x^2 + x^4] + 27 x^2 Sqrt[-71 - 96 x + 10 x^2 + x^4];"
        " Collect[#, Power[_, 1/2]]& //@ (Log[u]/8)",
        "1/8 Log[10001 + 3124 x^2 - 1408 x^3 + 54 x^4 - 128 x^5 + 20 x^6 + x^8"
        " + (781 - 528 x + 27 x^2 - 80 x^3 + 15 x^4 + x^6) Sqrt[-71 - 96 x + 10 x^2 + x^4]]",
        0);
    assert_eval_eq("Collect[a f[1] + b f[1] + c f[2] + d f[2], f[_]]",
                   "(a + b) f[1] + (c + d) f[2]", 0);
    /* Each match is its own variable, collected in canonical order. */
    assert_eval_eq("Collect[a f[x] + b f[x] + c f[y] + d f[y] + f[x] f[y], f[_]]",
                   "(c + d) f[y] + f[x] (a + b + f[y])", 0);
    assert_eval_eq("Collect[a x^2 + b x^2 + c x^3 + d x^3, x^_]",
                   "(a + b) x^2 + (c + d) x^3", 0);
    assert_eval_eq("Collect[a x^2 + b x^2 + c x + d x, x^_.]",
                   "(c + d) x + (a + b) x^2", 0);
    /* Nested matches are both variables. */
    assert_eval_eq("Collect[a Sqrt[1 + Sqrt[x]] + b Sqrt[1 + Sqrt[x]] + c Sqrt[x] + d Sqrt[x], Power[_, 1/2]]",
                   "(c + d) Sqrt[x] + (a + b) Sqrt[1 + Sqrt[x]]", 0);
    /* No match: unchanged, and the h wrapper still applies. */
    assert_eval_eq("Collect[c x + d x, x^_]", "c x + d x", 0);
    assert_eval_eq("Collect[c + d, f[_], h]", "h[c + d]", 0);
    assert_eval_eq("Collect[a Sqrt[y] + b Sqrt[y] + c, Sqrt[_], h]",
                   "h[c] + h[a + b] Sqrt[y]", 0);
    /* Alternatives, Condition, mixed pattern/concrete lists, equations. */
    assert_eval_eq("Collect[a f[x] + b f[x] + c, f[_] | g[_]]", "c + (a + b) f[x]", 0);
    assert_eval_eq("Collect[a f[1] + b f[1], f[n_] /; n > 0]", "(a + b) f[1]", 0);
    assert_eval_eq("Collect[a f[x] + b f[x] + c x + d x, {f[_], x}]",
                   "(c + d) x + (a + b) f[x]", 0);
    assert_eval_eq("Collect[a f[x] + b f[x] == c f[y] + d f[y], f[_]]",
                   "(a + b) f[x] == (c + d) f[y]", 0);
    /* A match buried in an opaque head does not collect. */
    assert_eval_eq("Collect[a f[g[x]] + b f[g[x]], g[_]]", "a f[g[x]] + b f[g[x]]", 0);
}

/* --- 14 user-supplied "passing" cases ------------------------------- */
static void test_pass_corpus(void) {
    /* P1: Multi-target Times grouping. */
    assert_eval_eq(
        "Collect[x^2*y^4 + z*(x*y^2)^2 + z + a*z, {x*y^2, z}]",
        "(1 + a) z + x^2 y^4 (1 + z)", 0);
    /* P2: Integer-power target. */
    assert_eval_eq("Collect[a^2*(a^2 + 1), a^2]",
                   "a^2 + a^4", 0);
    /* P3: Derivative as exact target. */
    assert_eval_eq("Collect[a*D[f[x], x] + b*D[f[x], x], D[f[x], x]]",
                   "(a + b) Derivative[1][f][x]", 0);
    /* P4: Derivative target with mixed terms. */
    assert_eval_eq(
        "Collect[f[x] + f[x]*D[f[x], x] + x*D[f[x], x]*f[x], D[f[x], x]]",
        "f[x] + Derivative[1][f][x] (x f[x] + f[x])", 0);
    /* P5: Negative-power target (1/f[x]) inside coefficients. */
    assert_eval_eq(
        "Collect[1/f[x] + 1/f[x]*D[f[x], x] + x*D[f[x], x]/f[x], D[f[x], x]]",
        "1/f[x] + Derivative[1][f][x] (1/f[x] + x/f[x])", 0);
    /* P6: Same after Expand. */
    assert_eval_eq(
        "Collect[Expand[1/f[x] + 1/f[x]*D[f[x], x] + x*D[f[x], x]/f[x]], "
        "D[f[x], x]]",
        "1/f[x] + Derivative[1][f][x] (1/f[x] + x/f[x])", 0);
    /* P7: Polynomial expansion. */
    assert_eval_eq("Collect[(x + a + 1)^3, x]",
                   "1 + 3 a + 3 a^2 + a^3 + (3 + 6 a + 3 a^2) x"
                   " + (3 + 3 a) x^2 + x^3", 0);
    /* P8: Simple polynomial. */
    assert_eval_eq(
        "Collect[a*x^4 + b*x^4 + 2*a^2*x - 3*b*x + x - 7, x]",
        "-7 + (1 + 2 a^2 - 3 b) x + (a + b) x^4", 0);
    /* P9: Mixed Sqrt[x] + x^(2/3) + x. */
    assert_eval_eq(
        "Collect[a*Sqrt[x] + Sqrt[x] + x^(2/3) - c*x + 3*x"
        " - 2*b*x^(2/3) + 5,  x]",
        "5 + (1 + a) Sqrt[x] + (1 - 2 b) x^(2/3) + (3 - c) x", 0);
    /* P10: Multivariate cubic. */
    assert_eval_eq(
        "Collect[(x*y + x*z + y*z + x + y)^3, {x, y}]",
        "y^3 (1 + 3 z + 3 z^2 + z^3)"
        " + x (y^2 (3 + 9 z + 9 z^2 + 3 z^3) + y^3 (3 + 6 z + 3 z^2))"
        " + x^2 (y (3 + 9 z + 9 z^2 + 3 z^3) + y^2 (6 + 12 z + 6 z^2)"
        " + y^3 (3 + 3 z))"
        " + x^3 (1 + y^3 + 3 z + 3 z^2 + z^3 + y (3 + 6 z + 3 z^2)"
        " + y^2 (3 + 3 z))", 0);
    /* P11: Larger multivariate cubic. */
    assert_eval_eq(
        "Collect[Expand[(y + 2*x + 3*x*y + 4*x*z + 5*y*z)^3], {x, y}]",
        "y^3 (1 + 15 z + 75 z^2 + 125 z^3)"
        " + x (y^2 (6 + 72 z + 270 z^2 + 300 z^3)"
        " + y^3 (9 + 90 z + 225 z^2))"
        " + x^2 (y (12 + 108 z + 288 z^2 + 240 z^3)"
        " + y^2 (36 + 252 z + 360 z^2)"
        " + y^3 (27 + 135 z))"
        " + x^3 (8 + 27 y^3 + 48 z + 96 z^2 + 64 z^3"
        " + y (36 + 144 z + 144 z^2)"
        " + y^2 (54 + 108 z))", 0);
    /* P13: Higher-order derivative as target. */
    assert_eval_eq(
        "Collect[a*D[f[x],{x,2}] + b*D[f[x],{x,2}], D[f[x],{x,2}]]",
        "(a + b) Derivative[2][f][x]", 0);
    /* P14: Nested radical pattern target. */
    assert_eval_eq(
        "Collect[(Sqrt[15 + 5*Sqrt[2]]*x + Sqrt[3 + Sqrt[2]]*y)*2, "
        "Sqrt[_]]",
        "2 Sqrt[15 + 5 Sqrt[2]] x + 2 Sqrt[3 + Sqrt[2]] y", 0);
}

/* P12 uses Series, which the parser/REPL accepts but the inputs in
 * test files cannot easily span the multi-statement form. Inline it
 * here using the post-Normal closed-form expression instead. */
static void test_pass_p12_normal_series(void) {
    /* The post-Normal[Series[...]] form can be written directly. */
    assert_eval_eq(
        "Collect[Normal[Series[Sin[a + b], {b, 0, 10}]],"
        " {Sin[a], Cos[a]}]",
        "Cos[a] (b - 1/6 b^3 + 1/120 b^5 - 1/5040 b^7 + 1/362880 b^9)"
        " + Sin[a] (1 - 1/2 b^2 + 1/24 b^4 - 1/720 b^6 + 1/40320 b^8"
        " - 1/3628800 b^10)", 0);
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_bug1_symbolic_exp_ratio);
    TEST(test_bug2_pattern_target);
    TEST(test_pass_corpus);
    TEST(test_pass_p12_normal_series);

    printf("All collect_corpus tests passed!\n");
    return 0;
}
