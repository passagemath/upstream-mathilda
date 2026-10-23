/*
 * test_subvalues.c -- SubValues: definitions whose left-hand side has a
 * compound head, f[x_][y_] := ..., keyed on the innermost head symbol f.
 *
 * Covers definition by SetDelayed / Set / TagSet, any depth, pattern and
 * condition forms, memoisation inside a right-hand side, head-first
 * evaluation (a DownValue of f[a] wins), SubValues[], Clear / ClearAll /
 * Unset, Protected tags, OptionsPattern, and that ordinary compound-head
 * calls (pure functions, Derivative) are unaffected.
 */
#include "test_utils.h"
#include "symtab.h"
#include "core.h"

#include <stdio.h>

static void test_basic_and_depth(void) {
    assert_eval_eq("ev[h_][y_] := {h, y}; {ev[1][2], ev[1]}", "{{1, 2}, ev[1]}", 0);
    assert_eval_eq("g3[a_][b_][c_] := {a, b, c}; {g3[1][2][3], g3[1][2]}", "{{1, 2, 3}, g3[1][2]}", 0);
    assert_eval_eq("SubValues[ev]", "{HoldPattern[ev[h_][y_]] :> {h, y}}", 0);
    assert_eval_eq("Attributes[SubValues]", "{HoldAll, Protected}", 0);
    assert_eval_eq("{DownValues[ev], SubValues[nosuchsv]}", "{{}, {}}", 0);
}

static void test_patterns_and_conditions(void) {
    assert_eval_eq("sp[x_, All][a__] := {\"all\", x, a}; sp[x_, n_Integer][a_] := {\"int\", x, n, a}; "
                   "{sp[q, All][1, 2], sp[q, 3][7], sp[q, 3][7, 8]}",
                   "{{\"all\", q, 1, 2}, {\"int\", q, 3, 7}, sp[q, 3][7, 8]}", 0);
    /* LHS guard, then a fallback with the same shape. */
    assert_eval_eq("h2[a_][b_] /; a > b := \"gt\"; h2[a_][b_] := \"other\"; {h2[5][1], h2[1][5]}",
                   "{\"gt\", \"other\"}", 0);
    /* RHS guard (hoisted). */
    assert_eval_eq("rg[a_][b_] := \"big\" /; a + b > 10; {rg[5][6], rg[1][1]}", "{\"big\", rg[1][1]}", 0);
    /* Scoped guard. */
    assert_eval_eq("mg[a_][b_] := Module[{s = a + b}, s /; s > 2]; mg[a_][b_] := 0; {mg[1][5], mg[0][1]}",
                   "{6, 0}", 0);
}

static void test_memoisation(void) {
    assert_eval_eq("mem[h_][y_, d_] := mem[h][y, d] = {\"computed\", h, y, d}; mem[k][2, 3]",
                   "{\"computed\", k, 2, 3}", 0);
    /* The second call uses the stored value; two SubValues now. */
    assert_eval_eq("{mem[k][2, 3], Length[SubValues[mem]]}", "{{\"computed\", k, 2, 3}, 2}", 0);
    assert_eval_eq("First[SubValues[mem]]", "HoldPattern[mem[k][2, 3]] :> {\"computed\", k, 2, 3}", 0);
    /* Immediate Set evaluates the arguments inside the head too. */
    assert_eval_eq("svx = 5; fq[svx][1] = 3; SubValues[fq]", "{HoldPattern[fq[5][1]] :> 3}", 0);
    assert_eval_eq("k2[x_][y_] = x + y; {k2[1][2], k2[a][b]}", "{3, a + b}", 0);
}

static void test_head_first(void) {
    /* The head is evaluated first: a DownValue of f[1] makes the call
     * "down"[3]; other heads reach the SubValue. */
    assert_eval_eq("f9[1] = \"down\"; f9[x_][y_] := \"sub\"; {f9[1], f9[2][3], f9[1][3]}",
                   "{\"down\", \"sub\", \"down\"[3]}", 0);
    assert_eval_eq("{DownValues[f9], SubValues[f9]}",
                   "{{HoldPattern[f9[1]] :> \"down\"}, {HoldPattern[f9[x_][y_]] :> \"sub\"}}", 0);
}

static void test_clear_unset_tags(void) {
    assert_eval_eq("ClearAll[ev]; {SubValues[ev], ev[1][2]}", "{{}, ev[1][2]}", 0);
    assert_eval_eq("Clear[mem]; SubValues[mem]", "{}", 0);
    assert_eval_eq("k2[x_][y_] =.; {SubValues[k2], k2[1][2]}", "{{}, k2[1][2]}", 0);
    /* Unset of a literal call that only a pattern rule covers: Unset::norep. */
    assert_eval_eq("fr[x_][y_] := {x, y}; Check[fr[1][2] =., \"msg\"]", "\"msg\"", 0);
    assert_eval_eq("fr[1][2]", "{1, 2}", 0);
    /* TagSet with the tag as the innermost head. */
    assert_eval_eq("fu /: fu[a_][b_] := \"tagsub\"; {fu[1][2], SubValues[fu]}",
                   "{\"tagsub\", {HoldPattern[fu[a_][b_]] :> \"tagsub\"}}", 0);
    /* A Protected innermost head cannot get a SubValue. */
    assert_eval_eq("Check[Sin[x_][y_] := 1, \"msg\"]", "\"msg\"", 0);
    assert_eval_eq("Sin[1][2]", "Sin[1][2]", 0);
    /* Remove drops them too. */
    assert_eval_eq("rmv[a_][b_] := 1; Remove[rmv]; rmv[1][2]", "rmv[1][2]", 0);
}

static void test_options_and_unwind(void) {
    assert_eval_eq("Options[osv] = {o -> 1}; osv[a_][x_, OptionsPattern[]] := {a, x, OptionValue[o]}; "
                   "{osv[q][1], osv[q][1, o -> 7]}", "{{q, 1, 1}, {q, 1, 7}}", 0);
    assert_eval_eq("thr[a_][b_] := Throw[{a, b}]; Catch[thr[1][2]; \"no\"]", "{1, 2}", 0);
}

static void test_unaffected_compound_heads(void) {
    assert_eval_eq("{(# + 1 &)[2], Function[x, x^2][3], Map[#^2 &, {1, 2}]}", "{3, 9, {1, 4}}", 0);
    assert_eval_eq("{D[Sin[x], x], Derivative[1][#^2 &][3]}", "{Cos[x], 6}", 0);
    assert_eval_eq("dsv[x_] := x^3; dsv'[2]", "12", 0);
    /* An UpValue of an argument still fires on a compound-head call. */
    assert_eval_eq("uv /: (hh[1])[uv] := \"up\"; hh[1][uv]", "\"up\"", 0);
}

int main(void) {
    symtab_init();
    core_init();
    TEST(test_basic_and_depth);
    TEST(test_patterns_and_conditions);
    TEST(test_memoisation);
    TEST(test_head_first);
    TEST(test_clear_unset_tags);
    TEST(test_options_and_unwind);
    TEST(test_unaffected_compound_heads);
    printf("All SubValues tests passed.\n");
    return 0;
}
