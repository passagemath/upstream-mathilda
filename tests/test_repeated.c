/* Repeated (p..) and RepeatedNull (p...) in the pattern matcher.
 *
 * Companion to the tokenizer cases in test_parse.c (issue #100: `0..` was a
 * parse error).  Every expected value below was produced by Mathematica 14 on
 * the same input -- none is Mathilda's own output read back -- so the suite
 * pins WL semantics, including the ones that are easy to get backwards:
 *   - an unconstrained `x:p..` followed by `y___` takes the SHORTEST run (it
 *     is not greedy; Longest[] makes it so);
 *   - a named pattern inside a repeat must bind the same value in every
 *     repetition ({1, 2} does not match {(x_ /; x > 0)..});
 *   - repeat bounds count UNITS, so Repeated[PatternSequence[1, 2], {2}]
 *     consumes four arguments;
 *   - a malformed spec emits Repeated::range / Repeated::order and matches
 *     nothing.
 *
 * Failures are counted rather than assert()ed: the unit-test build may define
 * NDEBUG, and a soft assert would let the summary line report success. */
#include "test_utils.h"
#include "expr.h"
#include "parse.h"
#include "eval.h"
#include "print.h"
#include "symtab.h"
#include "core.h"
#include "message.h"
#include <stdio.h>
#include <string.h>
#include <time.h>

static int g_failures = 0;
static int g_checks = 0;

/* Evaluate `input` and compare its InputForm-style rendering to `expected`. */
static void check(const char* input, const char* expected) {
    g_checks++;
    Expr* parsed = parse_expression(input);
    if (!parsed) {
        printf("FAIL: %s\n  did not parse\n", input);
        g_failures++;
        return;
    }
    Expr* r = evaluate(parsed);
    expr_free(parsed);
    char* s = expr_to_string(r);
    if (strcmp(s, expected) != 0) {
        printf("FAIL: %s\n  Expected: %s\n  Actual:   %s\n", input, expected, s);
        g_failures++;
    }
    free(s);
    expr_free(r);
}

/* As check(), and additionally require that evaluating `input` fired (or did
 * not fire) a message -- read off the Check[] counter, so this also proves the
 * diagnostic is routed through the message funnel. */
static void check_msg(const char* input, const char* expected, bool want_msg) {
    unsigned long before = mth_msg_fired_count();
    check(input, expected);
    bool fired = mth_msg_fired_count() != before;
    if (fired != want_msg) {
        printf("FAIL: %s\n  expected %s message\n", input, want_msg ? "a" : "no");
        g_failures++;
    }
}

/* The issue itself, and the plain single-argument repeat. */
static void test_basic(void) {
    check("MatchQ[{0, 0, 0}, {0..}]", "True");
    check("MatchQ[{}, {0..}]", "False");
    check("MatchQ[{}, {0...}]", "True");
    check("MatchQ[{0, 1}, {0..}]", "False");
    check("MatchQ[{1, 1, 2, 2}, {1.., 2..}]", "True");
    check("MatchQ[{1, 1, 2, 2}, {1.., 1.., 2..}]", "True");
    check("MatchQ[f[a, a], f[a..]]", "True");
    check("MatchQ[f[], f[a...]]", "True");
    check("MatchQ[f[], f[a..]]", "False");
    check("MatchQ[{1.5, 1.5}, {1.5..}]", "True");
    /* A repeat as the whole pattern: the subject is a one-argument sequence. */
    check("MatchQ[a, a..]", "True");
    check("MatchQ[a, a...]", "True");
    check("MatchQ[a, Repeated[a, {2}]]", "False");
}

static void test_element_patterns(void) {
    check("MatchQ[{1, 2, 3}, {_Integer..}]", "True");
    check("MatchQ[{1, 2., 3}, {_Integer..}]", "False");
    check("MatchQ[{1., 2.}, {_Real...}]", "True");
    check("MatchQ[{1, 3, 5}, {_?OddQ..}]", "True");
    check("MatchQ[{1, 2, 3}, {_?OddQ..}]", "False");
    check("MatchQ[{a, b, b}, {(a | b)..}]", "True");
    check("MatchQ[{{1, 2}, {3, 4}}, {{_, _}..}]", "True");
    check("MatchQ[{{1, 2}, {3}}, {{_, _}..}]", "False");
    /* Named sub-patterns bind once, across all repetitions. */
    check("MatchQ[{1, 1}, {(x_ /; x > 0)..}]", "True");
    check("MatchQ[{1, 2}, {(x_ /; x > 0)..}]", "False");
    check("MatchQ[{x, x, x}, {(a_)..}]", "True");
    check("MatchQ[{x, y}, {(a_)..}]", "False");
    check("f[1, 1, 1] /. f[Repeated[x_]] :> {x}", "{1}");
    check("MatchQ[{2, 2, 2}, {Repeated[x_, 3]}]", "True");
    check("MatchQ[{2, 2, 3}, {Repeated[x_, 3]}]", "False");
}

static void test_bounds(void) {
    check("MatchQ[{1, 1, 1}, {Repeated[1, {2, 3}]}]", "True");
    check("MatchQ[{1, 1, 1, 1}, {Repeated[1, {2, 3}]}]", "False");
    check("MatchQ[{1, 1}, {Repeated[1, 2]}]", "True");
    check("MatchQ[{1, 1, 1}, {Repeated[1, 2]}]", "False");
    check("MatchQ[{1, 1}, {Repeated[1, {2}]}]", "True");
    check("MatchQ[{1, 1, 1}, {Repeated[1, {2, Infinity}]}]", "True");
    check("MatchQ[{1}, {Repeated[1, {2, Infinity}]}]", "False");
    check("MatchQ[{}, {Repeated[1, {0}]}]", "True");
    check("MatchQ[{1, 1, 1}, {Repeated[1, {0, 1}], 1, 1}]", "True");
    check("MatchQ[{1, 1, 1}, {RepeatedNull[1, {1, 2}]}]", "False");
    check("MatchQ[{}, {RepeatedNull[1, 2]}]", "True");
    check("MatchQ[{}, {RepeatedNull[1, 0]}]", "True");
    check("MatchQ[{1}, {RepeatedNull[1, 0]}]", "False");
    check("MatchQ[{1, 2, 3}, {Repeated[_Integer, {0, 5}]}]", "True");
    /* A bound too large for any argument list is effectively unbounded. */
    check("MatchQ[{1, 1}, {Repeated[1, 10^30]}]", "True");
    check("MatchQ[{1, 1}, {Repeated[1, {0, 10^30}]}]", "True");
}

static void test_invalid_specs(void) {
    check_msg("MatchQ[{1, 1}, {Repeated[1, -1]}]", "False", true);
    check_msg("MatchQ[{1, 1}, {Repeated[1, {a, b}]}]", "False", true);
    check_msg("MatchQ[{1, 1}, {Repeated[1, {2.5}]}]", "False", true);
    check_msg("MatchQ[{1, 1}, {Repeated[1, {3, 2}]}]", "False", true);
    check_msg("MatchQ[{}, {Repeated[1, 0]}]", "False", true);
    /* Routed through the funnel: Quiet silences it, Check sees it. */
    check("Quiet[MatchQ[{1, 1}, {Repeated[1, -1]}]]", "False");
    check("Quiet[Check[MatchQ[{1, 1}, {Repeated[1, -1]}], err]]", "err");
    check_msg("MatchQ[{1, 1}, {Repeated[1, {1, 2}]}]", "True", false);
}

/* Binding a repeat to a name, and how a repeat shares a list with others. */
static void test_named_and_split(void) {
    check("{1, 1, 2} /. {x:1.., y___} :> {{x}, {y}}", "{{1}, {1, 2}}");
    check("{1, 1, 1, 2} /. {x:1..., y___} :> {{x}, {y}}", "{{}, {1, 1, 1, 2}}");
    check("{1, 1, 1, 2} /. {Longest[x:1...], y___} :> {{x}, {y}}", "{{1, 1, 1}, {2}}");
    check("f[1, 1, 1] /. f[x:1..] :> {x}", "{1, 1, 1}");
    check("{f[a, a, b]} /. f[x:a.., y_] :> g[{x}, y]", "{g[{a, a}, b]}");
    check("ReplaceList[{1, 1, 1}, {x:1..., y:1...} :> {{x}, {y}}]",
          "{{{}, {1, 1, 1}}, {{1}, {1, 1}}, {{1, 1}, {1}}, {{1, 1, 1}, {}}}");
    check("ReplaceList[{1, 1, 1}, {x__, y:1..} :> {{x}, {y}}]",
          "{{{1}, {1, 1}}, {{1, 1}, {1}}}");
    check("ReplaceList[{1, 1, 1, 1}, {Repeated[1, {1, 2}], y___} :> {y}]",
          "{{1, 1, 1}, {1, 1}}");
}

/* A repeat whose unit spans several arguments. */
static void test_sequence_units(void) {
    check("MatchQ[{a, b, a, b}, {PatternSequence[a, b]..}]", "True");
    check("MatchQ[{a, b, a}, {PatternSequence[a, b]..}]", "False");
    check("MatchQ[{1, 2, 1, 2, 1, 2}, {PatternSequence[1, 2]..}]", "True");
    check("MatchQ[{1, 2, 3, 1, 2, 3}, {PatternSequence[1, _, 3]..}]", "True");
    check("MatchQ[f[a, b, a, b], f[PatternSequence[a, b]..]]", "True");
    check("MatchQ[a, PatternSequence[a]..]", "True");
    /* Bounds count units. */
    check("MatchQ[{1, 2, 1, 2}, {Repeated[PatternSequence[1, 2], {2}]}]", "True");
    check("MatchQ[{1, 2, 1, 2, 1, 2}, {Repeated[PatternSequence[1, 2], {2}]}]", "False");
    check("MatchQ[{1, 2, 1, 2, 1, 2}, {Repeated[PatternSequence[1, 2], 2], ___}]", "True");
    check("{1, 2, 1, 2, 3} /. {x:Repeated[PatternSequence[1, 2], 1], y___} :> {{x}, {y}}",
          "{{1, 2}, {1, 2, 3}}");
    check("{a, b, a, b, c} /. {x:PatternSequence[a, b].., y___} :> {{x}, {y}}",
          "{{a, b}, {a, b, c}}");
    /* Empty subjects and empty units. */
    check("MatchQ[{}, {PatternSequence[a, b]...}]", "True");
    check("MatchQ[{}, {PatternSequence[]...}]", "True");
    check("MatchQ[{}, {PatternSequence[]..}]", "True");
    check("MatchQ[{1, 2, 3, 4}, {PatternSequence[_, _]...}]", "True");
    check("MatchQ[{1, 2, 3}, {PatternSequence[_, _]...}]", "False");
    /* Named sub-patterns are shared by every unit. */
    check("MatchQ[{1, 2, 1, 2}, {PatternSequence[x_, y_]..}]", "True");
    check("MatchQ[{1, 2, 3, 4}, {PatternSequence[x_, y_]..}]", "False");
    check("MatchQ[{a, b, a, b}, {x:PatternSequence[_, _].., x}]", "False");
    /* Variable-width units. */
    check("MatchQ[{1, 2, 3, 4, 5}, {PatternSequence[_, __]..}]", "True");
    check("{1, 2, 3, 4, 5} /. {x:PatternSequence[_, __]..} :> {x}", "{1, 2, 3, 4, 5}");
    check("ReplaceList[{1, 2, 3, 4}, {PatternSequence[x_, y__]..} :> {x, {y}}]",
          "{{1, {2, 3, 4}}}");
    check("MatchQ[{1, 2, 3}, {x__..}]", "True");
    check("MatchQ[{1, 2, 1, 2}, {x__..}]", "True");
    check("{1, 2, 1, 2} /. {x__..} :> {x}", "{1, 2}");
    /* Alternatives with a multi-argument branch, repeated or not. */
    check("MatchQ[{a, b, a, b}, {(PatternSequence[a, b] | c)..}]", "True");
    check("MatchQ[{c, a, b}, {(PatternSequence[a, b] | c)..}]", "True");
    check("MatchQ[{c, a, c}, {(PatternSequence[a, b] | c)..}]", "False");
    check("MatchQ[{a, b}, {PatternSequence[a, b] | c}]", "True");
    check("MatchQ[{a, b, d}, {PatternSequence[a, b] | c, d}]", "True");
    check("{a, b, d} /. {x:(PatternSequence[a, b] | c), y___} :> {{x}, {y}}",
          "{{a, b}, {d}}");
}

/* Repeats inside the pattern-driven list functions. */
static void test_consumers(void) {
    check("Cases[{{1, 1}, {1, 2}, {}}, {1..}]", "{{1, 1}}");
    check("Position[{1, {1, 1}, {2}}, {1..}]", "{{2}}");
    check("Count[{{1}, {1, 1}, {2}, {}}, {1...}]", "3");
    check("DeleteCases[{{1}, {1, 1}, {2}, {}}, {1..}]", "{{2}, {}}");
    check("FreeQ[{1, 2}, 1..]", "False");
    check("StringMatchQ[\"aaa\", \"a\"..]", "True");
    check("StringMatchQ[\"aab\", \"a\"..]", "False");
    check("StringMatchQ[\"\", \"a\"...]", "True");
    check("Attributes[Repeated]", "{Protected}");
    check("Attributes[RepeatedNull]", "{Protected}");
}

/* Long subjects stay linear: fixed-width units are checked in one pass and
 * the variable-width search only tries feasible unit widths (it used to be
 * quadratic on a failing match), without recursing per unit. */
static void test_scaling(void) {
    struct { const char* in; const char* out; } cases[] = {
        { "MatchQ[ConstantArray[1, 20000], {1..}]", "True" },
        { "MatchQ[Append[ConstantArray[1, 20000], 2], {1..}]", "False" },
        { "MatchQ[Range[20000], {_Integer..}]", "True" },
        { "MatchQ[Flatten[ConstantArray[{a, b}, 20000]], {PatternSequence[a, b]..}]", "True" },
        { "MatchQ[Flatten[ConstantArray[{a, b}, 20000]], {(PatternSequence[a, b] | c)..}]", "True" },
        { "MatchQ[Append[Flatten[ConstantArray[{a, b}, 20000]], a], {(PatternSequence[a, b] | c)..}]", "False" },
    };
    for (size_t i = 0; i < sizeof cases / sizeof cases[0]; i++) {
        clock_t t0 = clock();
        check(cases[i].in, cases[i].out);
        double secs = (double)(clock() - t0) / CLOCKS_PER_SEC;
        if (secs > 2.0) {
            printf("FAIL: %s took %.2fs (expected linear time)\n", cases[i].in, secs);
            g_failures++;
        }
    }
}

int main(void) {
    symtab_init();
    core_init();
    TEST(test_basic);
    TEST(test_element_patterns);
    TEST(test_bounds);
    TEST(test_invalid_specs);
    TEST(test_named_and_split);
    TEST(test_sequence_units);
    TEST(test_consumers);
    TEST(test_scaling);
    if (g_failures) {
        printf("\n%d of %d Repeated checks FAILED\n", g_failures, g_checks);
        return 1;
    }
    printf("\nAll %d Repeated tests passed!\n", g_checks);
    return 0;
}
