/*
 * test_setpredicates.c -- IntersectingQ, DuplicateFreeQ and Permute.
 *
 * Each head is checked on plain lists, on packed int64 lists (the machine-word
 * paths in src/list/setops.c), on visible NDArrays (which must delist rather
 * than stay unevaluated), on non-List heads, and on the message-and-unevaluated
 * contract for bad arguments.
 */
#include "test_utils.h"
#include "symtab.h"
#include "core.h"

#include <stdio.h>

static void test_intersectingq(void) {
    assert_eval_eq("{IntersectingQ[{1, 2}, {2, 3}], IntersectingQ[{1, 2}, {3}], IntersectingQ[{}, {}]}",
                   "{True, False, False}", 0);
    /* SameQ by default: 1. and 1 are different elements. */
    assert_eval_eq("IntersectingQ[{1., 2}, {1}]", "False", 0);
    assert_eval_eq("IntersectingQ[{1, 2}, {3, 4}, SameTest -> (Abs[#1 - #2] == 1 &)]", "True", 0);
    assert_eval_eq("IntersectingQ[{1, 2}, {5}, SameTest -> (Abs[#1 - #2] == 1 &)]", "False", 0);
    assert_eval_eq("IntersectingQ[{1, 2}, {1}, SameTest -> Automatic]", "True", 0);
    assert_eval_eq("IntersectingQ[f[1, 2], f[2]]", "True", 0);
    assert_eval_eq("IntersectingQ[{a, {b}}, {{b}}]", "True", 0);
    /* Packed int64 buffers, and a visible NDArray. */
    assert_eval_eq("{IntersectingQ[Range[10^5], Range[10^5, 2 10^5]], "
                   "IntersectingQ[Range[10^5], Range[10^5 + 1, 2 10^5]]}", "{True, False}", 0);
    assert_eval_eq("{IntersectingQ[Range[5], {}], IntersectingQ[Range[3], {3., x}]}", "{False, False}", 0);
    /* A visible NDArray: int64 elements are Integers, float64 elements Reals. */
    assert_eval_eq("{IntersectingQ[NDArray[{1, 2, 3}, DataType -> \"int64\"], {3, 4}], "
                   "IntersectingQ[NDArray[{1, 2, 3}], {3, 4}], IntersectingQ[NDArray[{1, 2, 3}], {3., 4}]}",
                   "{True, False, True}", 0);
    assert_eval_eq("IntersectingQ[N[Range[3]], {2.}]", "True", 0);
    /* Bad arguments: a message, unevaluated. */
    assert_eval_eq("Quiet[IntersectingQ[1, {1}]]", "IntersectingQ[1, {1}]", 0);
    assert_eval_eq("Check[IntersectingQ[f[1], {1}], \"msg\"]", "\"msg\"", 0);
    assert_eval_eq("Check[IntersectingQ[{1}], \"msg\"]", "\"msg\"", 0);
    assert_eval_eq("Attributes[IntersectingQ]", "{Protected}", 0);
    assert_eval_eq("Options[IntersectingQ]", "{SameTest -> Automatic}", 0);
}

static void test_duplicatefreeq(void) {
    assert_eval_eq("{DuplicateFreeQ[{1, 2, 3}], DuplicateFreeQ[{1, 2, 1}], DuplicateFreeQ[{}]}",
                   "{True, False, True}", 0);
    assert_eval_eq("{DuplicateFreeQ[{1, 1.}], DuplicateFreeQ[{\"a\", \"A\"}], DuplicateFreeQ[f[a, a]]}",
                   "{True, True, False}", 0);
    assert_eval_eq("DuplicateFreeQ[{1, 2, 3}, Abs[#1 - #2] == 1 &]", "False", 0);
    assert_eval_eq("DuplicateFreeQ[{1, 3, 5}, Abs[#1 - #2] == 1 &]", "True", 0);
    /* Packed int64 (bounded range and wide range), empty, and a visible NDArray. */
    assert_eval_eq("{DuplicateFreeQ[Range[10^5]], DuplicateFreeQ[Append[Range[10^5], 7]], "
                   "DuplicateFreeQ[10^12 Range[1000]], DuplicateFreeQ[Join[10^12 Range[1000], {10^12}]]}",
                   "{True, False, True, False}", 0);
    assert_eval_eq("{DuplicateFreeQ[NDArray[{1, 2, 3}, DataType -> \"int64\"]], DuplicateFreeQ[NDArray[{1, 2, 1}]], "
                   "DuplicateFreeQ[N[{1, 2, 1}]]}", "{True, False, False}", 0);
    assert_eval_eq("Quiet[DuplicateFreeQ[3]]", "DuplicateFreeQ[3]", 0);
    assert_eval_eq("Check[DuplicateFreeQ[3], \"msg\"]", "\"msg\"", 0);
    assert_eval_eq("Attributes[DuplicateFreeQ]", "{Protected}", 0);
}

static void test_permute(void) {
    /* Element i moves to position perm[[i]]. */
    assert_eval_eq("Permute[{a, b, c}, {3, 1, 2}]", "{b, c, a}", 0);
    assert_eval_eq("Permute[{a, b, c, d}, {2, 3, 1, 4}]", "{c, a, b, d}", 0);
    assert_eval_eq("Permute[{a, b, c}, {2, 1}]", "{b, a, c}", 0);
    /* Cycles: c1 -> c2 -> ... -> c1. */
    assert_eval_eq("Permute[{a, b, c, d}, Cycles[{{1, 3}}]]", "{c, b, a, d}", 0);
    assert_eval_eq("Permute[f[a, b, c], Cycles[{{1, 2, 3}}]]", "f[c, a, b]", 0);
    assert_eval_eq("Permute[{a, b, c, d}, Cycles[{{1, 2}, {3, 4}}]]", "{b, a, d, c}", 0);
    assert_eval_eq("Permute[{a, b, c, d}, Cycles[{}]]", "{a, b, c, d}", 0);
    /* The FeynCalc Tdec idiom: swap position 2 with each later one. */
    assert_eval_eq("Table[Permute[{p, q, r, s}, Cycles[{{2, i}}]], {i, 3, 4}]",
                   "{{p, r, q, s}, {p, s, r, q}}", 0);
    /* Packed input comes back permuted. */
    assert_eval_eq("Permute[Range[4], Cycles[{{1, 4}}]]", "{4, 2, 3, 1}", 0);
    assert_eval_eq("Permute[NDArray[{1, 2, 3}, DataType -> \"int64\"], {3, 1, 2}]", "{2, 3, 1}", 0);
    /* Invalid permutations: a message, unevaluated. */
    assert_eval_eq("Quiet[Permute[{a, b}, {1, 1}]]", "Permute[{a, b}, {1, 1}]", 0);
    assert_eval_eq("Check[Permute[{a, b}, {1, 1}], \"msg\"]", "\"msg\"", 0);
    assert_eval_eq("Quiet[Permute[{a, b}, Cycles[{{1, 5}}]]]", "Permute[{a, b}, Cycles[{{1, 5}}]]", 0);
    assert_eval_eq("Check[Permute[{a, b}, Cycles[{{1, 5}}]], \"msg\"]", "\"msg\"", 0);
    assert_eval_eq("Check[Permute[{a, b}, {3, 1, 2}], \"msg\"]", "\"msg\"", 0);
    assert_eval_eq("Permute[{a, b}, Cycles[{{1, 2}, {2, 1}}]]", "Permute[{a, b}, Cycles[{{1, 2}, {2, 1}}]]", 0);
    assert_eval_eq("Permute[x, {1}]", "Permute[x, {1}]", 0);
    assert_eval_eq("Attributes[Permute]", "{Protected}", 0);
}

int main(void) {
    symtab_init();
    core_init();
    TEST(test_intersectingq);
    TEST(test_duplicatefreeq);
    TEST(test_permute);
    printf("All set-predicate tests passed.\n");
    return 0;
}
