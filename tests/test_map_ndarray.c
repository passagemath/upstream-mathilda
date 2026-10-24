#include "eval.h"
#include "parse.h"
#include "expr.h"
#include "symtab.h"
#include "core.h"
#include "print.h"
#include "test_utils.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* NDArray support for the Map family (Map, MapIndexed, MapThread, MapAll).
 * An NDArray is an atomic value, so before these changes each of these was a
 * silent no-op over one. Assertions avoid bare machine reals (whose printed
 * form is formatting-sensitive) by mapping Head over elements, checking a Part
 * via equality, or counting _Real leaves. */
static void run(const char* input, const char* expected) {
    Expr* e = parse_expression(input);
    Expr* r = evaluate(e);
    char* s = expr_to_string(r);
    ASSERT_MSG(strcmp(s, expected) == 0, "%s: expected %s, got %s", input, expected, s);
    free(s);
    expr_free(e);
    expr_free(r);
}

/* ---------- Map ---------- */

static void test_map_ndarray_rank1(void) {
    /* f is applied to each scalar element (boxed as Real). */
    run("Map[Head, NDArray[Range[4]]]", "{Real, Real, Real, Real}");
    /* An explicit {1} level spec takes the same fast path. */
    run("Map[Head, NDArray[Range[4]], {1}]", "{Real, Real, Real, Real}");
}

static void test_map_ndarray_rank2_rows(void) {
    /* Level-1 parts of a rank-2 array are sub-NDArray rows. */
    run("Map[Head, NDArray[Table[i+j,{i,2},{j,3}]]]", "{NDArray, NDArray}");
}

static void test_map_ndarray_numeric_repacks(void) {
    /* A numeric result repacks into an NDArray (packed in -> packed out). */
    run("Head[Map[#^2 &, NDArray[Range[4]]]]", "NDArray");
    run("Map[#^2 &, NDArray[Range[4]]][[3]] == 9", "True");
}

static void test_map_ndarray_symbolic_is_list(void) {
    /* A symbolic result stays a List. */
    run("Head[Map[f, NDArray[Range[4]]]]", "List");
}

static void test_map_ndarray_deep_level(void) {
    /* A non-default level spec materializes and maps generically. */
    run("Map[Head, NDArray[Table[i+j,{i,2},{j,3}]], {2}]",
        "{{Real, Real, Real}, {Real, Real, Real}}");
}

/* ---------- MapIndexed ---------- */

static void test_mapindexed_ndarray_positions(void) {
    run("MapIndexed[#2 &, NDArray[Range[3]]]", "{{1}, {2}, {3}}");
    run("MapIndexed[#2 &, NDArray[Table[i+j,{i,2},{j,2}]]]", "{{1}, {2}}");
}

static void test_mapindexed_ndarray_values(void) {
    /* The value passed to f is the leading-axis part: scalar for rank 1, a
     * sub-NDArray row for rank 2. */
    run("MapIndexed[Head[#1] &, NDArray[Range[3]]]", "{Real, Real, Real}");
    run("MapIndexed[Head[#1] &, NDArray[Table[i+j,{i,2},{j,2}]]]", "{NDArray, NDArray}");
}

/* ---------- MapThread ---------- */

static void test_mapthread_ndarray_with_list(void) {
    /* An NDArray entry threads element-wise against a plain List. */
    run("MapThread[#2 &, {NDArray[Range[3]], {a, b, c}}]", "{a, b, c}");
    run("MapThread[#1 &, {NDArray[Range[3]], {a, b, c}}][[2]] == 2", "True");
}

static void test_mapthread_ndarray_pair(void) {
    run("MapThread[Plus, {NDArray[{1,2,3}], NDArray[{10,20,30}]}][[3]] == 33", "True");
}

/* ---------- MapAll ---------- */

static void test_mapall_ndarray_wraps_all_levels(void) {
    /* f wraps every level including the whole array at level 0 (head becomes f)
     * and each scalar leaf. */
    run("Head[MapAll[f, NDArray[Range[3]]]]", "f");
    run("Count[MapAll[f, NDArray[Range[3]]], _Real, Infinity]", "3");
    run("Count[MapAll[f, NDArray[Table[i+j,{i,2},{j,2}]]], _Real, Infinity]", "4");
}

/* ---------- Regression: plain lists still behave ---------- */

static void test_plain_list_unregressed(void) {
    run("Map[f, {a, b, c}]", "{f[a], f[b], f[c]}");
    run("Map[#^2 &, {1, 2, 3}]", "{1, 4, 9}");
    run("MapIndexed[f, {a, b, c}]", "{f[a, {1}], f[b, {2}], f[c, {3}]}");
    run("MapThread[f, {{1, 2}, {3, 4}}]", "{f[1, 3], f[2, 4]}");
    run("MapAll[f, {a, {b}}]", "f[{f[a], f[{f[b]}]}]");
}

/* ---------- Heads that walked args and saw an atom ----------
 * Apply, SortBy, Level, MaximalBy/MinimalBy, TakeLargestBy/TakeSmallestBy,
 * Thread and Distribute each treated a visible NDArray as an atom with no
 * elements. Several answered confidently and wrongly (SortBy unsorted,
 * Apply[Times, nd] returning nd, Level[nd, {2}] = {}). Each must now answer
 * as on the equivalent List, with a machine-array result presented as an
 * NDArray. An int64 buffer keeps the printed values exact. */
#define ND_M "NDArray[{{1, 2, 3}, {5, 6, 7}}, DataType -> \"int64\"]"
#define ND_V "NDArray[{3, 1, 2}, DataType -> \"int64\"]"

static void test_apply_ndarray(void) {
    run("Apply[f, " ND_M "]", "f[{1, 2, 3}, {5, 6, 7}]");
    run("Apply[f, " ND_M ", {1}]", "{f[1, 2, 3], f[5, 6, 7]}");
    run("f @@@ " ND_M, "{f[1, 2, 3], f[5, 6, 7]}");
    run("Apply[Plus, " ND_M ", {1}]", "NDArray[{6, 18}]");
    run("Apply[Times, " ND_M "]", "NDArray[{5, 12, 21}]");
    run("Apply[Plus, " ND_V "]", "6");
    /* Rationals have no machine buffer: an ordinary List. */
    run("Apply[#1/#2 &, " ND_M ", {1}]", "{1/2, 5/6}");
}

static void test_sortby_ndarray(void) {
    run("SortBy[" ND_M ", -Total[#] &]", "NDArray[{{5, 6, 7}, {1, 2, 3}}]");
    run("SortBy[" ND_V ", -# &]", "NDArray[{3, 2, 1}]");
    run("SortBy[" ND_V ", Minus, Greater]", "NDArray[{1, 2, 3}]");
    run("SortBy[" ND_V ", Identity]", "NDArray[{1, 2, 3}]");
}

static void test_maximalby_takeby_ndarray(void) {
    run("MaximalBy[" ND_M ", Total]", "NDArray[{{5, 6, 7}}]");
    run("MinimalBy[" ND_M ", Total]", "NDArray[{{1, 2, 3}}]");
    run("MaximalBy[" ND_V ", Minus]", "NDArray[{1}]");
    run("TakeLargestBy[" ND_M ", Total, 1]", "NDArray[{{5, 6, 7}}]");
    run("TakeSmallestBy[" ND_V ", Minus, 2]", "NDArray[{3, 2}]");
}

static void test_level_ndarray(void) {
    run("Level[" ND_M ", {2}]", "NDArray[{1, 2, 3, 5, 6, 7}]");
    run("Level[" ND_M ", {1}]", "NDArray[{{1, 2, 3}, {5, 6, 7}}]");
    run("Level[" ND_V ", {-1}, Heads -> True]", "{List, 3, 1, 2}");
}

static void test_thread_distribute_ndarray(void) {
    run("Thread[f[" ND_M "]]", "{f[{1, 2, 3}], f[{5, 6, 7}]}");
    run("Thread[f[" ND_V ", {a, b, c}]]", "{f[3, a], f[1, b], f[2, c]}");
    run("Thread[f[" ND_V "], List, 1]", "{f[3], f[1], f[2]}");
    run("Distribute[" ND_M ", List]",
        "NDArray[{{1, 5}, {1, 6}, {1, 7}, {2, 5}, {2, 6}, {2, 7}, {3, 5}, {3, 6}, {3, 7}}]");
    run("Distribute[f[" ND_V ", {a, b}], List]",
        "{f[3, a], f[3, b], f[1, a], f[1, b], f[2, a], f[2, b]}");
}

static void test_visible_delist_heads_match_list(void) {
    /* Same answer as the List, value for value, over a real buffer too. */
    run("r = {{1.5, 2.5}, {0.5, 4.}}; nd = NDArray[r]; "
        "{Normal[SortBy[nd, First]] == SortBy[r, First], "
        " Normal[Apply[Plus, nd, {1}]] == Apply[Plus, r, {1}], "
        " Normal[Level[nd, {2}]] == Level[r, {2}], "
        " Normal[MaximalBy[nd, Last]] == MaximalBy[r, Last], "
        " Normal[TakeLargestBy[nd, Total, 2]] == TakeLargestBy[r, Total, 2], "
        " Thread[f[nd, x]] == Thread[f[r, x]]}",
        "{True, True, True, True, True, True}");
    /* Symbolic and List arguments are untouched (an atom is SortBy::normal,
     * unevaluated, as in the reference). */
    run("Quiet[SortBy[x, f]]", "SortBy[x, f]");
    run("Apply[f, {1, 2}]", "f[1, 2]");
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_map_ndarray_rank1);
    TEST(test_map_ndarray_rank2_rows);
    TEST(test_map_ndarray_numeric_repacks);
    TEST(test_map_ndarray_symbolic_is_list);
    TEST(test_map_ndarray_deep_level);
    TEST(test_mapindexed_ndarray_positions);
    TEST(test_mapindexed_ndarray_values);
    TEST(test_mapthread_ndarray_with_list);
    TEST(test_mapthread_ndarray_pair);
    TEST(test_mapall_ndarray_wraps_all_levels);
    TEST(test_plain_list_unregressed);
    TEST(test_apply_ndarray);
    TEST(test_sortby_ndarray);
    TEST(test_maximalby_takeby_ndarray);
    TEST(test_level_ndarray);
    TEST(test_thread_distribute_ndarray);
    TEST(test_visible_delist_heads_match_list);

    printf("All Map-family NDArray tests passed!\n");
    return 0;
}
