/* test_sparse_func.c -- SparseArray through the functional and procedural
 * builtins:
 *
 *   - level specs: Map / MapIndexed / Apply / Scan at levels below 1, where
 *     the parts of a SparseArray's rows are reached (src/sparse/sparse_func.c);
 *   - iterators: Table / Do / Sum / Product over {i, s} (src/iter.c);
 *   - ordering, set and fold heads answered as on Normal[s] (Sort, Union,
 *     Split, Gather, Ordering, Fold, AnyTrue, MapAt, ...) and the atom
 *     ::normal message they give for a genuine atom;
 *   - mutation: AppendTo / PrependTo / AddTo / TimesBy on a sparse variable,
 *     s[[All, k]] = sparse;
 *   - threading and nesting: Thread, Flatten, Outer, Times with a dense List;
 *   - Clip / Chop and the background of a re-sparsified answer.
 *
 * Expected values were ground-truthed against the reference, except where
 * Mathilda deliberately differs: traversal visits elements in order (the
 * reference visits stored values first, then the background once), and
 * Sort / Union / Fold / ... answer rather than refuse with ::normal.
 *
 * Failures are counted rather than asserted: the CMake Release build passes
 * -DNDEBUG, which would no-op assert(). */

#include "test_utils.h"
#include "symtab.h"
#include "core.h"

static int failures = 0;

/* Evaluate `input`; fail unless it prints exactly `expected`. */
static void chk(const char* input, const char* expected) {
    struct Expr* p = parse_expression(input);
    if (!p) { fprintf(stderr, "FAIL(parse): %s\n", input); failures++; return; }
    struct Expr* e = evaluate(p);
    expr_free(p);
    char* s = expr_to_string(e);
    if (strcmp(s, expected) != 0) {
        fprintf(stderr, "FAIL: %s\n  Expected: %s\n  Actual:   %s\n", input, expected, s);
        failures++;
    }
    free(s);
    expr_free(e);
}

static void chk_true(const char* input) { chk(input, "True"); }

/* Run a statement for its side effect (definitions). */
static void run(const char* input) {
    struct Expr* p = parse_expression(input);
    if (!p) { fprintf(stderr, "FAIL(parse): %s\n", input); failures++; return; }
    struct Expr* e = evaluate(p);
    expr_free(p);
    expr_free(e);
}

/* `op` answers with head `head` and the Normal form `normal`. */
static void chk_hn(const char* op, const char* head, const char* normal) {
    char buf[4096];
    snprintf(buf, sizeof buf, "With[{r = %s}, {Head[r], Normal[r]}]", op);
    char want[4096];
    snprintf(want, sizeof want, "{%s, %s}", head, normal);
    chk(buf, want);
}

/* `op` agrees, after Normal, with the same operation on the Normal forms
 * (`ref`, which may spell it differently). */
static void chk_inv(const char* op, const char* ref) {
    char buf[4096];
    snprintf(buf, sizeof buf, "Normal[%s] === Normal[%s]", op, ref);
    chk_true(buf);
}

/* ------------------------------------------------------------------------ */

static void setup(void) {
    run("v = SparseArray[{2 -> 3, 5 -> -1}, 6]");
    run("m = SparseArray[{{1, 1} -> 1, {2, 3} -> 2.5, {3, 2} -> x}, {3, 4}]");
    run("mi = SparseArray[{{1, 1} -> 1, {2, 3} -> 2, {3, 2} -> -1}, {3, 4}]");
    run("t = SparseArray[{{1, 1, 1} -> 1, {2, 2, 2} -> 2}, {2, 2, 2}]");
    run("b = SparseArray[{1 -> 5}, 4, 7]");
    run("e = SparseArray[{}, {3}]");
    run("vn = Normal[v]; mn = Normal[m]; min = Normal[mi]; tn = Normal[t]; bn = Normal[b]");
}

static void test_map_levels(void) {
    /* leaf level: element map, sparse, background f[bg] */
    chk_hn("Map[f, m, {2}]", "SparseArray",
           "{{f[1], f[0], f[0], f[0]}, {f[0], f[0], f[2.5], f[0]}, {f[0], f[x], f[0], f[0]}}");
    chk("Map[f, m, {2}][\"Background\"]", "f[0]");
    chk_hn("Map[f, m, {-1}]", "SparseArray",
           "{{f[1], f[0], f[0], f[0]}, {f[0], f[0], f[2.5], f[0]}, {f[0], f[x], f[0], f[0]}}");
    chk_hn("Map[f, t, {3}]", "SparseArray",
           "{{{f[1], f[0]}, {f[0], f[0]}}, {{f[0], f[0]}, {f[0], f[2]}}}");
    chk_hn("Map[f, v, {1}]", "SparseArray", "{f[0], f[3], f[0], f[0], f[-1], f[0]}");
    chk_hn("Map[f, v, Infinity]", "SparseArray", "{f[0], f[3], f[0], f[0], f[-1], f[0]}");
    chk_hn("Map[f, b]", "SparseArray", "{f[5], f[7], f[7], f[7]}");
    chk_hn("Map[f, m, {2}, Heads -> True]", "SparseArray",
           "{{f[1], f[0], f[0], f[0]}, {f[0], f[0], f[2.5], f[0]}, {f[0], f[x], f[0], f[0]}}");
    chk_hn("Map[Sin, mi, {2}]", "SparseArray",
           "{{Sin[1], 0, 0, 0}, {0, 0, Sin[2], 0}, {0, -Sin[1], 0, 0}}");
    /* {-1} over a non-atomic element goes inside it, like the dense form */
    chk_hn("Map[f, SparseArray[{1 -> x^2}, 2], {-1}]", "SparseArray", "{f[x]^f[2], f[0]}");
    /* other levels: as on Normal, a dense answer */
    chk_hn("Map[f, m, {1, 2}]", "List",
           "{f[{f[1], f[0], f[0], f[0]}], f[{f[0], f[0], f[2.5], f[0]}], f[{f[0], f[x], f[0], f[0]}]}");
    chk("Map[f, m, All]",
        "f[{f[{f[1], f[0], f[0], f[0]}], f[{f[0], f[0], f[2.5], f[0]}], f[{f[0], f[x], f[0], f[0]}]}]");
    chk("Map[f, m, {0, 1}]", "f[{f[{1, 0, 0, 0}], f[{0, 0, 2.5, 0}], f[{0, x, 0, 0}]}]");
    chk_hn("Map[f, t, {2}]", "List", "{{f[{1, 0}], f[{0, 0}]}, {f[{0, 0}], f[{0, 2}]}}");
    chk_hn("Map[f, m, {-2}]", "List", "{f[{1, 0, 0, 0}], f[{0, 0, 2.5, 0}], f[{0, x, 0, 0}]}");
    chk("Normal[Map[f, m, 2]]",
        "{f[{f[1], f[0], f[0], f[0]}], f[{f[0], f[0], f[2.5], f[0]}], f[{f[0], f[x], f[0], f[0]}]}");
    /* level 1 of a matrix: the rows */
    chk("Normal[Map[f, m]]", "{f[{1, 0, 0, 0}], f[{0, 0, 2.5, 0}], f[{0, x, 0, 0}]}");
    chk("Normal[Map[f, m, {1}]]", "{f[{1, 0, 0, 0}], f[{0, 0, 2.5, 0}], f[{0, x, 0, 0}]}");
    /* the invariant over the level specs */
    const char* levels[] = { "{2}", "{-1}", "{1, 2}", "All", "{0, 1}", "{-2}", "2", "Infinity" };
    char a[256], r[256];
    for (size_t i = 0; i < sizeof levels / sizeof *levels; i++) {
        snprintf(a, sizeof a, "Map[g, mi, %s]", levels[i]);
        snprintf(r, sizeof r, "Map[g, min, %s]", levels[i]);
        chk_inv(a, r);
        snprintf(a, sizeof a, "Map[g, t, %s]", levels[i]);
        snprintf(r, sizeof r, "Map[g, tn, %s]", levels[i]);
        chk_inv(a, r);
    }
}

static void test_mapindexed_apply_scan(void) {
    chk_hn("MapIndexed[f, m, {2}]", "List",
           "{{f[1, {1, 1}], f[0, {1, 2}], f[0, {1, 3}], f[0, {1, 4}]}, "
           "{f[0, {2, 1}], f[0, {2, 2}], f[2.5, {2, 3}], f[0, {2, 4}]}, "
           "{f[0, {3, 1}], f[x, {3, 2}], f[0, {3, 3}], f[0, {3, 4}]}}");
    chk_hn("MapIndexed[f, v]", "List",
           "{f[0, {1}], f[3, {2}], f[0, {3}], f[0, {4}], f[-1, {5}], f[0, {6}]}");
    chk_inv("MapIndexed[f, mi, {1, 2}]", "MapIndexed[f, min, {1, 2}]");
    chk_inv("MapIndexed[f, t, {3}]", "MapIndexed[f, tn, {3}]");

    chk_hn("Apply[f, m, {1}]", "List", "{f[1, 0, 0, 0], f[0, 0, 2.5, 0], f[0, x, 0, 0]}");
    chk_hn("f @@@ m", "List", "{f[1, 0, 0, 0], f[0, 0, 2.5, 0], f[0, x, 0, 0]}");
    chk_hn("Apply[f, t, {2}]", "List", "{{f[1, 0], f[0, 0]}, {f[0, 0], f[0, 2]}}");
    chk("Apply[f, m, {0, 1}]", "f[f[1, 0, 0, 0], f[0, 0, 2.5, 0], f[0, x, 0, 0]]");
    chk("Normal[f @@ m]", "f[{1, 0, 0, 0}, {0, 0, 2.5, 0}, {0, x, 0, 0}]");
    chk("Plus @@ v", "2");
    chk("Normal[Plus @@ mi]", "{1, -1, 2, 0}");
    chk_hn("Apply[f, v, {1}]", "List", "{0, 3, 0, 0, -1, 0}");
    chk("Apply[List, m, {1}]", "{{1, 0, 0, 0}, {0, 0, 2.5, 0}, {0, x, 0, 0}}");

    /* Scan visits elements in order, the background included */
    chk("Reap[Scan[Sow, m, {2}]]", "{Null, {{1, 0, 0, 0, 0, 0, 2.5, 0, 0, x, 0, 0}}}");
    chk("Reap[Scan[Sow, v]]", "{Null, {{0, 3, 0, 0, -1, 0}}}");
    chk("Reap[Scan[Sow, t, {-1}]]", "{Null, {{1, 0, 0, 0, 0, 0, 0, 2}}}");
    chk("Catch[Scan[If[# != 0, Throw[#]] &, m, {2}]]", "1");
    chk("Catch[Scan[If[# < 0, Throw[#]] &, v]]", "-1");
    chk("Module[{s = 0}, Scan[(s += #) &, b]; s]", "26");
    chk("Normal[Reap[Scan[Sow, m]][[2, 1]]]", "{{1, 0, 0, 0}, {0, 0, 2.5, 0}, {0, x, 0, 0}}");
}

static void test_iterators(void) {
    chk("Table[i, {i, v}]", "{0, 3, 0, 0, -1, 0}");
    chk("Table[i + 1, {i, v}]", "{1, 4, 1, 1, 0, 1}");
    chk("Table[{i}, {i, v}]", "{{0}, {3}, {0}, {0}, {-1}, {0}}");
    chk("Table[i, {i, b}]", "{5, 7, 7, 7}");
    chk("Table[i, {i, e}]", "{0, 0, 0}");
    chk("Normal[Table[i^2, {i, m}]]", "{{1, 0, 0, 0}, {0, 0, 6.25, 0}, {0, x^2, 0, 0}}");
    chk("Table[i j, {i, v}, {j, 2}]", "{{0, 0}, {3, 6}, {0, 0}, {0, 0}, {-1, -2}, {0, 0}}");
    chk("Sum[i, {i, v}]", "2");
    chk("Sum[i^2, {i, v}]", "10");
    chk_hn("Sum[i, {i, m}]", "SparseArray", "{1, x, 2.5, 0}");
    chk("Sum[i, {i, e}]", "0");
    chk("Product[i + 1, {i, v}]", "0");
    chk("Product[i, {i, b}]", "1715");
    chk("Module[{l = {}}, Do[AppendTo[l, i], {i, v}]; l]", "{0, 3, 0, 0, -1, 0}");
    chk("Module[{c = 0}, Do[If[i == -1, Break[]]; c++, {i, v}]; c]", "4");
    chk("Module[{c = 0}, Do[If[i == 0, Continue[]]; c++, {i, v}]; c]", "2");
    chk("Module[{l = {}}, Do[AppendTo[l, Total[r]], {r, mi}]; l]", "{1, 2, -1}");
    chk("Catch[Do[If[i < 0, Throw[i]], {i, v}]]", "-1");
    chk("Reap[Do[Sow[i], {i, v}]]", "{Null, {{0, 3, 0, 0, -1, 0}}}");
    chk("f[] := Do[If[i == 3, Return[i]], {i, v}]; f[]", "3");
    chk("Table[i, {i, SparseArray[{}, {0}]}]", "{}");
}

static void test_ordering_sets_folds(void) {
    chk_hn("Sort[v]", "SparseArray", "{-1, 0, 0, 0, 0, 3}");
    chk_hn("Sort[b]", "SparseArray", "{5, 7, 7, 7}");
    chk_hn("SortBy[v, Minus]", "SparseArray", "{3, 0, 0, 0, 0, -1}");
    chk_hn("Sort[v, Greater]", "SparseArray", "{3, 0, 0, 0, 0, -1}");
    chk_hn("ReverseSort[v]", "SparseArray", "{3, 0, 0, 0, 0, -1}");
    chk("Union[v]", "{-1, 0, 3}");
    chk("Union[v, {7}]", "{-1, 0, 3, 7}");
    chk("DeleteDuplicates[v]", "{0, 3, -1}");
    chk("DeleteDuplicates[m]", "{{1, 0, 0, 0}, {0, 0, 2.5, 0}, {0, x, 0, 0}}");
    chk("Split[v]", "{{0}, {3}, {0, 0}, {-1}, {0}}");
    chk("Gather[v]", "{{0, 0, 0, 0}, {3}, {-1}}");
    chk("GatherBy[v, Sign]", "{{0, 0, 0, 0}, {3}, {-1}}");
    chk("Ordering[v]", "{5, 1, 3, 4, 6, 2}");
    chk("Counts[v]", "<|0 -> 4, 3 -> 1, -1 -> 1|>");
    chk("Nearest[v, 2]", "{3}");
    chk("Fold[Plus, 0, v]", "2");
    chk("FoldList[Plus, 0, v]", "{0, 0, 3, 3, 3, 2, 2}");
    chk("Fold[f, x, SparseArray[{2 -> 3}, 3]]", "f[f[f[x, 0], 3], 0]");
    chk("FoldList[f, SparseArray[{1 -> 1}, 2]]", "{1, f[1, 0]}");
    chk_hn("Fold[Plus, m]", "SparseArray", "{1, x, 2.5, 0}");
    /* a SparseArray initial value is a value, and stays sparse */
    chk_hn("Fold[Plus, v, {v, v}]", "SparseArray", "{0, 9, 0, 0, -3, 0}");
    chk("AnyTrue[v, Negative]", "True");
    chk("AllTrue[v, NonNegative]", "False");
    chk("AllTrue[e, # == 0 &]", "True");
    chk("NoneTrue[v, # > 5 &]", "True");
    chk("MapAt[f, v, 2]", "{0, f[3], 0, 0, -1, 0}");
    chk("MapAt[f, m, {2, 3}]", "{{1, 0, 0, 0}, {0, 0, f[2.5], 0}, {0, x, 0, 0}}");
    chk("MapAt[f, m, {All, 1}]", "{{f[1], 0, 0, 0}, {f[0], 0, 2.5, 0}, {f[0], x, 0, 0}}");
    chk("MapAt[f, m, 2]", "{{1, 0, 0, 0}, f[{0, 0, 2.5, 0}], {0, x, 0, 0}}");
    chk("Select[v, # != 0 &]", "{3, -1}");
    chk("Length[Subsets[SparseArray[{1 -> 1}, 3]]]", "8");
    /* a genuine atom: the reference's ::normal, unevaluated -- not the atom */
    chk("Quiet[Sort[x]]", "Sort[x]");
    chk("Check[Sort[x], failed]", "failed");
    chk("Check[Union[x], failed]", "failed");
    chk("Check[Split[x], failed]", "failed");
    chk("Check[DeleteDuplicates[x], failed]", "failed");
    chk("Check[SortBy[x, f], failed]", "failed");
    chk("Quiet[SortBy[5, f]]", "SortBy[5, f]");
}

static void test_mutation(void) {
    chk_hn("Module[{w = v}, AppendTo[w, 9]; w]", "SparseArray", "{0, 3, 0, 0, -1, 0, 9}");
    chk_hn("Module[{w = v}, PrependTo[w, 9]; w]", "SparseArray", "{9, 0, 3, 0, 0, -1, 0}");
    chk_hn("Module[{w = m}, AppendTo[w, {1, 2, 3, 4}]; w]", "SparseArray",
           "{{1, 0, 0, 0}, {0, 0, 2.5, 0}, {0, x, 0, 0}, {1, 2, 3, 4}}");
    chk("Module[{w = v}, AppendTo[w, {1}]; w]", "{0, 3, 0, 0, -1, 0, {1}}");
    chk("Module[{w = v}, AppendTo[w, 9]]", "SparseArray[<3>, {7}]");
    chk_hn("Module[{w = m}, w[[All, 1]] = SparseArray[{p, q, r}]; w]", "SparseArray",
           "{{p, 0, 0, 0}, {q, 0, 2.5, 0}, {r, x, 0, 0}}");
    chk("Check[Module[{w = m}, w[[All, 1]] = SparseArray[{p, q, r}]; Normal[w]], failed]",
        "{{p, 0, 0, 0}, {q, 0, 2.5, 0}, {r, x, 0, 0}}");
    chk_hn("Module[{w = m}, w[[2]] = SparseArray[{4 -> 9}, 4]; w]", "SparseArray",
           "{{1, 0, 0, 0}, {0, 0, 0, 9}, {0, x, 0, 0}}");
    chk_hn("Module[{w = v}, w[[2]] += 1; w]", "SparseArray", "{0, 4, 0, 0, -1, 0}");
    chk_hn("Module[{w = v}, w += 1; w]", "SparseArray", "{1, 4, 1, 1, 0, 1}");
    chk_hn("Module[{w = v}, w *= 2; w]", "SparseArray", "{0, 6, 0, 0, -2, 0}");
    chk("Module[{w = v, i}, For[i = 1, i <= 3, i++, w[[i]] = i]; Normal[w]]",
        "{1, 2, 3, 0, -1, 0}");
    chk("Module[{w = v, i = 1, s = 0}, While[i <= Length[w], s += w[[i]]; i++]; s]", "2");
    chk("Block[{w = v}, AppendTo[w, 1]; Length[w]]", "7");
    chk("With[{w = v}, Total[w]]", "2");
    /* an assignment that breaks the array leaves a List, as on the dense form */
    chk_hn("Module[{w = mi}, w[[1]] = x; w]", "List", "{x, {0, 0, 2, 0}, {0, -1, 0, 0}}");
    chk_hn("Module[{w = mi}, w[[All]] = 0; w]", "List", "{0, 0, 0}");
    chk_hn("Module[{w = mi}, w[[2]] = {1, 2}; w]", "List", "{{1, 0, 0, 0}, {1, 2}, {0, -1, 0, 0}}");
    chk_hn("Module[{w = mi}, w[[1]] = f[w[[1]]]; w]", "List",
           "{f[{1, 0, 0, 0}], {0, 0, 2, 0}, {0, -1, 0, 0}}");
    chk_hn("Module[{w = mi}, w[[2, 2]] = 5; w]", "SparseArray",
           "{{1, 0, 0, 0}, {0, 5, 2, 0}, {0, -1, 0, 0}}");
    /* AppendTo still refuses a genuine atom */
    chk("Module[{w = 5}, Head[Quiet[AppendTo[w, 1]]]]", "AppendTo");
}

static void test_thread_nesting(void) {
    chk("Thread[f[v, Range[6]]]", "{f[0, 1], f[3, 2], f[0, 3], f[0, 4], f[-1, 5], f[0, 6]}");
    chk("Thread[f[v]]", "{f[0], f[3], f[0], f[0], f[-1], f[0]}");
    chk("Thread[{v, Range[6]}]", "{{0, 1}, {3, 2}, {0, 3}, {0, 4}, {-1, 5}, {0, 6}}");
    chk("Thread[f[mi, {p, q, r}]]", "{f[{1, 0, 0, 0}, p], f[{0, 0, 2, 0}, q], f[{0, -1, 0, 0}, r]}");
    chk("Thread[v == Range[6]]", "False");   /* Equal evaluates first */
    chk_true("Quiet[Thread[f[v, {1, 2}]]] === Quiet[Thread[f[vn, {1, 2}]]]");
    chk_hn("Flatten[{v, v}]", "List", "{0, 3, 0, 0, -1, 0, 0, 3, 0, 0, -1, 0}");
    chk_hn("Flatten[{v, {1}}]", "List", "{0, 3, 0, 0, -1, 0, 1}");
    chk_hn("Flatten[{mi}, 1]", "List", "{{1, 0, 0, 0}, {0, 0, 2, 0}, {0, -1, 0, 0}}");
    chk_hn("Outer[f, {a, c}, v]", "List",
           "{{f[a, 0], f[a, 3], f[a, 0], f[a, 0], f[a, -1], f[a, 0]}, "
           "{f[c, 0], f[c, 3], f[c, 0], f[c, 0], f[c, -1], f[c, 0]}}");
    chk_hn("Outer[f, v, {a}]", "List", "{{f[0, a]}, {f[3, a]}, {f[0, a]}, {f[0, a]}, {f[-1, a]}, {f[0, a]}}");
    chk_hn("Outer[Times, v, {1, 2}]", "SparseArray",
           "{{0, 0}, {3, 6}, {0, 0}, {0, 0}, {-1, -2}, {0, 0}}");
    chk_hn("Outer[Times, {1, 2}, v]", "List", "{{0, 3, 0, 0, -1, 0}, {0, 6, 0, 0, -2, 0}}");
    chk("Head[Outer[Times, v, v]]", "SparseArray");
    chk_inv("Outer[Times, v, v]", "Outer[Times, vn, vn]");
    chk_hn("Outer[Plus, SparseArray[{1 -> 1}, 2], SparseArray[{2 -> 1}, 2]]", "List",
           "{{1, 2}, {0, 1}}");
    chk_hn("Outer[List, {a}, {v}]", "List", "{{{{a, 0}, {a, 3}, {a, 0}, {a, 0}, {a, -1}, {a, 0}}}}");
    /* Times with a dense array keeps a zero background's sparsity */
    chk_hn("v Range[6]", "SparseArray", "{0, 6, 0, 0, -5, 0}");
    chk_hn("Range[6] v", "SparseArray", "{0, 6, 0, 0, -5, 0}");
    chk_hn("m {1, 2, 3}", "SparseArray", "{{1, 0, 0, 0}, {0, 0, 5.0, 0}, {0, 3 x, 0, 0}}");
    chk_hn("2 v {p1, p2, p3, p4, p5, p6}", "SparseArray", "{0, 6 p2, 0, 0, -2 p5, 0}");
    chk_hn("v Range[6.]", "SparseArray", "{0, 6.0, 0, 0, -5.0, 0}");
    chk_hn("v + Range[6]", "List", "{1, 5, 3, 4, 4, 6}");
    chk_hn("b Range[4]", "List", "{5, 14, 21, 28}");
    chk_hn("v {1, 2, Infinity, 4, 5, 6}", "List", "{0, 6, Indeterminate, 0, -5, 0}");
    chk_hn("v x", "SparseArray", "{0, 3 x, 0, 0, -x, 0}");
    chk("Check[v Range[5], failed]", "failed");
    /* sparse operands of nested shapes broadcast as Listable does */
    run("v3 = SparseArray[{1 -> 2, 3 -> -1}, 3]; m3 = SparseArray[{{1, 1} -> 1, {2, 3} -> 2, {3, 2} -> -1}, {3, 3}]");
    chk_hn("m3 + v3", "SparseArray", "{{3, 2, 2}, {0, 0, 2}, {-1, -2, -1}}");
    chk_hn("v3 + m3", "SparseArray", "{{3, 2, 2}, {0, 0, 2}, {-1, -2, -1}}");
    chk_hn("m3 v3", "SparseArray", "{{2, 0, 0}, {0, 0, 0}, {0, 1, 0}}");
    chk_inv("t + SparseArray[{1 -> 5}, 2]", "tn + {5, 0}");
    chk_inv("f[m3, v3] /. f -> Plus", "Normal[m3] + Normal[v3]");
}

static void test_clip_background(void) {
    chk_hn("Clip[v - 1, {0, 10}]", "SparseArray", "{0, 2, 0, 0, 0, 0}");
    chk("Clip[v - 1, {0, 10}][\"Background\"]", "0");
    chk_hn("Clip[b, {0, 6}]", "SparseArray", "{5, 6, 6, 6}");
    chk("Clip[b, {0, 6}][\"Background\"]", "6");
    chk_hn("Clip[v, {0, 1}, {-9, 9}]", "SparseArray", "{0, 9, 0, 0, -9, 0}");
    chk_hn("Chop[SparseArray[{1 -> 1.*^-12, 2 -> 1.}, 3]]", "SparseArray", "{0, 1.0, 0}");
    chk_hn("FixedPoint[Clip[# - 1, {0, 10}] &, v]", "SparseArray", "{0, 0, 0, 0, 0, 0}");
    chk("Length[FixedPointList[Clip[# - 1, {0, 10}] &, v]]", "5");
    chk_true("Clip[Clip[v - 1, {0, 10}] - 1, {0, 10}] === Clip[Clip[v - 2, {0, 10}], {0, 10}]");
    /* a re-sparsified answer whose leading background vanished */
    chk("Rescale[v][\"Background\"] // Head", "Rational");
}

static void test_composition(void) {
    chk_hn("Nest[2 # &, v, 3]", "SparseArray", "{0, 24, 0, 0, -8, 0}");
    chk("Length[NestList[Rest, v, 3]]", "4");
    chk_hn("NestWhile[Rest, v, Length[#] > 3 &]", "SparseArray", "{0, -1, 0}");
    chk_hn("Composition[Reverse, Sort][v]", "SparseArray", "{3, 0, 0, 0, 0, -1}");
    chk("Through[{Total, Length}[v]]", "{2, 6}");
    chk_hn("Function[y, y^2, Listable][v]", "SparseArray", "{0, 9, 0, 0, 1, 0}");
    chk_hn("(#^2 + 1 &) /@ v", "SparseArray", "{1, 10, 1, 1, 2, 1}");
    chk("MapThread[f, {v, Range[6]}]", "{f[0, 1], f[3, 2], f[0, 3], f[0, 4], f[-1, 5], f[0, 6]}");
    chk("Inner[Times, v, Range[6], Plus]", "1");
    /* the pattern heads see an atom, as in the reference */
    chk("Count[m, 0, {2}]", "0");
    chk("Position[v, 3]", "{}");
    chk("Cases[m, _?Positive, {2}]", "{}");
    chk("MatchQ[v, _SparseArray]", "True");
}

/* Dense-side fixes the sparse comparisons exposed (the sparse answer must
 * equal the dense one, so the dense one has to be right). */
static void test_dense_fixes(void) {
    /* Rational and Complex are atoms to the level walkers */
    chk("Map[f, {1/2, 1 + 2 I, x/2}, {-1}]", "{f[1/2], f[1 + 2*I], f[1/2] f[x]}");
    chk("Map[f, {1/2}, Infinity]", "{f[1/2]}");
    chk("Apply[g, {1/2}, {1}]", "{1/2}");
    chk("Reap[Scan[Sow, {1/2, 2 + I}, {-1}]]", "{Null, {{1/2, 2 + I}}}");
    chk("MapIndexed[f, {1/2}, {-1}]", "{f[1/2, {1}]}");
    /* a packed list's [[All]] = scalar replaces the rows, as for any List */
    chk("Module[{w = {{1, 2}, {3, 4}}}, w[[All]] = 0; w]", "{0, 0}");
    chk("Module[{w = {{1, 2}, {3, 4}}}, w[[1 ;; 2]] = 0; w]", "{0, 0}");
    chk("Module[{w = {{1, 2}, {3, 4}}}, w[[1]] = 9; w]", "{9, {3, 4}}");
    chk("Module[{w = {{1, 2}, {3, 4}}}, w[[All, 1]] = 9; w]", "{{9, 2}, {9, 4}}");
    chk("Module[{w = {{1, 2}, {3, 4}}}, w[[1]] = {7, 8}; w]", "{{7, 8}, {3, 4}}");
    chk("Module[{w = Range[5]}, w[[2 ;; 3]] = 0; w]", "{1, 0, 0, 4, 5}");
    /* Map with All is {0, Infinity} */
    chk("Map[f, {{1}}, All]", "f[{f[{f[1]}]}]");
}

static void test_large(void) {
    run("big = SparseArray[{{1, 1} -> 1, {10^5, 10^5} -> 2}, {10^5, 10^5}]");
    run("bigv = SparseArray[{1 -> 1, 10^9 -> 2}, 10^9]");
    /* leaf-level Map never densifies */
    chk("With[{r = Map[f, big, {2}]}, {Head[r], r[\"Background\"], Length[r[\"NonzeroValues\"]]}]",
        "{SparseArray, f[0], 2}");
    chk("Map[Sin, big, {-1}][\"NonzeroValues\"]", "{Sin[1], Sin[2]}");
    chk("Map[g, bigv][[10^9]]", "g[2]");
    /* a level that needs every position is refused, once, not attempted */
    chk("Check[Scan[f, big, {2}], refused]", "refused");
    chk("Check[Table[i, {i, bigv}], refused]", "refused");
    chk("Check[Sort[bigv], refused]", "refused");
    /* rows of a large matrix iterate without densifying */
    chk("Module[{c = 0}, Do[c += Total[r], {r, Take[big, 3]}]; c]", "1");
    chk("Head[bigv Range[0]]", "Times");
}

int main(void) {
    symtab_init();
    core_init();
    setup();
    test_map_levels();
    test_mapindexed_apply_scan();
    test_iterators();
    test_ordering_sets_folds();
    test_mutation();
    test_thread_nesting();
    test_clip_background();
    test_composition();
    test_dense_fixes();
    test_large();
    if (failures) {
        fprintf(stderr, "\n%d sparse functional test(s) FAILED\n", failures);
        return 1;
    }
    printf("All sparse functional tests passed.\n");
    return 0;
}
