/* test_sparse_struct.c -- the native structural SparseArray kernels
 * (src/sparse/sparse_struct.c): Transpose, ConjugateTranspose, Flatten,
 * ArrayReshape, Partition, Join, Append/Prepend, Take/Drop/Rest/Most/First/
 * Last/Reverse, Diagonal, DiagonalMatrix, IdentityMatrix[n, SparseArray],
 * KroneckerProduct, ArrayFlatten, Upper/LowerTriangularize, ArrayPad/PadLeft/
 * PadRight and Map over a vector.
 *
 * Three kinds of check:
 *   - the exact canonical form (InputForm) of small answers, which pins the
 *     value kind, the background and the stored positions;
 *   - the KIND of the answer (sparse or dense);
 *   - the invariant Normal[op[s]] === op[Normal[s]] over random arrays, large
 *     enough to take the threaded paths.
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

/* The InputForm (canonical form) of a SparseArray-valued input. */
static void chk_if(const char* input, const char* expected) {
    char buf[4096];
    snprintf(buf, sizeof buf, "InputForm[%s]", input);
    chk(buf, expected);
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

/* op answers a SparseArray whose Normal is op applied to the Normal forms. */
static void chk_sparse_inv(const char* op) {
    char buf[4096];
    snprintf(buf, sizeof buf,
             "With[{r = %s}, Head[r] === SparseArray && Normal[r] === (%s /. a_SparseArray :> Normal[a])]",
             op, op);
    struct Expr* p = parse_expression(buf);
    struct Expr* e = evaluate(p);
    expr_free(p);
    char* s = expr_to_string(e);
    if (strcmp(s, "True") != 0) {
        fprintf(stderr, "FAIL(invariant): %s\n", op);
        failures++;
    }
    free(s);
    expr_free(e);
}

/* ------------------------------------------------------------------------ */

static void setup(void) {
    run("s = SparseArray[{{1,2}->3, {2,1}->4, {3,3}->5}, {3,4}]");
    run("v = SparseArray[{2->7, 4->1}, 5]");
    run("sx = SparseArray[{{1,2}->3, {2,1}->4}, {3,4}, x]");
    run("vx = SparseArray[{2->7}, 3, x]");
    run("s3 = SparseArray[{{1,2,3}->1, {2,1,1}->2}, {2,3,4}]");
    run("q = SparseArray[{{1,1}->1, {2,2}->2}, {2,2}]");
    run("s0 = SparseArray[{{1,2}->3., {2,1}->4.}, {2,2}, 0.]");
    run("sc = SparseArray[{{1,2}->1+2I, {2,1}->3.+I}, {2,2}]");
}

static void test_transpose(void) {
    chk_if("Transpose[s]",
           "SparseArray[Automatic, {4, 3}, 0, {1, {{0, 1, 2, 3, 3}, {{2}, {1}, {3}}}, {4, 3, 5}}]");
    chk_if("Transpose[sx]",
           "SparseArray[Automatic, {4, 3}, x, {1, {{0, 1, 2, 2, 2}, {{2}, {1}}}, {4, 3}}]");
    chk_if("Transpose[v]", "SparseArray[Automatic, {5}, 0, {1, {{0, 2}, {{2}, {4}}}, {7, 1}}]");
    chk("Normal[Transpose[s3, {2, 3, 1}]]",
        "{{{0, 0, 0}, {2, 0, 0}}, {{0, 0, 0}, {0, 0, 0}}, {{0, 1, 0}, {0, 0, 0}}, {{0, 0, 0}, {0, 0, 0}}}");
    chk("Normal[Transpose[s3, {2, 1}]]",
        "{{{0, 0, 0, 0}, {2, 0, 0, 0}}, {{0, 0, 1, 0}, {0, 0, 0, 0}}, {{0, 0, 0, 0}, {0, 0, 0, 0}}}");
    chk_true("Transpose[s3, {1, 2, 3}] === s3");
    chk("Normal[ConjugateTranspose[sc]]", "{{0, 3.0 - 1.0*I}, {1 - 2*I, 0}}");
    chk("ConjugateTranspose[SparseArray[{{1,2}->1.+2.I}, {2,2}, 1.+I]][\"Background\"]", "1.0 - 1.0*I");
    chk_true("Normal[ConjugateTranspose[s3, {2, 1, 3}]] === Transpose[Normal[s3], {2, 1, 3}]");
    /* repeated levels are the builtin's (a diagonal), not this path's */
    chk_sparse_inv("Transpose[SparseArray[{{1,1}->1, {2,2}->2}, {2,2}], {1, 1}]");
}

static void test_flatten_reshape(void) {
    chk_if("Flatten[s]", "SparseArray[Automatic, {12}, 0, {1, {{0, 3}, {{2}, {5}, {11}}}, {3, 4, 5}}]");
    chk("Normal[Flatten[s3, 1]]",
        "{{0, 0, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 0}, {2, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}}");
    chk_true("Flatten[s, 0] === s && Flatten[v] === v");
    chk("Normal[Flatten[s3, {{2}, {1, 3}}]]",
        "{{0, 0, 0, 0, 2, 0, 0, 0}, {0, 0, 1, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0}}");
    chk("Dimensions[Flatten[s3, {{3, 1}}]]", "{8, 3}");
    chk_true("Normal[Flatten[s, {{2}, {1}}]] === Transpose[Normal[s]]");
    chk_if("ArrayReshape[s, {2, 6}]",
           "SparseArray[Automatic, {2, 6}, 0, {1, {{0, 2, 3}, {{2}, {5}, {5}}}, {3, 4, 5}}]");
    chk("Normal[ArrayReshape[s, {5, 5}, 1]]",
        "{{0, 3, 0, 0, 4}, {0, 0, 0, 0, 0}, {5, 0, 1, 1, 1}, {1, 1, 1, 1, 1}, {1, 1, 1, 1, 1}}");
    chk_if("ArrayReshape[s, {2, 2}]", "SparseArray[Automatic, {2, 2}, 0, {1, {{0, 1, 1}, {{2}}}, {3}}]");
    chk("ArrayReshape[sx, {2, 2, 3}][\"Background\"]", "x");
    chk("Normal[Partition[v, 2, 1]]", "{{0, 7}, {7, 0}, {0, 1}, {1, 0}}");
    chk("Normal[Partition[s, 2, 1]]",
        "{{{0, 3, 0, 0}, {4, 0, 0, 0}}, {{4, 0, 0, 0}, {0, 0, 5, 0}}}");
    chk("Normal[Partition[v, 3, 2]]", "{{0, 7, 0}, {0, 1, 0}}");
    chk("Partition[v, 7]", "{}");
}

static void test_join_append(void) {
    chk_if("Join[v, v]",
           "SparseArray[Automatic, {10}, 0, {1, {{0, 4}, {{2}, {4}, {7}, {9}}}, {7, 1, 7, 1}}]");
    chk("Normal[Join[s, s, 2]]",
        "{{0, 3, 0, 0, 0, 3, 0, 0}, {4, 0, 0, 0, 4, 0, 0, 0}, {0, 0, 5, 0, 0, 0, 5, 0}}");
    chk_true("Head[Join[s, s]] === SparseArray && Normal[Join[s, s]] === Join[Normal[s], Normal[s]]");
    chk_true("Normal[Join[q, q, q]] === Join[Normal[q], Normal[q], Normal[q]]");
    /* backgrounds equal in value join; the first one is the answer's */
    chk("Join[v, SparseArray[{1->2.}, 2, 0.]][\"Background\"]", "0");
    chk("Normal[Join[q, s0, 2]]", "{{1, 0, 0, 3.0}, {0, 2, 4.0, 0}}");
    /* a dense operand, or arrays of other shapes, give a dense answer */
    chk("Head[Join[s, Normal[s]]]", "List");
    chk("Join[s, q] /. a_SparseArray :> Normal[a]",
        "{{0, 3, 0, 0}, {4, 0, 0, 0}, {0, 0, 5, 0}, {1, 0}, {0, 2}}");
    chk("Normal[Append[s, {1, 2, 3, 4}]]",
        "{{0, 3, 0, 0}, {4, 0, 0, 0}, {0, 0, 5, 0}, {1, 2, 3, 4}}");
    chk("Normal[Prepend[s, SparseArray[{1->9}, 4, x]]]",
        "{{9, x, x, x}, {0, 3, 0, 0}, {4, 0, 0, 0}, {0, 0, 5, 0}}");
    chk_if("Append[v, x]",
           "SparseArray[Automatic, {6}, 0, {1, {{0, 3}, {{2}, {4}, {6}}}, {7, 1, x}}]");
    chk_if("Prepend[vx, x]", "SparseArray[Automatic, {4}, x, {1, {{0, 1}, {{3}}}, {7}}]");
    chk("Append[s, {1, 2}] /. a_SparseArray :> Normal[a]",
        "{{0, 3, 0, 0}, {4, 0, 0, 0}, {0, 0, 5, 0}, {1, 2}}");
    chk("Normal[Append[s3, {{1, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 2}}]][[3]]",
        "{{1, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 2}}");
}

static void test_slices(void) {
    chk_if("Take[s, 2]", "SparseArray[Automatic, {2, 4}, 0, {1, {{0, 1, 2}, {{2}, {1}}}, {3, 4}}]");
    chk("Normal[Take[s, {2, 3}, {1, 2}]]", "{{4, 0}, {0, 0}}");
    chk("Normal[Take[v, {2, 4, 2}]]", "{7, 1}");
    chk("Normal[Take[v, {-1, 1, -2}]]", "{0, 0, 0}");
    chk("Normal[Take[s, All, {2}]]", "{{3}, {0}, {0}}");
    chk("Normal[Drop[s, 1, 1]]", "{{0, 0, 0}, {0, 5, 0}}");
    chk("Normal[Drop[s, {1, 3, 2}]]", "{{4, 0, 0, 0}}");
    chk("Normal[Drop[v, {2}]]", "{0, 0, 1, 0}");
    chk("Normal[Drop[s, None, {2}]]", "{{0, 0, 0}, {4, 0, 0}, {0, 5, 0}}");
    chk_if("Rest[v]", "SparseArray[Automatic, {4}, 0, {1, {{0, 2}, {{1}, {3}}}, {7, 1}}]");
    chk("Normal[Most[s]]", "{{0, 3, 0, 0}, {4, 0, 0, 0}}");
    chk_if("First[s]", "SparseArray[Automatic, {4}, 0, {1, {{0, 1}, {{2}}}, {3}}]");
    chk("Normal[Last[s3]]", "{{2, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}}");
    chk("First[v]", "0");
    chk("Last[vx]", "x");
    chk("Normal[Reverse[s]]", "{{0, 0, 5, 0}, {4, 0, 0, 0}, {0, 3, 0, 0}}");
    chk("Normal[Reverse[s, 2]]", "{{0, 0, 3, 0}, {0, 0, 0, 4}, {0, 5, 0, 0}}");
    chk("Normal[Reverse[s, {1, 2}]]", "{{0, 5, 0, 0}, {0, 0, 0, 4}, {0, 0, 3, 0}}");
    chk_true("Reverse[s, 3] === s");
    chk_true("Normal[Reverse[s3, {2, 3}]] === Reverse[Normal[s3], {2, 3}]");
    chk_true("Normal[Take[s3, All, 2, -1]] === Take[Normal[s3], All, 2, -1]");
    /* empty selections are the dense empty List */
    chk("Take[v, 0]", "{}");
    chk("Rest[SparseArray[{1 -> 1}, 1]]", "{}");
    chk("Drop[v, 5]", "{}");
}

static void test_diagonal(void) {
    chk_if("Diagonal[s]", "SparseArray[Automatic, {3}, 0, {1, {{0, 1}, {{3}}}, {5}}]");
    chk("Normal[Diagonal[s, 1]]", "{3, 0, 0}");
    chk("Normal[Diagonal[s, -1]]", "{4, 0}");
    chk("Diagonal[s, 5]", "{}");
    chk_if("Diagonal[sx]", "SparseArray[Automatic, {3}, x, {1, {{0, 0}, {}}, {}}]");
    chk("Normal[DiagonalMatrix[v, 1]][[2]]", "{0, 0, 7, 0, 0, 0}");
    chk("Normal[DiagonalMatrix[v, -1, 3]]", "{{0, 0, 0}, {0, 0, 0}, {0, 7, 0}}");
    chk("Normal[DiagonalMatrix[v, 0, {5, 7}]][[4]]", "{0, 0, 0, 1, 0, 0, 0}");
    /* the answer's background is 0: a non-zero vector background is stored */
    chk_if("DiagonalMatrix[vx]",
           "SparseArray[Automatic, {3, 3}, 0, {1, {{0, 1, 2, 3}, {{1}, {2}, {3}}}, {x, 7, x}}]");
    chk_if("DiagonalMatrix[SparseArray[{1->2}, 3, 1.]]",
           "SparseArray[Automatic, {3, 3}, 0.0, {1, {{0, 1, 2, 3}, {{1}, {2}, {3}}}, {2.0, 1.0, 1.0}}]");
    chk_if("IdentityMatrix[3, SparseArray]",
           "SparseArray[Automatic, {3, 3}, 0, {1, {{0, 1, 2, 3}, {{1}, {2}, {3}}}, {1, 1, 1}}]");
    chk("Normal[IdentityMatrix[{2, 3}, SparseArray]]", "{{1, 0, 0}, {0, 1, 0}}");
    chk("Head[IdentityMatrix[{3, 2}, SparseArray]]", "SparseArray");
}

static void test_kronecker_blocks(void) {
    chk_if("KroneckerProduct[q, q]",
           "SparseArray[Automatic, {4, 4}, 0, {1, {{0, 1, 2, 3, 4}, {{1}, {2}, {3}, {4}}}, {1, 2, 2, 4}}]");
    chk_true("Normal[KroneckerProduct[s, s]] === ArrayFlatten[Outer[Times, Normal[s], Normal[s]]]");
    chk_true("Head[KroneckerProduct[s, Normal[s]]] === SparseArray");
    chk_true("Normal[KroneckerProduct[q, q, q]] === Normal[KroneckerProduct[KroneckerProduct[q, q], q]]");
    /* a machine operand makes the product machine */
    chk_if("KroneckerProduct[q, s0]",
           "SparseArray[Automatic, {4, 4}, 0.0, {1, {{0, 1, 2, 3, 4}, {{2}, {1}, {4}, {3}}}, {3.0, 4.0, 6.0, 8.0}}]");
    chk("KroneckerProduct[SparseArray[{{1,1}->1, {2,2}->2.}, {2,2}], q][\"Background\"]", "0.0");
    /* vectors: a column when leading, a row after */
    chk("Dimensions[KroneckerProduct[SparseArray[{1->2}, 2], s]]", "{6, 4}");
    chk("Dimensions[KroneckerProduct[s, SparseArray[{1->2}, 2]]]", "{3, 8}");
    chk("Normal[KroneckerProduct[SparseArray[{1->2}, 2], v]]", "{{0, 14, 0, 2, 0}, {0, 0, 0, 0, 0}}");
    /* a non-zero background makes the product dense */
    chk("KroneckerProduct[SparseArray[{{1,1}->2}, {2,2}, 1], q]",
        "{{2, 0, 1, 0}, {0, 4, 0, 2}, {1, 0, 1, 0}, {0, 2, 0, 2}}");
    chk("Normal[ArrayFlatten[{{s, 0}, {0, s}}]][[4]]", "{0, 0, 0, 0, 0, 3, 0, 0}");
    chk_true("Normal[ArrayFlatten[{{q, q}, {q, q}}]] === ArrayFlatten[{{Normal[q], Normal[q]}, {Normal[q], Normal[q]}}]");
    chk("ArrayFlatten[{{sx, x}, {sx, sx}}][\"Background\"]", "x");
    /* a dense block, a scalar other than the background, or backgrounds that
     * differ: a dense answer */
    chk("Head[ArrayFlatten[{{q, Normal[q]}}]]", "List");
    chk("ArrayFlatten[{{q, 1}}]", "{{1, 0, 1}, {0, 2, 1}}");
    chk("Head[ArrayFlatten[{{q, s0}}]]", "List");
    chk("Head[ArrayFlatten[{{sx, 0}}]]", "List");
}

static void test_triangular_pad(void) {
    chk("Normal[UpperTriangularize[s]]", "{{0, 3, 0, 0}, {0, 0, 0, 0}, {0, 0, 5, 0}}");
    chk("Normal[UpperTriangularize[s, 1]]", "{{0, 3, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}}");
    chk("Normal[LowerTriangularize[s, -1]]", "{{0, 0, 0, 0}, {4, 0, 0, 0}, {0, 0, 0, 0}}");
    chk("UpperTriangularize[s0][\"Background\"]", "0.0");
    chk("UpperTriangularize[sx]", "{{x, 3, x, x}, {0, x, x, x}, {0, 0, x, x}}");
    chk("Normal[ArrayPad[s, 1]][[2]]", "{0, 0, 3, 0, 0, 0}");
    chk("Normal[ArrayPad[s, {{1, 0}, {0, 2}}]]",
        "{{0, 0, 0, 0, 0, 0}, {0, 3, 0, 0, 0, 0}, {4, 0, 0, 0, 0, 0}, {0, 0, 5, 0, 0, 0}}");
    chk("Normal[ArrayPad[s, 1, 7]][[1]]", "{7, 7, 7, 7, 7, 7}");
    chk("ArrayPad[s, 1, 7][\"Background\"]", "0");
    chk("Normal[ArrayPad[sx, 1]][[2]]", "{0, x, 3, x, x, 0}");
    chk("Normal[ArrayPad[s, -1]]", "{{0, 0}}");
    chk("Normal[ArrayPad[v, {1, 2}, x]]", "{x, 0, 7, 0, 1, 0, x, x}");
    chk("Normal[PadRight[s, {4, 5}, 9]][[4]]", "{9, 9, 9, 9, 9}");
    chk("Normal[PadLeft[s, {4, 5}]][[2]]", "{0, 0, 3, 0, 0}");
    chk("Normal[PadLeft[v, 3]]", "{0, 1, 0}");
    chk("Normal[PadRight[vx, 5, x]]", "{x, 7, x, x, x}");
    chk("Normal[PadRight[s, {2, 2}]]", "{{0, 3}, {4, 0}}");
    /* a scalar length on a matrix pads level 1 with scalars: dense */
    chk("PadRight[s, 4]", "{{0, 3, 0, 0}, {4, 0, 0, 0}, {0, 0, 5, 0}, 0}");
}

static void test_map(void) {
    chk_if("Map[f, v]",
           "SparseArray[Automatic, {5}, f[0], {1, {{0, 2}, {{2}, {4}}}, {f[7], f[1]}}]");
    chk_if("Map[Sin, SparseArray[{1->1.}, 3]]",
           "SparseArray[Automatic, {3}, 0, {1, {{0, 1}, {{1}}}, {0.8414709848078965}}]");
    chk("Map[f, vx][\"Background\"]", "f[x]");
    chk("Map[{#} &, v]", "{{0}, {7}, {0}, {1}, {0}}");
    chk("Map[f, s] /. a_SparseArray :> Normal[a]",
        "{f[{0, 3, 0, 0}], f[{4, 0, 0, 0}], f[{0, 0, 5, 0}]}");
}

/* Random arrays large enough for the threaded kernels (>= 2^17 entries). */
static void test_large(void) {
    run("SeedRandom[7]; big = SparseArray[RandomInteger[{1, 600}, {150000, 2}] -> RandomReal[{-1, 1}, 150000], {600, 700}]");
    run("bigI = SparseArray[RandomInteger[{1, 600}, {150000, 2}] -> RandomInteger[{-9, 9}, 150000], {600, 600}]");
    run("bigC = SparseArray[RandomInteger[{1, 400}, {150000, 2}] -> RandomComplex[{-1 - I, 1 + I}, 150000], {400, 500}]");
    chk_true("Normal[Transpose[big]] === Transpose[Normal[big]]");
    chk_true("Normal[Transpose[bigI]] === Transpose[Normal[bigI]]");
    chk_true("Normal[ConjugateTranspose[bigC]] === ConjugateTranspose[Normal[bigC]]");
    chk_true("Normal[Flatten[big]] === Flatten[Normal[big]]");
    chk_true("Normal[ArrayReshape[big, {700, 600}]] === ArrayReshape[Normal[big], {700, 600}]");
    chk_true("Normal[Join[big, big]] === Join[Normal[big], Normal[big]]");
    chk_true("Normal[Join[big, big, 2]] === Join[Normal[big], Normal[big], 2]");
    chk_true("Normal[Take[big, {100, 400, 3}, {-1, 1, -2}]] === Take[Normal[big], {100, 400, 3}, {-1, 1, -2}]");
    chk_true("Normal[Reverse[big, {1, 2}]] === Reverse[Normal[big], {1, 2}]");
    run("nI = Normal[bigI]");
    chk_true("Normal[UpperTriangularize[bigI, 3]] === Table[If[j - i >= 3, nI[[i, j]], 0], {i, 600}, {j, 600}]");
    chk_true("Normal[LowerTriangularize[bigI, -2]] === Table[If[j - i <= -2, nI[[i, j]], 0], {i, 600}, {j, 600}]");
    /* 1-D Laplacians: 3 * 600 entries each, so the product (3.2 * 10^6
     * entries) takes the threaded fill; checked on a block of rows */
    run("k1 = SparseArray[{{i_, i_} -> 2, {i_, j_} /; Abs[i - j] == 1 -> -1}, {600, 600}]");
    run("k2 = SparseArray[{{i_, i_} -> 2., {i_, j_} /; Abs[i - j] == 1 -> -1.}, {30, 30}]");
    run("kk = KroneckerProduct[k1, k1]");
    chk("Dimensions[kk]", "{360000, 360000}");
    chk("Table[Normal[kk[[1200 + j]]] === Flatten[Outer[Times, Normal[k1[[3]]], Normal[k1[[j]]]]], {j, 3}]",
        "{True, True, True}");
    chk_true("Normal[KroneckerProduct[k2, k2]] == ArrayFlatten[Outer[Times, Normal[k2], Normal[k2]]]");
    chk_true("Normal[ArrayFlatten[{{big, big}, {big, 0}}]] === ArrayFlatten[{{Normal[big], Normal[big]}, {Normal[big], 0}}]");
    chk("Total[IdentityMatrix[10^6, SparseArray][\"NonzeroValues\"]]", "1000000");
}

int main(void) {
    symtab_init();
    core_init();
    setup();
    test_transpose();
    test_flatten_reshape();
    test_join_append();
    test_slices();
    test_diagonal();
    test_kronecker_blocks();
    test_triangular_pad();
    test_map();
    test_large();
    if (failures) {
        fprintf(stderr, "\n%d sparse structural test(s) FAILED\n", failures);
        return 1;
    }
    printf("All sparse structural tests passed.\n");
    return 0;
}
