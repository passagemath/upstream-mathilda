/* Tests for Unitize.
 *
 * Cover: zero / non-zero across every number type (Integer, bigint, Rational,
 * Real incl. -0., arbitrary precision, Complex), the infinities and
 * Indeterminate, closed-form numeric constants decided by the zero test,
 * symbolic pass-through, the threshold form (machine, exact, closed-form and
 * complex x; the Unitize::post and ::argt paths), Listable threading, the
 * packed / NDArray / SparseArray surfaces agreeing bit for bit with the List
 * path, Compile[] at scalar and rank-1 shapes (incl. the run-time threshold
 * guard), the {Listable, Protected} attributes, and a memory-hygiene loop.
 */

#include "core.h"
#include "eval.h"
#include "expr.h"
#include "parse.h"
#include "print.h"
#include "symtab.h"
#include "test_utils.h"

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ------------------------------------------------------------------------
 *  One argument: numbers
 * ---------------------------------------------------------------------- */

static void test_unitize_numbers(void) {
    assert_eval_eq("Unitize[0]", "0", 0);
    assert_eval_eq("Unitize[5]", "1", 0);
    assert_eval_eq("Unitize[-7]", "1", 0);
    assert_eval_eq("Unitize[10^30]", "1", 0);
    assert_eval_eq("Unitize[1/2]", "1", 0);
    assert_eval_eq("Unitize[2.5]", "1", 0);
    assert_eval_eq("Unitize[0.]", "0", 0);
    assert_eval_eq("Unitize[-0.]", "0", 0);
    assert_eval_eq("Unitize[1.*^-300]", "1", 0);
    assert_eval_eq("Unitize[N[0, 20]]", "0", 0);
    assert_eval_eq("Unitize[N[Pi, 30] - N[Pi, 30]]", "0", 0);
    assert_eval_eq("Unitize[N[Pi, 30]]", "1", 0);
    assert_eval_eq("Unitize[1 + I]", "1", 0);
    assert_eval_eq("Unitize[1.5 I]", "1", 0);
    assert_eval_eq("Unitize[0. + 0. I]", "0", 0);
    /* Always the exact Integer, never a Real. */
    assert_eval_eq("Unitize[{1., 0.}]", "{1, 0}", 0);
    assert_eval_eq("Head[Unitize[2.5]]", "Integer", 0);
}

static void test_unitize_infinities(void) {
    assert_eval_eq("Unitize[Infinity]", "1", 0);
    assert_eval_eq("Unitize[-Infinity]", "1", 0);
    assert_eval_eq("Unitize[ComplexInfinity]", "1", 0);
    assert_eval_eq("Unitize[Indeterminate]", "Unitize[Indeterminate]", 0);
}

/* ------------------------------------------------------------------------
 *  One argument: closed-form numerics and symbols
 * ---------------------------------------------------------------------- */

static void test_unitize_numeric_constants(void) {
    assert_eval_eq("Unitize[Pi]", "1", 0);
    assert_eval_eq("Unitize[Sin[Pi]]", "0", 0);
    assert_eval_eq("Unitize[Sin[Pi/7]]", "1", 0);
    assert_eval_eq("Unitize[E - Exp[1]]", "0", 0);
    /* Not syntactically zero: the zero test proves it. */
    assert_eval_eq("Unitize[Sqrt[2] + Sqrt[3] - Sqrt[5 + 2 Sqrt[6]]]", "0", 0);
    assert_eval_eq("Unitize[Sin[1]^2 + Cos[1]^2 - 1]", "0", 0);
    assert_eval_eq("Unitize[GoldenRatio^2 - GoldenRatio - 1]", "0", 0);
}

static void test_unitize_symbolic(void) {
    assert_eval_eq("Unitize[x]", "Unitize[x]", 0);
    assert_eval_eq("Unitize[x - x]", "0", 0);
    assert_eval_eq("Unitize[True]", "Unitize[True]", 0);
    assert_eval_eq("Unitize[\"a\"]", "Unitize[\"a\"]", 0);
}

/* ------------------------------------------------------------------------
 *  Threshold form
 * ---------------------------------------------------------------------- */

static void test_unitize_threshold(void) {
    assert_eval_eq("Unitize[0.001, 0.01]", "0", 0);
    assert_eval_eq("Unitize[{0.5, 2.}, 1]", "{0, 1}", 0);
    assert_eval_eq("Unitize[1., 1.]", "1", 0);          /* |x| == dx -> 1 */
    assert_eval_eq("Unitize[-3, 3]", "1", 0);
    assert_eval_eq("Unitize[2, 3]", "0", 0);
    assert_eval_eq("Unitize[0.5, 1/2]", "1", 0);
    assert_eval_eq("Unitize[1/2, 0.5]", "1", 0);
    assert_eval_eq("Unitize[10^30, 10^31]", "0", 0);
    assert_eval_eq("Unitize[1 + I, 1]", "1", 0);
    assert_eval_eq("Unitize[1 + I, 2]", "0", 0);
    assert_eval_eq("Unitize[Pi, 4]", "0", 0);
    assert_eval_eq("Unitize[1, Pi]", "0", 0);
    assert_eval_eq("Unitize[Pi - 355/113, 1/10]", "0", 0);
    assert_eval_eq("Unitize[Pi - 355/113, 1/1000000000]", "1", 0);
    assert_eval_eq("Unitize[Infinity, 2]", "1", 0);
    assert_eval_eq("Unitize[{1, 2}, {1, 3}]", "{1, 0}", 0);
    /* -2^63 with an integer threshold: no Abs overflow. */
    assert_eval_eq("Unitize[-9223372036854775807 - 1, 5]", "1", 0);
}

static void test_unitize_threshold_unevaluated(void) {
    /* Non-positive threshold: Unitize::post, left unevaluated. */
    assert_eval_eq("Unitize[1, 0]", "Unitize[1, 0]", 0);
    assert_eval_eq("Unitize[1, -1]", "Unitize[1, -1]", 0);
    assert_eval_eq("Unitize[0, 0.]", "Unitize[0, 0.0]", 0);
    assert_eval_eq("Check[Unitize[1, 0], failed]", "failed", 0);
    assert_eval_eq("Quiet[Unitize[1, -1]]", "Unitize[1, -1]", 0);
    /* Symbolic, infinite or complex threshold, symbolic x: silent. */
    assert_eval_eq("Check[Unitize[1, x], failed]", "Unitize[1, x]", 0);
    assert_eval_eq("Unitize[2, Infinity]", "Unitize[2, Infinity]", 0);
    assert_eval_eq("Unitize[1, I]", "Unitize[1, I]", 0);
    assert_eval_eq("Unitize[x, 1]", "Unitize[x, 1]", 0);
    assert_eval_eq("Unitize[Indeterminate, 1]", "Unitize[Indeterminate, 1]", 0);
}

static void test_unitize_arity(void) {
    assert_eval_eq("Unitize[]", "Unitize[]", 0);
    assert_eval_eq("Unitize[1, 2, 3]", "Unitize[1, 2, 3]", 0);
    assert_eval_eq("Check[Unitize[], failed]", "failed", 0);
}

/* ------------------------------------------------------------------------
 *  Packed / NDArray / SparseArray surfaces
 * ---------------------------------------------------------------------- */

static void test_unitize_ndarray(void) {
    assert_eval_eq("Unitize[NDArray[{1., 0., -0., 3.}]]", "NDArray[{1, 0, 0, 1}]", 0);
    assert_eval_eq("Unitize[NDArray[{1, 0, -2}, DataType -> \"int64\"]]",
                   "NDArray[{1, 0, 1}]", 0);
    assert_eval_eq("Unitize[NDArray[{1., 0., 3.}], 2]", "NDArray[{0, 0, 1}]", 0);
}

static void test_unitize_packed(void) {
    /* float64 and int64 buffers map to a packed int64 buffer that bit-matches
     * the unpacked List path (machine scalar helper == kernel). */
    assert_eval_eq("v = RandomReal[{-1, 1}, 2000]; v[[5]] = 0.; v[[7]] = -0.; "
                   "Developer`PackedArrayQ[v]", "True", 0);
    assert_eval_eq("Developer`PackedArrayQ[Unitize[v]]", "True", 0);
    assert_eval_eq("Unitize[v] === Map[Unitize, Developer`FromPackedArray[v]]", "True", 0);
    assert_eval_eq("Developer`PackedArrayQ[Unitize[v, 0.5]]", "True", 0);
    assert_eval_eq("Unitize[v, 0.5] === Map[Unitize[#, 0.5] &, Developer`FromPackedArray[v]]",
                   "True", 0);
    assert_eval_eq("w = RandomInteger[{-3, 3}, 2000]; Developer`PackedArrayQ[Unitize[w]]",
                   "True", 0);
    assert_eval_eq("Unitize[w] === Map[Unitize, Developer`FromPackedArray[w]]", "True", 0);
    assert_eval_eq("Developer`PackedArrayQ[Unitize[w, 2]]", "True", 0);
    assert_eval_eq("Unitize[w, 2] === Map[Unitize[#, 2] &, Developer`FromPackedArray[w]]",
                   "True", 0);
    /* Real threshold over an int64 buffer: declines to the List path. */
    assert_eval_eq("Unitize[w, 2.] === Map[Unitize[#, 2.] &, Developer`FromPackedArray[w]]",
                   "True", 0);
    assert_eval_eq("Total[Unitize[Range[-100., 100.]]]", "200", 0);
    /* Non-positive threshold over a packed buffer: the kernel declines and the
     * List path leaves each element unevaluated (with Unitize::post). */
    assert_eval_eq("Quiet[Unitize[Developer`ToPackedArray[{1., 2.}], 0]]",
                   "{Unitize[1.0, 0], Unitize[2.0, 0]}", 0);
}

static void test_unitize_sparse(void) {
    assert_eval_eq("Normal[Unitize[SparseArray[{1 -> 2., 3 -> 0.}, {4}]]] === "
                   "Unitize[{2., 0, 0., 0}]", "True", 0);
    assert_eval_eq("Normal[Unitize[SparseArray[{1 -> 2., 3 -> 0.}, {4}]]]",
                   "{1, 0, 0, 0}", 0);
}

/* ------------------------------------------------------------------------
 *  Compile[]
 * ---------------------------------------------------------------------- */

static void test_unitize_compile(void) {
    assert_eval_eq("\"Compiled\" /. CompileDiagnostics[{{x, _Real}}, Unitize[x]]", "True", 0);
    assert_eval_eq("\"Compiled\" /. CompileDiagnostics[{{v, _Real, 1}}, Unitize[v]]", "True", 0);
    assert_eval_eq("\"Compiled\" /. CompileDiagnostics[{{n, _Integer}}, Unitize[n]]", "True", 0);
    assert_eval_eq("\"Compiled\" /. CompileDiagnostics[{{v, _Integer, 1}}, Unitize[v]]", "True", 0);
    assert_eval_eq("\"Compiled\" /. CompileDiagnostics[{{z, _Complex}}, Unitize[z]]", "True", 0);
    assert_eval_eq("\"Compiled\" /. CompileDiagnostics[{{x, _Real}, {d, _Real}}, Unitize[x, d]]",
                   "True", 0);
    assert_eval_eq("\"Compiled\" /. CompileDiagnostics[{{n, _Integer}, {d, _Integer}}, Unitize[n, d]]",
                   "True", 0);
    assert_eval_eq("\"Compiled\" /. CompileDiagnostics[{{v, _Real, 1}}, Unitize[v, 0.5]]",
                   "True", 0);

    assert_eval_eq("cf = Compile[{{x, _Real}}, Unitize[x]]; {cf[0.], cf[-0.], cf[2.5]}",
                   "{0, 0, 1}", 0);
    assert_eval_eq("cc = Compile[{{z, _Complex}}, Unitize[z]]; {cc[0. + 0. I], cc[1. I]}",
                   "{0, 1}", 0);
    assert_eval_eq("cd = Compile[{{x, _Real}, {d, _Real}}, Unitize[x, d]]; "
                   "{cd[0.001, 0.01], cd[0.5, 0.5], cd[-2., 1.]}", "{0, 1, 1}", 0);
    assert_eval_eq("ci = Compile[{{x, _Integer}, {d, _Integer}}, Unitize[x, d]]; "
                   "{ci[-3, 3], ci[2, 3], ci[-9223372036854775807 - 1, 5]}", "{1, 0, 1}", 0);
    /* A run-time threshold <= 0 hands the call back to the interpreter. */
    assert_eval_eq("Quiet[cd[1., -1.]]", "Unitize[1.0, -1.0]", 0);
    assert_eval_eq("Quiet[ci[1, 0]]", "Unitize[1, 0]", 0);
    /* Array surfaces bit-match the interpreter. */
    assert_eval_eq("u = RandomReal[{-1, 1}, 500]; u[[3]] = 0.; "
                   "Compile[{{v, _Real, 1}}, Unitize[v]][u] === Unitize[u]", "True", 0);
    assert_eval_eq("Compile[{{v, _Real, 1}}, Unitize[v, 0.5]][u] === Unitize[u, 0.5]", "True", 0);
    assert_eval_eq("k = RandomInteger[{-2, 2}, 500]; "
                   "Compile[{{v, _Integer, 1}}, Unitize[v]][k] === Unitize[k]", "True", 0);
}

/* ------------------------------------------------------------------------
 *  Attributes
 * ---------------------------------------------------------------------- */

static void test_unitize_attributes(void) {
    assert_eval_eq("Attributes[Unitize]", "{Listable, Protected}", 0);
}

/* ------------------------------------------------------------------------
 *  Memory hygiene
 * ---------------------------------------------------------------------- */

static void test_unitize_memory_loop(void) {
    const char* cases[] = {
        "Unitize[0]", "Unitize[2.5]", "Unitize[1/2]", "Unitize[1 + I]",
        "Unitize[Pi]", "Unitize[Sqrt[2] + Sqrt[3] - Sqrt[5 + 2 Sqrt[6]]]",
        "Unitize[x]", "Unitize[Indeterminate]", "Unitize[Infinity]",
        "Unitize[0.001, 0.01]", "Unitize[Pi, 4]", "Unitize[1 + I, 2]",
        "Quiet[Unitize[1, 0]]", "Unitize[1, x]", "Unitize[{0.5, 2.}, 1]",
        "Quiet[Unitize[]]", "Unitize[Range[-5., 5.]]", "Unitize[Range[-5, 5], 2]",
        NULL
    };
    for (int rep = 0; rep < 20; rep++) {
        for (int i = 0; cases[i]; i++) {
            Expr* p = parse_expression(cases[i]);
            ASSERT(p != NULL);
            Expr* v = evaluate(p);
            expr_free(p);
            expr_free(v);
        }
    }
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_unitize_numbers);
    TEST(test_unitize_infinities);
    TEST(test_unitize_numeric_constants);
    TEST(test_unitize_symbolic);
    TEST(test_unitize_threshold);
    TEST(test_unitize_threshold_unevaluated);
    TEST(test_unitize_arity);
    TEST(test_unitize_ndarray);
    TEST(test_unitize_packed);
    TEST(test_unitize_sparse);
    TEST(test_unitize_compile);
    TEST(test_unitize_attributes);
    TEST(test_unitize_memory_loop);

    printf("All Unitize tests passed.\n");
    return 0;
}
