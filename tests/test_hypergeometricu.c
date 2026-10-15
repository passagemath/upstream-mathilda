/* Tests for HypergeometricU[a, b, z], the Tricomi confluent hypergeometric U.
 *
 * Covers: the exact reductions U(0,b,z)=1 and U(a,a+1,z)=z^(-a); machine- and
 * arbitrary-precision (MPFR) real numerics across the integer-b (DLMF 13.2.9),
 * non-integer-b (Kummer connection) and large-z (asymptotic) regimes; complex
 * arguments; list threading; the NDArray / packed-array fast path (agreement with
 * the scalar path); Compile[]; the derivative rule; argument-count diagnostics;
 * attributes.
 *
 * Numeric reference values are from mpmath (mp.hyperu, 40 digits) -- an
 * independent implementation, not recalled constants.
 */

#include "eval.h"
#include "parse.h"
#include "expr.h"
#include "symtab.h"
#include "attr.h"
#include "core.h"
#include "print.h"
#include "test_utils.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ---- numeric helpers ------------------------------------------------- */

static double eval_real(const char* input) {
    Expr* e = parse_expression(input);
    ASSERT(e != NULL);
    Expr* r = evaluate(e);
    expr_free(e);
    ASSERT_MSG(r->type == EXPR_REAL, "%s: expected a Real result", input);
    double v = r->data.real;
    expr_free(r);
    return v;
}

static void assert_close(const char* input, double expected, double tol) {
    double v = eval_real(input);
    ASSERT_MSG(fabs(v - expected) <= tol,
               "%s: expected %.12g, got %.12g", input, expected, v);
}

/* ---- exact reductions and passthrough ------------------------------- */

void test_hyperu_exact() {
    assert_eval_eq("HypergeometricU[0, b, z]", "1", 0);            /* U(0,b,z) = 1 */
    assert_eval_eq("HypergeometricU[a, a + 1, z]", "z^(-a)", 0);   /* U(a,a+1,z) = z^-a */
    assert_eval_eq("HypergeometricU[1, 2, z]", "1/z", 0);
    assert_eval_eq("HypergeometricU[2, 3, 4]", "1/16", 0);
    assert_eval_eq("HypergeometricU[1, 2, 5]", "1/5", 0);
    /* All-symbolic / all-exact stays symbolic. */
    assert_eval_eq("HypergeometricU[3, 1, z]", "HypergeometricU[3, 1, z]", 0);
    assert_eval_eq("HypergeometricU[2, 5, 3]", "HypergeometricU[2, 5, 3]", 0);
}

/* ---- machine real numerics (all three regimes) ---------------------- */

void test_hyperu_machine_real() {
    /* integer b (DLMF 13.2.9 log-series). */
    assert_close("HypergeometricU[1.5, 7., 4.5]", 0.5549311339066449, 1e-10);
    assert_close("HypergeometricU[-1.5, 9., 1.7]", 373.9518871804518, 1e-8);
    assert_close("HypergeometricU[3., 2., 10.]", 0.0006199963612715087, 1e-12);
    /* non-integer b (Kummer connection). */
    assert_close("HypergeometricU[-0.278, 3.444, 1.0]", -0.4001016549990364, 1e-10);
    assert_close("HypergeometricU[2.3, 1.2, 0.7]", 0.2268899760107006, 1e-10);
    assert_close("HypergeometricU[0.5, 0.5, 2.0]", 0.5959060788258650, 1e-10);
    /* rational parameters with an inexact argument. */
    assert_close("HypergeometricU[0.6, 1/3, 1/7.]", 1.097260883254250, 1e-10);
    /* the log case b = 1. */
    assert_close("HypergeometricU[1., 1., 1.]", 0.5963473623231941, 1e-12);
}

/* ---- arbitrary precision (MPFR / FLINT) ----------------------------- */

void test_hyperu_arbitrary() {
    /* N[U(3/2,7,9/2),25]; first ~24 digits rounding-stable (mpmath). */
    assert_eval_startswith("N[HypergeometricU[3/2, 7, 9/2], 25]",
        "0.554931133906644949855594");
    /* precision tracks the input. */
    assert_eval_startswith("HypergeometricU[3/2, 7, 45/10`30]",
        "0.5549311339066449498555");
}

/* ---- complex argument ----------------------------------------------- */

void test_hyperu_complex() {
    /* U(3 - I, 2 I, 3.0) = -0.00526539 + 0.000793166 I (mpmath). */
    assert_close("Re[HypergeometricU[3 - I, 2 I, 3.0]]", -0.005265393320688201, 1e-12);
    assert_close("Im[HypergeometricU[3 - I, 2 I, 3.0]]", 0.0007931661643468159, 1e-12);
}

/* ---- list threading (Listable) -------------------------------------- */

void test_hyperu_listable() {
    assert_eval_eq("HypergeometricU[a, b, {p, q}]",
                   "{HypergeometricU[a, b, p], HypergeometricU[a, b, q]}", 0);
    assert_eval_eq("HypergeometricU[1, 2, {2., 4., 5.}]", "{0.5, 0.25, 0.2}", 0);
    assert_eval_eq("HypergeometricU[a, b, {}]", "{}", 0);
}

/* ---- NDArray / packed surface agreement ----------------------------- */

void test_hyperu_ndarray() {
    /* visible NDArray: element-wise, agreeing with the scalar path. */
    assert_close("HypergeometricU[1.5, 7., NDArray[{3.5, 4.5, 5.5}]][[2]]",
                 0.5549311339066449, 1e-10);
    /* a large numeric list auto-packs; the packed fast path shares the machine
     * helper, so elements match the scalar path bit-closely. */
    assert_close("HypergeometricU[1.5, 7., N[Range[200]/7]][[53]]",
                 eval_real("HypergeometricU[1.5, 7., N[53/7]]"), 1e-12);
    assert_close("HypergeometricU[2.3, 1.2, N[Range[200]/7]][[90]]",
                 eval_real("HypergeometricU[2.3, 1.2, N[90/7]]"), 1e-12);
}

/* ---- Compile[] ------------------------------------------------------ */

void test_hyperu_compile() {
    /* scalar form lowers (OP_KERNN via the N-ary machine kernel). */
    assert_eval_eq(
        "CompileDiagnostics[{{a, _Real}, {b, _Real}, {z, _Real}}, "
        "HypergeometricU[a, b, z]][[1]]", "\"Compiled\" -> True", 0);
    /* the compiled function runs and agrees with the scalar path. */
    assert_close(
        "Compile[{{a, _Real}, {b, _Real}, {z, _Real}}, HypergeometricU[a, b, z]][1.5, 7., 4.5]",
        0.5549311339066449, 1e-10);
}

/* ---- derivative ----------------------------------------------------- */

void test_hyperu_derivative() {
    assert_eval_eq("D[HypergeometricU[a, b, z], z]",
                   "-a HypergeometricU[1 + a, 1 + b, z]", 0);
    /* no elementary form for d/da, d/db -> inert Derivative. */
    assert_eval_eq("D[HypergeometricU[a, b, z], a]",
                   "Derivative[1, 0, 0][HypergeometricU][a, b, z]", 0);
    assert_eval_eq("D[HypergeometricU[a, b, z], b]",
                   "Derivative[0, 1, 0][HypergeometricU][a, b, z]", 0);
    /* no z-dependence -> 0. */
    assert_eval_eq("D[HypergeometricU[a, b, y], z]", "0", 0);
}

/* ---- argument count ------------------------------------------------- */

void test_hyperu_argcount() {
    assert_eval_eq("HypergeometricU[]", "HypergeometricU[]", 0);
    assert_eval_eq("HypergeometricU[1, 2]", "HypergeometricU[1, 2]", 0);
    assert_eval_eq("HypergeometricU[1, 2, 3, 4]", "HypergeometricU[1, 2, 3, 4]", 0);
}

/* ---- attributes ----------------------------------------------------- */

void test_hyperu_attributes() {
    SymbolDef* d = symtab_get_def("HypergeometricU");
    ASSERT_MSG(d != NULL, "HypergeometricU not registered");
    ASSERT_MSG((d->attributes & ATTR_LISTABLE) != 0, "HypergeometricU not Listable");
    ASSERT_MSG((d->attributes & ATTR_NUMERICFUNCTION) != 0, "HypergeometricU not NumericFunction");
    ASSERT_MSG((d->attributes & ATTR_PROTECTED) != 0, "HypergeometricU not Protected");
}

int main() {
    symtab_init();
    core_init();

    TEST(test_hyperu_exact);
    TEST(test_hyperu_machine_real);
    TEST(test_hyperu_arbitrary);
    TEST(test_hyperu_complex);
    TEST(test_hyperu_listable);
    TEST(test_hyperu_ndarray);
    TEST(test_hyperu_compile);
    TEST(test_hyperu_derivative);
    TEST(test_hyperu_argcount);
    TEST(test_hyperu_attributes);

    printf("All HypergeometricU tests passed.\n");
    return 0;
}
