/* Tests for ChebyshevU: the Chebyshev polynomial of the second kind U_n(x).
 *
 * Covers exact integer polynomials (including the U_10 / U_7 headline cases and
 * the negative-order reflection U_{-m} = -U_{m-2}), exact values at exact
 * arguments, the special values U_n(0) = Cos[n Pi/2] and U_n(1) = n+1, Listable
 * threading, machine and arbitrary-precision (MPFR) numerics (real and complex,
 * including the non-integer/complex order via the 2F1 and ArcCos closed forms and
 * the branch cut), the stable-recurrence path that avoids polynomial cancellation
 * for large integer order, the D[] derivative rule, symbolic fall-through,
 * argument-count diagnostics, and attributes. */

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

/* ---- exact integer polynomials -------------------------------------- */

void test_chebyshev_u_polynomials() {
    assert_eval_eq("Table[ChebyshevU[n, x], {n, 0, 4}]",
                   "{1, 2 x, -1 + 4 x^2, -4 x + 8 x^3, 1 - 12 x^2 + 16 x^4}", 0);
    assert_eval_eq("ChebyshevU[0, x]", "1", 0);
    assert_eval_eq("ChebyshevU[1, x]", "2 x", 0);
    assert_eval_eq("ChebyshevU[2, x]", "-1 + 4 x^2", 0);
    /* Headline examples from the Wolfram reference. */
    assert_eval_eq("ChebyshevU[10, x]",
                   "-1 + 60 x^2 - 560 x^4 + 1792 x^6 - 2304 x^8 + 1024 x^10", 0);
    assert_eval_eq("ChebyshevU[7, x]",
                   "-8 x + 80 x^3 - 192 x^5 + 128 x^7", 0);
}

void test_chebyshev_u_negative_order() {
    /* U_{-1} = 0; U_{-m} = -U_{m-2}, folded distributed into the coefficients. */
    assert_eval_eq("ChebyshevU[-1, x]", "0", 0);
    assert_eval_eq("ChebyshevU[-2, x]", "-1", 0);            /* -U_0 */
    assert_eval_eq("ChebyshevU[-3, x]", "-2 x", 0);          /* -U_1 */
    assert_eval_eq("ChebyshevU[-4, x]", "1 - 4 x^2", 0);     /* -U_2 */
    assert_eval_eq("ChebyshevU[-1, 2.5]", "0", 0);
}

void test_chebyshev_u_exact_values() {
    assert_eval_eq("ChebyshevU[2, 3]", "35", 0);
    assert_eval_eq("ChebyshevU[2, 5]", "99", 0);
    assert_eval_eq("ChebyshevU[0, 0]", "1", 0);
    assert_eval_eq("ChebyshevU[2, 1/3]", "-5/9", 0);
    /* U_0 is constant for any argument. */
    assert_eval_eq("ChebyshevU[0, Infinity]", "1", 0);
}

void test_chebyshev_u_special_values() {
    /* U_n(0) = Cos[n Pi/2] for symbolic order. */
    assert_eval_eq("ChebyshevU[n, 0]", "Cos[1/2 n Pi]", 0);
    /* U_n(1) = n + 1 for symbolic order. */
    assert_eval_eq("ChebyshevU[n, 1]", "1 + n", 0);
    /* Table at Pi/2 (integer order, symbolic argument). */
    assert_eval_eq("Table[ChebyshevU[n, Pi/2], {n, 0, 3}]",
                   "{1, Pi, -1 + Pi^2, -2 Pi + Pi^3}", 0);
}

/* ---- Listable threading --------------------------------------------- */

void test_chebyshev_u_listable() {
    assert_eval_eq("ChebyshevU[2, {1, 2, 3}]", "{3, 15, 35}", 0);
    assert_eval_eq("ChebyshevU[{0, 1, 2, 3}, x]",
                   "{1, 2 x, -1 + 4 x^2, -4 x + 8 x^3}", 0);
    assert_eval_eq("ChebyshevU[2, {}]", "{}", 0);
}

/* ---- machine-precision numerics ------------------------------------- */

void test_chebyshev_u_machine() {
    /* Integer order, inexact argument. */
    assert_close("ChebyshevU[12, 2.5]", 153027601.0, 1e3);
    /* Cancellation case: the stable recurrence gives ~1.15684, NOT 0. The
     * explicit degree-60 polynomial substituted at an inexact x loses all
     * significance here. */
    assert_close("ChebyshevU[60, 0.6666666666666666]", 1.1568372648685888, 1e-4);
    /* U_n(0) at inexact order evaluates through Cos[n Pi/2]. */
    assert_close("ChebyshevU[0.3, 0]", 0.8910065241883679, 1e-12);
    /* An inexact integer-valued order yields a clean real, not a complex with a
     * spurious tiny imaginary part. */
    assert_close("ChebyshevU[2., 3.]", 35.0, 1e-9);
    assert_close("ChebyshevU[12., 2.5]", 153027601.0, 1e3);
}

void test_chebyshev_u_complex() {
    /* Non-integer complex order at x = 3 (outside the 2F1 disk; the ArcCos form
     * answers): ChebyshevU[2-4I, 3] = 25.18033... - 24.31072... I. */
    assert_close("Re[ChebyshevU[2 - 4 I, 3.0]]", 25.1803326153886, 1e-6);
    assert_close("Im[ChebyshevU[2 - 4 I, 3.0]]", -24.3107208960466, 1e-6);
    /* Branch cut (-Infinity, -1): the ArcCos form takes Mathematica's branch,
     * Im ~ +0.597614 (the algebraic power form would give the wrong -0.383). */
    assert_close("Im[ChebyshevU[-1.5, -2.4]]", 0.5976143046, 1e-6);
}

/* ---- arbitrary-precision numerics ----------------------------------- */

void test_chebyshev_u_highprec() {
    /* Non-integer order, exact rational arguments, via N. */
    assert_eval_startswith(
        "N[ChebyshevU[1/7, 1/3], 50]",
        "1.0464308547237641232170382442932095700192279304889");
    /* Complex non-integer order at MPFR precision, boundary of the 2F1 disk
     * (the ArcCos form answers with the Mathematica branch). */
    assert_eval_startswith("ChebyshevU[2 - 4 I, 3`40]",
                           "25.1803326153885556188606621144528330");
}

/* ---- symbolic fall-through ------------------------------------------ */

void test_chebyshev_u_symbolic() {
    assert_eval_eq("ChebyshevU[n, x]", "ChebyshevU[n, x]", 0);
    /* Exact non-integer order with exact/symbolic argument stays symbolic
     * (half-integer closed forms are deferred). */
    assert_eval_eq("ChebyshevU[1/2, x]", "ChebyshevU[1/2, x]", 0);
    assert_eval_eq("ChebyshevU[3/2, 2]", "ChebyshevU[3/2, 2]", 0);
}

/* ---- derivative ----------------------------------------------------- */

void test_chebyshev_u_derivative() {
    assert_eval_eq(
        "D[ChebyshevU[n, x], x]",
        "((-1 - n) ChebyshevU[-1 + n, x] + n x ChebyshevU[n, x])/(-1 + x^2)", 0);
    /* The derivative of an explicit polynomial is the polynomial's derivative. */
    assert_eval_eq("D[ChebyshevU[2, x], x]", "8 x", 0);
}

/* ---- diagnostics & attributes --------------------------------------- */

void test_chebyshev_u_argcount() {
    /* Wrong arg counts stay unevaluated (an argrx/argx message goes to stderr). */
    assert_eval_eq("ChebyshevU[]", "ChebyshevU[]", 0);
    assert_eval_eq("ChebyshevU[3]", "ChebyshevU[3]", 0);
    assert_eval_eq("ChebyshevU[1, 2, 3]", "ChebyshevU[1, 2, 3]", 0);
}

void test_chebyshev_u_attributes() {
    SymbolDef* d = symtab_get_def("ChebyshevU");
    ASSERT_MSG(d != NULL, "ChebyshevU not registered");
    ASSERT_MSG((d->attributes & ATTR_LISTABLE) != 0, "ChebyshevU not Listable");
    ASSERT_MSG((d->attributes & ATTR_NUMERICFUNCTION) != 0,
               "ChebyshevU not NumericFunction");
    ASSERT_MSG((d->attributes & ATTR_PROTECTED) != 0, "ChebyshevU not Protected");
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_chebyshev_u_polynomials);
    TEST(test_chebyshev_u_negative_order);
    TEST(test_chebyshev_u_exact_values);
    TEST(test_chebyshev_u_special_values);
    TEST(test_chebyshev_u_listable);
    TEST(test_chebyshev_u_machine);
    TEST(test_chebyshev_u_complex);
    TEST(test_chebyshev_u_highprec);
    TEST(test_chebyshev_u_symbolic);
    TEST(test_chebyshev_u_derivative);
    TEST(test_chebyshev_u_argcount);
    TEST(test_chebyshev_u_attributes);

    printf("All ChebyshevU tests passed.\n");
    return 0;
}
