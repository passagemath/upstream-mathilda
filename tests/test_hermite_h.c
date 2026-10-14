/* Tests for HermiteH: the physicists' Hermite polynomial H_n(x).
 *
 * Covers exact integer polynomials (the H_10 / H_7 headline cases and the
 * Table[n,0..3] family), exact and fixed-point integer values, the special
 * value H_n(0) = 2^n Sqrt[Pi]/Gamma[(1-n)/2], Listable threading, machine and
 * arbitrary-precision (MPFR) numerics (real and complex, including the
 * non-integer/complex order via the Hypergeometric1F1 closed form), the stable
 * recurrence that avoids polynomial cancellation at large order
 * (HermiteH[100, N[5,16]]), the D[] derivative rule (2 n H_{n-1}), the finite
 * and at-Infinity Series expansions and SeriesCoefficient (concrete and
 * symbolic index), symbolic fall-through, argument-count diagnostics, and
 * attributes. */

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

void test_hermite_h_polynomials() {
    assert_eval_eq("Table[HermiteH[n, x], {n, 0, 3}]",
                   "{1, 2 x, -2 + 4 x^2, -12 x + 8 x^3}", 0);
    assert_eval_eq("HermiteH[0, x]", "1", 0);
    assert_eval_eq("HermiteH[1, x]", "2 x", 0);
    assert_eval_eq("HermiteH[2, x]", "-2 + 4 x^2", 0);
    /* Headline examples from the Wolfram reference. */
    assert_eval_eq("HermiteH[10, x]",
                   "-30240 + 302400 x^2 - 403200 x^4 + 161280 x^6 - 23040 x^8 + 1024 x^10", 0);
    assert_eval_eq("HermiteH[7, x]",
                   "-1680 x + 3360 x^3 - 1344 x^5 + 128 x^7", 0);
}

void test_hermite_h_exact_values() {
    assert_eval_eq("HermiteH[2, 3]", "34", 0);      /* 4*9 - 2 */
    assert_eval_eq("HermiteH[3, 5]", "940", 0);     /* 8*125 - 12*5 */
    assert_eval_eq("HermiteH[0, 0]", "1", 0);
    assert_eval_eq("HermiteH[1, 0]", "0", 0);
    /* Fixed-point integer table (Wolfram reference). */
    assert_eval_eq("Table[HermiteH[10, x], {x, 1, 5}]",
                   "{8224, 200416, -3093984, 125984224, 3275529760}", 0);
}

void test_hermite_h_special_values() {
    /* H_n(0) = 2^n Sqrt[Pi]/Gamma[(1-n)/2] for symbolic order. */
    assert_eval_eq("HermiteH[n, 0]", "(Sqrt[Pi] 2^n)/Gamma[1/2 (1 - n)]", 0);
    /* Half-integer order at 0 collapses to the Gamma form. */
    assert_eval_eq("HermiteH[1/2, 0]", "(Sqrt[2] Sqrt[Pi])/Gamma[1/4]", 0);
}

/* ---- Listable threading --------------------------------------------- */

void test_hermite_h_listable() {
    assert_eval_eq("HermiteH[3, {1, 2, 3}]", "{-4, 40, 180}", 0);
    assert_eval_eq("HermiteH[{0, 1, 2, 3}, x]",
                   "{1, 2 x, -2 + 4 x^2, -12 x + 8 x^3}", 0);
    assert_eval_eq("HermiteH[2, {}]", "{}", 0);
}

/* ---- machine-precision numerics ------------------------------------- */

void test_hermite_h_machine() {
    /* Integer order, inexact argument. */
    assert_close("HermiteH[3., 5]", 940.0, 1e-9);
    /* An inexact integer-valued order yields a clean real. */
    assert_close("HermiteH[2., 3.]", 34.0, 1e-9);
    /* Non-integer order, real argument (Hypergeometric1F1 form). */
    assert_close("HermiteH[0.71, 0.87]", 1.55872, 1e-4);
    /* Cancellation case: evaluating directly (the stable recurrence) gives the
     * right ~8.19e98, where substituting into the expanded polynomial at an
     * inexact x would lose all significance. */
    assert_eval_startswith("HermiteH[100, N[5, 16]]", "8.19319371305915");
}

void test_hermite_h_complex() {
    /* Complex non-integer order (the 1F1 closed form, with a cancellation guard
     * so the machine answer keeps its digits): HermiteH[3-4.583I, 5]. */
    assert_close("Re[HermiteH[3. - 4.583 I, 5]]", -731.082342, 1e-2);
    assert_close("Im[HermiteH[3. - 4.583 I, 5]]", 940.003462, 1e-2);
    /* Complex order and complex argument. */
    assert_close("Re[HermiteH[5.3 + I, .8 + I]]", 4.81841, 1e-3);
    assert_close("Im[HermiteH[5.3 + I, .8 + I]]", 157.49, 1e-2);
}

/* ---- arbitrary-precision numerics ----------------------------------- */

void test_hermite_h_highprec() {
    /* Non-integer order, exact rational arguments, via N (50 digits). */
    assert_eval_startswith(
        "N[HermiteH[2/3, 8/7], 50]",
        "1.795978863238586639435290581248609522279819602573");
    /* Complex non-integer order at MPFR precision tracks the input precision and
     * the cancellation guard keeps the digits. */
    assert_eval_startswith("Re[HermiteH[3.`30 - 4.583`30 I, 5]]",
                           "-731.08234218869585221");
}

/* ---- symbolic fall-through ------------------------------------------ */

void test_hermite_h_symbolic() {
    assert_eval_eq("HermiteH[n, x]", "HermiteH[n, x]", 0);
    /* Exact non-integer order with symbolic argument stays symbolic. */
    assert_eval_eq("HermiteH[1/2, x]", "HermiteH[1/2, x]", 0);
    /* Negative-integer order with an exact argument stays symbolic. */
    assert_eval_eq("HermiteH[-1, y]", "HermiteH[-1, y]", 0);
}

/* ---- derivative ----------------------------------------------------- */

void test_hermite_h_derivative() {
    /* d/dz H_n(z) = 2 n H_{n-1}(z). */
    assert_eval_eq("D[HermiteH[n, z], z]", "2 n HermiteH[-1 + n, z]", 0);
    /* Chain rule. */
    assert_eval_eq("D[HermiteH[n, x^2], x]", "4 n x HermiteH[-1 + n, x^2]", 0);
    /* The derivative of an explicit polynomial is its derivative:
     * d/dx (4 x^2 - 2) = 8 x. */
    assert_eval_eq("D[HermiteH[2, x], x]", "8 x", 0);
    assert_eval_eq("D[HermiteH[10, x], x]",
                   "604800 x - 1612800 x^3 + 967680 x^5 - 184320 x^7 + 10240 x^9", 0);
}

/* ---- Series and SeriesCoefficient ----------------------------------- */

void test_hermite_h_series() {
    /* Taylor at a generic point (free via the derivative rule). */
    assert_eval_eq("Series[HermiteH[n, x], {x, x0, 2}]",
        "HermiteH[n, x0] + 2 n HermiteH[-1 + n, x0] (x - x0) + "
        "2 n (-1 + n) HermiteH[-2 + n, x0] (x - x0)^2 + O[x - x0]^3", 0);
    /* Asymptotic expansion at Infinity (symbolic order). */
    assert_eval_eq("Series[HermiteH[n, x], {x, Infinity, 2}]",
                   "(2 x)^n (1 + 1/4 n (1 - n) (1/x)^2 + O[1/x]^4)", 0);
    /* An integer-order polynomial expands via the generic path. */
    assert_eval_eq("Normal[Series[HermiteH[2, x], {x, 0, 4}]]", "-2 + 4 x^2", 0);
}

void test_hermite_h_seriescoefficient() {
    /* Concrete index from the polynomial. */
    assert_eval_eq("SeriesCoefficient[HermiteH[10, x], {x, 0, 4}]", "-403200", 0);
    /* Symbolic order, concrete index (m = 2). */
    assert_eval_eq("SeriesCoefficient[HermiteH[n, x], {x, 0, 2}]",
                   "(2 Sqrt[Pi] n (-1 + n) 2^(-2 + n))/Gamma[1/2 (3 - n)]", 0);
    /* Symbolic index yields a Piecewise; check it reduces to the concrete value
     * at a sample. H_7 = -1680 x + 3360 x^3 - 1344 x^5 + 128 x^7, so the x^4
     * coefficient is 0 and the x^3 coefficient is 3360. */
    assert_eval_eq(
        "SeriesCoefficient[HermiteH[n, x], {x, 0, m}] /. {n -> 7, m -> 4}", "0", 0);
    assert_eval_eq(
        "SeriesCoefficient[HermiteH[n, x], {x, 0, m}] /. {n -> 7, m -> 3}", "3360", 0);
}

/* ---- diagnostics & attributes --------------------------------------- */

void test_hermite_h_argcount() {
    /* Wrong arg counts stay unevaluated (an argrx/argx message goes to stderr). */
    assert_eval_eq("HermiteH[]", "HermiteH[]", 0);
    assert_eval_eq("HermiteH[3]", "HermiteH[3]", 0);
    assert_eval_eq("HermiteH[1, 2, 3]", "HermiteH[1, 2, 3]", 0);
}

void test_hermite_h_attributes() {
    SymbolDef* d = symtab_get_def("HermiteH");
    ASSERT_MSG(d != NULL, "HermiteH not registered");
    ASSERT_MSG((d->attributes & ATTR_LISTABLE) != 0, "HermiteH not Listable");
    ASSERT_MSG((d->attributes & ATTR_NUMERICFUNCTION) != 0,
               "HermiteH not NumericFunction");
    ASSERT_MSG((d->attributes & ATTR_PROTECTED) != 0, "HermiteH not Protected");
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_hermite_h_polynomials);
    TEST(test_hermite_h_exact_values);
    TEST(test_hermite_h_special_values);
    TEST(test_hermite_h_listable);
    TEST(test_hermite_h_machine);
    TEST(test_hermite_h_complex);
    TEST(test_hermite_h_highprec);
    TEST(test_hermite_h_symbolic);
    TEST(test_hermite_h_derivative);
    TEST(test_hermite_h_series);
    TEST(test_hermite_h_seriescoefficient);
    TEST(test_hermite_h_argcount);
    TEST(test_hermite_h_attributes);

    printf("All HermiteH tests passed.\n");
    return 0;
}
