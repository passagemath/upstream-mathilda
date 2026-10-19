/* Tests for DirichletL[k, j, s].
 *
 * Reference values are Mathematica 14's (wolframscript, N[..., 20] or better);
 * the machine kernel is held to ~1e-13 relative, the MPFR path to its full
 * printed precision. Covers: the Zeta case k = 1, the exact closed forms
 * (non-positive integer s via generalized Bernoulli numbers; s = 1 for primitive
 * and imprimitive non-principal characters), machine real and complex s on both
 * sides of the critical line (Euler-Maclaurin and the functional equation),
 * real / complex / imprimitive characters, the s = 1 limit, arbitrary precision
 * and precision tracking, N keeping the labels exact, Interval enclosures,
 * Listable threading, the packed / NDArray / Compile surfaces (bit-for-bit
 * against the scalar path), messages and attributes. */

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

static void T(const char* in) { assert_eval_eq(in, "True", 0); }

static Expr* ev(const char* input) {
    Expr* e = parse_expression(input);
    ASSERT(e != NULL);
    Expr* r = evaluate(e);
    expr_free(e);
    return r;
}

static void assert_real_rel(const char* input, double expected, double rel) {
    Expr* r = ev(input);
    ASSERT_MSG(r->type == EXPR_REAL, "%s: expected a machine Real", input);
    double v = r->data.real;
    ASSERT_MSG(fabs(v - expected) <= rel * fabs(expected),
               "%s: expected %.17g, got %.17g", input, expected, v);
    expr_free(r);
}

static void assert_complex_rel(const char* input, double er, double ei, double rel) {
    Expr* r = ev(input);
    ASSERT_MSG(r->type == EXPR_FUNCTION && r->data.function.head->type == EXPR_SYMBOL
               && strcmp(r->data.function.head->data.symbol.name, "Complex") == 0
               && r->data.function.arg_count == 2
               && r->data.function.args[0]->type == EXPR_REAL
               && r->data.function.args[1]->type == EXPR_REAL,
               "%s: expected a machine Complex", input);
    double vr = r->data.function.args[0]->data.real;
    double vi = r->data.function.args[1]->data.real;
    double mag = hypot(er, ei);
    ASSERT_MSG(hypot(vr - er, vi - ei) <= rel * mag,
               "%s: expected %.17g %+.17g I, got %.17g %+.17g I", input, er, ei, vr, vi);
    expr_free(r);
}

/* ---- k = 1: Zeta ----------------------------------------------------- */

void test_dl_zeta() {
    T("DirichletL[1, 1, 2] === Pi^2/6");
    T("Table[DirichletL[1, 1, s], {s, -6, 6}] === {0, -1/252, 0, 1/120, 0, -1/12, -1/2, "
      "ComplexInfinity, Pi^2/6, Zeta[3], Pi^4/90, Zeta[5], Pi^6/945}");
    T("DirichletL[1, 1, {1, 2, 3, 4}] === {ComplexInfinity, Pi^2/6, Zeta[3], Pi^4/90}");
    T("DirichletL[1, 1, 2.] === Zeta[2.]");
    T("DirichletL[1, 1, 0.3 + 4 I] === Zeta[0.3 + 4 I]");
    assert_eval_startswith("DirichletL[1, 1, 2.`30]", "1.6449340668482264364724151666");
}

/* ---- exact closed forms ---------------------------------------------- */

void test_dl_exact_nonpositive() {
    /* Mathematica's table for k = 3: s = -4..0. */
    assert_eval_eq("Table[DirichletL[3, 2, s], {s, -4, 0}]", "{2/3, 0, -2/9, 0, 1/3}", 0);
    assert_eval_eq("Table[DirichletL[3, 1, s], {s, -4, 0}]", "{0, -13/60, 0, 1/6, 0}", 0);
    assert_eval_eq("Table[DirichletL[4, 2, s], {s, -4, 0}]", "{5/2, 0, -1/2, 0, 1/2}", 0);
    assert_eval_eq("Table[DirichletL[8, 4, s], {s, -4, 0}]", "{57, 0, -3, 0, 1}", 0);
    assert_eval_eq("Table[DirichletL[8, 2, s], {s, -4, 0}]", "{0, 11, 0, -1, 0}", 0);
    assert_eval_eq("Table[DirichletL[12, 3, s], {s, -4, 0}]", "{205, 0, -5, 0, 1}", 0);
    assert_eval_eq("DirichletL[12, 1, -3]", "91/60", 0);
    assert_eval_eq("DirichletL[7, 4, -4]", "32", 0);
    /* Complex characters: Gaussian rationals for order 4 ... */
    assert_eval_eq("DirichletL[5, 2, -4]", "148/25 + 86/25*I", 0);
    assert_eval_eq("DirichletL[5, 4, -2]", "-4/5 + 2/5*I", 0);
    assert_eval_eq("DirichletL[5, 2, 0]", "3/5 + 1/5*I", 0);
    /* ... and cyclotomic numbers otherwise (value checked numerically). */
    T("Abs[N[DirichletL[7, 2, -4]] - N[(-2*(-1)^(1/3)*(-34 - 79*(-1)^(1/3) + 67*(-1)^(2/3)))/7]] < 10^-12");
    T("Abs[N[DirichletL[7, 3, -3]] - N[-1/420*((-1)^(1/3)*(-599 + 1919*(-1)^(1/3) + 1321*(-1)^(2/3)))]] < 10^-12");
    /* s = -n against the numeric value from both sides. */
    T("Abs[N[DirichletL[45, 8, -3], 30] - DirichletL[45, 8, -3.`30]] < 10^-25");
    T("Abs[N[DirichletL[40, 6, -5], 30] - DirichletL[40, 6, -5.`30]] < 10^-20");
}

void test_dl_exact_one() {
    T("DirichletL[3, 2, 1] == Pi/(3 Sqrt[3])");
    T("DirichletL[4, 2, 1] == Pi/4");
    /* Mathematica's numeric values. */
    T("Abs[N[DirichletL[5, 3, 1], 20] - 0.43040894096400403888943323295060542542`20] < 10^-18");
    T("Abs[N[DirichletL[8, 2, 1], 20] - 0.62322524014023051339402008025056800265`20] < 10^-18");
    T("Abs[N[DirichletL[7, 3, 1], 20] - (0.53774738050490415906598555390898822376`20 - "
      "0.10529754563079555866579778611255217411`20 I)] < 10^-17");
    T("Abs[N[DirichletL[7, 2, 1], 20] - (0.80420572938678064453602680236665980362`20 + "
      "0.39866669881887309012652808390944874326`20 I)] < 10^-17");
    /* Imprimitive (Mathematica leaves these; the value is finite and exact). */
    T("Abs[N[DirichletL[24, 5, 1], 30] - DirichletL[24, 5, 1.`30]] < 10^-25");
    T("Abs[N[DirichletL[12, 2, 1], 30] - DirichletL[12, 2, 1.`30]] < 10^-25");
    /* Principal character mod k > 1 at s = 1: left alone, as Mathematica does. */
    assert_eval_eq("DirichletL[3, 1, 1]", "DirichletL[3, 1, 1]", 0);
    /* Positive integers s >= 2 and non-integers stay symbolic. */
    assert_eval_eq("DirichletL[3, 2, 2]", "DirichletL[3, 2, 2]", 0);
    assert_eval_eq("DirichletL[3, 2, 1/2]", "DirichletL[3, 2, 1/2]", 0);
    assert_eval_eq("DirichletL[3, 2, x]", "DirichletL[3, 2, x]", 0);
    assert_eval_eq("DirichletL[3, 2, 1/2 + I]", "DirichletL[3, 2, 1/2 + I]", 0);
}

/* ---- machine precision ----------------------------------------------- */

void test_dl_machine_real() {
    assert_real_rel("DirichletL[24, 5, 2.5]", 1.0094763148747815470, 1e-14);
    assert_real_rel("DirichletL[3, 2, 1.]", 0.60459978807807261686, 1e-14);      /* Pi/(3 Sqrt[3]) */
    assert_real_rel("DirichletL[24, 5, 1.0000001]", 1.0471975481551707048, 1e-13);  /* s -> 1, no cancellation */
    assert_real_rel("DirichletL[5, 3, 1.]", 0.43040894096400403889, 1e-14);
    assert_real_rel("DirichletL[3, 2, 0.5]", 0.48086755769682862618, 1e-14);
    assert_real_rel("DirichletL[5, 1, 1.23]", 4.2589889501084897507, 1e-13);
    /* Trivial zeros of the functional equation: exactly 0. */
    T("DirichletL[5, 3, -2.] === 0.");
    T("DirichletL[3, 2, -1.] === 0.");
}

void test_dl_machine_complex() {
    assert_complex_rel("DirichletL[2, 1, 1. + I]", 0.65435891740723345570, -0.38437635021479987609, 1e-14);
    assert_complex_rel("DirichletL[7, 2, 0.5]", 0.71394334376831949286, 0.47490218277139938264, 1e-14);
    /* Left half-plane via the functional equation. References are Mathematica's
     * N[..., 25] at the exact s: its own MACHINE values are off by up to 1e-9
     * relative here (7, 2, -7/2 + 10 I), so they are not used as references. */
    assert_complex_rel("DirichletL[7, 2, -3.5]", 4.2250481818457694345, 8.9783886305698590149, 1e-13);
    assert_complex_rel("DirichletL[7, 2, -3.5 + 10 I]", 15503.857766950871617, -6982.6912167172499759, 1e-13);
    assert_complex_rel("DirichletL[45, 8, 0.3 + 2 I]", -0.61080429689459664693, -0.76598560934285430571, 1e-13);
    assert_complex_rel("DirichletL[1000, 77, 1.5]", 0.84002610550110523317, 0.062162694867633627934, 1e-13);
    assert_complex_rel("DirichletL[5, 2, 1.5]", 0.92454210962186188202, 0.17701358631882574902, 1e-14);
    assert_complex_rel("DirichletL[5, 2, -2.5]", -0.83228142015092496780, -0.43943347553386544829, 1e-13);
    assert_complex_rel("DirichletL[45, 8, 0.7]", 1.2574782481722802360, -0.42842317814879006430, 1e-13);
    assert_complex_rel("DirichletL[45, 8, -0.7 + 3 I]", 28.190718984459987728, 14.358928715729273859, 1e-13);
    assert_complex_rel("DirichletL[20, 3, 0.5 + 14 I]", 1.4068306561625648022, 0.74744663475784097775, 1e-13);
    assert_complex_rel("DirichletL[13, 5, -5.5]", -3324.6073463682209524, -12222.291889806599664, 1e-13);
    assert_complex_rel("DirichletL[16, 2, 2.]", 0.97287220589824364879, -0.069139898897926274653, 1e-14);
    assert_complex_rel("DirichletL[16, 6, -1.5]", -2.2210543010682519025, -4.3617438609876288847, 1e-13);
    assert_complex_rel("DirichletL[100, 17, 0.75]", 0.77229889185147388591, -0.12130347010042804349, 1e-13);
    assert_complex_rel("DirichletL[9, 2, 1/3 + I/3 + 0.]", 1.1704822574537260181, 0.61779807092833671979, 1e-13);
    /* A real character at complex s, a complex character at real s. */
    T("Head[DirichletL[5, 3, 1.5 + 0.5 I]] === Complex");
    T("Head[DirichletL[5, 2, 1.5]] === Complex");
    T("Head[DirichletL[5, 3, 1.5]] === Real");
}

void test_dl_pole() {
    assert_eval_eq("DirichletL[5, 1, 1.]", "ComplexInfinity", 0);
    assert_eval_eq("DirichletL[1, 1, 1.]", "ComplexInfinity", 0);
}

/* ---- arbitrary precision --------------------------------------------- */

void test_dl_mpfr() {
    assert_eval_startswith("N[DirichletL[3, 2, 1/2], 50]",
                           "0.48086755769682862618122006323558987777682973083237");
    assert_eval_startswith("DirichletL[3, 2, 1.50000000000000000000000]",
                           "0.7039682448687332616676");
    T("Precision[DirichletL[3, 2, 1.50000000000000000000000]] > 20");
    T("Abs[DirichletL[7, 2, -3.5`30] - (4.22504818184576943452109639061240892151`17 + "
      "8.97838863056985901487075147069210985634`18 I)] < 10^-16");
    T("Abs[DirichletL[45, 8, 0.3`25 + 2 I] - (-0.61080429689459664693011518081554856897`22 - "
      "0.76598560934285430571199836367577561117`22 I)] < 10^-20");
    /* MPFR result agrees with the machine one. */
    T("Abs[DirichletL[100, 17, 0.75`40] - DirichletL[100, 17, 0.75]] < 10^-13");
    /* N leaves the labels exact. */
    T("FreeQ[Hold @@ {N[DirichletL[3, 2, s]]}, 3.]");
}

/* ---- Interval --------------------------------------------------------- */

void test_dl_interval() {
    /* Principal character mod 5 is decreasing on [1.23, 1.24]: the enclosure
     * is its endpoint values, rounded outward. */
    T("With[{iv = DirichletL[5, 1, Interval[{1.23, 1.24}]]}, "
      "Head[iv] === Interval && Min[iv] <= DirichletL[5, 1, 1.24] && "
      "Max[iv] >= DirichletL[5, 1, 1.23] && Max[iv] - Min[iv] < 0.1451]");
    /* A real non-principal character, on both sides of the critical line and
     * through a non-monotone stretch (left half-plane, functional equation):
     * the enclosure contains every sampled value and is tight. */
    T("With[{iv = DirichletL[5, 3, Interval[{1.1, 3.}]]}, "
      "And @@ Table[Min[iv] <= DirichletL[5, 3, x] <= Max[iv], {x, 1.1, 3., 0.05}] && "
      "Max[iv] - Min[iv] < 0.39]");
    T("With[{iv = DirichletL[5, 3, Interval[{0.1, 0.9}]]}, "
      "And @@ Table[Min[iv] <= DirichletL[5, 3, x] <= Max[iv], {x, 0.1, 0.9, 0.02}]]");
    T("With[{iv = DirichletL[5, 3, Interval[{-3., 0.1}]]}, "
      "And @@ Table[Min[iv] <= DirichletL[5, 3, x] <= Max[iv], {x, -3., 0.1, 0.01}] && "
      "Max[iv] - Min[iv] < 2.43]");
    /* A range through s = 1 cannot be enclosed by arb's sums: left alone. */
    T("Head[DirichletL[5, 3, Interval[{0.9, 1.1}]]] === DirichletL");
    T("Head[DirichletL[3, 1, Interval[{2., 3.}]]] === Interval");
    /* Complex characters have no real enclosure: left alone. */
    T("Head[DirichletL[5, 2, Interval[{1.23, 1.24}]]] === DirichletL");
}

/* ---- surfaces --------------------------------------------------------- */

void test_dl_listable_and_packed() {
    T("DirichletL[3, 2, {1, 2, -3}] === {Pi/(3 Sqrt[3]), DirichletL[3, 2, 2], 0}");
    T("PackedArrayQ[DirichletL[5, 3, Range[0.5, 5., 0.25]]]");
    T("PackedArrayQ[DirichletL[7, 2, Range[0.5, 5., 0.25]]]");
    /* a complex character on the real left half-line stays on the machine
     * kernel (the Gamma ratio of the functional equation is real there) */
    T("PackedArrayQ[DirichletL[7, 2, Range[-4., 0.25, 0.25]]]");
    /* (grids avoid the trivial zeros s = -1, -3, -5 of this odd character: a
     * packed complex 0 + 0 I materialises as Complex[0., 0.] while the scalar
     * canonicalises to 0., a core SameQ quirk unrelated to DirichletL) */
    T("With[{xs = Range[-3.9, 0.25, 0.125]}, Table[DirichletL[7, 2, x], {x, xs}] === DirichletL[7, 2, xs]]");
    T("Abs[DirichletL[7, 2, -3.3] - N[DirichletL[7, 2, -33/10], 30]] < 10^-13");
    /* bit-for-bit: scalar == unpacked List == packed == NDArray, real and complex characters */
    T("With[{xs = Range[-4.9, 10., 0.37]}, Table[DirichletL[5, 3, x], {x, xs}] === DirichletL[5, 3, xs] "
      "=== DirichletL[5, 3, FromPackedArray[xs]] === Normal[DirichletL[5, 3, NDArray[xs]]]]");
    T("With[{xs = Range[-4.9, 10., 0.37]}, Table[DirichletL[7, 2, x], {x, xs}] === DirichletL[7, 2, xs] "
      "=== DirichletL[7, 2, FromPackedArray[xs]] === Normal[DirichletL[7, 2, NDArray[xs]]]]");
    T("With[{zs = Table[x + 3.5 I, {x, -4.9, 10., 0.37}]}, Table[DirichletL[45, 8, z], {z, zs}] === DirichletL[45, 8, zs]]");
    /* a declined element (the pole at s = 1) degrades the buffer, never NaN */
    assert_eval_eq("DirichletL[3, 1, {1., 2.}][[1]]", "ComplexInfinity", 0);
}

void test_dl_compile() {
    T("(\"Compiled\" /. CompileDiagnostics[{{x, _Real}}, DirichletL[5, 3, x]])");
    T("(\"Compiled\" /. CompileDiagnostics[{{x, _Real}}, DirichletL[5, 2, x]])");
    T("(\"Compiled\" /. CompileDiagnostics[{{v, _Real, 1}}, DirichletL[5, 3, v]])");
    T("(\"Compiled\" /. CompileDiagnostics[{{v, _Real, 1}}, DirichletL[5, 2, v]])");
    T("(\"ResultType\" /. CompileDiagnostics[{{x, _Real}}, DirichletL[5, 3, x]]) === \"Real\"");
    T("(\"ResultType\" /. CompileDiagnostics[{{x, _Real}}, DirichletL[5, 2, x]]) === \"Complex\"");
    T("With[{xs = Range[-4.9, 10., 0.37], f = Compile[{{x, _Real}}, DirichletL[5, 2, x]]}, "
      "Table[f[x], {x, xs}] === Table[DirichletL[5, 2, x], {x, xs}]]");
    T("With[{xs = Range[-4.9, 10., 0.37], f = Compile[{{v, _Real, 1}}, DirichletL[5, 3, v]]}, "
      "f[xs] === DirichletL[5, 3, xs]]");
    T("With[{zs = Table[x + 3.5 I, {x, -4.9, 10., 0.37}], f = Compile[{{v, _Complex, 1}}, DirichletL[7, 2, v]]}, "
      "f[zs] === Table[DirichletL[7, 2, z], {z, zs}]]");
}

/* ---- messages / attributes ------------------------------------------- */

void test_dl_messages() {
    assert_eval_eq("DirichletL[]", "DirichletL[]", 0);
    assert_eval_eq("Check[DirichletL[], failed]", "failed", 0);
    assert_eval_eq("Check[DirichletL[7, 9, 2], failed]", "failed", 0);
    assert_eval_eq("Check[DirichletL[7, 0, 2], failed]", "failed", 0);
    assert_eval_eq("Check[DirichletL[3., 2, 0.5], failed]", "failed", 0);
    assert_eval_eq("DirichletL[0, 1, 2]", "DirichletL[0, 1, 2]", 0);
    T("Head[DirichletL[k, j, 2.]] === DirichletL");
}

void test_dl_attributes() {
    /* NHoldAll, not NumericFunction: NumericQ is False, N still evaluates. */
    assert_eval_eq("Attributes[DirichletL]", "{Listable, NHoldAll, Protected}", 0);
    T("!NumericQ[DirichletL[3, 2, 1/2]]");
    T("NumberQ[N[DirichletL[3, 2, 1/2]]]");
}

int main() {
    symtab_init();
    core_init();

    TEST(test_dl_zeta);
    TEST(test_dl_exact_nonpositive);
    TEST(test_dl_exact_one);
    TEST(test_dl_machine_real);
    TEST(test_dl_machine_complex);
    TEST(test_dl_pole);
    TEST(test_dl_mpfr);
    TEST(test_dl_interval);
    TEST(test_dl_listable_and_packed);
    TEST(test_dl_compile);
    TEST(test_dl_messages);
    TEST(test_dl_attributes);

    printf("All DirichletL tests passed.\n");
    return 0;
}
