/* test_divisibility_review.c -- the number-theory review, batch 1 (v0.471):
 * GCD / LCM / ExtendedGCD over the Gaussian integers and rationals, the
 * Mathematica messages these heads gained, and agreement across every machine
 * surface (scalar, unpacked List, packed List, NDArray, Compile) for EulerPhi,
 * MoebiusMu, JacobiSymbol, CoprimeQ and Divisible.
 *
 * Every literal expectation below was produced by Mathematica 14.0
 * (wolframscript), not recalled. */

#include "expr.h"
#include "eval.h"
#include "core.h"
#include "symtab.h"
#include "test_utils.h"
#include "parse.h"
#include "print.h"
#include <string.h>
#include <stdlib.h>

static void check(const char* input, const char* expected) {
    Expr* e = parse_expression(input);
    Expr* res = evaluate(e);
    char* res_str = expr_to_string(res);
    if (strcmp(res_str, expected) != 0) {
        printf("divisibility review test failed: %s\n  expected: %s\n  got:      %s\n",
               input, expected, res_str);
        ASSERT(0);
    }
    free(res_str);
    expr_free(e);
    expr_free(res);
}

/* ---- GCD / LCM in Z[i] and Q(i), first-quadrant normalised ---- */
static void test_gaussian_gcd_lcm(void) {
    check("GCD[3 + I, 2]", "1 + I");
    check("GCD[3 + I, 1 + 3 I]", "1 + I");
    check("GCD[6, 3 + 3 I]", "3 + 3*I");
    check("GCD[1/2 + I/3, 2]", "1/6");
    check("GCD[1/2 + I/3, 1/4]", "1/12");
    check("GCD[I, 0]", "1");
    check("GCD[0, 2 + 2 I]", "2 + 2*I");
    check("GCD[-1 - I]", "1 + I");
    check("GCD[2 I]", "2");
    check("GCD[2, 3 + I, x]", "GCD[2, 3 + I, x]");
    check("LCM[3 + I, 2]", "2 + 4*I");
    check("LCM[2 + I, 2 - I]", "5");
    check("LCM[1/2 + I/3, 2]", "6 + 4*I");
    check("LCM[1/2, 1 + I]", "1 + I");
    check("LCM[0, 1 + I]", "0");
    check("LCM[-1 - I]", "1 + I");
}

/* ---- ExtendedGCD in Z[i]: Mathematica's cofactors ---- */
static void test_gaussian_extendedgcd(void) {
    check("ExtendedGCD[3 + I, 2]", "{1 + I, {-I, 2*I}}");
    check("ExtendedGCD[2, 3 + I]", "{1 + I, {2*I, -I}}");
    check("ExtendedGCD[5 + I, 3 - 2 I]", "{2 + 3*I, {0, I}}");
    check("ExtendedGCD[1 + I]", "{1 + I, {1}}");
    check("ExtendedGCD[3 + I, 2, 5 + I]", "{1 + I, {-I, 2*I, 0}}");
    check("ExtendedGCD[0, 1 + I]", "{1 + I, {0, 1}}");
    check("Check[ExtendedGCD[1/2 + I, 2], m, ExtendedGCD::egcd]", "m");
    /* The Bezout identity holds on a random batch, whatever the cofactors. */
    check("With[{p = Table[{RandomInteger[{-50, 50}] + I RandomInteger[{-50, 50}], "
          "RandomInteger[{-50, 50}] + I RandomInteger[{-50, 50}]}, {300}]}, "
          "And @@ Map[With[{r = ExtendedGCD @@ #}, "
          "Expand[r[[2]] . #] == r[[1]] && r[[1]] == GCD @@ #] &, "
          "Select[p, # =!= {0, 0} &]]]", "True");
    /* int64 fast path == the GMP fold (a trailing 0 forces the fold). */
    check("With[{p = RandomInteger[{-10^15, 10^15}, {500, 2}]}, "
          "And @@ Map[ExtendedGCD[#[[1]], #[[2]]] === "
          "{ExtendedGCD[#[[1]], #[[2]], 0][[1]], Most[ExtendedGCD[#[[1]], #[[2]], 0][[2]]]} &, p]]",
          "True");
}

/* ---- Messages ---- */
static void test_messages(void) {
    check("Check[GCD[2.5, 5], m, GCD::exact]", "m");
    check("Check[LCM[2.5, 5], m, LCM::exact]", "m");
    check("Check[LCM[], m, LCM::argm]", "m");
    check("Quiet[LCM[]]", "LCM[]");
    check("Check[EulerPhi[1/2], m, EulerPhi::int]", "m");
    check("Check[EulerPhi[3 + I], m, EulerPhi::int]", "m");
    check("Check[EulerPhi[], m, EulerPhi::argx]", "m");
    check("Check[EulerPhi[x], m]", "EulerPhi[x]");
    check("Check[JacobiSymbol[2], m, JacobiSymbol::argr]", "m");
    check("Check[ChineseRemainder[{1, 2}, {3, 5}, d], m, ChineseRemainder::int]", "m");
    check("Check[ChineseRemainder[{1, 2}], m, ChineseRemainder::argr]", "m");
}

/* ---- Values that changed ---- */
static void test_values(void) {
    check("MoebiusMu[0]", "0");
    check("EulerPhi[0]", "0");
    check("{JacobiSymbol[2, 6], JacobiSymbol[2, -6], JacobiSymbol[0, 0], JacobiSymbol[3, 0]}",
          "{0, 0, 0, 0}");
    check("{CoprimeQ[1], CoprimeQ[-1], CoprimeQ[3], CoprimeQ[I], CoprimeQ[0], CoprimeQ[x]}",
          "{True, True, False, True, False, False}");
    check("Attributes[GCD]", "{Flat, Listable, OneIdentity, Orderless, Protected}");
    check("Attributes[LCM]", "{Flat, Listable, OneIdentity, Orderless, Protected}");
    check("Options[CoprimeQ]", "{}");
}

/* ---- Every machine surface agrees with the scalar builtin ---- */
static void test_surfaces(void) {
    const char* d = "Join[Range[-3000, -1], Range[0, 3000]]";
    char buf[512];
    const char* forms[] = {
        "EulerPhi[%s] === Map[EulerPhi, %s]",
        "MoebiusMu[%s] === Map[MoebiusMu, %s]",
        "JacobiSymbol[%s, 1009] === Map[JacobiSymbol[#, 1009] &, %s]",
        "JacobiSymbol[77, %s] === Map[JacobiSymbol[77, #] &, %s]",
        "CoprimeQ[%s, 30] === Map[CoprimeQ[#, 30] &, %s]",
        "Divisible[%s, 7] === Map[Divisible[#, 7] &, %s]",
        "Divisible[%s, -1] === Map[Divisible[#, -1] &, %s]",
    };
    for (size_t i = 0; i < sizeof forms / sizeof forms[0]; i++) {
        snprintf(buf, sizeof buf, forms[i], d, d);
        check(buf, "True");
    }
    /* Large enough to take the smallest-prime-factor sieve. */
    check("With[{b = Range[4*10^6, 4*10^6 + 50000]}, "
          "{EulerPhi[b] === Map[EulerPhi, b], MoebiusMu[b] === Map[MoebiusMu, b]}]",
          "{True, True}");
    /* Visible int64 NDArray. */
    check("Normal[EulerPhi[NDArray[Range[-50, 50], DataType -> \"int64\"]]] === "
          "EulerPhi[Range[-50, 50]]", "True");
    check("CoprimeQ[NDArray[Range[10], DataType -> \"int64\"], 6]",
          "{True, False, False, False, True, False, True, False, False, False}");
    /* A zero divisor degrades to the element path, with its message. */
    check("Quiet[Divisible[{1, 2}, 0]]", "{Divisible[1, 0], Divisible[2, 0]}");
    /* Compile lowers all three scalar heads, and agrees. */
    check("\"Compiled\" /. CompileDiagnostics[{{x, _Integer}}, EulerPhi[x] + MoebiusMu[x]]",
          "True");
    check("\"Compiled\" /. CompileDiagnostics[{{x, _Integer}, {y, _Integer}}, JacobiSymbol[x, y]]",
          "True");
    check("With[{f = Compile[{{x, _Integer}}, {EulerPhi[x], MoebiusMu[x], JacobiSymbol[x, 15]}]}, "
          "Map[f, Range[-200, 200]] === "
          "Table[{EulerPhi[x], MoebiusMu[x], JacobiSymbol[x, 15]}, {x, -200, 200}]]", "True");
}

int main() {
    symtab_init();
    core_init();

    TEST(test_gaussian_gcd_lcm);
    TEST(test_gaussian_extendedgcd);
    TEST(test_messages);
    TEST(test_values);
    TEST(test_surfaces);

    printf("All divisibility review tests passed!\n");
    return 0;
}
