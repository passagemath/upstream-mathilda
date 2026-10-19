/* test_congruence_review.c -- the number-theory review, batch 2 (v0.471):
 * Mod (zero modulus, infinities, Gaussian and exact-numeric operands, arity
 * messages), PowerMod (Gaussian base / modulus, ninv / root / divz / pmod,
 * negative modulus), MultiplicativeOrder (machine path, n <= 0, argtu) and
 * PrimitiveRootList (sign of n), and the machine surfaces they share.
 *
 * Every literal expectation below was verified against the reference
 * implementation, not recalled. */

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
        printf("congruence review test failed: %s\n  expected: %s\n  got:      %s\n",
               input, expected, res_str);
        ASSERT(0);
    }
    free(res_str);
    expr_free(e);
    expr_free(res);
}

/* ---- Mod ---- */
static void test_mod(void) {
    check("Quiet[Mod[7, 0]]", "Indeterminate");
    check("Check[Mod[7, 0], m, Mod::indet]", "m");
    check("Quiet[Mod[7, 0.]]", "Indeterminate");
    check("Quiet[Mod[Infinity, 3]]", "Indeterminate");
    check("Check[Mod[Infinity, 3], m, Infinity::indet]", "m");
    check("Mod[3, -Infinity]", "3");
    check("Mod[-3, Infinity]", "-3");
    check("Mod[Indeterminate, 3]", "Indeterminate");
    check("Mod[3 + 4 I, 2]", "-1");
    check("Mod[3 + 4 I, 1 + I]", "-1");
    check("Mod[5, 2 + I]", "0");
    check("Mod[1 + 2 I, 3 + 4 I, I]", "1 + 2*I");
    check("Mod[7, 3 + I, 1/2]", "I");
    check("Mod[Pi, 1]", "-3 + Pi");
    check("Mod[E, 2]", "-2 + E");
    check("Mod[Sqrt[2], 1]", "-1 + Sqrt[2]");
    check("Mod[x, 3]", "Mod[x, 3]");
    check("Check[Mod[7], m, Mod::argtu]", "m");
    check("Check[Mod[], m, Mod::argt]", "m");
    /* Gaussian Mod: machine and GMP routes agree with the nearest quotient. */
    check("With[{z = Table[RandomInteger[{-10^20, 10^20}] + I RandomInteger[{-10^20, 10^20}], {50}], "
          "w = Table[RandomInteger[{1, 10^6}] + I RandomInteger[{-10^6, 10^6}], {50}]}, "
          "And @@ MapThread[With[{r = Mod[#1, #2]}, IntegerQ[Re[(#1 - r)/#2]] && "
          "IntegerQ[Im[(#1 - r)/#2]] && Abs[Re[r/#2]] <= 1/2 && Abs[Im[r/#2]] <= 1/2] &, {z, w}]]",
          "True");
    /* Packed Mod with a constant divisor (reciprocal arm) == element path. */
    check("With[{e = Join[RandomInteger[{-10^15, 10^15}, 2000], {2^52 - 1, -2^52 + 1, 2^52, "
          "-2^52, 2^63 - 1, -2^63}]}, "
          "Table[Mod[e, k] === Map[Mod[#, k] &, e], {k, {2, 3, 7, 10^9 + 7, 2^40 + 1, 2^62}}]]",
          "{True, True, True, True, True, True}");
}

/* ---- PowerMod ---- */
static void test_powermod(void) {
    check("PowerMod[1 + I, 3, 5]", "-2 + 2*I");
    check("PowerMod[2, 10, 2 + I]", "-1");
    check("Check[PowerMod[3 + I, -1, 5], m, PowerMod::ninv]", "m");
    check("Check[PowerMod[2, -1, 4], m, PowerMod::ninv]", "m");
    check("Check[PowerMod[3, 1/2, 7], m, PowerMod::root]", "m");
    check("Check[PowerMod[2, 10, 0], m, PowerMod::divz]", "m");
    check("Check[PowerMod[2., 10, 7], m, PowerMod::pmod]", "m");
    check("Check[PowerMod[1/2, 3, 7], m, PowerMod::pmod]", "m");
    check("Check[PowerMod[2, 3], m, PowerMod::argrx]", "m");
    check("PowerMod[2, 10, -7]", "PowerMod[2, 10, -7]");
    check("PowerMod[x, 2, 7]", "PowerMod[x, 2, 7]");
    check("PowerMod[3, -1, 7]", "5");
    check("PowerMod[4, 1/2, 7]", "2");
    /* machine fast path == GMP (a bignum modulus forces GMP) */
    check("With[{p = RandomInteger[{-10^9, 10^9}, {300, 2}]}, And @@ Map[ "
          "PowerMod[#[[1]], Abs[#[[2]]], 4294967291] == "
          "Mod[PowerMod[#[[1]], Abs[#[[2]]], 4294967291 * 18446744073709551629], 4294967291] &, p]]",
          "True");
    /* Compile at both shapes */
    check("\"Compiled\" /. CompileDiagnostics[{{v, _Integer, 1}}, PowerMod[v, 5, 1009]]", "True");
    check("With[{f = Compile[{{v, _Integer, 1}}, PowerMod[v, 5, 1009]]}, "
          "f[Range[-50, 50]] === PowerMod[Range[-50, 50], 5, 1009]]", "True");
}

/* ---- MultiplicativeOrder, PrimitiveRootList ---- */
static void test_order_roots(void) {
    check("MultiplicativeOrder[3, 4294967291]", "2147483645");
    check("MultiplicativeOrder[3, 4294967311]", "4294967310");
    check("MultiplicativeOrder[2, 1000000007]", "500000003");
    check("MultiplicativeOrder[-5, 999999937]", "1338688");
    check("MultiplicativeOrder[10, 22]", "MultiplicativeOrder[10, 22]");
    check("MultiplicativeOrder[2, -7]", "MultiplicativeOrder[2, -7]");
    check("Check[MultiplicativeOrder[2], m, MultiplicativeOrder::argtu]", "m");
    check("PrimitiveRootList[-7]", "{3, 5}");
    check("PrimitiveRootList[-18]", "{5, 11}");
    check("PrimitiveRootList[0]", "{}");
}

int main() {
    symtab_init();
    core_init();

    TEST(test_mod);
    TEST(test_powermod);
    TEST(test_order_roots);

    printf("All congruence review tests passed!\n");
    return 0;
}
