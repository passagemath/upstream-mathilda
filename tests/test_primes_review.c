/* test_primes_review.c -- the number-theory review, batch 3 (v0.472):
 * the tuned Deleglise-Rivat behind PrimePi / Prime, Prime's sieve-the-gap
 * finish, NextPrime's signed-prime semantics, the PrimePi / PrimeQ / PrimeNu /
 * PrimeOmega messages and options, and PrimeQ's packed and Compile surfaces.
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
        printf("primes review test failed: %s\n  expected: %s\n  got:      %s\n",
               input, expected, res_str);
        ASSERT(0);
    }
    free(res_str);
    expr_free(e);
    expr_free(res);
}

/* ---- PrimePi: the tuned Deleglise-Rivat, and its front end ---- */
static void test_primepi(void) {
    check("PrimePi[10^9]", "50847534");
    check("PrimePi[10^12]", "37607912018");
    check("PrimePi[10^13]", "346065536839");
    check("PrimePi[5*10^13]", "1638923764567");
    check("PrimePi[10^10 + 1/2]", "455052511");
    check("PrimePi[Pi]", "2");
    check("PrimePi[E^10]", "2466");
    check("PrimePi[Infinity]", "Infinity");
    check("PrimePi[-Infinity]", "0");
    check("PrimePi[-3]", "0");
    check("PrimePi[2 + I]", "PrimePi[2 + I]");
    check("Check[PrimePi[], m, PrimePi::argx]", "m");
    check("Check[PrimePi[1, 2], m, PrimePi::nonopt]", "m");
    check("Check[PrimePi[10, Foo -> 1], m, PrimePi::optx]", "m");
    check("PrimePi[10, Method -> \"Foo\"]", "PrimePi[10, Method -> \"Foo\"]");
    check("PrimePi[100, ProgressReporting -> False]", "25");
    check("Options[PrimePi]", "{Method -> Automatic, ProgressReporting :> $ProgressReporting}");
    /* Every method agrees on random x (DR is Automatic; LMO and Meissel are
     * independent implementations). */
    check("And @@ Table[With[{x = RandomInteger[{10^6, 10^10}]}, "
          "PrimePi[x] == PrimePi[x, Method -> \"LMO\"] == PrimePi[x, Method -> \"Meissel\"]], {40}]",
          "True");
}

/* ---- Prime: R-estimate + one count + sieve the gap ---- */
static void test_prime(void) {
    check("Prime[123456789]", "2543568463");
    check("Prime[10^9]", "22801763489");
    check("Prime[10^10]", "252097800623");
    check("Prime[10^11]", "2760727302517");
    check("And @@ Table[With[{n = RandomInteger[{10^6, 10^8}]}, PrimePi[Prime[n]] == n && "
          "PrimeQ[Prime[n]]], {20}]", "True");
    check("Check[Prime[2.5], m, Prime::intpp]", "m");
}

/* ---- NextPrime: the negative primes count ---- */
static void test_nextprime(void) {
    check("NextPrime[-10]", "-7");
    check("NextPrime[3, -2]", "-2");
    check("NextPrime[2, -1]", "-2");
    check("NextPrime[-2]", "2");
    check("NextPrime[10, 0]", "NextPrime[10, 0]");
    check("NextPrime[7, 0]", "7");
    check("NextPrime[-7, 0]", "-7");
    check("NextPrime[-5/2, -2]", "-5");
    check("NextPrime[2.5, -1]", "2");
    check("NextPrime[0.5, -1]", "-2");
    check("NextPrime[{4, 10}, 2]", "{7, 13}");
    check("NextPrime[4, {1, 2}]", "{5, 7}");
    check("NextPrime[2^64, -1]", "18446744073709551557");
    check("NextPrime[10^20, 3]", "100000000000000000151");
    check("NextPrime[10, 1/2]", "NextPrime[10, 1/2]");
    check("Check[NextPrime[], m, NextPrime::argt]", "m");
    check("Attributes[NextPrime]", "{Protected}");
}

/* ---- PrimeQ / PrimeNu / PrimeOmega messages and options ---- */
static void test_messages_options(void) {
    check("Check[PrimeQ[], m, PrimeQ::argx]", "m");
    check("Check[PrimeQ[7, Foo -> 1], m, PrimeQ::optx]", "m");
    check("Check[PrimeQ[7, GaussianIntegers -> Automatic], m, PrimeQ::opttf]", "m");
    check("Check[PrimeQ[1, 2], m, PrimeQ::nonopt]", "m");
    check("PrimeQ[5, GaussianIntegers -> True]", "False");
    check("PrimeQ[3, GaussianIntegers -> True]", "True");
    check("Check[PrimeNu[2.], m, PrimeNu::exact]", "m");
    check("Check[PrimeOmega[2.], m, PrimeOmega::exact]", "m");
    check("Options[PrimeNu]", "{GaussianIntegers -> False}");
    check("Options[PrimeOmega]", "{GaussianIntegers -> False}");
    check("PrimeNu[5, GaussianIntegers -> True]", "2");
    check("PrimeOmega[-8]", "3");
    check("PrimeOmega[1/12]", "PrimeOmega[1/12]");
}

/* ---- PrimeQ machine surfaces ---- */
static void test_primeq_surfaces(void) {
    check("With[{d = Join[Range[-2000, 2000], RandomInteger[{-10^18, 10^18}, 2000], "
          "{2^63 - 1, -2^63 + 1}]}, PrimeQ[d] === Map[PrimeQ, d]]", "True");
    check("With[{r = Range[10^5]}, PrimeQ[r] === Map[PrimeQ, r]]", "True");
    check("Normal[PrimeQ[NDArray[Range[10], DataType -> \"int64\"]]]",
          "{False, True, True, False, True, False, True, False, False, False}");
    check("\"Compiled\" /. CompileDiagnostics[{{x, _Integer}}, PrimeQ[x]]", "True");
    check("With[{f = Compile[{{x, _Integer}}, PrimeQ[x]]}, "
          "Map[f, Range[-300, 300]] === Map[PrimeQ, Range[-300, 300]]]", "True");
}

int main() {
    symtab_init();
    core_init();

    TEST(test_primepi);
    TEST(test_prime);
    TEST(test_nextprime);
    TEST(test_messages_options);
    TEST(test_primeq_surfaces);

    printf("All primes review tests passed!\n");
    return 0;
}
