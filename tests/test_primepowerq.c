/* Unit tests for PrimePowerQ (src/facint.c), its packed / NDArray buffer path
 * (ndint_prime_power_predicate, src/ndinteger.c), its Compile[] lowering
 * (OP_PPQ_I) and the exported kernels mth_prime_power_base_u64 / _mpz.
 *
 * Every expected value in the generated groups below was produced by
 * Mathematica 14 (wolframscript), not recalled; the two printed-real cases
 * differ only in Mathilda's "2.0" spelling of WL's "2.". */
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <gmp.h>
#include "expr.h"
#include "eval.h"
#include "core.h"
#include "symtab.h"
#include "test_utils.h"
#include "parse.h"
#include "print.h"
#include "facint.h"

/* Shared driver: parse, evaluate, compare the printed form. */
static void check(const char* input, const char* expected) {
    Expr* e = parse_expression(input);
    Expr* res = evaluate(e);
    char* res_str = expr_to_string(res);
    if (strcmp(res_str, expected) != 0) {
        printf("PrimePowerQ test failed: %s\n  expected: %s\n  got:      %s\n",
               input, expected, res_str);
        ASSERT(0);
    }
    free(res_str);
    expr_free(e);
    expr_free(res);
}

/* ---- The examples from the Mathematica documentation (the request). */
void test_ppq_documentation() {
    check("PrimePowerQ[8]",
          "True");
    check("PrimePowerQ[6]",
          "False");
    check("PrimePowerQ[9]",
          "True");
    check("PrimePowerQ[2+I]",
          "True");
    check("PrimePowerQ[5, GaussianIntegers -> True]",
          "False");
    check("PrimePowerQ[1/25]",
          "True");
    check("PrimePowerQ[10^300+3]",
          "False");
    check("PrimePowerQ[{1, 2, 3, 4, 5, 6, 7}]",
          "{False, True, True, True, True, False, True}");
    check("Select[Range[100], PrimePowerQ[#] && !PrimeQ[#] &]",
          "{4, 8, 9, 16, 25, 27, 32, 49, 64, 81}");
}

/* ---- Small integers: units, primes, prime powers, and composites around the
 *      trial-division boundary (67^2 = 4489, 61*67 = 4087). */
void test_ppq_small() {
    check("PrimePowerQ[0]",
          "False");
    check("PrimePowerQ[1]",
          "False");
    check("PrimePowerQ[-1]",
          "False");
    check("PrimePowerQ[2]",
          "True");
    check("PrimePowerQ[3]",
          "True");
    check("PrimePowerQ[4]",
          "True");
    check("PrimePowerQ[12]",
          "False");
    check("PrimePowerQ[36]",
          "False");
    check("PrimePowerQ[49]",
          "True");
    check("PrimePowerQ[64]",
          "True");
    check("PrimePowerQ[100]",
          "False");
    check("PrimePowerQ[121]",
          "True");
    check("PrimePowerQ[4489]",
          "True");
    check("PrimePowerQ[4087]",
          "False");
    check("PrimePowerQ[4096]",
          "True");
    check("PrimePowerQ[1024]",
          "True");
    check("Select[Range[200], PrimePowerQ]",
          "{2, 3, 4, 5, 7, 8, 9, 11, 13, 16, 17, 19, 23, 25, 27, 29, 31, 32, 37, "
          "41, 43, 47, 49, 53, 59, 61, 64, 67, 71, 73, 79, 81, 83, 89, 97, 101, 1"
          "03, 107, 109, 113, 121, 125, 127, 128, 131, 137, 139, 149, 151, 157, 1"
          "63, 167, 169, 173, 179, 181, 191, 193, 197, 199}");
    check("Count[Range[10^5], _?PrimePowerQ]",
          "9700");
}

/* ---- Sign is ignored: -8 = -(2^3) is a prime power; INT64_MIN = -2^63 is
 *      answered through the unsigned magnitude. */
void test_ppq_negative() {
    check("PrimePowerQ[-2]",
          "True");
    check("PrimePowerQ[-8]",
          "True");
    check("PrimePowerQ[-6]",
          "False");
    check("PrimePowerQ[-27]",
          "True");
    check("PrimePowerQ[-1024]",
          "True");
    check("PrimePowerQ[-36]",
          "False");
    check("PrimePowerQ[-9223372036854775808]",
          "True");
    check("PrimePowerQ[-9223372036854775807]",
          "False");
}

/* ---- Around 2^63 / 2^64, and the u64 arm's exact k-th roots
 *      (k in {2,3,5,7}) for bases past the trial divisors. */
void test_ppq_machine_boundary() {
    check("PrimePowerQ[9223372036854775807]",
          "False");
    check("PrimePowerQ[2^62]",
          "True");
    check("PrimePowerQ[2^63]",
          "True");
    check("PrimePowerQ[2^64]",
          "True");
    check("PrimePowerQ[2^64 - 1]",
          "False");
    check("PrimePowerQ[2^64 + 1]",
          "False");
    check("PrimePowerQ[3^40]",
          "True");
    check("PrimePowerQ[3^41]",
          "True");
    check("PrimePowerQ[18446744073709551557]",
          "True");
    check("PrimePowerQ[4294967291^2]",
          "True");
    check("PrimePowerQ[4294967291*4294967279]",
          "False");
    check("PrimePowerQ[2642245^3]",
          "False");
    check("PrimePowerQ[2642243^3]",
          "False");
    check("PrimePowerQ[67^10]",
          "True");
    check("PrimePowerQ[71^10]",
          "True");
    check("PrimePowerQ[67^7]",
          "True");
    check("PrimePowerQ[67^9 * 71]",
          "False");
    check("PrimePowerQ[1031^6]",
          "True");
    check("PrimePowerQ[(67*71)^4]",
          "False");
}

/* ---- Bignums: trial-division hits (mpz_remove), primitive-root reduction
 *      through nested exponents, large semiprimes, strong pseudoprimes to small
 *      bases (3825123056546413051, 341550071728321), a 10^5-digit power. */
void test_ppq_bignum() {
    check("PrimePowerQ[2^1000]",
          "True");
    check("PrimePowerQ[2^1000 + 1]",
          "False");
    check("PrimePowerQ[3^1000]",
          "True");
    check("PrimePowerQ[6^1000]",
          "False");
    check("PrimePowerQ[(10^100 + 267)^7]",
          "True");
    check("PrimePowerQ[(10^100 + 267)^2 * 3]",
          "False");
    check("PrimePowerQ[(10^100 + 267)^2 * (10^100 + 949)]",
          "False");
    check("PrimePowerQ[(10^100 + 267)^210]",
          "True");
    check("PrimePowerQ[(10^100+267)^6 * (10^100+949)^6]",
          "False");
    check("PrimePowerQ[NextPrime[10^50]^(2*3*5*7)]",
          "True");
    check("PrimePowerQ[3^100000]",
          "True");
    check("PrimePowerQ[1031^5000]",
          "True");
    check("PrimePowerQ[1031^5000*1033]",
          "False");
    check("PrimePowerQ[(2^127 - 1)^3]",
          "True");
    check("PrimePowerQ[2^127 - 1]",
          "True");
    check("PrimePowerQ[2^128 + 1]",
          "False");
    check("PrimePowerQ[(2^61-1)*(2^89-1)]",
          "False");
    check("PrimePowerQ[10^300 + 7]",
          "False");
    check("PrimePowerQ[(10^20 + 39)^5]",
          "True");
    check("PrimePowerQ[((10^20+39)^2)^3]",
          "True");
    check("PrimePowerQ[-(10^20 + 39)^5]",
          "True");
    check("PrimePowerQ[3825123056546413051]",
          "False");
    check("PrimePowerQ[341550071728321]",
          "False");
    check("PrimePowerQ[561^2]",
          "False");
}

/* ---- A rational is a prime power iff its reciprocal is one. */
void test_ppq_rational() {
    check("PrimePowerQ[1/25]",
          "True");
    check("PrimePowerQ[-1/8]",
          "True");
    check("PrimePowerQ[2/3]",
          "False");
    check("PrimePowerQ[1/6]",
          "False");
    check("PrimePowerQ[8/9]",
          "False");
    check("PrimePowerQ[1/2]",
          "True");
    check("PrimePowerQ[-1/2]",
          "True");
    check("PrimePowerQ[3/2]",
          "False");
    check("PrimePowerQ[1/2^100]",
          "True");
    check("PrimePowerQ[1/(10^100+267)^3]",
          "True");
    check("PrimePowerQ[1/(6^50)]",
          "False");
}

/* ---- Gaussian integers: a unit times a power of ONE Gaussian prime (1+I;
 *      an inert q = 3 mod 4; or one split factor of p = 1 mod 4). */
void test_ppq_gaussian() {
    check("PrimePowerQ[2+I]",
          "True");
    check("PrimePowerQ[(2+I)^3]",
          "True");
    check("PrimePowerQ[5 (2+I)]",
          "False");
    check("PrimePowerQ[1+I]",
          "True");
    check("PrimePowerQ[(1+I)^5]",
          "True");
    check("PrimePowerQ[2 I]",
          "True");
    check("PrimePowerQ[3 I]",
          "True");
    check("PrimePowerQ[9 I]",
          "True");
    check("PrimePowerQ[5 I]",
          "False");
    check("PrimePowerQ[I]",
          "False");
    check("PrimePowerQ[-I]",
          "False");
    check("PrimePowerQ[2+2 I]",
          "True");
    check("PrimePowerQ[3+4 I]",
          "True");
    check("PrimePowerQ[-7-7 I]",
          "False");
    check("PrimePowerQ[7+7 I]",
          "False");
    check("PrimePowerQ[(2+I)^2 (2-I)]",
          "False");
    check("PrimePowerQ[(2+I)^2 I]",
          "True");
    check("PrimePowerQ[1+2 I]",
          "True");
    check("PrimePowerQ[3+2 I]",
          "True");
    check("PrimePowerQ[(3+2 I)^4]",
          "True");
    check("PrimePowerQ[(3+2 I)^4 (1+I)]",
          "False");
    check("PrimePowerQ[(1+I)^100]",
          "True");
    check("PrimePowerQ[(1+I)^101]",
          "True");
    check("PrimePowerQ[(10^20 + 39) I]",
          "True");
    check("PrimePowerQ[Expand[(1 + 10^10 I)^3]]",
          "False");
    check("PrimePowerQ[(10^50 + 151 + 2 I)^5]",
          "False");
    check("PrimePowerQ[0 + 3 I]",
          "True");
    check("PrimePowerQ[-3 - 0 I]",
          "True");
}

/* ---- Gaussian rationals: z is one iff 1/z is a Gaussian prime power. */
void test_ppq_gaussian_rational() {
    check("PrimePowerQ[1/(2+I)]",
          "True");
    check("PrimePowerQ[(1+I)/2]",
          "True");
    check("PrimePowerQ[I/2]",
          "True");
    check("PrimePowerQ[(1+I)^-3]",
          "True");
    check("PrimePowerQ[3/(2+I)]",
          "False");
    check("PrimePowerQ[(2+I)/(2-I)]",
          "False");
    check("PrimePowerQ[1/2 + I/3]",
          "False");
    check("PrimePowerQ[1/2 + I/2]",
          "True");
    check("PrimePowerQ[1/(3+4 I)]",
          "True");
    check("PrimePowerQ[1/(5 I)]",
          "False");
}

/* ---- GaussianIntegers -> True | False | Automatic; a well-formed :> value
 *      answers False, as Mathematica 14 does. */
void test_ppq_option() {
    check("PrimePowerQ[5, GaussianIntegers -> True]",
          "False");
    check("PrimePowerQ[2, GaussianIntegers -> True]",
          "True");
    check("PrimePowerQ[3, GaussianIntegers -> True]",
          "True");
    check("PrimePowerQ[9, GaussianIntegers -> True]",
          "True");
    check("PrimePowerQ[8, GaussianIntegers -> True]",
          "True");
    check("PrimePowerQ[-9, GaussianIntegers -> True]",
          "True");
    check("PrimePowerQ[25, GaussianIntegers -> True]",
          "False");
    check("PrimePowerQ[13, GaussianIntegers -> True]",
          "False");
    check("PrimePowerQ[0, GaussianIntegers -> True]",
          "False");
    check("PrimePowerQ[1, GaussianIntegers -> True]",
          "False");
    check("PrimePowerQ[1/9, GaussianIntegers -> True]",
          "True");
    check("PrimePowerQ[1/25, GaussianIntegers -> True]",
          "False");
    check("PrimePowerQ[1/2, GaussianIntegers -> True]",
          "True");
    check("PrimePowerQ[2+I, GaussianIntegers -> True]",
          "True");
    check("PrimePowerQ[5, GaussianIntegers -> False]",
          "True");
    check("PrimePowerQ[5, GaussianIntegers -> Automatic]",
          "True");
    check("PrimePowerQ[4, GaussianIntegers :> True]",
          "False");
    check("PrimePowerQ[5, GaussianIntegers :> True]",
          "False");
    check("PrimePowerQ[4, {GaussianIntegers -> True}]",
          "{True}");
    check("PrimePowerQ[Range[10], GaussianIntegers -> True]",
          "{False, True, True, True, False, False, True, True, True, False}");
    check("PrimePowerQ[2^1000, GaussianIntegers -> True]",
          "True");
    check("PrimePowerQ[3^1000, GaussianIntegers -> True]",
          "True");
    check("PrimePowerQ[(10^20+39)^3, GaussianIntegers -> True]",
          "True");
    check("PrimePowerQ[x, GaussianIntegers -> True]",
          "False");
    check("Options[PrimePowerQ]",
          "{GaussianIntegers -> Automatic}");
}

/* ---- Any other exact expression is False; Listable threads over lists. */
void test_ppq_nonnumeric() {
    check("PrimePowerQ[x]",
          "False");
    check("PrimePowerQ[\"s\"]",
          "False");
    check("PrimePowerQ[Pi]",
          "False");
    check("PrimePowerQ[Sqrt[2]]",
          "False");
    check("PrimePowerQ[f[2]]",
          "False");
    check("PrimePowerQ[2^x]",
          "False");
    check("PrimePowerQ[Infinity]",
          "False");
    check("PrimePowerQ[{}]",
          "{}");
    check("PrimePowerQ[{{4, 6}, {7, x}}]",
          "{{True, False}, {True, False}}");
}

/* ---- Messages route through the funnel: Quiet silences them, Check sees
 *      them, and the call stays unevaluated. */
void test_ppq_messages() {
    check("Quiet[PrimePowerQ[]]",
          "PrimePowerQ[]");
    check("Quiet[PrimePowerQ[3.2912]]",
          "PrimePowerQ[3.2912]");
    check("Quiet[PrimePowerQ[1, 2, 3]]",
          "PrimePowerQ[1, 2, 3]");
    check("Quiet[PrimePowerQ[4, 5]]",
          "PrimePowerQ[4, 5]");
    check("Quiet[PrimePowerQ[x, y]]",
          "PrimePowerQ[x, y]");
    check("Quiet[PrimePowerQ[4, foo -> 1]]",
          "PrimePowerQ[4, foo -> 1]");
    check("Quiet[PrimePowerQ[4, GaussianIntegers -> x]]",
          "PrimePowerQ[4, GaussianIntegers -> x]");
    check("Quiet[PrimePowerQ[2+I, GaussianIntegers -> False]]",
          "PrimePowerQ[2 + I, GaussianIntegers -> False]");
    check("Quiet[PrimePowerQ[2.5 + I]]",
          "PrimePowerQ[2.5 + 1.0*I]");
    check("Quiet[PrimePowerQ[{4, x, 2.}]]",
          "{True, False, PrimePowerQ[2.0]}");
    check("Check[PrimePowerQ[3.5], fail]",
          "fail");
    check("Check[PrimePowerQ[], fail]",
          "fail");
    check("Check[PrimePowerQ[4, 5], fail]",
          "fail");
    check("Check[PrimePowerQ[4, foo -> 1], fail]",
          "fail");
    check("Check[PrimePowerQ[4, GaussianIntegers -> x], fail]",
          "fail");
    check("Check[PrimePowerQ[2 + I, GaussianIntegers -> False], fail]",
          "fail");
    check("Check[PrimePowerQ[4], fail]",
          "True");
    check("Check[PrimePowerQ[x], fail]",
          "False");
    check("Check[PrimePowerQ[N[Pi, 30]], fail]",
          "fail");
}

/* ---- Cross-checks against FactorInteger and PrimeQ. */
void test_ppq_relations() {
    check("And @@ Table[PrimePowerQ[n] === (Abs[n] > 1 && Length[FactorInteger[Ab"
          "s[n]]] == 1), {n, -3000, 3000}]",
          "True");
    check("And @@ Flatten[Table[PrimePowerQ[p^k], {p, Prime[Range[25]]}, {k, 1, 1"
          "2}]]",
          "True");
    check("Select[Range[-30, 30], PrimePowerQ]",
          "{-29, -27, -25, -23, -19, -17, -16, -13, -11, -9, -8, -7, -5, -4, -3, "
          "-2, 2, 3, 4, 5, 7, 8, 9, 11, 13, 16, 17, 19, 23, 25, 27, 29}");
    check("PrimePowerQ[Range[2^40, 2^40 + 60]]",
          "{True, False, False, False, False, False, False, False, False, False, "
          "False, False, False, False, False, True, False, False, False, False, F"
          "alse, False, False, False, False, False, False, True, False, False, Fa"
          "lse, False, False, False, False, False, False, False, False, False, Fa"
          "lse, False, False, False, False, False, False, False, False, False, Fa"
          "lse, False, False, False, False, True, False, False, False, False, Fal"
          "se}");
}
/* ---- Attributes, registration, docstring ---- */
void test_ppq_registration() {
    check("Attributes[PrimePowerQ]", "{Listable, Protected}");
    check("Options[PrimePowerQ]", "{GaussianIntegers -> Automatic}");
    const char* doc = symtab_get_docstring("PrimePowerQ");
    ASSERT(doc && strstr(doc, "PrimePowerQ[expr]") && strstr(doc, "GaussianIntegers"));
    /* Protected: a definition is refused and the builtin still answers. */
    check("Quiet[PrimePowerQ[2] = 7]; PrimePowerQ[2]", "True");
}

/* ---- The packed and visible-NDArray surfaces agree with the List path ---- */
void test_ppq_packed() {
    /* Packed int64 in, packed bool out, element-for-element the boxed answer. */
    check("r = Range[-100, 5000]; {PackedArrayQ[r], PackedArrayQ[PrimePowerQ[r]]}",
          "{True, True}");
    check("PrimePowerQ[Range[-100, 5000]] === Map[PrimePowerQ, FromPackedArray[Range[-100, 5000]]]",
          "True");
    check("PrimePowerQ[Range[-100, 5000]] === Table[PrimePowerQ[k], {k, -100, 5000}]",
          "True");
    check("Count[PrimePowerQ[Range[10^6]], True]", "78734");
    /* INT64_MIN and INT64_MAX inside a buffer. */
    check("PrimePowerQ[ToPackedArray[{-9223372036854775807 - 1, 9223372036854775807, 1, 0, -1}]]",
          "{True, False, False, False, False}");
    /* Visible NDArray answers with the same List. */
    check("PrimePowerQ[NDArray[{1, 4, 6, 8, 9, 10, 2^62}, DataType -> \"int64\"]]",
          "{False, True, False, True, True, False, True}");
    /* A packed REAL buffer is not exact: per-element ::exact, unevaluated. */
    check("Quiet[PrimePowerQ[ToPackedArray[{2., 4.}]]]",
          "{PrimePowerQ[2.0], PrimePowerQ[4.0]}");
    check("Check[PrimePowerQ[ToPackedArray[{2., 4.}]], fail]", "fail");
    /* An option argument declines the buffer path and threads. */
    check("PrimePowerQ[Range[10], GaussianIntegers -> True] === "
          "Table[PrimePowerQ[k, GaussianIntegers -> True], {k, 10}]", "True");
    /* Rank 2 threads through the List path. */
    check("PrimePowerQ[ToPackedArray[{{4, 6}, {7, 10}}]]", "{{True, False}, {True, False}}");
}

/* ---- Compile[]: the scalar lowering is the interpreter's own kernel ---- */
void test_ppq_compile() {
    check("\"Compiled\" /. CompileDiagnostics[{{n, _Integer}}, PrimePowerQ[n]]", "True");
    check("cf = Compile[{{n, _Integer}}, PrimePowerQ[n]]; "
          "And @@ Table[cf[k] === PrimePowerQ[k], {k, -200, 5000}]", "True");
    check("cf = Compile[{{n, _Integer}}, PrimePowerQ[n]]; "
          "{cf[-9223372036854775807 - 1], cf[9223372036854775807], cf[2^62], cf[18446744073709551557 - 2^63]}",
          "{True, False, True, False}");
    check("g = Compile[{{n, _Integer}}, Module[{s = 0}, Do[If[PrimePowerQ[k], s++], {k, n}]; s]]; g[10^5]",
          "9700");
}

/* ---- The exported C kernels, including the prime base they report ---- */
void test_ppq_c_api() {
    uint64_t p = 0;
    ASSERT(!mth_prime_power_base_u64(0, &p));
    ASSERT(!mth_prime_power_base_u64(1, &p));
    ASSERT(mth_prime_power_base_u64(2, &p) && p == 2);
    ASSERT(mth_prime_power_base_u64(UINT64_C(1) << 63, &p) && p == 2);
    ASSERT(mth_prime_power_base_u64(UINT64_C(4294967291) * UINT64_C(4294967291), &p)
           && p == UINT64_C(4294967291));
    ASSERT(!mth_prime_power_base_u64(UINT64_C(4294967291) * UINT64_C(4294967279), &p));
    ASSERT(mth_prime_power_base_u64(UINT64_C(18446744073709551557), &p)
           && p == UINT64_C(18446744073709551557));
    ASSERT(!mth_prime_power_base_u64(UINT64_MAX, &p));        /* 3*5*17*257*641*65537*6700417 */
    /* 2642239 is the largest prime whose cube fits below 2^64 (k = 3 root). */
    ASSERT(mth_prime_power_base_u64(UINT64_C(2642239) * 2642239 * 2642239, &p) && p == 2642239);
    ASSERT(mth_prime_power_base_u64(UINT64_C(4489), &p) && p == 67);   /* 67^2 */
    ASSERT(mth_prime_power_base_u64(UINT64_C(1350125107), &p) && p == 67); /* 67^5 */
    ASSERT(mth_prime_power_base_u64(UINT64_C(4489), NULL));

    mpz_t n, b;
    mpz_init(n);
    mpz_init(b);
    /* (10^20 + 39)^6: primitive-root reduction through 2 and 3. */
    mpz_ui_pow_ui(n, 10, 20);
    mpz_add_ui(n, n, 39);
    mpz_pow_ui(n, n, 6);
    ASSERT(mth_prime_power_base_mpz(n, b));
    mpz_root(n, n, 6);
    ASSERT(mpz_cmp(n, b) == 0);
    /* 1031^500 * 1033: no trial divisor, not a perfect power. */
    mpz_ui_pow_ui(n, 1031, 500);
    mpz_mul_ui(n, n, 1033);
    ASSERT(!mth_prime_power_base_mpz(n, NULL));
    /* 1021^700: a trial-division hit, base reported. */
    mpz_ui_pow_ui(n, 1021, 700);
    ASSERT(mth_prime_power_base_mpz(n, b) && mpz_cmp_ui(b, 1021) == 0);
    /* Non-positive inputs are not prime powers (callers take |n|). */
    mpz_set_si(n, -8);
    ASSERT(!mth_prime_power_base_mpz(n, NULL));
    mpz_set_ui(n, 0);
    ASSERT(!mth_prime_power_base_mpz(n, NULL));
    mpz_clear(n);
    mpz_clear(b);
}

int main() {
    symtab_init();
    core_init();

    TEST(test_ppq_documentation);
    TEST(test_ppq_small);
    TEST(test_ppq_negative);
    TEST(test_ppq_machine_boundary);
    TEST(test_ppq_bignum);
    TEST(test_ppq_rational);
    TEST(test_ppq_gaussian);
    TEST(test_ppq_gaussian_rational);
    TEST(test_ppq_option);
    TEST(test_ppq_nonnumeric);
    TEST(test_ppq_messages);
    TEST(test_ppq_registration);
    TEST(test_ppq_packed);
    TEST(test_ppq_compile);
    TEST(test_ppq_c_api);

    printf("All PrimePowerQ tests passed!\n");
    return 0;
}
