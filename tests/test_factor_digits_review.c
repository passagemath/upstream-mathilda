/* test_factor_digits_review.c -- the number-theory review, batches 4-5
 * (v0.473): FactorInteger[n, k] (trial-then-cofactor grouping),
 * Divisors and SquareFreeQ options and messages (SquareFreeQ over GF(p)), and
 * the digit functions -- IntegerDigits, FromDigits (mantissa form, Roman
 * numerals, any base, Horner), DigitCount (reals, zero, digit lists),
 * IntegerExponent -- with their packed and Compile surfaces.
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
        printf("factor/digits review test failed: %s\n  expected: %s\n  got:      %s\n",
               input, expected, res_str);
        ASSERT(0);
    }
    free(res_str);
    expr_free(e);
    expr_free(res);
}

/* ---- FactorInteger[n, k] ---- */
static void test_factorinteger_k(void) {
    check("FactorInteger[360, 2]", "{{5, 1}, {72, 1}}");
    check("FactorInteger[30030, 3]", "{{11, 1}, {13, 1}, {210, 1}}");
    check("FactorInteger[100, 2]", "{{2, 2}, {5, 2}}");
    check("FactorInteger[-36, 2]", "{{-9, 1}, {2, 2}}");
    check("FactorInteger[-30030, 2]", "{{-2310, 1}, {13, 1}}");
    check("FactorInteger[-7/360, 2]", "{{-1/72, 1}, {5, -1}, {7, 1}}");
    check("FactorInteger[30030/17, 2]", "{{13, 1}, {17, -1}, {2310, 1}}");
    check("FactorInteger[10^30 + 1, 2]", "{{165573604901641, 1}, {6039610000603961, 1}}");
    check("FactorInteger[3*5*7*1000003*1000033, 3]", "{{7, 1}, {15, 1}, {1000036000099, 1}}");
    check("FactorInteger[-(3*5*7*1000003*1000033), 2]", "{{-105, 1}, {1000036000099, 1}}");
    check("FactorInteger[-2*(2^61 - 1), 2]", "{{-2305843009213693951, 1}, {2, 1}}");
    check("FactorInteger[1024, 1]", "{{1024, 1}}");
    check("FactorInteger[360, All]", "{{2, 3}, {3, 2}, {5, 1}}");
    check("FactorInteger[360, 0]", "FactorInteger[360, 0]");
    check("Check[FactorInteger[], m, FactorInteger::argt]", "m");
    check("Check[FactorInteger[360, 2, 3], m, FactorInteger::nonopt]", "m");
    check("FactorInteger[5, 1, GaussianIntegers -> True]", "{{5, 1}}");
}

/* ---- Divisors / SquareFreeQ ---- */
static void test_divisors_squarefree(void) {
    check("Check[Divisors[], m, Divisors::argx]", "m");
    check("Check[Divisors[1, 2], m, Divisors::nonopt]", "m");
    check("Check[Divisors[12, Foo -> 1], m, Divisors::optx]", "m");
    check("Check[Divisors[12, GaussianIntegers -> Automatic], m, Divisors::opttf]", "m");
    check("Divisors[{6, 8}]", "{{1, 2, 3, 6}, {1, 2, 4, 8}}");
    check("SquareFreeQ[(x - 1)^2, Modulus -> 2]", "False");
    check("SquareFreeQ[x^2 + 1, Modulus -> 3]", "True");
    check("SquareFreeQ[x^5 + 1, Modulus -> 5]", "False");
    check("SquareFreeQ[x^2 + 1, Modulus -> 0]", "True");
    check("SquareFreeQ[(x + y)^2, Modulus -> 5]", "SquareFreeQ[(x + y)^2, Modulus -> 5]");
    check("Check[SquareFreeQ[12, Modulus -> 5], m, SquareFreeQ::par]", "m");
    check("Check[SquareFreeQ[12, Foo -> 1], m, SquareFreeQ::par]", "m");
    check("SquareFreeQ[{12, 30}]", "SquareFreeQ[{12, 30}]");
    check("SquareFreeQ[9, GaussianIntegers -> True]", "False");
}

/* ---- IntegerDigits / IntegerExponent ---- */
static void test_integerdigits_exponent(void) {
    check("IntegerDigits[3^100, 60]",
          "{30, 12, 46, 34, 38, 45, 5, 20, 37, 10, 27, 49, 26, 25, 6, 58, 25, 8, 56, 57, "
          "19, 44, 23, 51, 7, 13, 21}");
    check("IntegerDigits[3^100, 7, 3]", "{3, 3, 4}");
    /* a bignum's digits come back packed */
    check("With[{d = IntegerDigits[3^2000]}, {PackedArrayQ[d], Length[d], Total[d], Head[d[[3]]]}]",
          "{True, 955, 4338, Integer}");
    check("IntegerDigits[3^100, 10, 0]", "{}");
    check("Take[IntegerDigits[3^100, 2, 300], 3]", "{0, 0, 0}");
    check("FromDigits[IntegerDigits[3^20000, 62], 62] == 3^20000", "True");
    check("Check[IntegerDigits[1234, 1/2], m, IntegerDigits::ibase]", "m");
    check("Check[IntegerDigits[1234, 10, -1], m, IntegerDigits::intnm]", "m");
    check("With[{r = Range[-3000, 3000]}, IntegerDigits[r] === Map[IntegerDigits, r]]", "True");
    check("Check[IntegerExponent[12, 2.], m, IntegerExponent::ibase]", "m");
    check("IntegerExponent[0]", "Infinity");
    check("With[{r = DeleteCases[Range[-3000, 3000], 0]}, "
          "IntegerExponent[r, 2] === Map[IntegerExponent[#, 2] &, r]]", "True");
    check("IntegerExponent[{0, 8, 12}, 2]", "{Infinity, 3, 2}");
    check("\"Compiled\" /. CompileDiagnostics[{{v, _Integer, 1}}, IntegerExponent[v, 2]]", "True");
}

/* ---- FromDigits ---- */
static void test_fromdigits(void) {
    check("FromDigits[{{1, 2}, -1}]", "3/250");
    check("FromDigits[{{1, 2}, 2}]", "12");
    check("FromDigits[{{1, 2}, 0}, 2]", "1");
    check("FromDigits[{{1, {3}}, 1}]", "4/3");
    check("FromDigits[{{1, 2, {3, 4}}, 2}]", "1222/99");
    check("FromDigits[{{{3}}, 0}]", "1/3");
    check("FromDigits[{{}, 2}]", "0");
    check("FromDigits[{{1, 2}, {3, 4}}]", "{120, 1200}");
    check("FromDigits[{{1, {3}}, {1, 2}}]", "{4/3, 40/3}");
    /* packed digit buffers: read in place, same answer as the List path */
    check("With[{r = RandomInteger[9, 2000]}, FromDigits[r] == Fold[10 #1 + #2 &, 0, r]]", "True");
    check("With[{r = RandomInteger[15, 2000]}, FromDigits[r, 16] == Fold[16 #1 + #2 &, 0, r]]", "True");
    check("FromDigits[Range[20], x] === FromDigits[Table[i, {i, 20}], x]", "True");
    check("FromDigits[{1., 2., 3.}] === 123.", "True");
    check("FromDigits[NDArray[{1, 2, 3}, DataType -> \"int64\"], 2]", "11");
    check("FromDigits[{1, 2, 3}, 1]", "6");
    check("FromDigits[{1, 2, 3}, -10]", "83");
    check("FromDigits[{1, 2, 3}, 1/2]", "17/4");
    /* the reference prints c + 10 (b + 10 a); compare the tree, not the term order */
    check("FromDigits[{a, b, c}] === c + 10 (b + 10 a)", "True");
    check("FromDigits[\"MCMXC\", \"Roman\"]", "1990");
    check("FromDigits[\"iv\", \"Roman\"]", "4");
    check("FromDigits[\"IIII\", \"Roman\"]", "4");
    check("Quiet[FromDigits[\"ABC\", \"Roman\"]]", "$Failed");
    check("FromDigits[\"12a\"]", "130");
    check("Check[FromDigits[], m, FromDigits::argt]", "m");
    check("FromDigits[StringJoin[Table[\"7\", 20000]]] == (10^20000 - 1) 7/9", "True");
}

/* ---- DigitCount ---- */
static void test_digitcount(void) {
    check("DigitCount[0]", "{0, 0, 0, 0, 0, 0, 0, 0, 0, 1}");
    check("DigitCount[0, 2]", "{0, 0}");
    check("DigitCount[0, 3]", "{0, 0, 1}");
    check("DigitCount[1234, 10, {1, 2}]", "{1, 1}");
    check("DigitCount[2.5, 2]", "{2, 51}");
    check("DigitCount[1.5]", "{1, 0, 0, 0, 1, 0, 0, 0, 0, 14}");
    check("Check[DigitCount[12, 10, 10], m, DigitCount::digp]", "m");
    check("Check[DigitCount[12, 1], m, DigitCount::base]", "m");
    check("Check[DigitCount[1/2], m, DigitCount::int]", "m");
    /* threads over the first argument only (not Listable) */
    check("DigitCount[{1, 0, 1}, 2]", "{{1, 0}, {0, 0}, {1, 0}}");
    check("DigitCount[{12, 30}, 10, {1, 3}]", "{{1, 0}, {0, 1}}");
    check("DigitCount[{12.5, 3}]", "{{1, 1, 0, 0, 1, 0, 0, 0, 0, 13}, {0, 0, 1, 0, 0, 0, 0, 0, 0, 0}}");
    check("DigitCount[{}]", "{}");
    check("Check[DigitCount[{12, x}, 10, 1], m, DigitCount::int]", "m");
    check("Check[DigitCount[12, {2, 3}], m, DigitCount::base]", "m");
    check("DigitSum[1234]", "10");
}

int main() {
    symtab_init();
    core_init();

    TEST(test_factorinteger_k);
    TEST(test_divisors_squarefree);
    TEST(test_integerdigits_exponent);
    TEST(test_fromdigits);
    TEST(test_digitcount);

    printf("All factor/digits review tests passed!\n");
    return 0;
}
