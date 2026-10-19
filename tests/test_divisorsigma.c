#include "expr.h"
#include "eval.h"
#include "core.h"
#include "symtab.h"
#include "test_utils.h"
#include "parse.h"
#include "print.h"
#include <string.h>
#include <stdlib.h>

/* Shared driver: parse, evaluate, compare the printed form. */
static void check(const char* input, const char* expected) {
    Expr* e = parse_expression(input);
    Expr* res = evaluate(e);
    char* res_str = expr_to_string(res);
    if (strcmp(res_str, expected) != 0) {
        printf("DivisorSigma test failed: %s\n  expected: %s\n  got:      %s\n",
               input, expected, res_str);
        ASSERT(0);
    }
    free(res_str);
    expr_free(e);
    expr_free(res);
}

/* ---- Basic integer sigma_k ---- */
void test_divisorsigma_basic() {
    check("DivisorSigma[1, 20]", "42");
    check("DivisorSigma[2, 20]", "546");
    check("DivisorSigma[2, 6]", "50");
    check("DivisorSigma[1, 6]", "12");
    check("DivisorSigma[1, 1]", "1");        /* empty product */
    check("DivisorSigma[1, 28]", "56");      /* perfect number */
    check("DivisorSigma[3, 4]", "73");       /* 1 + 8 + 64 */
    check("DivisorSigma[1, 7]", "8");        /* prime */
}

/* ---- sigma_0 is the divisor count ---- */
void test_divisorsigma_count() {
    check("DivisorSigma[0, 12]", "6");
    check("DivisorSigma[0, 1]", "1");
    check("DivisorSigma[0, 7]", "2");
    check("DivisorSigma[0, 60]", "12");
    /* Consistency with Divisors. */
    check("DivisorSigma[0, 360] == Length[Divisors[360]]", "True");
    check("DivisorSigma[1, 360] == Total[Divisors[360]]", "True");
    check("DivisorSigma[2, 84] == Total[Divisors[84]^2]", "True");
}

/* ---- Negative powers give rationals ---- */
void test_divisorsigma_negative() {
    check("DivisorSigma[-2, 10]", "13/10");
    check("DivisorSigma[-1, 6]", "2");       /* 1 + 1/2 + 1/3 + 1/6 */
    check("DivisorSigma[-1, 28] == 2", "True");  /* perfect: sigma_-1 == 2 */
}

/* ---- Sign of n is ignored ---- */
void test_divisorsigma_sign() {
    check("DivisorSigma[1, -20]", "42");
    check("DivisorSigma[2, -6]", "50");
}

/* ---- Big integers ---- */
void test_divisorsigma_bignum() {
    check("DivisorSigma[2, 10^30]",
          "1388888888888888888587721618754026546109040824596043713112451");
    check("DivisorSigma[1, 2^60] == 2^61 - 1", "True");  /* sigma_1(2^60) */
}

/* ---- Listable threading ---- */
void test_divisorsigma_listable() {
    check("DivisorSigma[2, {1, 2, 3, 4, 5}]", "{1, 5, 10, 21, 26}");
    check("DivisorSigma[1, {6, 28, 12}]", "{12, 56, 28}");
    check("MemberQ[Attributes[DivisorSigma], Listable]", "True");
}

/* ---- Attributes ---- */
void test_divisorsigma_attributes() {
    check("Attributes[DivisorSigma]", "{Listable, NHoldAll, Protected}");
    check("MemberQ[Attributes[DivisorSigma], Protected]", "True");
    check("MemberQ[Attributes[DivisorSigma], NHoldAll]", "True");
}

/* ---- Symbolic exponent: multiplicative formula ---- */
void test_divisorsigma_symbolic() {
    check("DivisorSigma[k, 2]", "(-1 + 2^(2 k))/(-1 + 2^k)");
    check("DivisorSigma[k, 1]", "1");
    /* Substitution back to a concrete exponent recovers the numeric value. */
    check("DivisorSigma[k, 20] /. k -> 2", "546");
    check("DivisorSigma[k, 20] /. k -> 1", "42");
    check("DivisorSigma[k, 6] /. k -> 2", "50");
}

/* ---- Rational exponent: matches the sum of rational powers of divisors ---- */
void test_divisorsigma_rational() {
    check("DivisorSigma[1/2, 12]",
          "(2 (-1 + 2 Sqrt[2]))/((-1 + Sqrt[2]) (-1 + Sqrt[3]))");
    /* Algebraically equal to Total[Divisors[12]^(1/2)]. */
    check("PossibleZeroQ[DivisorSigma[1/2, 12] - Total[Divisors[12]^(1/2)]]",
          "True");
    check("PossibleZeroQ[DivisorSigma[1/3, 18] - Total[Divisors[18]^(1/3)]]",
          "True");
}

/* ---- Gaussian integers ---- */
void test_divisorsigma_gaussian() {
    /* Non-real input auto-enables Z[i]. */
    check("DivisorSigma[1, 3 + I]", "2 + 6*I");
    /* GaussianIntegers -> True uses the multiplicative formula over Z[i]. */
    check("DivisorSigma[2, 6, GaussianIntegers -> True]", "-30 + 20*I");
    check("DivisorSigma[1, 6, GaussianIntegers -> True]", "8 + 12*I");
    check("DivisorSigma[2, 100, GaussianIntegers -> True]", "6479 - 6018*I");
    /* sigma_0 over Z[i] counts the divisor classes. */
    check("DivisorSigma[0, 5, GaussianIntegers -> True]", "4");
    check("DivisorSigma[0, 5, GaussianIntegers -> True] == "
          "Length[Divisors[5, GaussianIntegers -> True]]", "True");
}

/* ---- Unevaluated cases ---- */
void test_divisorsigma_unevaluated() {
    check("DivisorSigma[2, x]", "DivisorSigma[2, x]");
    check("DivisorSigma[k, n]", "DivisorSigma[k, n]");
    check("DivisorSigma[2, 0]", "DivisorSigma[2, 0]");
    check("DivisorSigma[2, 5/2]", "DivisorSigma[2, 5/2]");
    check("DivisorSigma[2, 2.5]", "DivisorSigma[2, 2.5]");
    /* Wrong argument counts stay unevaluated (argr / argrx / nonopt). */
    check("DivisorSigma[]", "DivisorSigma[]");
    check("DivisorSigma[5]", "DivisorSigma[5]");
    check("DivisorSigma[1, 2, 3]", "DivisorSigma[1, 2, 3]");
}

/* ---- Options and messages (expectations verified against the reference) ---- */
void test_divisorsigma_options() {
    check("Options[DivisorSigma]", "{GaussianIntegers -> False}");
    check("DivisorSigma[1, 6, GaussianIntegers :> True]", "8 + 12*I");
    /* Of repeated settings the FIRST wins. */
    check("DivisorSigma[1, 6, GaussianIntegers -> True, GaussianIntegers -> False]",
          "8 + 12*I");
    check("DivisorSigma[1, 6, GaussianIntegers -> False]", "12");
    /* Each error is a real message: Check sees it by name, Quiet silences it,
     * and the call stays unevaluated. */
    check("Check[DivisorSigma[1], m, DivisorSigma::argr]", "m");
    check("Check[DivisorSigma[], m, DivisorSigma::argrx]", "m");
    check("Check[DivisorSigma[1, 2, 3], m, DivisorSigma::nonopt]", "m");
    check("Check[DivisorSigma[1, 6, Foo -> 1], m, DivisorSigma::optx]", "m");
    check("Check[DivisorSigma[1, 6, GaussianIntegers -> 3], m, DivisorSigma::opttf]", "m");
    check("Check[DivisorSigma[1, 6, GaussianIntegers -> Automatic], m, DivisorSigma::opttf]",
          "m");
    check("Quiet[DivisorSigma[1, 6, Foo -> 1]]", "DivisorSigma[1, 6, Foo -> 1]");
    check("Quiet[DivisorSigma[1, 6, GaussianIntegers -> Automatic]]",
          "DivisorSigma[1, 6, GaussianIntegers -> Automatic]");
    check("Check[DivisorSigma[1, 6], m]", "12");
}

/* ---- Exact machine-k fast path: int64, GMP and Rational branches ---- */
void test_divisorsigma_fastpath() {
    check("DivisorSigma[30, 6]", "221073919926625563736250");          /* int64 overflow */
    check("DivisorSigma[25, 10]", "10000000298023223910507558");
    check("DivisorSigma[-30, 6]",
          "110536959963312781868125/110536959860366678949888");
    check("DivisorSigma[-1, -12]", "7/3");
    check("DivisorSigma[-2, 2^64]",
          "453709822561251284617832809909024281941/340282366920938463463374607431768211456");
    check("DivisorSigma[2, -2^63]", "113427455640312821154458202477256070485");
    check("DivisorSigma[1, -9223372036854775807]", "10994507040830097408");
    check("DivisorSigma[1, 2^64 + 1]", "18446811354131136516");
    /* Either side of the scalar trial-division ceiling (2^20). */
    check("DivisorSigma[2, 1048573]", "1099505336330");
    check("DivisorSigma[3, 1048583]", "1152944594505171288");
    check("DivisorSigma[10^20, 1]", "1");
}

/* ---- Machine-real k: Horner in double ---- */
void test_divisorsigma_real() {
    check("DivisorSigma[2., 6]", "50.0");
    check("DivisorSigma[-2., 10]", "1.3");
    /* n = 1 is the empty product, an exact 1. */
    check("DivisorSigma[2., {1, 2, 3}]", "{1, 5.0, 10.0}");
    check("Abs[DivisorSigma[1.5, 6] - 23.72151800415189] < 10^-13", "True");
    /* The closed form is 0/0 here (Indeterminate); Horner gives the
     * divisor count, which is the value of sigma_0. A deliberate difference. */
    check("DivisorSigma[0., 6]", "4.0");
    /* Arbitrary-precision k keeps the expression path. */
    check("Precision[DivisorSigma[2.5`30, 6]] > 25", "True");
}

/* ---- Every machine surface gives the same answer ---- */
void test_divisorsigma_surfaces() {
    /* Packed buffer (per-element kernel) vs element-by-element scalar. */
    check("With[{d = DeleteCases[Range[-3000, 3000], 0]}, "
          "DivisorSigma[3, d] === Map[DivisorSigma[3, #] &, d]]", "True");
    /* Large enough to take the smallest-prime-factor sieve. */
    check("With[{d = Range[5*10^6, 5*10^6 + 20000]}, "
          "DivisorSigma[2, d] === Map[DivisorSigma[2, #] &, d]]", "True");
    check("DivisorSigma[2, Range[5*10^6, 5*10^6 + 5]]",
          "{34720102945736, 27790363736200, 31976112505000, 25005921785000, "
          "36863484909342, 26002601040104}");
    check("With[{d = Join[Range[-40000, -1], Range[1, 40000]]}, "
          "DivisorSigma[1, d] === Map[DivisorSigma[1, #] &, d]]", "True");
    /* An element with no int64 answer abandons the buffer, never wraps. */
    check("DivisorSigma[1, {1, 2, 0, 4}]", "{1, 3, DivisorSigma[1, 0], 7}");
    check("DivisorSigma[25, Range[10]] === Map[DivisorSigma[25, #] &, Range[10]]", "True");
    check("DivisorSigma[-1, Range[6]]", "{1, 3/2, 4/3, 7/4, 6/5, 2}");
    /* Visible int64 NDArray. */
    check("Normal[DivisorSigma[3, NDArray[Range[200], DataType -> \"int64\"]]] === "
          "DivisorSigma[3, Range[200]]", "True");
    /* Compile: scalar and array forms both lower, and agree. */
    check("\"Compiled\" /. CompileDiagnostics[{{k, _Integer}, {x, _Integer}}, DivisorSigma[k, x]]",
          "True");
    check("\"Compiled\" /. CompileDiagnostics[{{v, _Integer, 1}}, DivisorSigma[1, v]]", "True");
    check("With[{f = Compile[{{x, _Integer}}, DivisorSigma[3, x]]}, "
          "Map[f, Range[-500, -1]] === DivisorSigma[3, Range[-500, -1]]]", "True");
    check("With[{f = Compile[{{v, _Integer, 1}}, DivisorSigma[3, v]]}, "
          "f[Range[500]] === DivisorSigma[3, Range[500]]]", "True");
    /* A compiled call with no machine answer hands back to the interpreter. */
    check("With[{f = Compile[{{k, _Integer}, {x, _Integer}}, DivisorSigma[k, x]]}, "
          "{f[1, 12], f[-1, 6], f[40, 6] == DivisorSigma[40, 6]}]", "{28, 2, True}");
}

int main() {
    symtab_init();
    core_init();

    TEST(test_divisorsigma_basic);
    TEST(test_divisorsigma_count);
    TEST(test_divisorsigma_negative);
    TEST(test_divisorsigma_sign);
    TEST(test_divisorsigma_bignum);
    TEST(test_divisorsigma_listable);
    TEST(test_divisorsigma_attributes);
    TEST(test_divisorsigma_symbolic);
    TEST(test_divisorsigma_rational);
    TEST(test_divisorsigma_gaussian);
    TEST(test_divisorsigma_unevaluated);
    TEST(test_divisorsigma_options);
    TEST(test_divisorsigma_fastpath);
    TEST(test_divisorsigma_real);
    TEST(test_divisorsigma_surfaces);

    printf("All DivisorSigma tests passed!\n");
    return 0;
}
