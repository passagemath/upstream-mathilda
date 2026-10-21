/* Unit tests for the bit-level integer builtins (src/bitwise/).
 *
 *   BitLength[n] -- number of binary bits needed to represent the integer n.
 *
 * BitXor[n1, n2, ...] -- bitwise XOR (two's complement), Flat/Orderless/
 * OneIdentity/Listable: documented examples, bignums across the int64 boundary,
 * pair cancellation / integer folding / BitNot extraction, properties against
 * an IntegerDigits oracle, packed + NDArray (binary and n-ary), cross-surface
 * bit-equality incl. Compile, CompileDiagnostics, ::int + Quiet/Check.
 *
 * BitLength coverage:
 *   - Small positives, powers of two (2^k -> k+1, 2^k-1 -> k), zero.
 *   - Negatives: BitLength[n] = BitLength[BitNot[n]], BitNot[n] = -n-1, so
 *     BitLength[-1]=0, BitLength[-2]=1, BitLength[-2^k]=k.
 *   - int64 boundary incl. INT64_MIN (-2^63 -> 63), where the int64 kernel's
 *     complement path is exact and IntegerLength's negate path would overflow.
 *   - Bignums: 2^100 -> 101, 2^100-1 -> 100, -2^100 -> 100, 100! -> 525.
 *   - Listable threading (flat + nested) and the packed / visible-NDArray fast
 *     paths (single-headed int64 result, packed == visible == List).
 *   - Compile[] round-trip at scalar and rank-1 array shapes.
 *   - Property tests: BitLength[n] == IntegerLength[n, 2] and
 *     == Floor[Log[2, n]] + 1 for n > 0.
 *   - Error paths: real / rational / complex arg (::int), wrong arity (::argx),
 *     symbolic arg (silent, unevaluated).
 *   - Attribute, docstring, interned-symbol introspection.
 *   - Repeated-evaluation stress loop to catch double-frees under valgrind.
 */

#include "expr.h"
#include "eval.h"
#include "core.h"
#include "symtab.h"
#include "test_utils.h"
#include "parse.h"
#include "print.h"
#include "attr.h"
#include "sym_names.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>

/* Capture stderr while `input` is parsed + evaluated.  Returns the collected
 * stderr text as a heap string (caller frees) and writes the printed result
 * into *out_result_str (also heap-allocated).  Fixed temp path; safe because
 * tests run serially.  The original stderr fd is saved and restored with
 * dup/dup2 (as test_digit_sum.c does) rather than reopened as /dev/tty, which
 * does not exist in a headless run and left stderr closed for the next
 * capture. */
static char* eval_capturing_stderr(const char* input, char** out_result_str) {
    const char* path = "/tmp/mathilda_bitwise_stderr.log";
    fflush(stderr);
    int saved_stderr = dup(STDERR_FILENO);
    if (!freopen(path, "w+", stderr)) {
        if (saved_stderr != -1) close(saved_stderr);
        if (out_result_str) *out_result_str = NULL;
        return NULL;
    }

    Expr* p = parse_expression(input);
    Expr* e = evaluate(p);
    if (out_result_str) *out_result_str = expr_to_string(e);
    expr_free(p);
    expr_free(e);

    fflush(stderr);
    if (saved_stderr != -1) {
        dup2(saved_stderr, STDERR_FILENO);
        close(saved_stderr);
    }

    FILE* f = fopen(path, "r");
    if (!f) return NULL;
    fseek(f, 0, SEEK_END);
    long n = ftell(f);
    if (n < 0) { fclose(f); return NULL; }
    fseek(f, 0, SEEK_SET);
    char* buf = malloc((size_t)n + 1);
    if (!buf) { fclose(f); return NULL; }
    size_t got = fread(buf, 1, (size_t)n, f);
    buf[got] = '\0';
    fclose(f);
    remove(path);
    return buf;
}

/* --- Small positives ------------------------------------------------- */

static void test_zero(void) {
    assert_eval_eq("BitLength[0]", "0", 0);
}

static void test_small_positives(void) {
    assert_eval_eq("BitLength[1]", "1", 0);
    assert_eval_eq("BitLength[2]", "2", 0);
    assert_eval_eq("BitLength[3]", "2", 0);
    assert_eval_eq("BitLength[4]", "3", 0);
    assert_eval_eq("BitLength[5]", "3", 0);
    assert_eval_eq("BitLength[7]", "3", 0);
    assert_eval_eq("BitLength[8]", "4", 0);
}

static void test_byte_boundaries(void) {
    assert_eval_eq("BitLength[255]", "8", 0);
    assert_eval_eq("BitLength[256]", "9", 0);
    assert_eval_eq("BitLength[1023]", "10", 0);
    assert_eval_eq("BitLength[1024]", "11", 0);
}

static void test_powers_of_two(void) {
    /* 2^k needs k+1 bits; 2^k - 1 needs k bits. */
    assert_eval_eq("BitLength[2^10]", "11", 0);
    assert_eval_eq("BitLength[2^10 - 1]", "10", 0);
    assert_eval_eq("BitLength[2^30]", "31", 0);
    assert_eval_eq("BitLength[2^30 - 1]", "30", 0);
}

/* --- Negatives (BitNot semantics) ------------------------------------ */

static void test_negatives(void) {
    /* BitLength[n] = BitLength[BitNot[n]], BitNot[n] = -n-1. */
    assert_eval_eq("BitLength[-1]", "0", 0);   /* BitNot[-1] = 0  */
    assert_eval_eq("BitLength[-2]", "1", 0);   /* BitNot[-2] = 1  */
    assert_eval_eq("BitLength[-3]", "2", 0);   /* BitNot[-3] = 2  */
    assert_eval_eq("BitLength[-4]", "2", 0);   /* BitNot[-4] = 3, BitLength[3] = 2 (WL: 2) */
    assert_eval_eq("BitLength[-8]", "3", 0);   /* BitNot[-8] = 7  */
    assert_eval_eq("BitLength[-256]", "8", 0); /* BitNot[-256] = 255 */
}

static void test_negative_powers_of_two(void) {
    /* BitLength[-2^k] = k, since BitNot[-2^k] = 2^k - 1. */
    assert_eval_eq("BitLength[-2^10]", "10", 0);
    assert_eval_eq("BitLength[-2^30]", "30", 0);
}

/* --- int64 boundary, including INT64_MIN ----------------------------- */

static void test_int64_boundary(void) {
    assert_eval_eq("BitLength[2^62]", "63", 0);
    assert_eval_eq("BitLength[2^63]", "64", 0);      /* just past int64 max */
    assert_eval_eq("BitLength[2^63 - 1]", "63", 0);  /* int64 max          */
}

static void test_int64_min(void) {
    /* -2^63 is INT64_MIN.  BitNot[-2^63] = 2^63 - 1, which has 63 bits.  The
     * int64 kernel complements the unsigned magnitude, so this is exact where
     * a negate-based kernel would overflow. */
    assert_eval_eq("BitLength[-9223372036854775808]", "63", 0);
    assert_eval_eq("BitLength[-2^63]", "63", 0);
}

/* --- Bignums --------------------------------------------------------- */

static void test_bignums(void) {
    assert_eval_eq("BitLength[2^100]", "101", 0);
    assert_eval_eq("BitLength[2^100 - 1]", "100", 0);
    assert_eval_eq("BitLength[-2^100]", "100", 0);
    /* Documented: IntegerLength[100!, 2] == 525, so BitLength[100!] == 525. */
    assert_eval_eq("BitLength[100!]", "525", 0);
}

/* --- Listable threading ---------------------------------------------- */

static void test_listable_flat(void) {
    assert_eval_eq("BitLength[{0, 1, 2, 7, 8, 255, 256}]",
                   "{0, 1, 2, 3, 4, 8, 9}", 0);
}

static void test_listable_negatives(void) {
    assert_eval_eq("BitLength[{-1, -2, -8, -256}]", "{0, 1, 3, 8}", 0);
}

static void test_listable_nested(void) {
    assert_eval_eq("BitLength[{{1, 10}, {100, 1000}}]",
                   "{{1, 4}, {7, 10}}", 0);
}

/* --- Packed / NDArray fast paths ------------------------------------- */

static void test_packed_range(void) {
    assert_eval_eq("BitLength[Range[8]]", "{1, 2, 2, 3, 3, 3, 3, 4}", 0);
}

static void test_visible_ndarray_int64(void) {
    /* A visible int64 NDArray must produce a single-headed int64 array whose
     * elements agree with the List path, INT64_MIN included. */
    assert_eval_eq(
        "BitLength[NDArray[{-9223372036854775808, -1, 0, 1, 255}, "
        "DataType -> \"int64\"]]",
        "NDArray[{63, 0, 0, 1, 8}]", 0);
}

static void test_packed_agrees_with_list(void) {
    /* The packed producer (Range) and the plain List must give the same
     * answer -- the surfaces are not allowed to diverge. */
    assert_eval_eq("BitLength[Range[16]] == BitLength[Table[n, {n, 1, 16}]]",
                   "True", 0);
}

/* --- Compile[] round-trip -------------------------------------------- */

static void test_compile_scalar(void) {
    assert_eval_eq("Compile[{{x, _Integer}}, BitLength[x]][13]", "4", 0);
    assert_eval_eq("Compile[{{x, _Integer}}, BitLength[x]][-8]", "3", 0);
    assert_eval_eq("Compile[{{x, _Integer}}, BitLength[x]][0]", "0", 0);
}

static void test_compile_array(void) {
    assert_eval_eq(
        "Compile[{{v, _Integer, 1}}, BitLength[v]][{1, 7, 8, 255}]",
        "{1, 3, 4, 8}", 0);
}

/* --- Property tests -------------------------------------------------- */

static void test_matches_integerlength_base2(void) {
    assert_eval_eq(
        "Table[BitLength[n] == IntegerLength[n, 2], {n, 1, 64}] === "
        "Table[True, {64}]",
        "True", 0);
}

static void test_floor_log_identity(void) {
    assert_eval_eq(
        "Table[BitLength[n] == Floor[Log[2, n]] + 1, {n, 1, 40}] === "
        "Table[True, {40}]",
        "True", 0);
}

/* --- Error paths ----------------------------------------------------- */

static void test_real_int_diagnostic(void) {
    char* result = NULL;
    char* err = eval_capturing_stderr("BitLength[3.5]", &result);
    ASSERT(result != NULL);
    ASSERT_STR_EQ(result, "BitLength[3.5]");
    ASSERT(err != NULL);
    ASSERT_MSG(strstr(err, "BitLength::int") != NULL,
               "expected int diagnostic, got: %s", err);
    ASSERT(strstr(err, "Integer expected at position 1") != NULL);
    free(result);
    free(err);
}

static void test_rational_int_diagnostic(void) {
    char* result = NULL;
    char* err = eval_capturing_stderr("BitLength[3/2]", &result);
    ASSERT(result != NULL);
    ASSERT(err != NULL);
    ASSERT(strstr(err, "BitLength::int") != NULL);
    ASSERT(strstr(err, "Integer expected at position 1") != NULL);
    free(result);
    free(err);
}

static void test_complex_int_diagnostic(void) {
    char* result = NULL;
    char* err = eval_capturing_stderr("BitLength[1 + 2 I]", &result);
    ASSERT(result != NULL);
    ASSERT(err != NULL);
    ASSERT_MSG(strstr(err, "BitLength::int") != NULL,
               "expected int diagnostic, got: %s", err);
    free(result);
    free(err);
}

static void test_zero_args_argx(void) {
    char* result = NULL;
    char* err = eval_capturing_stderr("BitLength[]", &result);
    ASSERT(result != NULL);
    ASSERT_STR_EQ(result, "BitLength[]");
    ASSERT(err != NULL);
    ASSERT_MSG(strstr(err, "BitLength::argx") != NULL,
               "expected argx diagnostic, got: %s", err);
    ASSERT(strstr(err, "called with 0 arguments") != NULL);
    ASSERT(strstr(err, "1 argument is expected") != NULL);
    free(result);
    free(err);
}

static void test_two_args_argx(void) {
    char* result = NULL;
    char* err = eval_capturing_stderr("BitLength[1, 2]", &result);
    ASSERT(result != NULL);
    ASSERT_STR_EQ(result, "BitLength[1, 2]");
    ASSERT(err != NULL);
    ASSERT_MSG(strstr(err, "BitLength::argx") != NULL,
               "expected argx diagnostic, got: %s", err);
    ASSERT(strstr(err, "called with 2 arguments") != NULL);
    free(result);
    free(err);
}

static void test_symbolic_left_unevaluated(void) {
    char* result = NULL;
    char* err = eval_capturing_stderr("BitLength[x]", &result);
    ASSERT(result != NULL);
    ASSERT(strstr(result, "BitLength[x]") != NULL);
    ASSERT_MSG(err == NULL || strstr(err, "BitLength::") == NULL,
               "expected no diagnostic for symbolic input, got: %s",
               err ? err : "(null)");
    free(result);
    free(err);
}

/* --- Introspection --------------------------------------------------- */

static void test_attributes(void) {
    SymbolDef* def = symtab_get_def("BitLength");
    ASSERT(def != NULL);
    uint32_t a = get_attributes("BitLength");
    ASSERT((a & ATTR_PROTECTED) != 0);
    ASSERT((a & ATTR_LISTABLE) != 0);
}

static void test_docstring_set(void) {
    SymbolDef* def = symtab_get_def("BitLength");
    ASSERT(def != NULL);
    ASSERT(def->docstring != NULL);
    ASSERT(strstr(def->docstring, "binary bits") != NULL);
}

static void test_sym_pointer_interned(void) {
    ASSERT(SYM_BitLength != NULL);
    ASSERT(strcmp(SYM_BitLength, "BitLength") == 0);
}

/* --- Memory-safety stress loop --------------------------------------- */

static void test_repeated_evaluation_does_not_corrupt(void) {
    /* Mix machine-int / bignum / negative / INT64_MIN / listable / packed /
     * error paths so any mishandled mpz_t or args[] slot surfaces under
     * valgrind. */
    for (int k = 0; k < 50; k++) {
        assert_eval_eq("BitLength[255]", "8", 0);
        assert_eval_eq("BitLength[-8]", "3", 0);
        assert_eval_eq("BitLength[0]", "0", 0);
        assert_eval_eq("BitLength[-2^63]", "63", 0);
        assert_eval_eq("BitLength[2^100]", "101", 0);
        assert_eval_eq("BitLength[{1, 10, 100, 1000}]", "{1, 4, 7, 10}", 0);
        assert_eval_eq("BitLength[Range[8]]", "{1, 2, 2, 3, 3, 3, 3, 4}", 0);
    }
}


/* ===================================================================== *
 *  BitXor[n1, n2, ...]                                                   *
 * ===================================================================== */

/* --- Documented examples --------------------------------------------- */

static void test_bitxor_basic(void) {
    assert_eval_eq("BitXor[61, 15]", "50", 0);
    assert_eval_eq("BitXor[3333, 5555, 7777, 9999]", "8664", 0);
    assert_eval_eq("BitXor[-2, -3]", "3", 0);
    assert_eval_eq("BitXor[-6, 3]", "-7", 0);
    assert_eval_eq("BitXor[]", "0", 0);
    assert_eval_eq("BitXor[5]", "5", 0);
    assert_eval_eq("BitXor[0, 0]", "0", 0);
    assert_eval_eq("BitXor[7, 7]", "0", 0);
    assert_eval_eq("BitXor[-1, -1]", "0", 0);
}

static void test_bitxor_range_threading(void) {
    assert_eval_eq("BitXor[Range[120], 5]",
        "{4, 7, 6, 1, 0, 3, 2, 13, 12, 15, 14, 9, 8, 11, 10, 21, 20, 23, 22, "
        "17, 16, 19, 18, 29, 28, 31, 30, 25, 24, 27, 26, 37, 36, 39, 38, 33, "
        "32, 35, 34, 45, 44, 47, 46, 41, 40, 43, 42, 53, 52, 55, 54, 49, 48, "
        "51, 50, 61, 60, 63, 62, 57, 56, 59, 58, 69, 68, 71, 70, 65, 64, 67, "
        "66, 77, 76, 79, 78, 73, 72, 75, 74, 85, 84, 87, 86, 81, 80, 83, 82, "
        "93, 92, 95, 94, 89, 88, 91, 90, 101, 100, 103, 102, 97, 96, 99, 98, "
        "109, 108, 111, 110, 105, 104, 107, 106, 117, 116, 119, 118, 113, "
        "112, 115, 114, 125}", 0);
}

static void test_bitxor_truth_table(void) {
    assert_eval_eq("Outer[BitXor, {1, 0}, {1, 0}]", "{{0, 1}, {1, 0}}", 0);
}

static void test_bitxor_empty_lists(void) {
    assert_eval_eq("BitXor[{}]", "{}", 0);
    assert_eval_eq("BitXor[{}, {}]", "{}", 0);
}

static void test_bitxor_listable_shapes(void) {
    assert_eval_eq("BitXor[{1, 2}, 3, {4, 5}]", "{6, 4}", 0);
    assert_eval_eq("BitXor[{{1, 2}, {3, 4}}, 1]", "{{0, 3}, {2, 5}}", 0);
    assert_eval_eq("BitXor[{1, 2}, {3, 4}]", "{2, 6}", 0);
}

/* --- Bignums --------------------------------------------------------- */

static void test_bitxor_bignums(void) {
    assert_eval_eq("BitXor[3^100, 5^100]",
        "7888609052210118054116793203790437104101594041461795136781437785275168", 0);
    /* Two's complement crosses the int64 boundary both ways. */
    assert_eval_eq("BitXor[2^63, -1]", "-9223372036854775809", 0);
    assert_eval_eq("BitXor[-2^63, 1]", "-9223372036854775807", 0);
    assert_eval_eq("BitXor[-2^63, -1]", "9223372036854775807", 0);
    assert_eval_eq("BitXor[2^63 - 1, -2^63]", "-1", 0);
    /* A bignum result that fits normalises back to a machine Integer. */
    assert_eval_eq("BitXor[2^70, 2^70 + 5]", "5", 0);
    assert_eval_eq("MachineIntegerQ[BitXor[2^70, 2^70 + 5]]", "True", 0);
    assert_eval_eq("BitXor[2^100, -2^100]", "-2535301200456458802993406410752", 0);
    assert_eval_eq("BitXor[-3^50, 7, -5^30]",
        "718817681789702120559719", 0);
    assert_eval_eq("BitXor[{2^64, 1}, 2^64]", "{0, 18446744073709551617}", 0);
}

static void test_bitxor_bignum_operands_untouched(void) {
    /* The bignum path XORs in place into a private result node: the operands
     * themselves -- bound, shared, or repeated -- must come out unchanged. */
    assert_eval_eq("p = 3^200; q = 5^150; r = BitXor[p, q, p, 7]; "
                   "{p == 3^200, q == 5^150, r == BitXor[q, 7]}",
                   "{True, True, True}", 0);
    assert_eval_eq("p = 2^100; {BitXor[p, p], p == 2^100}", "{0, True}", 0);
    assert_eval_eq("l = {2^90, 2^91}; BitXor[l, 2^90]; l == {2^90, 2^91}", "True", 0);
    /* machine constants folded into a bignum once, mixed signs */
    assert_eval_eq("BitXor[2^100, 1, 2, 4, -8]", "-1267650600228229401496703205377", 0);
}

/* --- Symbolic simplification ----------------------------------------- */

static void test_bitxor_pair_cancellation(void) {
    assert_eval_eq("BitXor[x, x]", "0", 0);
    assert_eval_eq("BitXor[x, x, x]", "x", 0);
    assert_eval_eq("BitXor[x, x, x, x]", "0", 0);
    assert_eval_eq("BitXor[x, x, x, x, x]", "x", 0);
    assert_eval_eq("BitXor[x, y, x, y]", "0", 0);
    assert_eval_eq("BitXor[x, y, x, y, x]", "x", 0);
    assert_eval_eq("BitXor[x, y, z, x]", "BitXor[y, z]", 0);
    assert_eval_eq("BitXor[f[a], f[a], 3]", "3", 0);
    assert_eval_eq("BitXor[a, 3, a]", "3", 0);
    /* A cancelled non-integer pair is gone, so it is never an error. */
    assert_eval_eq("BitXor[x, 2.5, 2.5]", "x", 0);
    assert_eval_eq("BitXor[2.5, 2.5]", "0", 0);
}

static void test_bitxor_integer_folding(void) {
    assert_eval_eq("BitXor[1, 2, x]", "BitXor[3, x]", 0);
    assert_eval_eq("BitXor[0, x]", "x", 0);
    assert_eval_eq("BitXor[x, 0, 0]", "x", 0);
    assert_eval_eq("BitXor[2, x, 2]", "x", 0);
    assert_eval_eq("BitXor[1, 2, x, x, y]", "BitXor[3, y]", 0);
    assert_eval_eq("BitXor[-5, x]", "BitXor[-5, x]", 0);
    assert_eval_eq("BitXor[x, 2^70, 2^70 + 1]", "BitXor[1, x]", 0);
}

static void test_bitxor_bitnot(void) {
    assert_eval_eq("BitXor[-1, x]", "BitNot[x]", 0);
    assert_eval_eq("BitXor[x, y, -1]", "BitNot[BitXor[x, y]]", 0);
    /* -1 is pulled out before the integers combine ... */
    assert_eval_eq("BitXor[1, -1, x]", "BitNot[BitXor[1, x]]", 0);
    assert_eval_eq("BitXor[3, -1, x]", "BitNot[BitXor[3, x]]", 0);
    assert_eval_eq("BitXor[-1, x, y, y, 3]", "BitNot[BitXor[3, x]]", 0);
    assert_eval_eq("BitXor[-1, 2^70, x]",
                   "BitNot[BitXor[1180591620717411303424, x]]", 0);
    /* ... and a combined constant of -1 is one too. */
    assert_eval_eq("BitXor[-2, 1, x]", "BitNot[x]", 0);
    assert_eval_eq("BitXor[-1, 0, x]", "BitNot[x]", 0);
    /* An even count cancels; everything else gone leaves the integer. */
    assert_eval_eq("BitXor[x, -1, -1]", "x", 0);
    assert_eval_eq("BitXor[-1, -1, -1, x]", "BitNot[x]", 0);
    assert_eval_eq("BitXor[-1, x, x]", "-1", 0);
    assert_eval_eq("BitXor[x, -1, y, y]", "BitNot[x]", 0);
    assert_eval_eq("BitXor[x, BitNot[x]]", "BitXor[x, BitNot[x]]", 0);
}

static void test_bitxor_orderless_flat(void) {
    assert_eval_eq("BitXor[y, x]", "BitXor[x, y]", 0);
    assert_eval_eq("BitXor[x]", "x", 0);
    assert_eval_eq("BitXor[BitXor[x, 1], BitXor[2, y]]", "BitXor[3, x, y]", 0);
    assert_eval_eq("BitXor[BitXor[x, y], x]", "y", 0);
}

/* --- Properties ------------------------------------------------------ */

static void test_bitxor_properties(void) {
    /* Self-inverse and associativity over a grid that includes negatives. */
    assert_eval_eq(
        "And @@ Flatten[Table[BitXor[BitXor[a, b], b] == a, "
        "{a, -20, 20}, {b, -20, 20}]]", "True", 0);
    assert_eval_eq(
        "And @@ Flatten[Table[BitXor[BitXor[a, b], c] == BitXor[a, BitXor[b, c]], "
        "{a, -6, 6}, {b, -6, 6}, {c, -6, 6}]]", "True", 0);
    /* Independent oracle: bitwise parity from the binary digits. */
    assert_eval_eq(
        "And @@ Flatten[Table[BitXor[a, b] == "
        "FromDigits[Mod[IntegerDigits[a, 2, 8] + IntegerDigits[b, 2, 8], 2], 2], "
        "{a, 0, 63}, {b, 0, 63}]]", "True", 0);
    /* Two's complement: BitXor[n, -1] == -n - 1 (BitNot). */
    assert_eval_eq(
        "And @@ Table[BitXor[n, -1] == -n - 1, {n, {-2^80, -7, 0, 9, 2^63, 3^60}}]",
        "True", 0);
}

/* --- Packed / NDArray surfaces --------------------------------------- */

static void test_bitxor_packed_and_ndarray(void) {
    assert_eval_eq("PackedArrayQ[BitXor[Range[10], 3]]", "True", 0);
    assert_eval_eq("BitXor[NDArray[{1, 2, 3}, DataType -> \"int64\"], 6]",
                   "NDArray[{7, 4, 5}]", 0);
    assert_eval_eq(
        "BitXor[NDArray[{-9223372036854775808, -1, 0}, DataType -> \"int64\"], "
        "NDArray[{1, -1, -1}, DataType -> \"int64\"]]",
        "NDArray[{-9223372036854775807, 0, -1}]", 0);
    /* n-ary: folded on the buffer, scalars combined first. */
    assert_eval_eq("BitXor[Range[10], 1, 2]",
                   "{2, 1, 0, 7, 6, 5, 4, 11, 10, 9}", 0);
    assert_eval_eq("PackedArrayQ[BitXor[Range[10], 1, 2]]", "True", 0);
    assert_eval_eq("BitXor[Range[5], Range[5], Reverse[Range[5]]]",
                   "{5, 4, 3, 2, 1}", 0);
    assert_eval_eq("BitXor[NDArray[{1, 2, 3}, DataType -> \"int64\"], 1, 2]",
                   "NDArray[{2, 1, 0}]", 0);
    assert_eval_eq("BitXor[Range[4], 1, 1]", "{1, 2, 3, 4}", 0);
    /* Off-buffer operands take the exact List path. */
    assert_eval_eq("BitXor[Range[3], 2^70]",
        "{1180591620717411303425, 1180591620717411303426, 1180591620717411303427}", 0);
    assert_eval_eq("BitXor[Range[3], x, x]", "{1, 2, 3}", 0);
    assert_eval_eq("BitXor[Range[3], -1]", "{-2, -3, -4}", 0);
}

static void test_bitxor_surfaces_agree(void) {
    assert_eval_eq(
        "SeedRandom[3]; a = Join[RandomInteger[{-2^62, 2^62}, 300], {-2^63, 2^63 - 1, 0, -1}]; "
        "b = Join[RandomInteger[{-2^62, 2^62}, 300], {-1, -2^63, 5, -2^63}]; "
        "ref = Table[BitXor[a[[i]], b[[i]]], {i, Length[a]}]; "
        "{BitXor[a, b] === ref, "
        " Normal[BitXor[NDArray[a, DataType -> \"int64\"], NDArray[b, DataType -> \"int64\"]]] === ref, "
        " Compile[{{v, _Integer, 1}, {w, _Integer, 1}}, BitXor[v, w]][a, b] === ref, "
        " Table[Compile[{{x, _Integer}, {y, _Integer}}, BitXor[x, y]][a[[i]], b[[i]]], {i, Length[a]}] === ref, "
        " BitXor[a, b, 99] === Table[BitXor[a[[i]], b[[i]], 99], {i, Length[a]}]}",
        "{True, True, True, True, True}", 0);
}

/* --- Compile[] ------------------------------------------------------- */

static void test_bitxor_compile(void) {
    assert_eval_eq("Compile[{{a, _Integer}, {b, _Integer}}, BitXor[a, b]][61, 15]", "50", 0);
    assert_eval_eq("Compile[{{a, _Integer}, {b, _Integer}}, BitXor[a, b]][-2, -3]", "3", 0);
    assert_eval_eq(
        "Compile[{{a, _Integer}, {b, _Integer}, {c, _Integer}, {d, _Integer}}, "
        "BitXor[a, b, c, d]][3333, 5555, 7777, 9999]", "8664", 0);
    assert_eval_eq("Compile[{{v, _Integer, 1}}, BitXor[v, 5]][Range[8]]",
                   "{4, 7, 6, 1, 0, 3, 2, 13}", 0);
    assert_eval_eq("Compile[{{v, _Integer, 1}}, BitXor[6, v, 4]][Range[4]]",
                   "{3, 0, 1, 6}", 0);
    assert_eval_eq(
        "Compile[{{v, _Integer, 1}, {w, _Integer, 1}}, BitXor[v, w]][{1, 2, 3}, {3, 2, 1}]",
        "{2, 0, 2}", 0);
}

static void test_bitxor_compile_diagnostics(void) {
    assert_eval_eq("\"Compiled\" /. CompileDiagnostics[{{a, _Integer}, {b, _Integer}}, BitXor[a, b]]", "True", 0);
    assert_eval_eq("\"Compiled\" /. CompileDiagnostics[{{a, _Integer}, {b, _Integer}, {c, _Integer}}, BitXor[a, b, c]]", "True", 0);
    assert_eval_eq("\"Compiled\" /. CompileDiagnostics[{{v, _Integer, 1}}, BitXor[v, 5]]", "True", 0);
    assert_eval_eq("\"Compiled\" /. CompileDiagnostics[{{v, _Integer, 1}, {w, _Integer, 1}}, BitXor[v, w]]", "True", 0);
    assert_eval_eq("\"Compiled\" /. CompileDiagnostics[{{v, _Integer, 1}}, BitXor[v, 6, 4]]", "True", 0);
}

/* --- Error paths ----------------------------------------------------- */

static void test_bitxor_real_int_diagnostic(void) {
    char* result = NULL;
    char* err = eval_capturing_stderr("BitXor[2.5, 3]", &result);
    ASSERT(result != NULL);
    ASSERT_STR_EQ(result, "BitXor[2.5, 3]");
    ASSERT(err != NULL);
    ASSERT_MSG(strstr(err, "BitXor::int") != NULL,
               "expected int diagnostic, got: %s", err);
    ASSERT(strstr(err, "Integer expected at position 1 in BitXor[2.5, 3]") != NULL);
    free(result);
    free(err);
}

static void test_bitxor_rational_complex_diagnostic(void) {
    const char* in[] = { "BitXor[1/2, x]", "BitXor[I, 3]", "BitXor[-1, x, 2.5]" };
    for (size_t i = 0; i < sizeof in / sizeof in[0]; i++) {
        char* result = NULL;
        char* err = eval_capturing_stderr(in[i], &result);
        ASSERT(result != NULL);
        ASSERT(strstr(result, "BitXor[") != NULL);    /* left unevaluated */
        ASSERT_MSG(err && strstr(err, "BitXor::int") != NULL,
                   "expected int diagnostic for %s, got: %s", in[i], err ? err : "(null)");
        free(result);
        free(err);
    }
}

static void test_bitxor_quiet_check(void) {
    assert_eval_eq("Check[BitXor[2.5, 3], failed]", "failed", 0);
    assert_eval_eq("Check[BitXor[x, x, 7], failed]", "7", 0);
    char* result = NULL;
    char* err = eval_capturing_stderr("Quiet[BitXor[2.5, 3]]", &result);
    ASSERT(result != NULL);
    ASSERT_STR_EQ(result, "BitXor[2.5, 3]");
    ASSERT_MSG(err == NULL || strstr(err, "BitXor::") == NULL,
               "Quiet leaked: %s", err ? err : "(null)");
    free(result);
    free(err);
}

static void test_bitxor_symbolic_silent(void) {
    char* result = NULL;
    char* err = eval_capturing_stderr("BitXor[x, y]", &result);
    ASSERT(result != NULL);
    ASSERT_STR_EQ(result, "BitXor[x, y]");
    ASSERT_MSG(err == NULL || strstr(err, "BitXor::") == NULL,
               "expected no diagnostic for symbolic input, got: %s", err ? err : "(null)");
    free(result);
    free(err);
}

/* --- Introspection --------------------------------------------------- */

static void test_bitxor_attributes(void) {
    assert_eval_eq("Attributes[BitXor]",
                   "{Flat, Listable, OneIdentity, Orderless, Protected}", 0);
}

static void test_bitxor_docstring_and_sym(void) {
    SymbolDef* def = symtab_get_def("BitXor");
    ASSERT(def != NULL);
    ASSERT(def->docstring != NULL);
    ASSERT(strstr(def->docstring, "bitwise XOR") != NULL);
    ASSERT(SYM_BitXor != NULL && strcmp(SYM_BitXor, "BitXor") == 0);
    ASSERT(SYM_BitNot != NULL && strcmp(SYM_BitNot, "BitNot") == 0);
}

static void test_bitxor_stress(void) {
    for (int k = 0; k < 50; k++) {
        assert_eval_eq("BitXor[61, 15]", "50", 0);
        assert_eval_eq("BitXor[3^100, 5^100, -7]", "-7888609052210118054116793203790437104101594041461795136781437785275175", 0);
        assert_eval_eq("BitXor[1, -1, x]", "BitNot[BitXor[1, x]]", 0);
        assert_eval_eq("BitXor[x, y, x, 2^70, 3]", "BitXor[1180591620717411303427, y]", 0);
        assert_eval_eq("BitXor[Range[6], 1, 2]", "{2, 1, 0, 7, 6, 5}", 0);
        assert_eval_eq("BitXor[Range[3], 2^70]",
            "{1180591620717411303425, 1180591620717411303426, 1180591620717411303427}", 0);
    }
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_zero);
    TEST(test_small_positives);
    TEST(test_byte_boundaries);
    TEST(test_powers_of_two);

    TEST(test_negatives);
    TEST(test_negative_powers_of_two);

    TEST(test_int64_boundary);
    TEST(test_int64_min);

    TEST(test_bignums);

    TEST(test_listable_flat);
    TEST(test_listable_negatives);
    TEST(test_listable_nested);

    TEST(test_packed_range);
    TEST(test_visible_ndarray_int64);
    TEST(test_packed_agrees_with_list);

    TEST(test_compile_scalar);
    TEST(test_compile_array);

    TEST(test_matches_integerlength_base2);
    TEST(test_floor_log_identity);

    TEST(test_real_int_diagnostic);
    TEST(test_rational_int_diagnostic);
    TEST(test_complex_int_diagnostic);
    TEST(test_zero_args_argx);
    TEST(test_two_args_argx);
    TEST(test_symbolic_left_unevaluated);

    TEST(test_attributes);
    TEST(test_docstring_set);
    TEST(test_sym_pointer_interned);

    TEST(test_repeated_evaluation_does_not_corrupt);

    TEST(test_bitxor_basic);
    TEST(test_bitxor_range_threading);
    TEST(test_bitxor_truth_table);
    TEST(test_bitxor_empty_lists);
    TEST(test_bitxor_listable_shapes);
    TEST(test_bitxor_bignums);
    TEST(test_bitxor_bignum_operands_untouched);
    TEST(test_bitxor_pair_cancellation);
    TEST(test_bitxor_integer_folding);
    TEST(test_bitxor_bitnot);
    TEST(test_bitxor_orderless_flat);
    TEST(test_bitxor_properties);
    TEST(test_bitxor_packed_and_ndarray);
    TEST(test_bitxor_surfaces_agree);
    TEST(test_bitxor_compile);
    TEST(test_bitxor_compile_diagnostics);
    TEST(test_bitxor_real_int_diagnostic);
    TEST(test_bitxor_rational_complex_diagnostic);
    TEST(test_bitxor_quiet_check);
    TEST(test_bitxor_symbolic_silent);
    TEST(test_bitxor_attributes);
    TEST(test_bitxor_docstring_and_sym);
    TEST(test_bitxor_stress);

    printf("All BitLength / BitXor tests passed!\n");
    return 0;
}
