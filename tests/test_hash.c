/* Unit tests for the Hash builtin (the default "Expression" hash type).
 *
 * The "Expression" hash is implementation-defined -- its numeric value may
 * change between builds (if expr_hash changes) and between systems -- so the
 * tests are written as BUILD-INDEPENDENT invariants (determinism, argument
 * evaluation, format round-trips, structural sensitivity) plus a handful of
 * direct-C cross-checks that recompute expr_hash() locally and format the
 * result exactly as the builtin should, pinning every string format (notably
 * the little-endian byte order) regardless of the actual hash value. */

#include "expr.h"
#include "eval.h"
#include "core.h"
#include "symtab.h"
#include "test_utils.h"
#include "parse.h"
#include "print.h"
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <inttypes.h>

/* ---- shared drivers ---------------------------------------------------- */

/* Parse, evaluate, compare the printed (InputForm) result. */
static void check(const char* input, const char* expected) {
    Expr* e = parse_expression(input);
    Expr* res = evaluate(e);
    char* res_str = expr_to_string(res);
    if (strcmp(res_str, expected) != 0) {
        printf("Hash test failed: %s\n  expected: %s\n  got:      %s\n",
               input, expected, res_str);
        ASSERT(0);
    }
    free(res_str);
    expr_free(e);
    expr_free(res);
}

/* Evaluate a string, returning the result Expr (caller frees). */
static Expr* eval_str(const char* input) {
    Expr* e = parse_expression(input);
    Expr* res = evaluate(e);
    expr_free(e);
    return res;
}

/* The Expression hash of the evaluated form of `input`. */
static uint64_t hash_of(const char* input) {
    Expr* r = eval_str(input);
    uint64_t h = expr_hash(r);
    expr_free(r);
    return h;
}

/* Reference base-36 render (mirrors core.c hash_format). out must hold >= 14. */
static void u64_to_base36(uint64_t h, char* out) {
    static const char digits[] = "0123456789abcdefghijklmnopqrstuvwxyz";
    char tmp[14];
    int pos = 14;
    tmp[--pos] = '\0';
    if (h == 0) {
        tmp[--pos] = '0';
    } else {
        while (h > 0) { tmp[--pos] = digits[h % 36]; h /= 36; }
    }
    memcpy(out, tmp + pos, (size_t)(14 - pos));
}

/* Assert that `input` evaluates to a String equal to `expected`. */
static void expect_string(const char* input, const char* expected) {
    Expr* r = eval_str(input);
    ASSERT_MSG(r->type == EXPR_STRING, "expected a String result from %s", input);
    ASSERT_STR_EQ(r->data.string, expected);
    expr_free(r);
}

/* ---- core behaviour ---------------------------------------------------- */

void test_hash_basic_integer() {
    /* Hash[expr] returns a non-negative Integer. */
    check("Head[Hash[5]]", "Integer");
    check("Hash[5] >= 0", "True");
    check("Head[Hash[\"abc\"]]", "Integer");
    check("Head[Hash[1.5]]", "Integer");
    check("Head[Hash[3/2]]", "Integer");
    check("Head[Hash[f[x, y]]]", "Integer");
    check("Head[Hash[{1, 2, 3}]]", "Integer");
    check("Head[Hash[2^100]]", "Integer");   /* bigint argument */
}

void test_hash_determinism_and_eval() {
    /* Same expression -> same hash. */
    check("Hash[1] == Hash[1]", "True");
    check("Hash[\"abc\"] == Hash[\"abc\"]", "True");
    /* Distinct expressions -> distinct hash (collision astronomically unlikely). */
    check("Hash[1] == Hash[2]", "False");
    /* The argument is evaluated first (Hash is not HoldAll). */
    check("Hash[1 + 1] == Hash[2]", "True");
    check("Hash[Range[3]] == Hash[{1, 2, 3}]", "True");
}

void test_hash_structural_sensitivity() {
    /* Order matters; different heads differ. */
    check("Hash[{1, 2}] == Hash[{2, 1}]", "False");
    check("Hash[f[x]] == Hash[g[x]]", "False");
    check("Hash[f[x]] == Hash[f[y]]", "False");
}

void test_hash_type_argument() {
    /* "Expression" is accepted and equals the bare form. */
    check("Hash[5, \"Expression\"] == Hash[5]", "True");
    check("Hash[42, \"Expression\"] == Hash[42]", "True");
}

void test_hash_format_equivalence() {
    /* The "Integer" format is the default. */
    check("Hash[5, \"Expression\", \"Integer\"] == Hash[5]", "True");
}

void test_hash_string_formats_roundtrip() {
    /* Build-independent: the string formats reproduce the integer hash. */
    check("Hash[5, \"Expression\", \"DecimalString\"] == ToString[Hash[5]]", "True");
    check("FromDigits[Hash[5, \"Expression\", \"DecimalString\"]] == Hash[5]", "True");
    check("FromDigits[Hash[5, \"Expression\", \"HexString\"], 16] == Hash[5]", "True");
    check("FromDigits[Hash[5, \"Expression\", \"Base36String\"], 36] == Hash[5]", "True");
    /* HexString is zero-padded to the full 64-bit width. */
    check("StringLength[Hash[5, \"Expression\", \"HexString\"]]", "16");
    check("StringLength[Hash[5, \"Expression\", \"HexStringLittleEndian\"]]", "16");
    check("Head[Hash[5, \"Expression\", \"HexStringLittleEndian\"]]", "String");
    check("Head[Hash[5, \"Expression\", \"Base36String\"]]", "String");
}

/* Direct-C exact cross-check: recompute the hash and format it locally, then
 * require the builtin to produce exactly that -- pins every format including
 * the little-endian byte order, independent of the hash value. */
void test_hash_format_exact() {
    const char* inputs[] = { "5", "\"abc\"", "f[x, y]", "{1, 2, 3}", "1.5" };
    for (size_t k = 0; k < sizeof(inputs) / sizeof(inputs[0]); k++) {
        uint64_t h = hash_of(inputs[k]);

        char dec[21], hex[17], hexle[17], b36[14];
        snprintf(dec, sizeof dec, "%" PRIu64, h);
        snprintf(hex, sizeof hex, "%016" PRIx64, h);
        for (int i = 0; i < 8; i++)
            snprintf(hexle + 2 * i, 3, "%02x",
                     (unsigned)((h >> (8 * i)) & 0xFFu));
        u64_to_base36(h, b36);

        char expr[128];
        snprintf(expr, sizeof expr, "Hash[%s, \"Expression\", \"DecimalString\"]", inputs[k]);
        expect_string(expr, dec);
        snprintf(expr, sizeof expr, "Hash[%s, \"Expression\", \"HexString\"]", inputs[k]);
        expect_string(expr, hex);
        snprintf(expr, sizeof expr, "Hash[%s, \"Expression\", \"HexStringLittleEndian\"]", inputs[k]);
        expect_string(expr, hexle);
        snprintf(expr, sizeof expr, "Hash[%s, \"Expression\", \"Base36String\"]", inputs[k]);
        expect_string(expr, b36);
    }
}

/* Exercise the bigint return path: find an input whose 64-bit hash has its
 * high bit set (so it cannot fit a signed int64) and confirm Hash returns an
 * EXPR_BIGINT whose decimal digits equal the unsigned hash value. */
void test_hash_bigint_return() {
    const char* candidates[] = {
        "1", "2", "3", "100", "\"x\"", "\"abc\"", "foo",
        "f[x, y]", "{1, 2, 3}", "3/2", "1.5", "g[a, b, c]"
    };
    int found = 0;
    for (size_t k = 0; k < sizeof(candidates) / sizeof(candidates[0]); k++) {
        uint64_t h = hash_of(candidates[k]);
        if (h <= (uint64_t)INT64_MAX) continue;   /* fits int64; not this case */
        found = 1;

        char buf[64];
        snprintf(buf, sizeof buf, "Hash[%s]", candidates[k]);
        Expr* r = eval_str(buf);
        ASSERT_MSG(r->type == EXPR_BIGINT,
                   "expected a bigint hash for %s", candidates[k]);
        char dec[21];
        snprintf(dec, sizeof dec, "%" PRIu64, h);
        char* s = expr_to_string(r);
        ASSERT_STR_EQ(s, dec);                     /* unsigned decimal, no sign */
        free(s);
        expr_free(r);
        break;
    }
    ASSERT_MSG(found, "no candidate hashed above INT64_MAX to exercise bigint path");
}

void test_hash_unsupported_stays_unevaluated() {
    /* Unsupported hash types stay symbolic (no message). */
    check("Head[Hash[5, \"SHA256\"]]", "Hash");
    check("Head[Hash[5, \"MD5\"]]", "Hash");
    /* Unsupported formats of the Expression hash stay symbolic. */
    check("Head[Hash[5, \"Expression\", \"ByteArray\"]]", "Hash");
    check("Head[Hash[5, \"Expression\", \"Base64Encoding\"]]", "Hash");
    check("Head[Hash[5, \"Expression\", \"Nonsense\"]]", "Hash");
    /* Non-string type / format arguments stay symbolic. */
    check("Head[Hash[5, 7]]", "Hash");
    check("Head[Hash[5, \"Expression\", 7]]", "Hash");
    /* Wrong arg counts stay symbolic. */
    check("Head[Hash[]]", "Hash");
    check("Head[Hash[1, 2, 3, 4]]", "Hash");
}

void test_hash_symbol_binding() {
    /* A bound symbol hashes as its value (argument is evaluated). */
    check("x = 5; Hash[x] == Hash[5]", "True");
}

void test_hash_attributes() {
    /* Hash is Protected (and only Protected -- not Listable). */
    check("MemberQ[Attributes[Hash], Protected]", "True");
    check("MemberQ[Attributes[Hash], Listable]", "False");
    /* Not Listable: threads over the whole list, not element-wise. */
    check("Head[Hash[{1, 2, 3}]]", "Integer");
}

int main() {
    symtab_init();
    core_init();

    TEST(test_hash_basic_integer);
    TEST(test_hash_determinism_and_eval);
    TEST(test_hash_structural_sensitivity);
    TEST(test_hash_type_argument);
    TEST(test_hash_format_equivalence);
    TEST(test_hash_string_formats_roundtrip);
    TEST(test_hash_format_exact);
    TEST(test_hash_bigint_return);
    TEST(test_hash_unsupported_stays_unevaluated);
    TEST(test_hash_symbol_binding);
    TEST(test_hash_attributes);

    printf("All Hash tests passed!\n");
    symtab_clear();   /* test_hash_symbol_binding binds x */
    return 0;
}
