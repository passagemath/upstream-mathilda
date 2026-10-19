/* test_partitions_review.c -- the number-theory review, batch 6 (v0.474):
 * IntegerPartitions (its four argument messages, a negative length step, a
 * number as the part spec, the None / Take-sequence forms of m with ::take
 * clamping, bignum n, unsorted part lists through the int64 count engine),
 * and PartitionsP / PartitionsQ (the startup tables, the per-term-precision
 * q series, and agreement across scalar, packed, NDArray and Compile).
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
        printf("partitions review test failed: %s\n  expected: %s\n  got:      %s\n",
               input, expected, res_str);
        ASSERT(0);
    }
    free(res_str);
    expr_free(e);
    expr_free(res);
}

/* ---- IntegerPartitions: argument messages ---- */
static void test_ip_messages(void) {
    check("Check[IntegerPartitions[x], m, IntegerPartitions::rat]", "m");
    check("Check[IntegerPartitions[2.5], m, IntegerPartitions::rat]", "m");
    check("Check[IntegerPartitions[{5, 6}], m, IntegerPartitions::rat]", "m");
    check("Check[IntegerPartitions[5, -1], m, IntegerPartitions::nninfseq]", "m");
    check("Check[IntegerPartitions[5, {1, 2, 0}], m, IntegerPartitions::nninfseq]", "m");
    check("Check[IntegerPartitions[8, {Infinity}], m, IntegerPartitions::nninfseq]", "m");
    check("Check[IntegerPartitions[5, 2.], m, IntegerPartitions::nninfseq]", "m");
    check("Check[IntegerPartitions[6, x, y], m, IntegerPartitions::nninfseq]", "m");
    check("Check[IntegerPartitions[5, All, x], m, IntegerPartitions::rseq]", "m");
    check("Check[IntegerPartitions[6, All, {1, 2.}], m, IntegerPartitions::rseq]", "m");
    check("Check[IntegerPartitions[5, All, All, x], m, IntegerPartitions::seq]", "m");
    check("Check[IntegerPartitions[6, All, All, {1, 5, 0}], m, IntegerPartitions::seq]", "m");
    check("Check[IntegerPartitions[5, All, {0, 1}], m, IntegerPartitions::undef]", "m");
    check("Check[IntegerPartitions[], m, IntegerPartitions::argb]", "m");
}

/* ---- IntegerPartitions: length, part and m specifications ---- */
static void test_ip_specs(void) {
    check("IntegerPartitions[8, {6, 2, -2}]",
          "{{7, 1}, {6, 2}, {5, 3}, {5, 1, 1, 1}, {4, 4}, {4, 2, 1, 1}, {3, 3, 1, 1}, "
          "{3, 2, 2, 1}, {3, 1, 1, 1, 1, 1}, {2, 2, 2, 2}, {2, 2, 1, 1, 1, 1}}");
    check("IntegerPartitions[8, {2, 4, -1}]", "{}");
    check("IntegerPartitions[8, {2, Infinity, -1}]", "{}");
    check("IntegerPartitions[8, {3, 3, -1}]", "{{6, 1, 1}, {5, 2, 1}, {4, 3, 1}, {4, 2, 2}, {3, 3, 2}}");
    check("IntegerPartitions[8, {2, Infinity, 3}]",
          "{{7, 1}, {6, 2}, {5, 3}, {4, 4}, {4, 1, 1, 1, 1}, {3, 2, 1, 1, 1}, {2, 2, 2, 1, 1}, "
          "{1, 1, 1, 1, 1, 1, 1, 1}}");
    check("IntegerPartitions[6, All, 2]", "{{2, 2, 2}, {2, 2, 1, 1}, {2, 1, 1, 1, 1}, {1, 1, 1, 1, 1, 1}}");
    check("IntegerPartitions[6, All, 5/2]", "{{2, 2, 2}, {2, 2, 1, 1}, {2, 1, 1, 1, 1}, {1, 1, 1, 1, 1, 1}}");
    check("IntegerPartitions[6, All, -3]", "{}");
    check("IntegerPartitions[6, All, {}]", "{}");
    check("IntegerPartitions[0, All, {}]", "{{}}");
    check("IntegerPartitions[6, 10^30] === IntegerPartitions[6]", "True");
    check("IntegerPartitions[6, All, All, None]", "{}");
    check("IntegerPartitions[6, All, All, {3}]", "{{4, 2}}");
    check("IntegerPartitions[6, All, All, {2, 4}]", "{{5, 1}, {4, 2}, {4, 1, 1}}");
    check("IntegerPartitions[6, All, All, {1, -1, 3}]", "{{6}, {4, 1, 1}, {3, 1, 1, 1}, {2, 1, 1, 1, 1}}");
    check("IntegerPartitions[6, All, All, {-2}]", "{{2, 1, 1, 1, 1}}");
    check("IntegerPartitions[6, All, All, {-1, 1, -1}] === Reverse[IntegerPartitions[6]]", "True");
    check("IntegerPartitions[6, All, All, {4, 2}]", "{}");
    check("IntegerPartitions[6, All, All, {0}]", "{}");
    /* a sequence past the end warns and answers what is there */
    check("Quiet[IntegerPartitions[6, All, All, {20}]]", "{}");
    check("Quiet[IntegerPartitions[6, All, All, {2, 30}]] === Rest[IntegerPartitions[6]]", "True");
    check("Quiet[IntegerPartitions[6, All, All, -20]] === IntegerPartitions[6]", "True");
    check("Check[IntegerPartitions[6, All, All, {2, 30}], m, IntegerPartitions::take]", "m");
    /* bignum n with a length bound */
    check("IntegerPartitions[10^20, 1]", "{{100000000000000000000}}");
    check("IntegerPartitions[10^20, {3}, All, 2]",
          "{{99999999999999999998, 1, 1}, {99999999999999999997, 2, 1}}");
    check("IntegerPartitions[10^20, All, All, 3]",
          "{{100000000000000000000}, {99999999999999999999, 1}, {99999999999999999998, 2}}");
    check("IntegerPartitions[10^20, {2}, {1, 10^20 - 1}]", "{{99999999999999999999, 1}}");
    /* a huge n with a small length bound needs no part table */
    check("Length[IntegerPartitions[10^7, {2}]]", "5000000");
}

/* ---- IntegerPartitions: unsorted and repeated part lists ---- */
static void test_ip_unsorted(void) {
    check("IntegerPartitions[6, All, {3, 1, 2}]",
          "{{2, 2, 2}, {2, 2, 1, 1}, {2, 1, 1, 1, 1}, {2, 1, 3}, {1, 1, 1, 1, 1, 1}, {1, 1, 1, 3}, {3, 3}}");
    check("IntegerPartitions[6, All, {2, 2, 1}]",
          "{{1, 1, 1, 1, 1, 1}, {1, 1, 1, 1, 2}, {1, 1, 1, 1, 2}, {1, 1, 2, 2}, {1, 1, 2, 2}, "
          "{1, 1, 2, 2}, {2, 2, 2}, {2, 2, 2}, {2, 2, 2}, {2, 2, 2}}");
    check("IntegerPartitions[10, All, {7, 2, 5, 2}, {2, 4}]",
          "{{2, 2, 2, 2, 2}, {2, 2, 2, 2, 2}, {2, 2, 2, 2, 2}}");
    check("Length[IntegerPartitions[60, All, {11, 7, 5, 3, 2, 1}]]", "8021");
    check("Total[Length /@ IntegerPartitions[60, All, {11, 7, 5, 3, 2, 1}]]", "199462");
    check("IntegerPartitions[5/2, All, {1/2, 1}]", "{{1, 1, 1/2}, {1, 1/2, 1/2, 1/2}, {1/2, 1/2, 1/2, 1/2, 1/2}}");
}

/* ---- PartitionsP / PartitionsQ ---- */
static void test_partitions_counts(void) {
    check("PartitionsP[{-3, 0, 10}]", "{0, 1, 42}");
    check("PartitionsP[405]", "9147679068859117602");
    check("PartitionsP[406]", "9725512513742021729");
    check("PartitionsQ[769]", "9031248611609428978");
    check("PartitionsQ[770]", "9322334643320220726");
    check("PartitionsQ[999]", "8397642283676830881826");
    check("PartitionsQ[1000]", "8635565795744155161506");
    check("Mod[PartitionsP[10^5], 10^9 + 7]", "49037875");
    check("Mod[PartitionsQ[10^5], 10^9 + 7]", "119726462");
    check("PartitionsQ[{1, 2, 3}]", "{1, 1, 2}");
    check("PartitionsP[2.5]", "PartitionsP[2.5]");
    check("Check[PartitionsP[], m, PartitionsP::argx]", "m");
    check("Check[PartitionsQ[1, 2], m, PartitionsQ::argx]", "m");
    /* every surface answers from the same int64 table */
    check("With[{r = Range[-5, 405]}, {PartitionsP[r] === Map[PartitionsP, r], "
          "PackedArrayQ[PartitionsP[r]], "
          "Normal[PartitionsP[NDArray[r, DataType -> \"int64\"]]] === PartitionsP[r], "
          "Compile[{{v, _Integer, 1}}, PartitionsP[v]][r] === PartitionsP[r]}]",
          "{True, True, True, True}");
    check("With[{r = Range[-5, 769]}, {PartitionsQ[r] === Map[PartitionsQ, r], "
          "Compile[{{x, _Integer}}, PartitionsQ[x]] /@ r === PartitionsQ[r]}]",
          "{True, True}");
    /* past the int64 range the packed answer falls back to bignums */
    check("PartitionsP[Range[400, 410]] === Map[PartitionsP, Range[400, 410]]", "True");
    check("\"Compiled\" /. CompileDiagnostics[{{v, _Integer, 1}}, PartitionsQ[v]]", "True");
}

int main() {
    symtab_init();
    core_init();

    TEST(test_ip_messages);
    TEST(test_ip_specs);
    TEST(test_ip_unsorted);
    TEST(test_partitions_counts);

    printf("All partitions review tests passed!\n");
    return 0;
}
