/* Tests for DirichletCharacter[k, j, n].
 *
 * The labelling must be Mathematica's: every table below was produced by
 * Mathematica 14 (wolframscript), including the non-cyclic 2-adic moduli whose
 * digit order (the (-1) component before the <5> component) the spec does not
 * state. Covers prime and prime-power moduli, composite moduli (mixed-radix
 * index), real characters, the multiplicative / periodic structure, large
 * moduli (20!), the exact normal form of the values, list threading, the packed
 * int64 surface, the messages (intp / invl / int / argrx) and the attributes. */

#include "eval.h"
#include "parse.h"
#include "expr.h"
#include "symtab.h"
#include "attr.h"
#include "core.h"
#include "print.h"
#include "test_utils.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static void T(const char* in) { assert_eval_eq(in, "True", 0); }

void test_dc_basic() {
    assert_eval_eq("Table[DirichletCharacter[2, 1, n], {n, 10}]",
                   "{1, 0, 1, 0, 1, 0, 1, 0, 1, 0}", 0);
    assert_eval_eq("Table[DirichletCharacter[1, 1, n], {n, -5, 5}]",
                   "{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}", 0);
    assert_eval_eq("Table[DirichletCharacter[5, 1, n], {n, EulerPhi[5]}]", "{1, 1, 1, 1}", 0);
    assert_eval_eq("Table[DirichletCharacter[9, 1, n], {n, EulerPhi[9]}]",
                   "{1, 1, 0, 1, 1, 0}", 0);
}

void test_dc_mod7_table() {
    /* All six characters mod 7, n = 0..6 (Mathematica's grid). */
    T("Table[DirichletCharacter[7, j, n], {j, 1, 6}, {n, 0, 6}] == "
      "{{0, 1, 1, 1, 1, 1, 1},"
      " {0, 1, E^(2 I Pi/3), E^(I Pi/3), E^(-2 I Pi/3), E^(-I Pi/3), -1},"
      " {0, 1, E^(-2 I Pi/3), E^(2 I Pi/3), E^(2 I Pi/3), E^(-2 I Pi/3), 1},"
      " {0, 1, 1, -1, 1, -1, -1},"
      " {0, 1, E^(2 I Pi/3), E^(-2 I Pi/3), E^(-2 I Pi/3), E^(2 I Pi/3), 1},"
      " {0, 1, E^(-2 I Pi/3), E^(-I Pi/3), E^(2 I Pi/3), E^(I Pi/3), -1}}");
    /* At the primitive root 3: exp(2 Pi I (j-1)/6). */
    T("Table[DirichletCharacter[7, j, 3], {j, 6}] == "
      "{1, E^(I Pi/3), E^(2 I Pi/3), -1, E^(-2 I Pi/3), E^(-I Pi/3)}");
    T("Table[DirichletCharacter[7, 2, 3^n] == DirichletCharacter[7, 2, 3]^n, {n, 6}] == "
      "{True, True, True, True, True, True}");
}

void test_dc_real_characters() {
    assert_eval_eq("Table[DirichletCharacter[k, 1, n], {k, 3, 10}, {n, 0, k - 1}]",
                   "{{0, 1, 1}, {0, 1, 0, 1}, {0, 1, 1, 1, 1}, {0, 1, 0, 0, 0, 1}, "
                   "{0, 1, 1, 1, 1, 1, 1}, {0, 1, 0, 1, 0, 1, 0, 1}, "
                   "{0, 1, 1, 0, 1, 1, 0, 1, 1}, {0, 1, 0, 1, 0, 0, 0, 1, 0, 1}}", 0);
    assert_eval_eq("Table[DirichletCharacter[k, EulerPhi[k]/2 + 1, n], {k, 3, 10}, {n, 0, k - 1}]",
                   "{{0, 1, -1}, {0, 1, 0, -1}, {0, 1, -1, -1, 1}, {0, 1, 0, 0, 0, -1}, "
                   "{0, 1, 1, -1, 1, -1, -1}, {0, 1, 0, -1, 0, 1, 0, -1}, "
                   "{0, 1, -1, 0, 1, -1, 0, 1, -1}, {0, 1, 0, -1, 0, 0, 0, -1, 0, 1}}", 0);
    /* JacobiSymbol agreements from the spec. */
    T("Table[DirichletCharacter[5, 3, n], {n, 0, 4}] == Table[JacobiSymbol[n, 5], {n, 0, 4}]");
    assert_eval_eq("Table[DirichletCharacter[7, 4, n], {n, 0, 6}]", "{0, 1, 1, -1, 1, -1, -1}", 0);
    assert_eval_eq("DirichletCharacter[12, 3, {1, 5, 7, 11}]", "{1, 1, -1, -1}", 0);
}

void test_dc_two_adic() {
    /* Mathematica's order mod 8 and mod 16: the (-1) digit is most significant. */
    assert_eval_eq("Table[DirichletCharacter[8, j, n], {j, 4}, {n, 0, 7}]",
                   "{{0, 1, 0, 1, 0, 1, 0, 1}, {0, 1, 0, -1, 0, -1, 0, 1}, "
                   "{0, 1, 0, -1, 0, 1, 0, -1}, {0, 1, 0, 1, 0, -1, 0, -1}}", 0);
    T("Table[DirichletCharacter[16, j, n], {j, 8}, {n, 1, 15, 2}] == "
      "{{1, 1, 1, 1, 1, 1, 1, 1}, {1, -I, I, -1, -1, I, -I, 1}, "
      "{1, -1, -1, 1, 1, -1, -1, 1}, {1, I, -I, -1, -1, -I, I, 1}, "
      "{1, -1, 1, -1, 1, -1, 1, -1}, {1, I, I, 1, -1, -I, -I, -1}, "
      "{1, 1, -1, -1, 1, 1, -1, -1}, {1, -I, -I, 1, -1, I, I, -1}}");
    /* mod 32, j = 10: a unit times an exponential, uncombined (-E^(3 I Pi/4)). */
    T("Table[DirichletCharacter[32, 10, n], {n, 1, 7, 2}] == "
      "{1, -E^(3 I Pi/4), E^(I Pi/4), -I}");
    /* mod 24 and mod 40: the 2-adic part leads the mixed radix. */
    assert_eval_eq("Table[DirichletCharacter[24, j, 5], {j, 8}]",
                   "{1, -1, -1, 1, 1, -1, -1, 1}", 0);
    assert_eval_eq("Table[DirichletCharacter[40, j, 3], {j, 4}]", "{1, -I, -1, I}", 0);
}

void test_dc_composite() {
    T("DirichletCharacter[3^2 5, 8, 7] == -I E^(-2 I Pi/3)");
    T("DirichletCharacter[3^2 5, 8, 7] == DirichletCharacter[3^2, 2, 7] DirichletCharacter[5, 4, 2]");
    /* Lift of the character mod 3 to 3^3 and 3^4 (index (j-1) phi(mu)/phi(k) + 1). */
    assert_eval_eq("Table[DirichletCharacter[3^3, 10, n], {n, {1, 2, 4, 7, 8, 11, 13, 14}}]",
                   "{1, -1, 1, 1, -1, -1, 1, -1}", 0);
    assert_eval_eq("Table[DirichletCharacter[3^4, 28, n], {n, {1, 2, 4, 7, 8, 11, 13, 14}}]",
                   "{1, -1, 1, 1, -1, -1, 1, -1}", 0);
    /* The full decomposition check over the units mod 45. */
    T("With[{c = Select[Range[45], CoprimeQ[#, 45] &]}, "
      "Table[DirichletCharacter[45, 8, n], {n, c}] == "
      "Table[DirichletCharacter[9, 2, Mod[n, 9]] DirichletCharacter[5, 4, Mod[n, 5]], {n, c}]]");
}

void test_dc_structure() {
    T("Table[DirichletCharacter[7, 3, n] == DirichletCharacter[7, 3, 7 + n], {n, 7}] == "
      "{True, True, True, True, True, True, True}");
    T("DirichletCharacter[7, 3, 4 5] == DirichletCharacter[7, 3, 4] DirichletCharacter[7, 3, 5]");
    assert_eval_eq("Table[DirichletCharacter[7, 3, 7 n], {n, 6}]", "{0, 0, 0, 0, 0, 0}", 0);
    T("Chop[Abs[N[Table[DirichletCharacter[7, 5, n], {n, 0, 4}]]] - {0, 1, 1, 1, 1}] == {0, 0, 0, 0, 0}");
    /* Negative and big n reduce mod k. */
    T("DirichletCharacter[7, 2, -3] == E^(-2 I Pi/3)");
    T("DirichletCharacter[7, 2, 10^30 + 3] == DirichletCharacter[7, 2, Mod[10^30 + 3, 7]]");
    /* Completely multiplicative over a composite modulus. */
    T("And @@ Flatten[Table[DirichletCharacter[60, 11, a b] == "
      "DirichletCharacter[60, 11, a] DirichletCharacter[60, 11, b], {a, 1, 30}, {b, 1, 30}]]");
}

void test_dc_large_modulus() {
    T("DirichletCharacter[20!, 300, 23] == E^(I Pi/9)");
    T("DirichletCharacter[20!, 1, 23] == 1");
    assert_eval_eq("DirichletCharacter[20!, 300, 22]", "0", 0);
    /* A prime modulus near 2^40: Pohlig-Hellman discrete logarithms. */
    T("With[{p = NextPrime[2^40]}, DirichletCharacter[p, 2, PrimitiveRoot[p]] == E^(2 Pi I/(p - 1))]");
    T("With[{p = NextPrime[2^40]}, DirichletCharacter[p, 3, PrimitiveRoot[p]^5] == E^(20 Pi I/(p - 1))]");
}

void test_dc_threading_and_packed() {
    assert_eval_eq("DirichletCharacter[3, 2, {1, 2, 3, 4, 5}]", "{1, -1, 0, 1, -1}", 0);
    assert_eval_eq("DirichletCharacter[{3, 4}, 2, {1, 2}]", "{1, 0}", 0);
    /* A packed int64 n with a real character answers with a packed array. */
    assert_eval_eq("PackedArrayQ[DirichletCharacter[12, 3, Range[-5, 30]]]", "True", 0);
    T("DirichletCharacter[12, 3, Range[-5, 30]] == Table[DirichletCharacter[12, 3, n], {n, -5, 30}]");
    /* A complex character answers with exact roots of unity. */
    T("DirichletCharacter[7, 2, Range[0, 20]] == Table[DirichletCharacter[7, 2, n], {n, 0, 20}]");
    T("Normal[DirichletCharacter[8, 3, NDArray[Range[0, 15], DataType -> \"int64\"]]] == Table[DirichletCharacter[8, 3, n], {n, 0, 15}]");
}

void test_dc_messages() {
    assert_eval_eq("DirichletCharacter[7, 9, 2]", "DirichletCharacter[7, 9, 2]", 0);
    assert_eval_eq("DirichletCharacter[7, 0, 2]", "DirichletCharacter[7, 0, 2]", 0);
    assert_eval_eq("DirichletCharacter[7, 2, 2.5]", "DirichletCharacter[7, 2, 2.5]", 0);
    assert_eval_eq("DirichletCharacter[0, 1, 2]", "DirichletCharacter[0, 1, 2]", 0);
    assert_eval_eq("DirichletCharacter[-7, 2, 3]", "DirichletCharacter[-7, 2, 3]", 0);
    assert_eval_eq("DirichletCharacter[k, j, n]", "DirichletCharacter[k, j, n]", 0);
    assert_eval_eq("DirichletCharacter[]", "DirichletCharacter[]", 0);
    /* Messages route through the funnel: Check sees them. */
    assert_eval_eq("Check[DirichletCharacter[7, 9, 2], failed]", "failed", 0);
    assert_eval_eq("Check[DirichletCharacter[7, 2, 2.5], failed]", "failed", 0);
    assert_eval_eq("Check[DirichletCharacter[], failed]", "failed", 0);
    assert_eval_eq("Quiet[DirichletCharacter[0, 1, 2]]", "DirichletCharacter[0, 1, 2]", 0);
}

void test_dc_attributes() {
    assert_eval_eq("Attributes[DirichletCharacter]", "{Protected}", 0);
}

int main() {
    symtab_init();
    core_init();

    TEST(test_dc_basic);
    TEST(test_dc_mod7_table);
    TEST(test_dc_real_characters);
    TEST(test_dc_two_adic);
    TEST(test_dc_composite);
    TEST(test_dc_structure);
    TEST(test_dc_large_modulus);
    TEST(test_dc_threading_and_packed);
    TEST(test_dc_messages);
    TEST(test_dc_attributes);

    printf("All DirichletCharacter tests passed.\n");
    return 0;
}
