#include "test_utils.h"
#include "symtab.h"
#include "core.h"

/*
 * Assumed-INTEGER Simplify stress corpus. The integer-assumption sibling of
 * the assumed-trig (test_logexp_simplify.c / test_invtrig_simplify.c) and
 * assumed-hyperbolic (test_hyp_assume_simplify.c) corpora. Each case is a
 * (numerically verified) identity that must collapse to 0 under
 * Element[_, Integers], or a soundness decline that must NOT collapse.
 *
 * Grown phase by phase (see tasks/todo.md / the plan). This file currently
 * covers:
 *   Phase 1 -- root-of-unity powers ((-1)^P, I^P, (-I)^P, E^(I Pi X)),
 *              reduced by the Boolean-ring mod-2 engine in
 *              src/simp/simp_assume_rewrite.c (simp_reduce_mod2 /
 *              canon_neg1_pow), and Pi-periodicity with a compound integer
 *              coefficient (trig_pi_shift_symbolic de-restriction).
 *
 * Two corpus cases are NOT unconditional identities and are asserted to stay
 * non-zero (Mathilda is already correct on them):
 *   #29  Exp[I Pi(2n+1/2)] - I(-1)^n   (= 2I at n=3, not 0)
 *   #88  Sum[(-1)^k Binomial[n,k],{k,0,n}]  (= 1 at n=0; it is 0^n)
 */

/* ---- Phase 1 Group A: root-of-unity powers ---- */
void test_int_a_neg1_even(void) {
    assert_eval_eq("Simplify[(-1)^(2 n) - 1, Element[n,Integers]]", "0", 0);
}
void test_int_a_neg1_odd(void) {
    assert_eval_eq("Simplify[(-1)^(2 n + 1) + 1, Element[n,Integers]]", "0", 0);
}
void test_int_a_neg1_4n(void) {
    assert_eval_eq("Simplify[(-1)^(4 n) - 1, Element[n,Integers]]", "0", 0);
}
void test_int_a_neg1_3n(void) {
    assert_eval_eq("Simplify[(-1)^(3 n) - (-1)^n, Element[n,Integers]]", "0", 0);
}
void test_int_a_neg1_nsq(void) {
    assert_eval_eq("Simplify[(-1)^(n^2) - (-1)^n, Element[n,Integers]]", "0", 0);
}
void test_int_a_neg1_nnp1(void) {
    assert_eval_eq("Simplify[(-1)^(n (n + 1)) - 1, Element[n,Integers]]", "0", 0);
}
void test_int_a_neg1_nnm1(void) {
    assert_eval_eq("Simplify[(-1)^(n (n - 1)) - 1, Element[n,Integers]]", "0", 0);
}
void test_int_a_i_4n(void) {
    assert_eval_eq("Simplify[I^(4 n) - 1, Element[n,Integers]]", "0", 0);
}
void test_int_a_i_4np2(void) {
    assert_eval_eq("Simplify[I^(4 n + 2) + 1, Element[n,Integers]]", "0", 0);
}
void test_int_a_i_2n(void) {
    assert_eval_eq("Simplify[I^(2 n) - (-1)^n, Element[n,Integers]]", "0", 0);
}
void test_int_a_negi_2n(void) {
    assert_eval_eq("Simplify[(-I)^(2 n) - (-1)^n, Element[n,Integers]]", "0", 0);
}
void test_int_a_exp_ipin(void) {
    assert_eval_eq("Simplify[E^(I Pi n) - (-1)^n, Element[n,Integers]]", "0", 0);
}
void test_int_a_neg1_modparity(void) {
    assert_eval_eq("Simplify[(-1)^Mod[n, 2] - (-1)^n, Element[n,Integers]]", "0", 0);
}
void test_int_a_neg1_sumdiff(void) {
    assert_eval_eq("Simplify[(-1)^(n + m) - (-1)^n (-1)^m, Element[{n,m},Integers]]", "0", 0);
}
void test_int_a_neg1_minus(void) {
    assert_eval_eq("Simplify[(-1)^(n - m) - (-1)^n / (-1)^m, Element[{n,m},Integers]]", "0", 0);
}
void test_int_a_exp_ipi_nm(void) {
    assert_eval_eq("Simplify[Exp[I Pi (n - m)] - Exp[I Pi (n + m)], Element[{n,m},Integers]]", "0", 0);
}

/* ---- Phase 1 Group B: Pi-periodicity, compound integer coefficient ---- */
void test_int_b_sin_sum(void) {
    assert_eval_eq("Simplify[Sin[(n + m) Pi], Element[{n,m},Integers]]", "0", 0);
}
void test_int_b_cos_sum(void) {
    assert_eval_eq("Simplify[Cos[(n + m) Pi] - (-1)^(n + m), Element[{n,m},Integers]]", "0", 0);
}

/* ---- Phase 1: the #94 Arg correctness fix (was returning -Pi/2) ---- */
void test_int_arg_exp(void) {
    assert_eval_eq("Simplify[Arg[Exp[I Pi (2 n + 1/2)]] - Pi/2, Element[n,Integers]]", "0", 0);
}

/* ---- Soundness: these must DECLINE off the integer assumption ---- */
void test_int_soundness_declines(void) {
    /* Parity folding needs integrality: n real => (-1)^(2n) is not 1. */
    assert_eval_eq("Simplify[(-1)^(2 n) - 1, Element[n,Reals]] === 0", "False", 0);
    assert_eval_eq("Simplify[(-1)^(2 n + 1) + 1, Element[n,Reals]] === 0", "False", 0);
    assert_eval_eq("Simplify[(-1)^(3 n) - (-1)^n, Element[n,Reals]] === 0", "False", 0);
    assert_eval_eq("Simplify[I^(2 n) - (-1)^n, Element[n,Reals]] === 0", "False", 0);
    /* Pi-periodicity needs integrality. */
    assert_eval_eq("Simplify[Sin[(n + m) Pi], Element[{n,m},Reals]] === 0", "False", 0);
    /* The Arg fold must not fire for a non-integer phase. */
    assert_eval_eq("Simplify[Arg[Exp[I Pi (2 n + 1/2)]] - Pi/2, Element[n,Reals]] === 0", "False", 0);
}

/* ---- Non-identities: Mathilda is correct that these do NOT vanish ---- */
void test_int_nonidentities(void) {
    /* #29: Exp[I Pi(2n+1/2)] = I, but I(-1)^n alternates -> = 2I at odd n. */
    assert_eval_eq("Simplify[Exp[I Pi (2 n + 1/2)] - I (-1)^n, Element[n,Integers]] === 0", "False", 0);
    /* #88: Sum = 0^n, which is 1 at n=0, not unconditionally 0. */
    assert_eval_eq("Simplify[Sum[(-1)^k Binomial[n, k], {k, 0, n}], Element[n,Integers]] === 0", "False", 0);
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_int_a_neg1_even);
    TEST(test_int_a_neg1_odd);
    TEST(test_int_a_neg1_4n);
    TEST(test_int_a_neg1_3n);
    TEST(test_int_a_neg1_nsq);
    TEST(test_int_a_neg1_nnp1);
    TEST(test_int_a_neg1_nnm1);
    TEST(test_int_a_i_4n);
    TEST(test_int_a_i_4np2);
    TEST(test_int_a_i_2n);
    TEST(test_int_a_negi_2n);
    TEST(test_int_a_exp_ipin);
    TEST(test_int_a_neg1_modparity);
    TEST(test_int_a_neg1_sumdiff);
    TEST(test_int_a_neg1_minus);
    TEST(test_int_a_exp_ipi_nm);

    TEST(test_int_b_sin_sum);
    TEST(test_int_b_cos_sum);

    TEST(test_int_arg_exp);

    TEST(test_int_soundness_declines);
    TEST(test_int_nonidentities);

    printf("All assumed-integer Simplify tests passed!\n");
    return 0;
}
