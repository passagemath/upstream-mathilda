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

/* ---- Phase 2 Group C: rounding / modular of integers ---- */
void test_int_c_fracpart(void) {
    assert_eval_eq("Simplify[FractionalPart[n], Element[n,Integers]]", "0", 0);
}
void test_int_c_floor(void) {
    assert_eval_eq("Simplify[Floor[n] - n, Element[n,Integers]]", "0", 0);
}
void test_int_c_ceiling(void) {
    assert_eval_eq("Simplify[Ceiling[n] - n, Element[n,Integers]]", "0", 0);
}
void test_int_c_round(void) {
    assert_eval_eq("Simplify[Round[n] - n, Element[n,Integers]]", "0", 0);
}
void test_int_c_mod1(void) {
    assert_eval_eq("Simplify[Mod[n, 1], Element[n,Integers]]", "0", 0);
}
void test_int_c_mod2n(void) {
    assert_eval_eq("Simplify[Mod[2 n, 2], Element[n,Integers]]", "0", 0);
}
void test_int_c_floorceil_half(void) {
    assert_eval_eq("Simplify[Floor[n/2] + Ceiling[n/2] - n, Element[n,Integers]]", "0", 0);
}
void test_int_c_mod_sq(void) {
    assert_eval_eq("Simplify[Mod[n^2, 2] - Mod[n, 2], Element[n,Integers]]", "0", 0);
}
void test_int_c_floor_shift(void) {
    assert_eval_eq("Simplify[Floor[n + x] - (n + Floor[x]), Element[n,Integers]]", "0", 0);
}
void test_int_c_ceil_shift(void) {
    assert_eval_eq("Simplify[Ceiling[n + x] - (n + Ceiling[x]), Element[n,Integers]]", "0", 0);
}
void test_int_c_fracpart_shift(void) {
    assert_eval_eq("Simplify[FractionalPart[n + x] - FractionalPart[x], Element[n,Integers]]", "0", 0);
}
void test_int_c_mod1_sum(void) {
    assert_eval_eq("Simplify[Mod[n + m, 1] - Mod[m, 1], Element[{n,m},Integers]]", "0", 0);
}
void test_int_c_mod2_nested(void) {
    assert_eval_eq("Simplify[Mod[n + m, 2] - Mod[Mod[n, 2] + Mod[m, 2], 2], Element[{n,m},Integers]]", "0", 0);
}
void test_int_c_mod_divis(void) {
    assert_eval_eq("Simplify[Mod[n m, m], Element[{n,m},Integers]]", "0", 0);
}

/* ---- Phase 3 Group D: Gamma / factorial / combinatorics ---- */
void test_int_d_gamma_fact(void) {
    assert_eval_eq("Simplify[Gamma[n + 1] - n!, Element[n,Integers]]", "0", 0);
}
void test_int_d_binom_symmetry(void) {
    assert_eval_eq("Simplify[Binomial[n, k] - Binomial[n, n - k], Element[{n,k},Integers]]", "0", 0);
}
void test_int_d_pascal(void) {
    assert_eval_eq("Simplify[Binomial[n + 1, k] - (Binomial[n, k] + Binomial[n, k - 1]), Element[{n,k},Integers]]", "0", 0);
}
void test_int_d_pochhammer(void) {
    assert_eval_eq("Simplify[Pochhammer[1, n] - n!, Element[n,Integers]]", "0", 0);
}
void test_int_d_gamma_ratio2(void) {
    assert_eval_eq("Simplify[Gamma[n + 2]/Gamma[n] - n (n + 1), Element[n,Integers]]", "0", 0);
}
void test_int_d_gamma_ratio1(void) {
    assert_eval_eq("Simplify[Gamma[n]/Gamma[n - 1] - (n - 1), Element[n,Integers]]", "0", 0);
}
void test_int_d_factorial2(void) {
    assert_eval_eq("Simplify[Factorial2[2 n] - 2^n n!, Element[n,Integers]]", "0", 0);
}

/* ---- Phase 4 Group E/F: orthogonal-poly value + symbolic GCD/LCM ---- */
void test_int_ef_legendre_m1(void) {
    assert_eval_eq("Simplify[LegendreP[n, -1] - (-1)^n, Element[n,Integers]]", "0", 0);
}
void test_int_ef_gcd_self(void) {
    assert_eval_eq("Simplify[GCD[n, n] - Abs[n], Element[n,Integers]]", "0", 0);
}
void test_int_ef_gcd_one(void) {
    assert_eval_eq("Simplify[GCD[n, 1] - 1, Element[n,Integers]]", "0", 0);
}
void test_int_ef_lcm_one(void) {
    assert_eval_eq("Simplify[LCM[n, 1] - Abs[n], Element[n,Integers]]", "0", 0);
}
void test_int_ef_gcd_consecutive(void) {
    assert_eval_eq("Simplify[GCD[n, n + 1] - 1, Element[n,Integers]]", "0", 0);
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
    /* ... nor its three ingredients. E^(2 I n Pi) is 1 only for integer n (an
     * even coefficient must not make 2 n vanish mod 2 before n is shown to be
     * an integer); it is not positive for real n (a positive base needs a REAL
     * exponent); and Arg[E^(I x)] is not 0 for real x (Reduce must not decide
     * E^(I x) > 0 by a real sign diagram). */
    assert_eval_eq("Simplify[E^(2 I n Pi), Element[n,Reals]] === E^(2 I n Pi)", "True", 0);
    assert_eval_eq("Simplify[E^(2 I n Pi) > 0, Element[n,Reals]] === True", "False", 0);
    assert_eval_eq("Refine[Arg[E^(I x)], Element[x,Reals]] === Arg[E^(I x)]", "True", 0);
    assert_eval_eq("Refine[Arg[I E^(I x)], Element[x,Reals]] === Arg[I E^(I x)]", "True", 0);
    /* ... while the integer case still folds. */
    assert_eval_eq("Simplify[E^(2 I n Pi), Element[n,Integers]]", "1", 0);
    assert_eval_eq("Simplify[E^(I n Pi), Element[n,Integers]]", "(-1)^n", 0);
    assert_eval_eq("Simplify[E^(2 I n Pi + x), Element[n,Integers]]", "E^x", 0);
    /* Rounding/mod rules need integrality: Floor[n] != n, Floor[n+x] shift,
     * and Mod[n,1] != 0 (it is FractionalPart[n]) for real n. */
    assert_eval_eq("Simplify[Floor[n] - n, Element[n,Reals]] === 0", "False", 0);
    assert_eval_eq("Simplify[Floor[n + x] - (n + Floor[x]), Element[x,Reals]] === 0", "False", 0);
    assert_eval_eq("Simplify[Mod[n, 1], Element[n,Reals]] === 0", "False", 0);
    assert_eval_eq("Simplify[Floor[n/2] + Ceiling[n/2] - n, Element[n,Reals]] === 0", "False", 0);
    /* Symbolic GCD/LCM rules need integrality. */
    assert_eval_eq("Simplify[GCD[n, n] - Abs[n]] === 0", "False", 0);
    assert_eval_eq("Simplify[GCD[n, n + 1] - 1] === 0", "False", 0);
}

/*
 * Documented residual (true identity, not yet collapsed):
 *   #54  Mod[n,2] + Mod[n+1,2] - 1   (Element[n,Integers])
 * Needs the parity-complement identity Mod[n+1,2] == 1 - Mod[n,2]; the
 * mod-2 argument reducer keeps both as Mod[.,2] and they do not cancel. Not
 * asserted (it would pin the current non-reduced form).
 */

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

    TEST(test_int_c_fracpart);
    TEST(test_int_c_floor);
    TEST(test_int_c_ceiling);
    TEST(test_int_c_round);
    TEST(test_int_c_mod1);
    TEST(test_int_c_mod2n);
    TEST(test_int_c_floorceil_half);
    TEST(test_int_c_mod_sq);
    TEST(test_int_c_floor_shift);
    TEST(test_int_c_ceil_shift);
    TEST(test_int_c_fracpart_shift);
    TEST(test_int_c_mod1_sum);
    TEST(test_int_c_mod2_nested);
    TEST(test_int_c_mod_divis);

    TEST(test_int_d_gamma_fact);
    TEST(test_int_d_binom_symmetry);
    TEST(test_int_d_pascal);
    TEST(test_int_d_pochhammer);
    TEST(test_int_d_gamma_ratio2);
    TEST(test_int_d_gamma_ratio1);
    TEST(test_int_d_factorial2);

    TEST(test_int_ef_legendre_m1);
    TEST(test_int_ef_gcd_self);
    TEST(test_int_ef_gcd_one);
    TEST(test_int_ef_lcm_one);
    TEST(test_int_ef_gcd_consecutive);

    TEST(test_int_soundness_declines);
    TEST(test_int_nonidentities);

    printf("All assumed-integer Simplify tests passed!\n");
    return 0;
}
