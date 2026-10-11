#include "test_utils.h"
#include "symtab.h"
#include "core.h"

/*
 * Assumed-hyperbolic Simplify stress corpus (v0.346) -- the hyperbolic sibling
 * of the assumed-trig corpus guarded in test_logexp_simplify.c /
 * test_invtrig_simplify.c. For each of the 50 trig cases there is a genuine
 * (numerically verified) hyperbolic counterpart that must collapse to 0 under a
 * domain assumption. Five assumption buckets, mirroring the trig corpus:
 *
 *   Group 1 (1-10)   complex components     Element[x,Reals] && Element[y,Reals]
 *   Group 2 (11-20)  principal sqrt / Abs   x > 0
 *   Group 3 (21-30)  inverse-hyperbolic     0 < x < 1/2
 *   Group 4 (31-40)  imaginary-period       Element[n,Integers] && Element[x,Reals]
 *   Group 5 (41-50)  pathological sign flip x < 0
 *
 * 48/50 collapse. Two are documented residuals (see the final block): both are
 * true identities, blocked by pre-existing general radical / inverse-hyperbolic
 * constant zero-testing gaps that are independent of the hyperbolic work.
 *
 * Mechanisms exercised (all in src/simp/, extended from trig to hyperbolic):
 *   M1 sign oracle (simp_assume.c): sign(Sinh/Tanh/Coth/Csch[u]) = sign(u),
 *      Sech>0; plus a prov_re fix for a rational coefficient (prov_re[x/2]).
 *   M4 radicand prep (simp_logexp_abs.c): 1 +/- Cosh[u] -> 2 {Cosh,Sinh}[u/2]^2.
 *   M5 inverse-of-forward (simp_assume_rewrite.c): ArcSech/ArcCsch/ArcCoth.
 *   Imaginary-period shift (simp_assume_rewrite.c): H[t + m Pi I].
 * The Group-1 complex-component rewrites and the forward-of-inverse (Group 3)
 * identities ride on pre-existing machinery.
 */

/* ---- Group 1: complex components, x,y in Reals ---- */
void test_hyp_g1_conjugate_sinh(void) {
    assert_eval_eq("Simplify[Conjugate[Sinh[x + I y]] - Sinh[x - I y], Element[x,Reals] && Element[y,Reals]]", "0", 0);
}
void test_hyp_g1_conjugate_cosh(void) {
    assert_eval_eq("Simplify[Conjugate[Cosh[x + I y]] - Cosh[x - I y], Element[x,Reals] && Element[y,Reals]]", "0", 0);
}
void test_hyp_g1_conjugate_tanh(void) {
    assert_eval_eq("Simplify[Conjugate[Tanh[x + I y]] - Tanh[x - I y], Element[x,Reals] && Element[y,Reals]]", "0", 0);
}
void test_hyp_g1_re_sinh(void) {
    assert_eval_eq("Simplify[Re[Sinh[x + I y]] - Sinh[x] Cos[y], Element[x,Reals] && Element[y,Reals]]", "0", 0);
}
void test_hyp_g1_im_sinh(void) {
    assert_eval_eq("Simplify[Im[Sinh[x + I y]] - Cosh[x] Sin[y], Element[x,Reals] && Element[y,Reals]]", "0", 0);
}
void test_hyp_g1_re_cosh(void) {
    assert_eval_eq("Simplify[Re[Cosh[x + I y]] - Cosh[x] Cos[y], Element[x,Reals] && Element[y,Reals]]", "0", 0);
}
void test_hyp_g1_im_cosh(void) {
    assert_eval_eq("Simplify[Im[Cosh[x + I y]] - Sinh[x] Sin[y], Element[x,Reals] && Element[y,Reals]]", "0", 0);
}
void test_hyp_g1_abs_sinh_sq(void) {
    assert_eval_eq("Simplify[Abs[Sinh[x + I y]]^2 - (Sinh[x]^2 + Sin[y]^2), Element[x,Reals] && Element[y,Reals]]", "0", 0);
}
void test_hyp_g1_abs_cosh_sq(void) {
    assert_eval_eq("Simplify[Abs[Cosh[x + I y]]^2 - (Sinh[x]^2 + Cos[y]^2), Element[x,Reals] && Element[y,Reals]]", "0", 0);
}
void test_hyp_g1_im_tanh_coth(void) {
    assert_eval_eq("Simplify[Im[Tanh[x]] + Im[Coth[x]], Element[x,Reals] && Element[y,Reals]]", "0", 0);
}

/* ---- Group 2: principal sqrt-of-square / Abs, x > 0 ---- */
void test_hyp_g2_sqrt_cosh2m1(void) {
    assert_eval_eq("Simplify[Sqrt[Cosh[x]^2 - 1] - Sinh[x], x > 0]", "0", 0);
}
void test_hyp_g2_sqrt_sinh2p1(void) {
    assert_eval_eq("Simplify[Sqrt[Sinh[x]^2 + 1] - Cosh[x], x > 0]", "0", 0);
}
void test_hyp_g2_sqrt_cosh2x_m1(void) {
    assert_eval_eq("Simplify[Sqrt[Cosh[2 x] - 1] - Sqrt[2] Sinh[x], x > 0]", "0", 0);
}
void test_hyp_g2_sqrt_cosh2x_p1(void) {
    assert_eval_eq("Simplify[Sqrt[Cosh[2 x] + 1] - Sqrt[2] Cosh[x], x > 0]", "0", 0);
}
void test_hyp_g2_half_sinh(void) {
    assert_eval_eq("Simplify[Sqrt[(Cosh[x] - 1)/2] - Sinh[x/2], x > 0]", "0", 0);
}
void test_hyp_g2_half_cosh(void) {
    assert_eval_eq("Simplify[Sqrt[(Cosh[x] + 1)/2] - Cosh[x/2], x > 0]", "0", 0);
}
void test_hyp_g2_sqrt_1msech2(void) {
    assert_eval_eq("Simplify[Sqrt[1 - Sech[x]^2] - Tanh[x], x > 0]", "0", 0);
}
void test_hyp_g2_sqrt_coth2m1(void) {
    assert_eval_eq("Simplify[Sqrt[Coth[x]^2 - 1] - Csch[x], x > 0]", "0", 0);
}
void test_hyp_g2_abs_sinh_cosh(void) {
    assert_eval_eq("Simplify[Abs[Sinh[x] + Cosh[x]] - (Sinh[x] + Cosh[x]), x > 0]", "0", 0);
}
void test_hyp_g2_log_cosh_sinh(void) {
    assert_eval_eq("Simplify[Log[Cosh[x] + Sinh[x]] - x, x > 0]", "0", 0);
}

/* ---- Group 3: inverse-hyperbolic compositions, 0 < x < 1/2 ---- */
void test_hyp_g3_cosh_arcsinh(void) {
    assert_eval_eq("Simplify[Cosh[ArcSinh[x]] - Sqrt[1 + x^2], 0 < x < 1/2]", "0", 0);
}
void test_hyp_g3_sinh_arctanh(void) {
    assert_eval_eq("Simplify[Sinh[ArcTanh[x]] - x/Sqrt[1 - x^2], 0 < x < 1/2]", "0", 0);
}
void test_hyp_g3_cosh_arctanh(void) {
    assert_eval_eq("Simplify[Cosh[ArcTanh[x]] - 1/Sqrt[1 - x^2], 0 < x < 1/2]", "0", 0);
}
void test_hyp_g3_tanh_arcsinh(void) {
    assert_eval_eq("Simplify[Tanh[ArcSinh[x]] - x/Sqrt[1 + x^2], 0 < x < 1/2]", "0", 0);
}
void test_hyp_g3_coth_arcsinh(void) {
    assert_eval_eq("Simplify[Coth[ArcSinh[x]] - Sqrt[1 + x^2]/x, 0 < x < 1/2]", "0", 0);
}
void test_hyp_g3_arctanh_log(void) {
    assert_eval_eq("Simplify[ArcTanh[x] - (1/2) Log[(1 + x)/(1 - x)], 0 < x < 1/2]", "0", 0);
}
void test_hyp_g3_arcsinh_double(void) {
    assert_eval_eq("Simplify[ArcSinh[2 x Sqrt[1 + x^2]] - 2 ArcSinh[x], 0 < x < 1/2]", "0", 0);
}
void test_hyp_g3_arctanh_double(void) {
    assert_eval_eq("Simplify[ArcTanh[2 x/(1 + x^2)] - 2 ArcTanh[x], 0 < x < 1/2]", "0", 0);
}

/* ---- Group 4: imaginary-period / phase, n in Integers ---- */
void test_hyp_g4_sinh_npi(void) {
    assert_eval_eq("Simplify[Sinh[n Pi I], Element[n,Integers] && Element[x,Reals]]", "0", 0);
}
void test_hyp_g4_cosh_npi(void) {
    assert_eval_eq("Simplify[Cosh[n Pi I] - (-1)^n, Element[n,Integers] && Element[x,Reals]]", "0", 0);
}
void test_hyp_g4_tanh_npi(void) {
    assert_eval_eq("Simplify[Tanh[n Pi I], Element[n,Integers] && Element[x,Reals]]", "0", 0);
}
void test_hyp_g4_sinh_shift(void) {
    assert_eval_eq("Simplify[Sinh[x + 2 n Pi I] - Sinh[x], Element[n,Integers] && Element[x,Reals]]", "0", 0);
}
void test_hyp_g4_cosh_shift(void) {
    assert_eval_eq("Simplify[Cosh[x + 2 n Pi I] - Cosh[x], Element[n,Integers] && Element[x,Reals]]", "0", 0);
}
void test_hyp_g4_tanh_shift(void) {
    assert_eval_eq("Simplify[Tanh[x + n Pi I] - Tanh[x], Element[n,Integers] && Element[x,Reals]]", "0", 0);
}
void test_hyp_g4_coth_shift(void) {
    assert_eval_eq("Simplify[Coth[x + n Pi I] - Coth[x], Element[n,Integers] && Element[x,Reals]]", "0", 0);
}
void test_hyp_g4_sinh_halfshift(void) {
    assert_eval_eq("Simplify[Sinh[x + (2 n + 1) Pi I] + Sinh[x], Element[n,Integers] && Element[x,Reals]]", "0", 0);
}
void test_hyp_g4_cosh_halfshift(void) {
    assert_eval_eq("Simplify[Cosh[x + (2 n + 1) Pi I] + Cosh[x], Element[n,Integers] && Element[x,Reals]]", "0", 0);
}
void test_hyp_g4_exp_euler(void) {
    assert_eval_eq("Simplify[Exp[I n Pi] - Cosh[n Pi I], Element[n,Integers] && Element[x,Reals]]", "0", 0);
}

/* ---- Group 5: pathological sign flips, x < 0 ---- */
void test_hyp_g5_sqrt_sinh2(void) {
    assert_eval_eq("Simplify[Sqrt[Sinh[x]^2] + Sinh[x], x < 0]", "0", 0);
}
void test_hyp_g5_sqrt_cosh2m1(void) {
    assert_eval_eq("Simplify[Sqrt[Cosh[x]^2 - 1] + Sinh[x], x < 0]", "0", 0);
}
void test_hyp_g5_abs_sinh(void) {
    assert_eval_eq("Simplify[Abs[Sinh[x]] + Sinh[x], x < 0]", "0", 0);
}
void test_hyp_g5_abs_tanh(void) {
    assert_eval_eq("Simplify[Abs[Tanh[x]] + Tanh[x], x < 0]", "0", 0);
}
void test_hyp_g5_sqrt_1msech2(void) {
    assert_eval_eq("Simplify[Sqrt[1 - Sech[x]^2] + Tanh[x], x < 0]", "0", 0);
}
void test_hyp_g5_sqrt_cosh2x_p1(void) {
    assert_eval_eq("Simplify[Sqrt[Cosh[2 x] + 1] - Sqrt[2] Cosh[x], x < 0]", "0", 0);
}
void test_hyp_g5_sqrt_cosh2x_m1(void) {
    assert_eval_eq("Simplify[Sqrt[Cosh[2 x] - 1] + Sqrt[2] Sinh[x], x < 0]", "0", 0);
}
void test_hyp_g5_arccosh_cosh(void) {
    assert_eval_eq("Simplify[ArcCosh[Cosh[x]] + x, x < 0]", "0", 0);
}
void test_hyp_g5_arcsinh_sinh(void) {
    assert_eval_eq("Simplify[ArcSinh[Sinh[x]] - x, x < 0]", "0", 0);
}
void test_hyp_g5_arcsech_sech(void) {
    assert_eval_eq("Simplify[ArcSech[Sech[x]] + x, x < 0]", "0", 0);
}

/* ---- Soundness: each new mechanism must DECLINE off-identity ---- */
void test_hyp_soundness_declines(void) {
    /* Sign unknown => Sqrt[Sinh^2] stays Abs, must not assume Sinh>=0. */
    assert_eval_eq("Simplify[Sqrt[Sinh[x]^2] - Sinh[x], Element[x,Reals]] === 0", "False", 0);
    /* Wrong region: on x<0, Sqrt[Cosh^2-1] = -Sinh, so the x>0 identity fails. */
    assert_eval_eq("Simplify[Sqrt[Cosh[x]^2 - 1] - Sinh[x], x < 0] === 0", "False", 0);
    /* No integrality: the imaginary-period shift must not fire for free n. */
    assert_eval_eq("Simplify[Sinh[x + 2 n Pi I] - Sinh[x], Element[x,Reals]] === 0", "False", 0);
    /* A REAL multiple of Pi is NOT a hyperbolic period: must not reduce. */
    assert_eval_eq("Simplify[Sinh[x + 2 Pi] - Sinh[x], Element[x,Reals] && Element[n,Integers]] === 0", "False", 0);
    /* ArcSech[Sech[x]] = |x|; with the sign of x unknown it must not collapse. */
    assert_eval_eq("Simplify[ArcSech[Sech[x]] + x, Element[x,Reals]] === 0", "False", 0);
}

/*
 * Documented residuals (true identities, not yet collapsed by Simplify):
 *   #26  ArcSinh[x/Sqrt[1-x^2]] - ArcTanh[x]        (0 < x < 1/2)
 *   #29  ArcCosh[1 + 2 x^2]     - 2 ArcSinh[x]      (0 < x < 1/2)
 * Both are verified identities (the corpus's numeric pass). They reach the
 * inverse-function region prover (simp_region_identity_is_zero), where:
 *   - #26 needs a radical-product / rational-in-Sqrt normalization of the
 *     gradient that Simplify lacks (Sqrt[1-x^2] Sqrt[1+x^2/(1-x^2)] -> 1);
 *   - #29 needs a k*ArcSinh[rational] + Sum Log[algebraic] = 0 constant
 *     zero-test at the sample point.
 * Neither gap is hyperbolic-specific (e.g. Simplify[Sqrt[x^2+x^4] - x Sqrt[1+x^2],
 * x>0] likewise does not vanish). They are left for dedicated radical /
 * inverse-function zero-testing work, as the trig corpus left its one-off #20.
 * A guard here would only pin the current non-reduced printed form, which is
 * brittle, so they are intentionally not asserted.
 */

int main(void) {
    symtab_init();
    core_init();

    /* Group 1 */
    TEST(test_hyp_g1_conjugate_sinh);
    TEST(test_hyp_g1_conjugate_cosh);
    TEST(test_hyp_g1_conjugate_tanh);
    TEST(test_hyp_g1_re_sinh);
    TEST(test_hyp_g1_im_sinh);
    TEST(test_hyp_g1_re_cosh);
    TEST(test_hyp_g1_im_cosh);
    TEST(test_hyp_g1_abs_sinh_sq);
    TEST(test_hyp_g1_abs_cosh_sq);
    TEST(test_hyp_g1_im_tanh_coth);

    /* Group 2 */
    TEST(test_hyp_g2_sqrt_cosh2m1);
    TEST(test_hyp_g2_sqrt_sinh2p1);
    TEST(test_hyp_g2_sqrt_cosh2x_m1);
    TEST(test_hyp_g2_sqrt_cosh2x_p1);
    TEST(test_hyp_g2_half_sinh);
    TEST(test_hyp_g2_half_cosh);
    TEST(test_hyp_g2_sqrt_1msech2);
    TEST(test_hyp_g2_sqrt_coth2m1);
    TEST(test_hyp_g2_abs_sinh_cosh);
    TEST(test_hyp_g2_log_cosh_sinh);

    /* Group 3 (26 and 29 are documented residuals) */
    TEST(test_hyp_g3_cosh_arcsinh);
    TEST(test_hyp_g3_sinh_arctanh);
    TEST(test_hyp_g3_cosh_arctanh);
    TEST(test_hyp_g3_tanh_arcsinh);
    TEST(test_hyp_g3_coth_arcsinh);
    TEST(test_hyp_g3_arctanh_log);
    TEST(test_hyp_g3_arcsinh_double);
    TEST(test_hyp_g3_arctanh_double);

    /* Group 4 */
    TEST(test_hyp_g4_sinh_npi);
    TEST(test_hyp_g4_cosh_npi);
    TEST(test_hyp_g4_tanh_npi);
    TEST(test_hyp_g4_sinh_shift);
    TEST(test_hyp_g4_cosh_shift);
    TEST(test_hyp_g4_tanh_shift);
    TEST(test_hyp_g4_coth_shift);
    TEST(test_hyp_g4_sinh_halfshift);
    TEST(test_hyp_g4_cosh_halfshift);
    TEST(test_hyp_g4_exp_euler);

    /* Group 5 */
    TEST(test_hyp_g5_sqrt_sinh2);
    TEST(test_hyp_g5_sqrt_cosh2m1);
    TEST(test_hyp_g5_abs_sinh);
    TEST(test_hyp_g5_abs_tanh);
    TEST(test_hyp_g5_sqrt_1msech2);
    TEST(test_hyp_g5_sqrt_cosh2x_p1);
    TEST(test_hyp_g5_sqrt_cosh2x_m1);
    TEST(test_hyp_g5_arccosh_cosh);
    TEST(test_hyp_g5_arcsinh_sinh);
    TEST(test_hyp_g5_arcsech_sech);

    /* Soundness */
    TEST(test_hyp_soundness_declines);

    printf("All assumed-hyperbolic Simplify tests passed!\n");
    return 0;
}
