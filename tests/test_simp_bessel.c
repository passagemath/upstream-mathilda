#include "test_utils.h"
#include "symtab.h"
#include "core.h"

/*
 * Simplify: Bessel functions of contiguous order (src/simp/simp_bessel.c).
 *
 * Every expectation here is DERIVED, not transcribed from another system:
 *
 *   recurrences (DLMF 10.6.1, 10.29.1)
 *     J, Y:  C[m-1] + C[m+1] =  (2m/z) C[m]
 *     I:     I[m-1] - I[m+1] =  (2m/z) I[m]
 *     K:     K[m-1] - K[m+1] = -(2m/z) K[m]
 *   derivatives (DLMF 10.6.2, 10.29.2)
 *     J, Y:  C' = C[m-1] - (m/z) C = -C[m+1] + (m/z) C
 *     I:     I' = I[m-1] - (m/z) I =  I[m+1] + (m/z) I
 *     K:     K' = -K[m-1] - (m/z) K = -K[m+1] + (m/z) K
 *   the Bessel ODEs  z^2 w'' + z w' + (z^2 -+ m^2) w = 0,
 *   d/dx[x^m C[m]] and d/dx[x^-m C[m]] (DLMF 10.6.6, 10.29.4),
 *   and Integrate/D round trips (D[Integrate[f]] == f).
 *
 * and each family is cross-checked NUMERICALLY at two sample points to 30
 * digits (bessel_numeric_zero), so a wrongly-remembered identity fails the
 * test rather than being enshrined by it.  The soundness block checks that
 * expressions with no applicable identity are left unchanged and that
 * every Simplify result equals its input numerically.
 */

/* Assert Simplify[e] == 0 exactly, and that e is numerically 0 (so the
 * identity itself, not just the simplifier, is verified). */
static void bessel_zero(const char* e) {
    char buf[4096];
    snprintf(buf, sizeof buf, "Simplify[%s]", e);
    assert_eval_eq(buf, "0", 0);
    snprintf(buf, sizeof buf,
             "Chop[N[{(%s) /. {x -> 17/10, n -> 3/7, a -> 6/5}, "
             "(%s) /. {x -> 37/10, n -> 5/9, a -> 2/3}}, 30], 10^-20]", e, e);
    assert_eval_eq(buf, "{0, 0}", 0);
}

/* Assert Simplify[e] equals `expected` (printed form) and that both agree
 * numerically with e. */
static void bessel_to(const char* e, const char* expected) {
    char buf[4096];
    snprintf(buf, sizeof buf, "Simplify[%s]", e);
    assert_eval_eq(buf, expected, 0);
    snprintf(buf, sizeof buf,
             "Chop[N[{((%s) - (%s)) /. {x -> 17/10, n -> 3/7}, "
             "((%s) - (%s)) /. {x -> 37/10, n -> 5/9}}, 30], 10^-20]",
             e, expected, e, expected);
    assert_eval_eq(buf, "{0, 0}", 0);
}

/* Soundness: Simplify[e] agrees with e numerically (form unconstrained). */
static void bessel_sound(const char* e) {
    char buf[4096];
    snprintf(buf, sizeof buf,
             "Chop[N[{(Simplify[%s] - (%s)) /. {x -> 17/10, y -> 29/10, n -> 3/7}, "
             "(Simplify[%s] - (%s)) /. {x -> 37/10, y -> 11/10, n -> 5/9}}, 30], 10^-20]",
             e, e, e, e);
    assert_eval_eq(buf, "{0, 0}", 0);
}

/* No identity applies: Simplify must leave e as it is. */
static void bessel_unchanged(const char* e) {
    char buf[4096];
    snprintf(buf, sizeof buf, "Simplify[%s] === (%s)", e, e);
    assert_eval_eq(buf, "True", 0);
}

/* ---- A. three-term recurrences ---------------------------------------- */

void test_recurrence_J(void) {
    bessel_zero("BesselJ[-2/3, x] + BesselJ[4/3, x] - 2/(3 x) BesselJ[1/3, x]");
    bessel_zero("BesselJ[1/4, x] + BesselJ[9/4, x] - 5/(2 x) BesselJ[5/4, x]");
    bessel_zero("BesselJ[-7/4, x] + BesselJ[1/4, x] + 3/(2 x) BesselJ[-3/4, x]");
    bessel_zero("BesselJ[0, x] + BesselJ[2, x] - 2/x BesselJ[1, x]");
    bessel_zero("BesselJ[2, x] + BesselJ[4, x] - 6/x BesselJ[3, x]");
    bessel_zero("BesselJ[n - 1, x] + BesselJ[n + 1, x] - (2 n)/x BesselJ[n, x]");
}

void test_recurrence_Y(void) {
    bessel_zero("BesselY[-2/3, x] + BesselY[4/3, x] - 2/(3 x) BesselY[1/3, x]");
    bessel_zero("BesselY[-3/5, x] + BesselY[7/5, x] - 4/(5 x) BesselY[2/5, x]");
    bessel_zero("BesselY[1, x] + BesselY[3, x] - 4/x BesselY[2, x]");
    bessel_zero("BesselY[n - 1, x] + BesselY[n + 1, x] - (2 n)/x BesselY[n, x]");
}

void test_recurrence_I(void) {
    bessel_zero("BesselI[-2/3, x] - BesselI[4/3, x] - 2/(3 x) BesselI[1/3, x]");
    bessel_zero("BesselI[2/5, x] - BesselI[12/5, x] - 14/(5 x) BesselI[7/5, x]");
    bessel_zero("BesselI[0, x] - BesselI[2, x] - 2/x BesselI[1, x]");
    bessel_zero("BesselI[n - 1, x] - BesselI[n + 1, x] - (2 n)/x BesselI[n, x]");
}

void test_recurrence_K(void) {
    bessel_zero("BesselK[-2/3, x] - BesselK[4/3, x] + 2/(3 x) BesselK[1/3, x]");
    bessel_zero("BesselK[-5/7, x] - BesselK[9/7, x] + 4/(7 x) BesselK[2/7, x]");
    bessel_zero("BesselK[1, x] - BesselK[3, x] + 4/x BesselK[2, x]");
    bessel_zero("BesselK[n - 1, x] - BesselK[n + 1, x] + (2 n)/x BesselK[n, x]");
}

void test_recurrence_scaled_arguments(void) {
    bessel_zero("BesselJ[-2/3, 2 x] + BesselJ[4/3, 2 x] - 2/(3 (2 x)) BesselJ[1/3, 2 x]");
    bessel_zero("BesselJ[n - 1, x^2] + BesselJ[n + 1, x^2] - (2 n)/x^2 BesselJ[n, x^2]");
    bessel_zero("BesselY[-1/3, a x] + BesselY[5/3, a x] - 4/(3 a x) BesselY[2/3, a x]");
    bessel_zero("BesselI[-1/4, 3 x] - BesselI[7/4, 3 x] - 3/(2 (3 x)) BesselI[3/4, 3 x]");
    bessel_zero("BesselK[n - 1, x/2] - BesselK[n + 1, x/2] + (2 n)/(x/2) BesselK[n, x/2]");
}

void test_recurrence_half_integer(void) {
    /* Elementary orders: still an identity (not a complete zero test there). */
    bessel_zero("BesselJ[1/2, x] + BesselJ[5/2, x] - 3/x BesselJ[3/2, x]");
    bessel_zero("BesselI[-1/2, x] - BesselI[3/2, x] - 1/x BesselI[1/2, x]");
}

void test_recurrence_two_step(void) {
    /* C[m+2] = (4m(m+1)/z^2 - 1) C[m] - (2(m+1)/z) C[m-1]  (J, Y) */
    bessel_zero("BesselJ[7/3, x] - ((4 (1/3) (4/3)/x^2 - 1) BesselJ[1/3, x]"
                " - (2 (4/3)/x) BesselJ[-2/3, x])");
    bessel_zero("BesselJ[n + 2, x] - ((4 n (n + 1)/x^2 - 1) BesselJ[n, x]"
                " - (2 (n + 1)/x) BesselJ[n - 1, x])");
    bessel_zero("BesselY[13/5, x] - ((4 (3/5) (8/5)/x^2 - 1) BesselY[3/5, x]"
                " - (2 (8/5)/x) BesselY[-2/5, x])");
}

void test_recurrence_polynomial_combinations(void) {
    bessel_zero("Expand[x (BesselJ[-2/3, x] + BesselJ[4/3, x] - 2/(3 x) BesselJ[1/3, x])]");
    bessel_zero("Expand[(x^2 + 1) (BesselK[-3/4, x] - BesselK[5/4, x] + 1/(2 x) BesselK[1/4, x])]");
    bessel_zero("Expand[(BesselJ[n - 1, x] + BesselJ[n + 1, x])^2 - 4 n^2/x^2 BesselJ[n, x]^2]");
    bessel_zero("Expand[(BesselI[-1/3, x] - BesselI[5/3, x]) BesselI[1/3, x]"
                " - 4/(3 x) BesselI[2/3, x] BesselI[1/3, x]]");
    /* two families (J and Y) in one expression */
    bessel_zero("Expand[BesselJ[1/3, x] (BesselY[2/3, x] + BesselY[8/3, x] - 10/(3 x) BesselY[5/3, x])"
                " + BesselY[2/3, x] (BesselJ[-2/3, x] + BesselJ[4/3, x] - 2/(3 x) BesselJ[1/3, x])]");
}

/* ---- B. derivative identities ----------------------------------------- */

void test_derivative_forms(void) {
    bessel_zero("D[BesselJ[1/3, x], x] - (BesselJ[-2/3, x] - (1/3)/x BesselJ[1/3, x])");
    bessel_zero("D[BesselJ[1/3, x], x] - (-BesselJ[4/3, x] + (1/3)/x BesselJ[1/3, x])");
    bessel_zero("D[BesselY[n, x], x] - (BesselY[n - 1, x] - n/x BesselY[n, x])");
    bessel_zero("D[BesselY[n, x], x] - (-BesselY[n + 1, x] + n/x BesselY[n, x])");
    bessel_zero("D[BesselI[2/5, x], x] - (BesselI[-3/5, x] - (2/5)/x BesselI[2/5, x])");
    bessel_zero("D[BesselI[n, x], x] - (BesselI[n + 1, x] + n/x BesselI[n, x])");
    bessel_zero("D[BesselK[1/3, x], x] - (-BesselK[-2/3, x] - (1/3)/x BesselK[1/3, x])");
    bessel_zero("D[BesselK[n, x], x] - (-BesselK[n + 1, x] + n/x BesselK[n, x])");
}

void test_derivative_power_weighted(void) {
    bessel_zero("D[x^(1/3) BesselJ[1/3, x], x] - x^(1/3) BesselJ[-2/3, x]");
    bessel_zero("D[x^(-1/3) BesselJ[1/3, x], x] + x^(-1/3) BesselJ[4/3, x]");
    bessel_zero("D[x^n BesselY[n, x], x] - x^n BesselY[n - 1, x]");
    bessel_zero("D[x^(2/5) BesselI[2/5, x], x] - x^(2/5) BesselI[-3/5, x]");
    bessel_zero("D[x^(-n) BesselI[n, x], x] - x^(-n) BesselI[n + 1, x]");
    bessel_zero("D[x^(1/4) BesselK[1/4, x], x] + x^(1/4) BesselK[-3/4, x]");
    bessel_zero("D[x^(-n) BesselK[n, x], x] + x^(-n) BesselK[n + 1, x]");
}

void test_bessel_odes(void) {
    bessel_zero("x^2 D[BesselJ[1/3, x], {x, 2}] + x D[BesselJ[1/3, x], x] + (x^2 - 1/9) BesselJ[1/3, x]");
    bessel_zero("x^2 D[BesselY[n, x], {x, 2}] + x D[BesselY[n, x], x] + (x^2 - n^2) BesselY[n, x]");
    bessel_zero("x^2 D[BesselI[2/3, x], {x, 2}] + x D[BesselI[2/3, x], x] - (x^2 + 4/9) BesselI[2/3, x]");
    bessel_zero("x^2 D[BesselK[n, x], {x, 2}] + x D[BesselK[n, x], x] - (x^2 + n^2) BesselK[n, x]");
}

void test_derivative_products(void) {
    bessel_zero("D[BesselJ[1/3, x]^3, x] - 3 BesselJ[1/3, x]^2 (1/(3 x) BesselJ[1/3, x] - BesselJ[4/3, x])");
    bessel_zero("D[BesselJ[n, x] BesselY[n, x], x] - (BesselJ[n, x] (n/x BesselY[n, x] - BesselY[n + 1, x])"
                " + BesselY[n, x] (BesselJ[n - 1, x] - n/x BesselJ[n, x]))");
}

/* ---- C. Integrate / D round trips ------------------------------------- */

void test_integrate_roundtrips(void) {
    /* The reported case: D[Integrate[f]] came back as a 4-order polynomial in
     * BesselJ that Simplify could not collapse (and Factor leaked indet). */
    bessel_zero("D[Integrate[x^2 BesselJ[1/3, x]^3, x], x] - x^2 BesselJ[1/3, x]^3");
    bessel_to("D[Integrate[x^2 BesselJ[1/3, x]^3, x], x]", "x^2 BesselJ[1/3, x]^3");
    bessel_zero("D[Integrate[x BesselJ[0, x]^2, x], x] - x BesselJ[0, x]^2");
    bessel_zero("D[Integrate[x BesselJ[1, x]^2, x], x] - x BesselJ[1, x]^2");
    bessel_zero("D[Integrate[x BesselI[0, x]^2, x], x] - x BesselI[0, x]^2");
    bessel_zero("D[Integrate[x^(4/3) BesselJ[1/3, x], x], x] - x^(4/3) BesselJ[1/3, x]");
    bessel_zero("D[Integrate[x^(-1/3) BesselJ[4/3, x], x], x] - x^(-1/3) BesselJ[4/3, x]");
    bessel_zero("D[Integrate[x BesselJ[1/3, x]^2, x], x] - x BesselJ[1/3, x]^2");
    bessel_zero("D[Integrate[x BesselK[1/3, x]^2, x], x] - x BesselK[1/3, x]^2");
    bessel_zero("D[Integrate[x BesselJ[n, x]^2, x], x] - x BesselJ[n, x]^2");
    bessel_zero("D[Integrate[x^3 BesselJ[0, x], x], x] - x^3 BesselJ[0, x]");
    bessel_zero("D[Integrate[x BesselY[2/5, x]^2, x], x] - x BesselY[2/5, x]^2");
    bessel_zero("D[Integrate[x BesselJ[0, x] BesselY[0, x], x], x] - x BesselJ[0, x] BesselY[0, x]");
    bessel_zero("D[Integrate[x^2 BesselI[1/3, x]^3, x], x] - x^2 BesselI[1/3, x]^3");
}

void test_integrate_roundtrip_no_messages(void) {
    /* Factor (inside Simplify) substituted x -> 0 into BesselJ[-2/3, x] and
     * printed Infinity::indet; Check must now return the value. */
    assert_eval_eq("Check[Simplify[D[Integrate[x^2 BesselJ[1/3, x]^3, x], x]], $Failed]",
                   "x^2 BesselJ[1/3, x]^3", 0);
}

/* ---- D. collapses to a non-zero form ---------------------------------- */

void test_collapse_nonzero(void) {
    bessel_to("BesselJ[n - 1, x] + BesselJ[n + 1, x]", "(2 n BesselJ[n, x])/x");
    bessel_to("BesselY[n - 1, x] + BesselY[n + 1, x]", "(2 n BesselY[n, x])/x");
    bessel_to("BesselI[n - 1, x] - BesselI[n + 1, x]", "(2 n BesselI[n, x])/x");
    bessel_to("BesselK[n + 1, x] - BesselK[n - 1, x]", "(2 n BesselK[n, x])/x");
    /* x (J[-2/3] + J[4/3]) = 2/3 J[1/3] */
    bessel_to("x BesselJ[-2/3, x] + x BesselJ[4/3, x] + BesselJ[1/3, x]", "5/3 BesselJ[1/3, x]");
}

void test_collapse_derivatives_sound(void) {
    bessel_sound("D[x^2 BesselK[1/3, x]^2 BesselI[2/3, x], x]");
    bessel_sound("D[x BesselY[2/5, x]^2 BesselJ[2/5, x], x]");
    bessel_sound("D[BesselJ[n, x]^2 + BesselJ[n + 1, x]^2, x]");
    bessel_sound("BesselJ[-2/3, x] + BesselJ[4/3, x]");
}

/* ---- E. soundness: no identity applies -------------------------------- */

void test_unchanged(void) {
    bessel_unchanged("BesselJ[1/3, x] + BesselJ[4/3, x]");       /* adjacent basis pair */
    bessel_unchanged("BesselJ[0, x] + BesselJ[2, x] + BesselJ[4, x]");
    bessel_unchanged("BesselJ[1/3, x] + BesselJ[1/3, 2 x]");     /* different arguments */
    bessel_unchanged("BesselJ[1/3, x] - BesselJ[1/3, y]");
    bessel_unchanged("BesselJ[1/3, x] + BesselJ[2/3, x]");       /* orders not congruent mod 1 */
    /* (the J/Y Casoratian is NOT unchanged: it is -2/(Pi x), see
     * test_simp_sf_bessel.c) */
    /* J[1/3] + J[7/3] = (8/(3x)) J[4/3] (recurrence at m = 4/3): one kernel
     * fewer, accepted by the kernel-count tiebreak */
    bessel_to("BesselJ[1/3, x] + BesselJ[7/3, x] + BesselY[7/3, x]",
              "(8 BesselJ[4/3, x])/(3 x) + BesselY[7/3, x]");
    bessel_unchanged("BesselJ[n, x] + x^2 BesselJ[n + 3, x]");
    bessel_unchanged("BesselJ[1/3, x]^2 + BesselJ[4/3, x]^2 + BesselJ[7/3, x]^2");
    bessel_unchanged("BesselJ[1/3, x] + BesselJ[301/3, x]");      /* wide span */
    bessel_unchanged("BesselJ[1/3, x] + BesselI[7/3, x]");        /* different kinds */
}

void test_sound_all_unchanged(void) {
    bessel_sound("BesselJ[1/3, x] BesselY[4/3, x] - BesselJ[4/3, x] BesselY[1/3, x]");
    bessel_sound("BesselJ[1/3, x]^2 + BesselJ[4/3, x]^2 + BesselJ[7/3, x]^2");
    bessel_sound("BesselJ[n, x] + x^2 BesselJ[n + 3, x]");
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_recurrence_J);
    TEST(test_recurrence_Y);
    TEST(test_recurrence_I);
    TEST(test_recurrence_K);
    TEST(test_recurrence_scaled_arguments);
    TEST(test_recurrence_half_integer);
    TEST(test_recurrence_two_step);
    TEST(test_recurrence_polynomial_combinations);
    TEST(test_derivative_forms);
    TEST(test_derivative_power_weighted);
    TEST(test_bessel_odes);
    TEST(test_derivative_products);
    TEST(test_integrate_roundtrips);
    TEST(test_integrate_roundtrip_no_messages);
    TEST(test_collapse_nonzero);
    TEST(test_collapse_derivatives_sound);
    TEST(test_unchanged);
    TEST(test_sound_all_unchanged);

    printf("All Bessel Simplify tests passed!\n");
    return 0;
}
