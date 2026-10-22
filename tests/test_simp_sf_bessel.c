#include "simp_sf_test_util.h"
#include "symtab.h"
#include "core.h"

/*
 * Simplify, special-function engine, Bessel module (src/simp/simp_bessel.c on
 * src/simp/simp_sf_engine.c): the relations ADDED by the engine port --
 * Wronskian/Casoratian pairs, BesselK order evenness, integer-order
 * reflection under assumptions.  The contiguous-order recurrences themselves
 * are covered by tests/test_simp_bessel.c.
 *
 * Derivations (docs/design/simp_special_functions.md, section Bessel):
 *
 *   W{J_nu, Y_nu}(z) = 2/(Pi z)  (Abel: W ~ 1/z for z^2 w'' + z w' + ...;
 *   constant from the large-z asymptotics), and J' = J[nu-1] - (nu/z) J
 *   turns it into the Casoratian
 *       J[nu+1] Y[nu] - J[nu] Y[nu+1] = 2/(Pi z).
 *   Likewise W{I_nu, K_nu} = -1/z with I' = I[nu+1] + (nu/z) I,
 *   K' = -K[nu+1] + (nu/z) K gives
 *       I[nu+1] K[nu] + I[nu] K[nu+1] = 1/z.
 *   The k-step Casoratian: with C[n+k] = A_k C[n] + B_k C[n+1] (shared by J
 *   and Y), J[n+k] Y[n] - J[n] Y[n+k] = B_k (J[n+1] Y[n] - J[n] Y[n+1]),
 *   B_2 = 2(n+1)/z, B_3 = 4(n+1)(n+2)/z^2 - 1.
 *
 *   K_nu(z) = Int_0^Inf exp(-z cosh t) cosh(nu t) dt is even in nu.
 *   For integer n the generating function exp(z/2 (t - 1/t)) = Sum J_n t^n
 *   is invariant under t -> -1/t, so J_{-n} = (-1)^n J_n; Y_{-n} = (-1)^n
 *   Y_n follows from Y's limit definition; I_{-n} = I_n likewise.
 */

static const char* INT_PTS[] = {
    "{x -> 17/10, y -> 29/10, n -> 3}",
    "{x -> 37/10, y -> 11/10, n -> 2}",
    "{x -> 9/10 + 2/5 I, y -> 7/5 - 1/3 I, n -> 4}",
};

/* ---- A. J/Y Casoratian ---------------------------------------------- */

static void test_jy_casoratian(void) {
    sf_family("BesselJ");
    sf_to("BesselJ[n + 1, x] BesselY[n, x] - BesselJ[n, x] BesselY[n + 1, x]", "2/(Pi x)");
    sf_to("BesselJ[1/3, x] BesselY[4/3, x] - BesselJ[4/3, x] BesselY[1/3, x]", "-2/(Pi x)");
    sf_to("BesselJ[2, x] BesselY[1, x] - BesselJ[1, x] BesselY[2, x]", "2/(Pi x)");
    sf_to("BesselJ[0, x] BesselY[1, x] - BesselJ[1, x] BesselY[0, x]", "-2/(Pi x)");
    sf_to("BesselJ[n - 1, x] BesselY[n, x] - BesselJ[n, x] BesselY[n - 1, x]", "-2/(Pi x)");
    sf_to("BesselJ[n + 1, 2 x] BesselY[n, 2 x] - BesselJ[n, 2 x] BesselY[n + 1, 2 x]", "1/(Pi x)");
    sf_to("BesselJ[5/4, y] BesselY[1/4, y] - BesselJ[1/4, y] BesselY[5/4, y]", "2/(Pi y)");
    sf_zero("Pi x (BesselJ[n + 1, x] BesselY[n, x] - BesselJ[n, x] BesselY[n + 1, x]) - 2");
    sf_zero("x^2 (BesselJ[n + 1, x] BesselY[n, x] - BesselJ[n, x] BesselY[n + 1, x])^2 - 4/Pi^2");
    sf_zero("BesselJ[n + 2, x] BesselY[n, x] - BesselJ[n, x] BesselY[n + 2, x] - 4 (n + 1)/(Pi x^2)");
    sf_zero("BesselJ[n + 3, x] BesselY[n, x] - BesselJ[n, x] BesselY[n + 3, x]"
            " - 2/(Pi x) (4 (n + 1) (n + 2)/x^2 - 1)");
    sf_zero("BesselJ[7/3, x] BesselY[1/3, x] - BesselJ[1/3, x] BesselY[7/3, x] - 16/(3 Pi x^2)");
    sf_zero("BesselJ[n + 1, x] BesselY[n, x] - BesselJ[n, x] BesselY[n + 1, x]"
            " + BesselJ[n, x] - 2/(Pi x) - BesselJ[n, x]");
    sf_identity("BesselJ[n + 1, x] BesselY[n, x] - BesselJ[n, x] BesselY[n + 1, x] - 2/(Pi x)");
    sf_identity("BesselJ[n + 3, x] BesselY[n, x] - BesselJ[n, x] BesselY[n + 3, x]"
                " - 2/(Pi x) (4 (n + 1) (n + 2)/x^2 - 1)");
}

static void test_jy_wronskian_derivatives(void) {
    sf_family("BesselY");
    sf_to("x (BesselJ[n, x] D[BesselY[n, x], x] - D[BesselJ[n, x], x] BesselY[n, x])", "2/Pi");
    sf_zero("BesselJ[1/3, x] D[BesselY[1/3, x], x] - D[BesselJ[1/3, x], x] BesselY[1/3, x] - 2/(Pi x)");
    sf_zero("BesselJ[2, x] D[BesselY[2, x], x] - D[BesselJ[2, x], x] BesselY[2, x] - 2/(Pi x)");
    sf_zero("BesselJ[n, y] D[BesselY[n, y], y] - D[BesselJ[n, y], y] BesselY[n, y] - 2/(Pi y)");
    sf_zero("Pi x^2 (BesselJ[n, x] D[BesselY[n, x], x] - D[BesselJ[n, x], x] BesselY[n, x])^2 - 4/Pi");
}

/* ---- B. I/K pair ----------------------------------------------------- */

static void test_ik_pair(void) {
    sf_family("BesselI");
    sf_to("BesselI[n + 1, x] BesselK[n, x] + BesselI[n, x] BesselK[n + 1, x]", "1/x");
    sf_to("BesselI[1/4, x] BesselK[5/4, x] + BesselI[5/4, x] BesselK[1/4, x]", "1/x");
    sf_to("BesselI[1, x] BesselK[0, x] + BesselI[0, x] BesselK[1, x]", "1/x");
    sf_to("BesselI[n, x] BesselK[n - 1, x] + BesselI[n - 1, x] BesselK[n, x]", "1/x");
    sf_to("x (BesselI[n, x] D[BesselK[n, x], x] - D[BesselI[n, x], x] BesselK[n, x])", "-1");
    sf_zero("x (BesselI[n + 1, x] BesselK[n, x] + BesselI[n, x] BesselK[n + 1, x]) - 1");
    sf_zero("BesselI[n + 1, 3 x] BesselK[n, 3 x] + BesselI[n, 3 x] BesselK[n + 1, 3 x] - 1/(3 x)");
    sf_zero("BesselI[2/3, x] D[BesselK[2/3, x], x] - D[BesselI[2/3, x], x] BesselK[2/3, x] + 1/x");
    /* two-step: pure recurrence, I and K recurrences differ in sign */
    sf_zero("BesselI[n + 2, x] BesselK[n, x] + BesselI[n, x] BesselK[n + 2, x] - 2 BesselI[n, x] BesselK[n, x]"
            " - 2 (n + 1)/x (BesselI[n, x] BesselK[n + 1, x] - BesselI[n + 1, x] BesselK[n, x])");
    sf_identity("BesselI[n + 1, x] BesselK[n, x] + BesselI[n, x] BesselK[n + 1, x] - 1/x");
}

/* ---- C. BesselK evenness --------------------------------------------- */

static void test_k_even(void) {
    sf_family("BesselK");
    sf_zero("BesselK[-n, x] - BesselK[n, x]");
    sf_zero("BesselK[1 - n, x] - BesselK[n - 1, x]");
    sf_zero("BesselK[-n - 1, x] - BesselK[n + 1, x]");
    sf_zero("BesselK[-2/3, x] - BesselK[2/3, x]");
    sf_zero("BesselK[-7/4, x] - BesselK[7/4, x]");
    sf_zero("BesselK[-n, x]^2 - BesselK[n, x]^2");
    sf_zero("BesselK[-n, 2 x] - BesselK[n, 2 x]");
    /* evenness feeding the recurrence: K[2/3] = K[-2/3] joins class 1/3 */
    sf_zero("BesselK[2/3, x] - BesselK[4/3, x] + 2/(3 x) BesselK[1/3, x]");
    sf_zero("BesselK[1/3, x] - BesselK[5/3, x] + 4/(3 x) BesselK[2/3, x]");
    sf_zero("BesselK[n + 1, x] - BesselK[1 - n, x] - 2 n/x BesselK[n, x]");
    sf_zero("x BesselK[-n, x] - x BesselK[n, x] + BesselJ[n, x] - BesselJ[n, x]");
    sf_identity("BesselK[-n, x] - BesselK[n, x]");
    sf_identity("BesselK[2/3, x] - BesselK[4/3, x] + 2/(3 x) BesselK[1/3, x]");
}

/* ---- D. integer-order reflection (assumption-gated) ------------------- */

static void test_integer_reflection(void) {
    sf_family("BesselJ");
    const char* A = "Element[n, Integers]";
    sf_zero_assuming("BesselJ[-n, x] - (-1)^n BesselJ[n, x]", A, INT_PTS, 3);
    sf_zero_assuming("BesselY[-n, x] - (-1)^n BesselY[n, x]", A, INT_PTS, 3);
    sf_zero_assuming("BesselI[-n, x] - BesselI[n, x]", A, INT_PTS, 3);
    sf_zero_assuming("BesselJ[1 - n, x] + (-1)^n BesselJ[n - 1, x]", A, INT_PTS, 3);
    sf_zero_assuming("BesselJ[-n, x] BesselJ[n, x] - (-1)^n BesselJ[n, x]^2", A, INT_PTS, 3);
    sf_zero_assuming("BesselI[-n, y] - BesselI[n, y]", A, INT_PTS, 3);
    sf_zero_assuming("BesselJ[-n - 1, x] + (-1)^n BesselJ[n + 1, x]", A, INT_PTS, 3);
    /* not identities for non-integer n: must stay */
    sf_unchanged("BesselJ[-n, x] - (-1)^n BesselJ[n, x]");
    sf_unchanged("BesselI[-n, x] - BesselI[n, x]");
    sf_unchanged("BesselY[-n, x] - (-1)^n BesselY[n, x]");
}

/* ---- E. soundness / no-ops --------------------------------------------- */

static void test_unchanged(void) {
    sf_family("BesselY");
    sf_unchanged("BesselJ[n, x] BesselY[n + 1, y]");                  /* different arguments */
    sf_unchanged("BesselJ[n, x] BesselY[m, x]");                      /* unrelated orders */
    sf_unchanged("BesselJ[1/3, x] BesselY[2/3, x]");                  /* classes differ */
    sf_unchanged("BesselK[n, x] + BesselK[n, y]");
    sf_unchanged("BesselK[2/3, x]");
    sf_unchanged("BesselJ[n, x] + BesselY[n, x]");
    sf_unchanged("BesselJ[1/3, x] BesselY[4/3, x]");
    sf_unchanged("BesselI[n, x] BesselK[n + 1, x]");
    sf_unchanged("x^2 + 1");                                           /* no special function */
    sf_sound("BesselJ[n, x] BesselY[n + 1, x] + BesselJ[n + 1, x] BesselY[n, x]");
    sf_sound("BesselJ[n + 2, x] BesselY[n + 1, x] - BesselJ[n, x] BesselY[n + 3, x]");
    sf_sound("BesselI[n + 2, x] BesselK[n, x] - BesselI[n, x] BesselK[n + 2, x]");
    sf_sound("BesselK[2/3, x] BesselK[-1/3, x] - BesselK[1/3, x]^2");
}

/* ---- F. performance guards ---------------------------------------------- */

static void test_fast(void) {
    sf_family("BesselI");
    sf_fast("BesselJ[n + 6, x] BesselY[n, x] - BesselJ[n, x] BesselY[n + 6, x]", 5.0);
    sf_fast("Sum[BesselK[k/3 - 2, x], {k, 0, 12}]", 5.0);
    sf_fast("(x^2 + 1)^4 + Sin[x]^2 + Log[x]", 2.0);       /* no registered head: prefilter */
}

/* ---- G. generated per-kind identities ----------------------------------- */

static char g_buf[8192];
#define FMT(...) (snprintf(g_buf, sizeof g_buf, __VA_ARGS__), g_buf)

/* order shapes (classes 0, 1/3, 2/5, symbolic) and argument shapes */
static const char* M_SHAPES[] = { "n", "n + 1", "1/3", "4/3", "2/5", "n/2", "2 n", "a + b" };
#define N_M (sizeof M_SHAPES / sizeof M_SHAPES[0])
static const char* Z_SHAPES[] = { "x", "2 x", "x/3" };

/* For each kind: the three-term recurrence (DLMF 10.6.1 / 10.29.1) at every
 * order and argument shape, the derivative rule (10.6.2 / 10.29.2), the ODE
 * residual z^2 w'' + z w' + s1 (z^2 + s2 m^2) w with (s1, s2) = ode_sign
 * ("+-" for J, Y: + (z^2 - m^2); "-+" for I, K: - (z^2 + m^2)), and the pair
 * relation. */
static void test_generated_kind(const char* H, const char* rec_fmt, const char* der_fmt,
                                const char* ode_sign, const char* pair_fmt) {
    sf_family(H);
    for (size_t i = 0; i < N_M; i++)
        for (size_t j = 0; j < 3; j++) {
            const char* m = M_SHAPES[i];
            const char* z = Z_SHAPES[j];
            sf_zero(FMT(rec_fmt, H, m, z, H, m, z, m, z, H, m, z));
        }
    for (size_t i = 0; i < N_M; i++) {
        const char* m = M_SHAPES[i];
        sf_zero(FMT(der_fmt, H, m, H, m, m, H, m));
        sf_zero(FMT("x^2 D[%s[%s, x], {x, 2}] + x D[%s[%s, x], x] %s (x^2 %s (%s)^2) %s[%s, x]",
                    H, m, H, m, ode_sign[0] == '+' ? "+" : "-", ode_sign[1] == '+' ? "+" : "-", m, H, m));
        if (pair_fmt) sf_zero(FMT(pair_fmt, m, m, m, m));
    }
}

static void test_generated(void) {
    /* J, Y: C[m-1] + C[m+1] - (2m/z) C[m];  C' = C[m-1] - (m/z) C */
    test_generated_kind("BesselJ",
        "%s[(%s) - 1, %s] + %s[(%s) + 1, %s] - 2 (%s)/(%s) %s[%s, %s]",
        "D[%s[%s, x], x] - (%s[(%s) - 1, x] - (%s)/x %s[%s, x])", "+-",
        "BesselJ[(%s) + 1, x] BesselY[%s, x] - BesselJ[%s, x] BesselY[(%s) + 1, x] - 2/(Pi x)");
    test_generated_kind("BesselY",
        "%s[(%s) - 1, %s] + %s[(%s) + 1, %s] - 2 (%s)/(%s) %s[%s, %s]",
        "D[%s[%s, x], x] - (%s[(%s) - 1, x] - (%s)/x %s[%s, x])", "+-",
        NULL);
    /* I: I[m-1] - I[m+1] - (2m/z) I[m];  I' = I[m-1] - (m/z) I */
    test_generated_kind("BesselI",
        "%s[(%s) - 1, %s] - %s[(%s) + 1, %s] - 2 (%s)/(%s) %s[%s, %s]",
        "D[%s[%s, x], x] - (%s[(%s) - 1, x] - (%s)/x %s[%s, x])", "-+",
        "BesselI[(%s) + 1, x] BesselK[%s, x] + BesselI[%s, x] BesselK[(%s) + 1, x] - 1/x");
    /* K: K[m-1] - K[m+1] + (2m/z) K[m];  K' = -K[m-1] - (m/z) K */
    test_generated_kind("BesselK",
        "%s[(%s) - 1, %s] - %s[(%s) + 1, %s] + 2 (%s)/(%s) %s[%s, %s]",
        "D[%s[%s, x], x] - (-%s[(%s) - 1, x] - (%s)/x %s[%s, x])", "-+",
        NULL);
}

int main(void) {
    symtab_init();
    core_init();
    TEST(test_jy_casoratian);
    TEST(test_jy_wronskian_derivatives);
    TEST(test_ik_pair);
    TEST(test_k_even);
    TEST(test_integer_reflection);
    TEST(test_unchanged);
    TEST(test_fast);
    TEST(test_generated);
    return sf_finish("simp_sf_bessel", 900);
}
