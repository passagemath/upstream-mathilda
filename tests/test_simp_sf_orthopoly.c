#include "simp_sf_test_util.h"
#include "symtab.h"
#include "core.h"

/*
 * Simplify, special-function engine, orthogonal-polynomial module
 * (src/simp/simp_sf_orthopoly.c).  Every identity is derived
 * (docs/design/simp_special_functions.md, section 8) and was checked
 * numerically at NON-INTEGER and complex degree before it was encoded:
 *
 *   (n+1) P[n+1] = (2n+1) x P[n] - n P[n-1]           (P and Q, all n)
 *   P[n+1] Q[n] - P[n] Q[n+1] = 1/(n+1)               (Casoratian)
 *   P[-n-1] = P[n]                                     (ODE in n(n+1))
 *   C[n+1] = 2x C[n] - C[n-1]                          (T and U)
 *   T[-n] = T[n];  U[-n] = -U[n-2];  T[n] = (U[n] - U[n-2])/2
 *   H[n+1] = 2x H[n] - 2n H[n-1]
 *   (n+1) L^a[n+1] = (2n+1+a-x) L^a[n] - (n+a) L^a[n-1]
 *   x L^(a+1)[n] = (n+a+1) L^a[n] - (n+1) L^a[n+1];  L^a[n] = L^(a+1)[n] - L^(a+1)[n-1]
 *   ODEs: (1-x^2) P'' - 2x P' + n(n+1) P = 0,  (1-x^2) T'' - x T' + n^2 T = 0,
 *         (1-x^2) U'' - 3x U' + n(n+2) U = 0,  H'' - 2x H' + 2n H = 0,
 *         x L'' + (a+1-x) L' + n L = 0
 *   parity C[n, -x] = (-1)^n C[n, x]   (P, T, U, H; integer n)
 */

static char g_buf[8192];
#define FMT(...) (snprintf(g_buf, sizeof g_buf, __VA_ARGS__), g_buf)

static const char* N[] = { "n", "n + 1", "1/3", "4/3", "2/5", "n/2", "2 n", "a + b" };
#define N_N (sizeof N / sizeof N[0])
static const char* X[] = { "x", "x/2", "2 x/3" };

/* Legendre P/Q: inside (-1, 1) and off the real axis, where Q's Casoratian
 * was verified */
static const char* LEG_PTS[] = {
    "{x -> 13/20, n -> 3/7, a -> 6/5, b -> 7/3}",
    "{x -> -3/10, n -> 5/9, a -> 2/3, b -> 13/5}",
    "{x -> 1/4 + I/5, n -> 2/7 + I/9, a -> 4/3 + I/4, b -> 5/2 - I/5}",
};
static const char* INT_PTS[] = { "{x -> 13/20, n -> 3}", "{x -> -3/10, n -> 4}", "{x -> 1/4 + I/5, n -> 5}" };

/* ---- Legendre P ------------------------------------------------------------ */

static void test_legendreP(void) {
    sf_family("LegendreP");
    for (size_t i = 0; i < N_N; i++)
        for (size_t j = 0; j < 3; j++)
            sf_zero_at(FMT("((%s) + 1) LegendreP[(%s) + 1, %s] - (2 (%s) + 1) (%s) LegendreP[%s, %s] + (%s) LegendreP[(%s) - 1, %s]",
                           N[i], N[i], X[j], N[i], X[j], N[i], X[j], N[i], N[i], X[j]), LEG_PTS, 3);
    for (size_t i = 0; i < N_N; i++) {
        const char* n = N[i];
        sf_zero_at(FMT("(1 - x^2) D[LegendreP[%s, x], {x, 2}] - 2 x D[LegendreP[%s, x], x] + (%s) ((%s) + 1) LegendreP[%s, x]",
                       n, n, n, n, n), LEG_PTS, 3);
        sf_zero_at(FMT("LegendreP[-(%s) - 1, x] - LegendreP[%s, x]", n, n), LEG_PTS, 3);
        sf_zero_at(FMT("x (((%s) + 1) LegendreP[(%s) + 1, x] - (2 (%s) + 1) x LegendreP[%s, x] + (%s) LegendreP[(%s) - 1, x])",
                       n, n, n, n, n, n), LEG_PTS, 3);
        sf_zero_at(FMT("LegendreP[(%s) + 1, x] LegendreQ[%s, x] - LegendreP[%s, x] LegendreQ[(%s) + 1, x] - 1/((%s) + 1)",
                       n, n, n, n, n), LEG_PTS, 3);
    }
    sf_to("LegendreP[n + 1, x] LegendreQ[n, x] - LegendreP[n, x] LegendreQ[n + 1, x]", "1/(1 + n)");
    sf_unchanged("LegendreP[n, x] + LegendreP[n, y]");
    sf_zero_assuming("LegendreP[n, -x] - (-1)^n LegendreP[n, x]", "Element[n, Integers]", INT_PTS, 3);
}

/* ---- Legendre Q ------------------------------------------------------------ */

static void test_legendreQ(void) {
    sf_family("LegendreQ");
    for (size_t i = 0; i < N_N; i++)
        for (size_t j = 0; j < 3; j++)
            sf_zero_at(FMT("((%s) + 1) LegendreQ[(%s) + 1, %s] - (2 (%s) + 1) (%s) LegendreQ[%s, %s] + (%s) LegendreQ[(%s) - 1, %s]",
                           N[i], N[i], X[j], N[i], X[j], N[i], X[j], N[i], N[i], X[j]), LEG_PTS, 3);
    for (size_t i = 0; i < N_N; i++) {
        const char* n = N[i];
        sf_zero_at(FMT("(1 - x^2) D[LegendreQ[%s, x], {x, 2}] - 2 x D[LegendreQ[%s, x], x] + (%s) ((%s) + 1) LegendreQ[%s, x]",
                       n, n, n, n, n), LEG_PTS, 3);
        sf_zero_at(FMT("x (((%s) + 1) LegendreQ[(%s) + 1, x] - (2 (%s) + 1) x LegendreQ[%s, x] + (%s) LegendreQ[(%s) - 1, x])",
                       n, n, n, n, n, n), LEG_PTS, 3);
        sf_zero_at(FMT("LegendreQ[(%s) + 2, x] - ((2 (%s) + 3) x LegendreQ[(%s) + 1, x] - ((%s) + 1) LegendreQ[%s, x])/((%s) + 2)",
                       n, n, n, n, n, n), LEG_PTS, 3);
    }
    sf_unchanged("LegendreQ[n, x] + LegendreQ[n + 1/2, x]");
    sf_unchanged("LegendreQ[-n - 1, x] - LegendreQ[n, x]");   /* no reflection for Q */
}

/* ---- Chebyshev T and U --------------------------------------------------- */

static void test_chebyshevT(void) {
    sf_family("ChebyshevT");
    for (size_t i = 0; i < N_N; i++)
        for (size_t j = 0; j < 3; j++)
            sf_zero(FMT("ChebyshevT[(%s) + 1, %s] - 2 (%s) ChebyshevT[%s, %s] + ChebyshevT[(%s) - 1, %s]",
                        N[i], X[j], X[j], N[i], X[j], N[i], X[j]));
    for (size_t i = 0; i < N_N; i++) {
        const char* n = N[i];
        sf_zero(FMT("(1 - x^2) D[ChebyshevT[%s, x], {x, 2}] - x D[ChebyshevT[%s, x], x] + (%s)^2 ChebyshevT[%s, x]", n, n, n, n));
        sf_zero(FMT("ChebyshevT[-(%s), x] - ChebyshevT[%s, x]", n, n));
        sf_zero(FMT("ChebyshevT[%s, x] - (ChebyshevU[%s, x] - ChebyshevU[(%s) - 2, x])/2", n, n, n));
    }
    for (int k = 0; k < 4; k++)
        sf_zero_assuming(FMT("ChebyshevT[n + %d, -x] - (-1)^(n + %d) ChebyshevT[n + %d, x]", k, k, k),
                         "Element[n, Integers]", INT_PTS, 3);
    sf_to("ChebyshevT[n + 1, x] + ChebyshevT[n - 1, x]", "2 x ChebyshevT[n, x]");
    sf_unchanged("ChebyshevT[n, x] + ChebyshevT[n, y]");
}

static void test_chebyshevU(void) {
    sf_family("ChebyshevU");
    for (size_t i = 0; i < N_N; i++)
        for (size_t j = 0; j < 3; j++)
            sf_zero(FMT("ChebyshevU[(%s) + 1, %s] - 2 (%s) ChebyshevU[%s, %s] + ChebyshevU[(%s) - 1, %s]",
                        N[i], X[j], X[j], N[i], X[j], N[i], X[j]));
    for (size_t i = 0; i < N_N; i++) {
        const char* n = N[i];
        sf_zero(FMT("(1 - x^2) D[ChebyshevU[%s, x], {x, 2}] - 3 x D[ChebyshevU[%s, x], x] + (%s) ((%s) + 2) ChebyshevU[%s, x]", n, n, n, n, n));
        sf_zero(FMT("ChebyshevU[-(%s), x] + ChebyshevU[(%s) - 2, x]", n, n));
        sf_zero(FMT("x (ChebyshevU[(%s) + 1, x] - 2 x ChebyshevU[%s, x] + ChebyshevU[(%s) - 1, x])", n, n, n));
    }
    for (int k = 0; k < 4; k++)
        sf_zero(FMT("(ChebyshevU[n + %d, x] - 2 x ChebyshevU[n + %d, x] + ChebyshevU[n + %d, x]) ChebyshevU[n, x]", k + 1, k, k - 1));
    sf_unchanged("ChebyshevU[n, x] + ChebyshevU[n + 1/2, x]");
}

/* ---- Hermite --------------------------------------------------------------- */

static void test_hermite(void) {
    sf_family("HermiteH");
    for (size_t i = 0; i < N_N; i++)
        for (size_t j = 0; j < 3; j++)
            sf_zero(FMT("HermiteH[(%s) + 1, %s] - 2 (%s) HermiteH[%s, %s] + 2 (%s) HermiteH[(%s) - 1, %s]",
                        N[i], X[j], X[j], N[i], X[j], N[i], N[i], X[j]));
    for (size_t i = 0; i < N_N; i++) {
        const char* n = N[i];
        sf_zero(FMT("D[HermiteH[%s, x], {x, 2}] - 2 x D[HermiteH[%s, x], x] + 2 (%s) HermiteH[%s, x]", n, n, n, n));
        sf_zero(FMT("x (HermiteH[(%s) + 1, x] - 2 x HermiteH[%s, x] + 2 (%s) HermiteH[(%s) - 1, x])", n, n, n, n));
        sf_zero(FMT("HermiteH[(%s) + 2, x] - 2 x HermiteH[(%s) + 1, x] + 2 ((%s) + 1) HermiteH[%s, x]", n, n, n, n));
    }
    for (int k = 0; k < 4; k++)
        sf_zero_assuming(FMT("HermiteH[n + %d, -x] - (-1)^(n + %d) HermiteH[n + %d, x]", k, k, k),
                         "Element[n, Integers]", INT_PTS, 3);
    sf_unchanged("HermiteH[n, x] + HermiteH[n, y]");
    sf_unchanged("HermiteH[n, -x] - (-1)^n HermiteH[n, x]");     /* integer n only */
}

/* ---- Laguerre -------------------------------------------------------------- */

static void test_laguerre(void) {
    sf_family("LaguerreL");
    for (size_t i = 0; i < N_N; i++)
        for (size_t j = 0; j < 3; j++)
            sf_zero(FMT("((%s) + 1) LaguerreL[(%s) + 1, %s] - (2 (%s) + 1 - (%s)) LaguerreL[%s, %s] + (%s) LaguerreL[(%s) - 1, %s]",
                        N[i], N[i], X[j], N[i], X[j], N[i], X[j], N[i], N[i], X[j]));
    static const char* A[] = { "a", "a + 1/2" };
    for (size_t i = 0; i < N_N; i++)
        for (size_t j = 0; j < 2; j++)
            sf_zero(FMT("((%s) + 1) LaguerreL[(%s) + 1, %s, x] - (2 (%s) + 1 + (%s) - x) LaguerreL[%s, %s, x] + ((%s) + (%s)) LaguerreL[(%s) - 1, %s, x]",
                        N[i], N[i], A[j], N[i], A[j], N[i], A[j], N[i], A[j], N[i], A[j]));
    for (size_t i = 0; i < N_N; i++)
        sf_zero(FMT("x LaguerreL[%s, a + 1, x] - ((%s) + a + 1) LaguerreL[%s, a, x] + ((%s) + 1) LaguerreL[(%s) + 1, a, x]",
                    N[i], N[i], N[i], N[i], N[i]));
    for (size_t i = 0; i < 4; i++) {
        const char* n = N[i];
        sf_zero(FMT("x D[LaguerreL[%s, x], {x, 2}] + (1 - x) D[LaguerreL[%s, x], x] + (%s) LaguerreL[%s, x]", n, n, n, n));
        sf_zero(FMT("x D[LaguerreL[%s, a, x], {x, 2}] + (a + 1 - x) D[LaguerreL[%s, a, x], x] + (%s) LaguerreL[%s, a, x]", n, n, n, n));
        /* L^a[n+1] = L^(a+1)[n+1] - L^(a+1)[n] (from the a-shift and the
         * n-recurrence); at symbolic a -- numeric LaguerreL does not yet
         * evaluate at a = -1, so the a = -1 instance cannot be checked */
        sf_zero(FMT("LaguerreL[(%s) + 1, a, x] - LaguerreL[(%s) + 1, a + 1, x] + LaguerreL[%s, a + 1, x]", n, n, n));
    }
    sf_unchanged("LaguerreL[n, x] + LaguerreL[n, y]");
    sf_unchanged("LaguerreL[n, a, x] + LaguerreL[n, b, x]");
}

int main(void) {
    symtab_init();
    core_init();
    TEST(test_legendreP);
    TEST(test_legendreQ);
    TEST(test_chebyshevT);
    TEST(test_chebyshevU);
    TEST(test_hermite);
    TEST(test_laguerre);
    return sf_finish("simp_sf_orthopoly", 1000);
}
