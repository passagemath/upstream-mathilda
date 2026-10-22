#include "simp_sf_test_util.h"
#include "symtab.h"
#include "core.h"

/*
 * Simplify, special-function engine, elliptic / ProductLog / QPochhammer
 * module (src/simp/simp_sf_misc.c).  Derived identities
 * (docs/design/simp_special_functions.md, section 11), each checked
 * numerically:
 *
 *   E(m) K(1-m) + E(1-m) K(m) - K(m) K(1-m) = Pi/2       (Legendre; every m)
 *   W E^W = z,  E^(k W) = (z/W)^k,  Log[W] = Log[z] - W    (W = ProductLog[z])
 *   ProductLog[u E^u] = u                                 (u real, u >= -1)
 *   (a; q)_{n+1} = (a; q)_n (1 - a q^n),  (a; q)_{n+k} = (a; q)_n (a q^n; q)_k
 */

static char g_buf[8192];
#define FMT(...) (snprintf(g_buf, sizeof g_buf, __VA_ARGS__), g_buf)

static const char* POS_PTS[] = { "{x -> 17/10, a -> 6/5}", "{x -> 3/10, a -> 2/3}", "{x -> 4, a -> 1/9}" };

/* ---- elliptic --------------------------------------------------------------- */

static void test_elliptic(void) {
    sf_family("EllipticK/E");
    static const char* M[] = { "m", "m/2", "m^2", "1/3", "2/5", "a", "m + 1/10", "1 - m" };
    for (size_t i = 0; i < 8; i++) {
        const char* p = M[i];
        char R[1024];
        snprintf(R, sizeof R,
                 "(EllipticE[%s] EllipticK[1 - (%s)] + EllipticE[1 - (%s)] EllipticK[%s] - EllipticK[%s] EllipticK[1 - (%s)] - Pi/2)",
                 p, p, p, p, p, p);
        sf_zero(R);
        sf_zero(FMT("x %s", R));
        sf_zero(FMT("%s EllipticK[%s]", R, p));
        sf_zero(FMT("%s^2", R));
        sf_zero(FMT("2/Pi (EllipticE[%s] EllipticK[1 - (%s)] + EllipticE[1 - (%s)] EllipticK[%s] - EllipticK[%s] EllipticK[1 - (%s)]) - 1",
                    p, p, p, p, p, p));
        sf_zero(FMT("%s + EllipticE[%s] - EllipticE[%s]", R, p, p));
    }
    sf_zero("D[EllipticK[m], m] - (EllipticE[m] - (1 - m) EllipticK[m])/(2 m (1 - m))");
    sf_zero("D[EllipticE[m], m] - (EllipticE[m] - EllipticK[m])/(2 m)");
    sf_to("EllipticE[m] EllipticK[1 - m] + EllipticE[1 - m] EllipticK[m] - EllipticK[m] EllipticK[1 - m]", "1/2 Pi");
    sf_unchanged("EllipticE[m] EllipticK[1 - y] + EllipticE[1 - m] EllipticK[m]");
}

/* ---- ProductLog --------------------------------------------------------------- */

static void test_productlog(void) {
    sf_family("ProductLog");
    static const char* Z[] = { "z", "2 z", "z^2", "x", "z + 1/3", "a b", "3/z", "-z" };
    for (size_t i = 0; i < 8; i++) {
        const char* z = Z[i];
        sf_zero(FMT("ProductLog[%s] E^ProductLog[%s] - (%s)", z, z, z));
        sf_zero(FMT("E^ProductLog[%s] - (%s)/ProductLog[%s]", z, z, z));
        sf_zero(FMT("E^(2 ProductLog[%s]) - ((%s)/ProductLog[%s])^2", z, z, z));
        sf_zero(FMT("Log[ProductLog[%s]] - Log[%s] + ProductLog[%s]", z, z, z));
        sf_zero(FMT("ProductLog[%s] E^ProductLog[%s]/(%s) - 1", z, z, z));
        sf_zero(FMT("x (ProductLog[%s] E^ProductLog[%s] - (%s))", z, z, z));
    }
    static const char* U[] = { "x", "2 x", "x/3", "x + a" };
    for (size_t i = 0; i < 4; i++)
        sf_zero_assuming(FMT("ProductLog[(%s) E^(%s)] - (%s)", U[i], U[i], U[i]), "x > 0 && a > 0", POS_PTS, 3);
    sf_unchanged("ProductLog[x E^x]");                   /* fails for x < -1 */
    sf_to("ProductLog[z] E^ProductLog[z]", "z");
}

/* ---- QPochhammer ---------------------------------------------------------------- */

static void test_qpochhammer(void) {
    sf_family("QPochhammer");
    static const char* A[] = { "a", "a/2", "a^2", "2 a" };
    static const char* N[] = { "n", "n + 1", "2 n" };
    for (size_t i = 0; i < 4; i++)
        for (size_t j = 0; j < 3; j++)
            for (int k = 1; k <= 4; k++)
                sf_zero(FMT("QPochhammer[%s, q, (%s) + %d] - QPochhammer[%s, q, %s] QPochhammer[(%s) q^(%s), q, %d]",
                            A[i], N[j], k, A[i], N[j], A[i], N[j], k));
    sf_to("QPochhammer[a, q, n + 1]/QPochhammer[a, q, n]", "1 - a q^n");
    sf_zero("QPochhammer[a, q, n + 1] - (1 - a q^n) QPochhammer[a, q, n]");
    sf_zero("QPochhammer[a, q, n] - (1 - a q^(n - 1)) QPochhammer[a, q, n - 1]");
    sf_unchanged("QPochhammer[a, q, n] + QPochhammer[b, q, n]");
}

int main(void) {
    symtab_init();
    core_init();
    TEST(test_elliptic);
    TEST(test_productlog);
    TEST(test_qpochhammer);
    return sf_finish("simp_sf_misc", 500);
}
