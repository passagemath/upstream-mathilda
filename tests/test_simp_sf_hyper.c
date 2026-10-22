#include "simp_sf_test_util.h"
#include "symtab.h"
#include "core.h"

/*
 * Simplify, special-function engine, hypergeometric module
 * (src/simp/simp_sf_hyper.c).  Every identity is derived
 * (docs/design/simp_special_functions.md, section 10) and checked
 * numerically (|z| < 1 at every sample point):
 *
 *   2F1:  a F(a+1) = (theta + a) F,  b F(b+1) = (theta + b) F,
 *         (c-1) F(c-1) = (theta + c - 1) F                     (the series)
 *         c(c-1)(z-1) F(c-1) + c(c-1-(2c-a-b-1)z) F + (c-a)(c-b) z F(c+1) = 0
 *         (c-a) F(a-1) + (2a - c + (b-a) z) F + a (z-1) F(a+1) = 0   (Gauss)
 *         z(1-z) F'' + (c - (a+b+1) z) F' - a b F = 0
 *         F = (1-z)^(c-a-b) F(c-a, c-b; c; z)                  (Euler; z < 1)
 *   1F1:  a M(a+1) = (theta + a) M,  b(b-1) M(b-1) + b(1-b-z) M + z(b-a) M(b+1) = 0,
 *         (b-a) M(a-1) + (2a - b + z) M - a M(a+1) = 0,
 *         z M'' + (b - z) M' - a M = 0,  M = E^z M(b-a; b; -z)  (Kummer)
 *   0F1:  G(b+1) = (b/z) theta G,  (b-1) G(b-1) - (b-1) G - (z/b) G(b+1) = 0,
 *         z G'' + b G' - G = 0
 *   U:    U(a-1) + (b-2a-z) U(a) + a(a-b+1) U(a+1) = 0,
 *         U(a, b, z) = z^(1-b) U(a-b+1, 2-b, z)
 *   PFQ:  the parameter lists are unordered.
 */

static char g_buf[8192];
#define FMT(...) (snprintf(g_buf, sizeof g_buf, __VA_ARGS__), g_buf)

/* real points with |z| < 1, and z < 1 for the gated Euler transformation */
static const char* LT1_PTS[] = {
    "{a -> 6/5, b -> 7/3, c -> 19/4, z -> 13/20}",
    "{a -> 2/3, b -> 13/5, c -> 31/7, z -> -3/10}",
    "{a -> 1/5, b -> 9/4, c -> 7/2, z -> 1/3}",
};

/* ---- 2F1 -------------------------------------------------------------------- */

static void test_2f1(void) {
    sf_family("Hypergeometric2F1");
    static const char* S[][3] = { { "a", "b", "c" }, { "a + 1/2", "b", "c" }, { "2 a", "b", "c" },
                                  { "a", "b", "c + 1/3" }, { "a", "a + b", "c" } };
    for (size_t i = 0; i < 5; i++) {
        const char *A = S[i][0], *B = S[i][1], *C = S[i][2];
        char F[256];
        snprintf(F, sizeof F, "Hypergeometric2F1[%s, %s, %s, z]", A, B, C);
        sf_zero(FMT("(%s) Hypergeometric2F1[(%s) + 1, %s, %s, z] - z D[%s, z] - (%s) %s", A, A, B, C, F, A, F));
        sf_zero(FMT("(%s) Hypergeometric2F1[%s, (%s) + 1, %s, z] - z D[%s, z] - (%s) %s", B, A, B, C, F, B, F));
        sf_zero(FMT("((%s) - 1) Hypergeometric2F1[%s, %s, (%s) - 1, z] - z D[%s, z] - ((%s) - 1) %s", C, A, B, C, F, C, F));
        char R4[1024];
        snprintf(R4, sizeof R4,
                 "((%s) ((%s) - 1) (z - 1) Hypergeometric2F1[%s, %s, (%s) - 1, z] + (%s) ((%s) - 1 - (2 (%s) - (%s) - (%s) - 1) z) %s "
                 "+ ((%s) - (%s)) ((%s) - (%s)) z Hypergeometric2F1[%s, %s, (%s) + 1, z])",
                 C, C, A, B, C, C, C, C, A, B, F, C, A, C, B, A, B, C);
        sf_zero(R4);
        char R5[1024];
        snprintf(R5, sizeof R5,
                 "(((%s) - (%s)) Hypergeometric2F1[(%s) - 1, %s, %s, z] + (2 (%s) - (%s) + ((%s) - (%s)) z) %s "
                 "+ (%s) (z - 1) Hypergeometric2F1[(%s) + 1, %s, %s, z])",
                 C, A, A, B, C, A, C, B, A, F, A, A, B, C);
        sf_zero(R5);
        sf_zero(FMT("z (1 - z) D[%s, {z, 2}] + ((%s) - ((%s) + (%s) + 1) z) D[%s, z] - (%s) (%s) %s", F, C, A, B, F, A, B, F));
        sf_zero(FMT("%s - Hypergeometric2F1[%s, %s, %s, z]", F, B, A, C));
        sf_zero(FMT("D[%s, z] - (%s) (%s)/(%s) Hypergeometric2F1[(%s) + 1, (%s) + 1, (%s) + 1, z]", F, A, B, C, A, B, C));
        sf_zero(FMT("x %s", R4));
        sf_zero(FMT("%s %s", R5, F));
    }
    for (size_t i = 0; i < 4; i++) {
        const char *A = S[i][0], *B = S[i][1], *C = S[i][2];
        sf_zero_assuming(FMT("Hypergeometric2F1[%s, %s, %s, z] - (1 - z)^((%s) - (%s) - (%s)) Hypergeometric2F1[(%s) - (%s), (%s) - (%s), %s, z]",
                             A, B, C, C, A, B, C, A, C, B, C), "z < 1", LT1_PTS, 3);
    }
    sf_unchanged("Hypergeometric2F1[a, b, c, z] + Hypergeometric2F1[a, b, c, w]");
    sf_unchanged("Hypergeometric2F1[a, b, c, z] - (1 - z)^(c - a - b) Hypergeometric2F1[c - a, c - b, c, z]");
}

/* ---- 1F1 -------------------------------------------------------------------- */

static void test_1f1(void) {
    sf_family("Hypergeometric1F1");
    static const char* S[][2] = { { "a", "b" }, { "a + 1/2", "b" }, { "2 a", "b" }, { "a", "b + 1/3" },
                                  { "a", "a + b" }, { "1/3", "b" } };
    for (size_t i = 0; i < 6; i++) {
        const char *A = S[i][0], *B = S[i][1];
        char M[256];
        snprintf(M, sizeof M, "Hypergeometric1F1[%s, %s, z]", A, B);
        sf_zero(FMT("(%s) Hypergeometric1F1[(%s) + 1, %s, z] - z D[%s, z] - (%s) %s", A, A, B, M, A, M));
        sf_zero(FMT("((%s) - 1) Hypergeometric1F1[%s, (%s) - 1, z] - z D[%s, z] - ((%s) - 1) %s", B, A, B, M, B, M));
        char R3[1024];
        snprintf(R3, sizeof R3,
                 "((%s) ((%s) - 1) Hypergeometric1F1[%s, (%s) - 1, z] + (%s) (1 - (%s) - z) %s + z ((%s) - (%s)) Hypergeometric1F1[%s, (%s) + 1, z])",
                 B, B, A, B, B, B, M, B, A, A, B);
        sf_zero(R3);
        sf_zero(FMT("((%s) - (%s)) Hypergeometric1F1[(%s) - 1, %s, z] + (2 (%s) - (%s) + z) %s - (%s) Hypergeometric1F1[(%s) + 1, %s, z]",
                    B, A, A, B, A, B, M, A, A, B));
        sf_zero(FMT("z D[%s, {z, 2}] + ((%s) - z) D[%s, z] - (%s) %s", M, B, M, A, M));
        sf_zero(FMT("%s - E^z Hypergeometric1F1[(%s) - (%s), %s, -z]", M, B, A, B));
        sf_zero(FMT("x %s", R3));
        sf_zero(FMT("(%s - E^z Hypergeometric1F1[(%s) - (%s), %s, -z]) %s", M, B, A, B, M));
        sf_zero(FMT("D[%s, z] - (%s)/(%s) Hypergeometric1F1[(%s) + 1, (%s) + 1, z]", M, A, B, A, B));
    }
    sf_to("Hypergeometric1F1[a, b, z] - E^z Hypergeometric1F1[b - a, b, -z]", "0");
    sf_unchanged("Hypergeometric1F1[a, b, z] + Hypergeometric1F1[a, b, w]");
}

/* ---- 0F1 -------------------------------------------------------------------- */

static void test_0f1(void) {
    sf_family("Hypergeometric0F1");
    /* (no half-integer b: there 0F1 evaluates to elementary functions) */
    static const char* B[] = { "b", "b + 1/2", "2 b", "b/3", "a + b", "1/3", "7/3", "b + 1" };
    static const char* Z[] = { "z", "2 z" };
    for (size_t i = 0; i < 8; i++)
        for (size_t j = 0; j < 2; j++) {
            const char *b = B[i], *z = Z[j];
            char G[256];
            snprintf(G, sizeof G, "Hypergeometric0F1[%s, %s]", b, z);
            sf_zero(FMT("Hypergeometric0F1[(%s) + 1, %s] - (%s)/(%s) (%s) D[%s, z]/D[%s, z]", b, z, b, z, z, G, z));
            sf_zero(FMT("((%s) - 1) Hypergeometric0F1[(%s) - 1, %s] - ((%s) - 1) %s - (%s)/(%s) Hypergeometric0F1[(%s) + 1, %s]",
                        b, b, z, b, G, z, b, b, z));
            if (j == 0) {
                sf_zero(FMT("z D[%s, {z, 2}] + (%s) D[%s, z] - %s", G, b, G, G));
                sf_zero(FMT("x (((%s) - 1) Hypergeometric0F1[(%s) - 1, z] - ((%s) - 1) %s - z/(%s) Hypergeometric0F1[(%s) + 1, z])",
                            b, b, b, G, b, b));
            } else {
                sf_zero(FMT("(((%s) - 1) Hypergeometric0F1[(%s) - 1, %s] - ((%s) - 1) %s - (%s)/(%s) Hypergeometric0F1[(%s) + 1, %s]) %s",
                            b, b, z, b, G, z, b, b, z, G));
                /* the three-term relation at b + 1: G(b+2) = b (b+1) (G(b) - G(b+1))/Z */
                sf_zero(FMT("Hypergeometric0F1[(%s) + 2, %s] - (%s) ((%s) + 1) (%s - Hypergeometric0F1[(%s) + 1, %s])/(%s)",
                            b, z, b, b, G, b, z, z));
            }
        }
    sf_unchanged("Hypergeometric0F1[b, z] + Hypergeometric0F1[b, w]");
}

/* ---- U -------------------------------------------------------------------- */

static void test_u(void) {
    sf_family("HypergeometricU");
    static const char* S[][2] = { { "a", "b" }, { "a + 1/2", "b" }, { "2 a", "b" }, { "a", "b + 1/3" },
                                  { "a", "a + b" }, { "1/3", "b" }, { "a", "2 b" }, { "a + 1", "b" } };
    static const char* Z[] = { "z", "2 z" };
    for (size_t i = 0; i < 8; i++) {
        const char *A = S[i][0], *B = S[i][1];
        for (size_t j = 0; j < 2; j++)
            sf_zero(FMT("HypergeometricU[(%s) - 1, %s, %s] + ((%s) - 2 (%s) - (%s)) HypergeometricU[%s, %s, %s] + (%s) ((%s) - (%s) + 1) HypergeometricU[(%s) + 1, %s, %s]",
                        A, B, Z[j], B, A, Z[j], A, B, Z[j], A, A, B, A, B, Z[j]));
        sf_zero(FMT("HypergeometricU[%s, %s, z] - z^(1 - (%s)) HypergeometricU[(%s) - (%s) + 1, 2 - (%s), z]", A, B, B, A, B, B));
        sf_zero(FMT("x (HypergeometricU[(%s) - 1, %s, z] + ((%s) - 2 (%s) - z) HypergeometricU[%s, %s, z] + (%s) ((%s) - (%s) + 1) HypergeometricU[(%s) + 1, %s, z])",
                    A, B, B, A, A, B, A, A, B, A, B));
        sf_zero(FMT("(HypergeometricU[(%s) - 1, %s, z] + ((%s) - 2 (%s) - z) HypergeometricU[%s, %s, z] + (%s) ((%s) - (%s) + 1) HypergeometricU[(%s) + 1, %s, z]) HypergeometricU[%s, %s, z]",
                    A, B, B, A, A, B, A, A, B, A, B, A, B));
        sf_zero(FMT("D[HypergeometricU[%s, %s, z], z] + (%s) HypergeometricU[(%s) + 1, (%s) + 1, z]", A, B, A, A, B));
        sf_zero(FMT("HypergeometricU[%s, %s, z] + ((%s) - 2 ((%s) + 1) - z) HypergeometricU[(%s) + 1, %s, z] + ((%s) + 1) ((%s) - (%s) + 2) HypergeometricU[(%s) + 2, %s, z]",
                    A, B, B, A, A, B, A, A, B, A, B));
    }
    sf_unchanged("HypergeometricU[a, b, z] + HypergeometricU[a, b, w]");
}

/* ---- PFQ parameter order ------------------------------------------------- */

static void test_pfq(void) {
    sf_family("HypergeometricPFQ");
    static const char* P3[][3] = { { "a", "b", "c" }, { "a + 1", "b", "c/2" }, { "k", "q", "s" } };
    static const int perm[6][3] = { { 0, 1, 2 }, { 0, 2, 1 }, { 1, 0, 2 }, { 1, 2, 0 }, { 2, 0, 1 }, { 2, 1, 0 } };
    for (size_t i = 0; i < 3; i++)
        for (int p = 0; p < 6; p++)
            for (int l = 0; l < 2; l++) {
                const char* u[3] = { P3[i][perm[p][0]], P3[i][perm[p][1]], P3[i][perm[p][2]] };
                const char* lo = l ? "{nu, m}" : "{m, nu}";
                sf_zero(FMT("HypergeometricPFQ[{%s, %s, %s}, %s, z/2] - HypergeometricPFQ[{%s, %s, %s}, {m, nu}, z/2]",
                            u[0], u[1], u[2], lo, P3[i][0], P3[i][1], P3[i][2]));
            }
    static const char* S2[][2] = { { "a + 1", "b" }, { "b + 1/2", "a" }, { "2 a", "b" }, { "a", "a + 1" },
                                   { "k + 1", "q" }, { "q", "k + 1/3" }, { "s", "a" }, { "b", "nu" } };
    for (size_t i = 0; i < 8; i++)
        sf_zero(FMT("HypergeometricPFQ[{%s, %s}, {c}, z] - HypergeometricPFQ[{%s, %s}, {c}, z]",
                    S2[i][0], S2[i][1], S2[i][1], S2[i][0]));
    sf_zero("HypergeometricPFQ[{a, b}, {a}, z] - (1 - z)^(-b)");
    sf_zero("HypergeometricPFQ[{b, a}, {c}, z] - Hypergeometric2F1[a, b, c, z]");
    sf_zero("HypergeometricPFQ[{b + 1, a}, {c}, z] - Hypergeometric2F1[a, b + 1, c, z]");
    sf_zero("HypergeometricPFQ[{q, k, s, a}, {nu, m, c}, z/3] - HypergeometricPFQ[{a, k, q, s}, {c, m, nu}, z/3]");
    sf_zero("x (HypergeometricPFQ[{s, k}, {m}, z] - HypergeometricPFQ[{k, s}, {m}, z])");
    sf_unchanged("HypergeometricPFQ[{a, b, c}, {m, nu}, z] + HypergeometricPFQ[{a, b, c}, {m, nu}, w]");
}

int main(void) {
    symtab_init();
    core_init();
    TEST(test_2f1);
    TEST(test_1f1);
    TEST(test_0f1);
    TEST(test_u);
    TEST(test_pfq);
    return sf_finish("simp_sf_hyper", 1000);
}
