#include "simp_sf_test_util.h"
#include "symtab.h"
#include "core.h"

/*
 * Simplify, special-function engine, error / exponential-integral / Fresnel
 * module (src/simp/simp_sf_errexp.c).
 *
 * Derived identities (docs/design/simp_special_functions.md, section 6),
 * each checked numerically; the domain of every one was probed in all four
 * quadrants and on the cuts, and the gated ones are tested only under the
 * assumptions that make them true:
 *
 *   Erfc = 1 - Erf;  Erf[z0, z1] = Erf[z1] - Erf[z0]          (definitions)
 *   Erfi[z] = -I Erf[I z]                    (t -> i t in the integral)
 *   InverseErfc[x] = InverseErf[1 - x];  Erf[InverseErf[x]] = x;
 *   Erfc[InverseErfc[x]] = x;  InverseErf[Erf[x]] = x (x real)
 *   n E_{n+1}(z) = e^-z - z E_n(z)           (E_n integral, by parts)
 *   E_n(z) = z^(n-1) Gamma[1 - n, z]         (z > 0; t -> t/z)
 *   LogIntegral[x] = ExpIntegralEi[Log[x]]   (li = Ei o log, all x)
 *   Ci[-x] = Ci[x] + I Pi, Chi[-x] = Chi[x] + I Pi   (x > 0)
 *   Si[I z] = I Shi[z], Shi[I z] = I Si[z]   (t -> i t)
 *   FresnelS[I z] = -I FresnelS[z], FresnelC[I z] = I FresnelC[z]
 */

static char g_buf[8192];
#define FMT(...) (snprintf(g_buf, sizeof g_buf, __VA_ARGS__), g_buf)

static const char* Z[] = { "x", "2 x", "x + 1/3", "a + b", "x/2", "-x", "x^2", "Sqrt[x]" };
#define N_Z (sizeof Z / sizeof Z[0])

/* real points inside (0, 1) for the inverse functions */
static const char* UNIT_PTS[] = { "{x -> 3/10}", "{x -> 7/10}", "{x -> 9/20}" };
/* positive reals for the gated reflections and the E_n / Gamma relation */
static const char* POS_PTS[] = { "{x -> 17/10, a -> 6/5, b -> 7/3, z -> 13/20, n -> 3/7, m -> 2/5}",
                                 "{x -> 37/10, a -> 2/3, b -> 13/5, z -> 3/10, n -> 5/9, m -> 1/7}",
                                 "{x -> 1/5, a -> 1/9, b -> 4, z -> 5/2, n -> 11/7, m -> 9/4}" };

/* ---- Erf / Erfc / Erf[z0, z1] ------------------------------------------- */

static void test_erf(void) {
    sf_family("Erf");
    for (size_t i = 0; i < N_Z; i++) {
        const char* z = Z[i];
        sf_zero(FMT("Erf[%s] + Erfc[%s] - 1", z, z));
        sf_zero(FMT("Erfc[%s] + Erfc[-(%s)] - 2", z, z));
        sf_zero(FMT("Erfc[-(%s)] - 1 - Erf[%s]", z, z));
        sf_zero(FMT("x (Erf[%s] + Erfc[%s]) - x", z, z));
        sf_zero(FMT("(Erf[%s] + Erfc[%s])^2 - 1", z, z));
        sf_zero(FMT("Erf[%s] - (1 - Erfc[%s])", z, z));
        sf_zero(FMT("Erf[%s, 2 (%s)] - Erf[2 (%s)] + Erf[%s]", z, z, z, z));
        sf_zero(FMT("Erfc[%s] Erf[%s] - Erf[%s] + Erf[%s]^2", z, z, z, z));
    }
    sf_to("Erf[x] + Erfc[x]", "1");
    sf_to("Erfc[x] + Erfc[-x]", "2");
    sf_to("Erfc[-x] - Erfc[x]", "2 Erf[x]");
    sf_to("Erf[x, y] - Erf[y] + Erf[x]", "0");
    sf_to("Erfc[x] (1 + Erf[x]) + Erf[x]^2 - Erf[x]", "Erfc[x]");
    sf_unchanged("Erf[x] + Erfc[y]");
    sf_unchanged("Erfc[x]");
    sf_identity("Erf[x] + Erfc[x] - 1");
}

/* ---- Erfi --------------------------------------------------------------- */

static void test_erfi(void) {
    sf_family("Erfi");
    for (size_t i = 0; i < N_Z; i++) {
        const char* z = Z[i];
        sf_zero(FMT("Erfi[%s] + I Erf[I (%s)]", z, z));
        sf_zero(FMT("Erfi[%s] + I - I Erfc[I (%s)]", z, z));
        sf_zero(FMT("x (Erfi[%s] + I Erf[I (%s)])", z, z));
        sf_zero(FMT("2 Erfi[%s] + 2 I Erf[I (%s)]", z, z));
        sf_zero(FMT("Erfi[%s]^2 + Erf[I (%s)]^2", z, z));
        sf_zero(FMT("Erfi[%s] Erf[I (%s)] + I Erf[I (%s)]^2", z, z, z));
    }
    sf_to("Erfi[x] + I Erf[I x]", "0");
    sf_unchanged("Erfi[x]");
    sf_unchanged("Erfi[x] + Erf[y]");
    sf_identity("Erfi[x] + I Erf[I x]");
}

/* ---- InverseErf / InverseErfc -------------------------------------------- */

static void test_inverse_erf(void) {
    sf_family("InverseErf");
    static const char* S[] = { "x", "x/2", "(1 + x)/2", "x^2", "1 - x", "2 x/3", "x/3 + 1/2", "Sqrt[x]" };
    for (size_t i = 0; i < 8; i++) {
        const char* s = S[i];
        sf_zero_assuming(FMT("InverseErfc[%s] - InverseErf[1 - (%s)]", s, s), "0 < x < 1", UNIT_PTS, 3);
        sf_zero_assuming(FMT("Erf[InverseErf[%s]] - (%s)", s, s), "0 < x < 1", UNIT_PTS, 3);
        sf_zero_assuming(FMT("Erfc[InverseErfc[%s]] - (%s)", s, s), "0 < x < 1", UNIT_PTS, 3);
        sf_zero_assuming(FMT("InverseErf[Erf[%s]] - (%s)", s, s), "0 < x < 1", UNIT_PTS, 3);
        sf_zero_assuming(FMT("Erfc[InverseErf[%s]] - 1 + (%s)", s, s), "0 < x < 1", UNIT_PTS, 3);
        sf_zero_assuming(FMT("Erf[InverseErfc[%s]] - 1 + (%s)", s, s), "0 < x < 1", UNIT_PTS, 3);
        sf_zero_assuming(FMT("x (InverseErfc[%s] - InverseErf[1 - (%s)])", s, s), "0 < x < 1", UNIT_PTS, 3);
    }
    sf_to_assuming("Erf[InverseErf[x]]", "0 < x < 1", "x", UNIT_PTS, 3);
    sf_to_assuming("Erfc[InverseErfc[x]]", "0 < x < 1", "x", UNIT_PTS, 3);
    sf_unchanged("InverseErf[Erf[x]]");          /* not an identity off the reals */
    sf_identity_real("InverseErfc[x/2] - InverseErf[1 - x/2]");
}

/* ---- ExpIntegralE ----------------------------------------------------------- */

static void test_expintegrale(void) {
    sf_family("ExpIntegralE");
    static const char* N[] = { "n", "n + 1/2", "2 n", "1/3", "m", "n + 1" };
    static const char* Zs[] = { "z", "2 z", "w" };
    for (size_t i = 0; i < 6; i++)
        for (size_t j = 0; j < 3; j++) {
            const char* n = N[i];
            const char* z = Zs[j];
            for (int k = 0; k <= 2; k++)
                sf_zero(FMT("(%s + %d) ExpIntegralE[%s + %d, %s] + (%s) ExpIntegralE[%s + %d, %s] - E^-(%s)",
                            n, k, n, k + 1, z, z, n, k, z, z));
            /* two steps: E_{n+2} = (e^-z - z (e^-z - z E_n)/n)/(n + 1) */
            sf_zero(FMT("ExpIntegralE[%s + 2, %s] - (E^-(%s) - (%s) (E^-(%s) - (%s) ExpIntegralE[%s, %s])/(%s))/((%s) + 1)",
                        n, z, z, z, z, z, n, z, n, n));
        }
    for (size_t i = 0; i < 6; i++)
        sf_zero_assuming(FMT("ExpIntegralE[%s, z] - z^((%s) - 1) Gamma[1 - (%s), z]", N[i], N[i], N[i]),
                         "z > 0", POS_PTS, 3);
    sf_to("ExpIntegralE[n + 1, z] n + z ExpIntegralE[n, z]", "E^(-z)");
    sf_unchanged("ExpIntegralE[n, z] + ExpIntegralE[n, w]");
    sf_unchanged("ExpIntegralE[n, z] - z^(n - 1) Gamma[1 - n, z]");   /* needs z > 0 */
    sf_identity("n ExpIntegralE[n + 1, z] + z ExpIntegralE[n, z] - E^-z");
}

/* ---- LogIntegral / ExpIntegralEi ----------------------------------------- */

static void test_logintegral(void) {
    sf_family("LogIntegral");
    for (size_t i = 0; i < N_Z; i++) {
        const char* z = Z[i];
        sf_zero(FMT("LogIntegral[%s] - ExpIntegralEi[Log[%s]]", z, z));
        sf_zero(FMT("x (LogIntegral[%s] - ExpIntegralEi[Log[%s]])", z, z));
        sf_zero(FMT("(LogIntegral[%s] - ExpIntegralEi[Log[%s]])^2", z, z));
        sf_zero(FMT("(LogIntegral[%s] - ExpIntegralEi[Log[%s]]) Erf[x]", z, z));
        sf_zero(FMT("LogIntegral[%s] + 1 - ExpIntegralEi[Log[%s]] - 1", z, z));
        sf_zero(FMT("2 LogIntegral[%s] - ExpIntegralEi[Log[%s]] - LogIntegral[%s]", z, z, z));
    }
    sf_to("LogIntegral[x] - ExpIntegralEi[Log[x]]", "0");
    sf_unchanged("LogIntegral[x]");
    sf_unchanged("LogIntegral[x] + ExpIntegralEi[x]");
    sf_identity("LogIntegral[x] - ExpIntegralEi[Log[x]]");
}

/* ---- Sin/Cos/Sinh/Cosh integrals ----------------------------------------- */

static void test_trig_integrals(void) {
    sf_family("SinCosIntegral");
    static const char* P[] = { "x", "2 x", "x + 1/3", "x^2", "Sqrt[x]", "x/2", "a + b", "3 x" };
    const char* A = "x > 0 && a > 0 && b > 0";
    for (size_t i = 0; i < 8; i++) {
        const char* p = P[i];
        sf_zero_assuming(FMT("CosIntegral[-(%s)] - CosIntegral[%s] - I Pi", p, p), A, POS_PTS, 3);
        sf_zero_assuming(FMT("CoshIntegral[-(%s)] - CoshIntegral[%s] - I Pi", p, p), A, POS_PTS, 3);
    }
    for (size_t i = 0; i < N_Z; i++) {
        const char* z = Z[i];
        sf_zero(FMT("SinIntegral[I (%s)] - I SinhIntegral[%s]", z, z));
        sf_zero(FMT("SinhIntegral[I (%s)] - I SinIntegral[%s]", z, z));
        sf_zero(FMT("x (SinIntegral[I (%s)] - I SinhIntegral[%s])", z, z));
        sf_zero(FMT("SinIntegral[-(%s)] + SinIntegral[%s]", z, z));
    }
    sf_to_assuming("CosIntegral[-x] - CosIntegral[x]", "x > 0", "I Pi", POS_PTS, 3);
    sf_unchanged("CosIntegral[-x] - CosIntegral[x]");        /* fails off the positive reals */
    sf_identity("SinIntegral[I x] - I SinhIntegral[x]");
}

/* ---- Fresnel ---------------------------------------------------------------- */

static void test_fresnel(void) {
    sf_family("Fresnel");
    for (size_t i = 0; i < N_Z; i++) {
        const char* z = Z[i];
        sf_zero(FMT("FresnelS[I (%s)] + I FresnelS[%s]", z, z));
        sf_zero(FMT("FresnelC[I (%s)] - I FresnelC[%s]", z, z));
        sf_zero(FMT("FresnelS[-I (%s)] - I FresnelS[%s]", z, z));
        sf_zero(FMT("FresnelC[-I (%s)] + I FresnelC[%s]", z, z));
        sf_zero(FMT("FresnelS[I (%s)]^2 + FresnelS[%s]^2", z, z));
        sf_zero(FMT("FresnelC[I (%s)] FresnelS[I (%s)] - FresnelC[%s] FresnelS[%s]", z, z, z, z));
    }
    /* a lone FresnelS[I x] stays: -I FresnelS[x] scores 7 against 6 */
    sf_unchanged("FresnelS[I x]");
    sf_unchanged("FresnelS[x] + FresnelC[x]");
    sf_identity("FresnelC[I x] - I FresnelC[x]");
}

int main(void) {
    symtab_init();
    core_init();
    TEST(test_erf);
    TEST(test_erfi);
    TEST(test_inverse_erf);
    TEST(test_expintegrale);
    TEST(test_logintegral);
    TEST(test_trig_integrals);
    TEST(test_fresnel);
    return sf_finish("simp_sf_errexp", 1200);
}
