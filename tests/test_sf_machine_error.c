/* Tests for the error-function family's machine kernels
 * (src/special_functions/sf_machine_error.c) and their wiring:
 *
 *   - kernel spot values against mpmath (40 digits), placed on both sides of
 *     every regime boundary (Maclaurin / Faddeeva table / Salzer series for
 *     complex erf; Maclaurin / auxiliary tables / asymptotics for real
 *     Fresnel; the w-segments and the Halley tail of InverseErf(c); ...),
 *     each within 4 ulps (and, next to an axis, each COMPONENT within 4 ulps
 *     of itself);
 *   - the declines (overflow, poles, subnormal InverseErfc) return false so
 *     the builtin falls through to its MPFR path, and that path still answers;
 *   - the five surfaces -- scalar == unpacked List == packed List == NDArray
 *     == Compile (scalar and rank-1) -- agree BIT FOR BIT for every head, for
 *     machine-real and machine-complex arguments;
 *   - Compile[] lowers every head at the scalar and rank-1 shapes, real and
 *     complex.
 *
 * Reference values: mpmath 1.3 (erf/erfc/erfi/fresnelc/fresnels/sin at 40
 * digits; erfinv/erfcinv by Newton on mpmath's erf/erfc at 40 digits, see
 * tools/sf_fit/fitlib.py).
 */
#include "eval.h"
#include "parse.h"
#include "expr.h"
#include "symtab.h"
#include "core.h"
#include "print.h"
#include "test_utils.h"
#include "special_functions/sf_machine_error.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct { const char* k; double x, y, re, im; } CRef;
typedef struct { const char* k; double x, v; } RRef;

/* ---- mpmath references (generated; see the header comment) ---------- */
static const CRef C_REFS[] = {
    { "sf_erf_c", 0.5, 0.6, 0.70492131565630478, 0.55930009597086028 },
    { "sf_erf_c", 0.55, 0.59, 0.74728769900470893, 0.51373135813113802 },
    { "sf_erf_c", 1.2, -0.7, 1.066370012980338, -0.12024336401600993 },
    { "sf_erf_c", 3.9, 3.9, 0.90039595948500351, -0.0229973246349751 },
    { "sf_erf_c", -2.5, 1.5, -1.0004844145745748, 0.0034035003087279406 },
    { "sf_erf_c", 4.1, 0.3, 1.0000000059853129, 4.2040903575910231e-09 },
    { "sf_erf_c", 0.2, 4.2, 6077869.6225590566, -355618.41659534985 },
    { "sf_erf_c", 1.0, 10.0, 4.8386522261160068e+41, 2.7770225151412483e+41 },
    { "sf_erf_c", 0.0, 3.0, 0, 1629.9946226015657 },
    { "sf_erf_c", 0.0, -0.3, 0, -0.34894933875893619 },
    { "sf_erf_c", 5.5, -6.5, 10736.796316979702, 1132.1988540287109 },
    { "sf_erfc_c", 0.3, 0.2, 0.65876251852786138, -0.20852883788276888 },
    { "sf_erfc_c", 1.0, 1.0, -0.31615128169794765, -0.19045346923783468 },
    { "sf_erfc_c", 3.9, -2.0, -1.6139534354328792e-06, -5.8233709909980613e-07 },
    { "sf_erfc_c", 5.0, 1.0, -2.9597765469100242e-12, 2.8460183820855939e-12 },
    { "sf_erfc_c", -1.5, 0.7, 2.0404046154368713, -0.033625498125576171 },
    { "sf_erfc_c", 0.45, -0.3, 0.48570541938551659, 0.28144890596834476 },
    { "sf_erfc_c", 12.0, 3.0, -1.0666962867763503e-60, -1.4047587472713185e-62 },
    { "sf_erfi_c", 0.7, 0.4, 0.7544497636553148, 0.66479328447962172 },
    { "sf_erfi_c", 2.0, -1.5, 0.69951168616312442, 0.10504928977401753 },
    { "sf_erfi_c", -0.3, 3.0, -2.2855283370730961e-05, 1.000007577302835 },
    { "sf_erfi_c", 1.5, 0.0, 4.5847332572844266, 0 },
    { "sf_fresnel_c_c", 0.9, 0.85, 1.6661184635299107, 1.3878808770221371 },
    { "sf_fresnel_s_c", 0.9, 0.85, -0.90411920054662676, 1.1778560388238348 },
    { "sf_fresnel_c_c", 0.95, 0.92, 2.0462507683067148, 1.8118368549191395 },
    { "sf_fresnel_s_c", 0.95, 0.92, -1.3223025944180999, 1.5547939060147387 },
    { "sf_fresnel_c_c", 1.5, 0.5, 0.29034452989915344, -1.0795825961178092 },
    { "sf_fresnel_s_c", 1.5, 0.5, 1.5958995966193614, -0.21646093372162054 },
    { "sf_fresnel_c_c", -2.0, 1.7, -1800.3439463669149, -2120.4568049309264 },
    { "sf_fresnel_s_c", -2.0, 1.7, -2120.956803561895, 1799.8439440718855 },
    { "sf_fresnel_c_c", 3.0, 2.0, 5815899.102940469, -3788100.2001828956 },
    { "sf_fresnel_s_c", 3.0, 2.0, 3788100.7001828952, 5815898.602940469 },
    { "sf_fresnel_c_c", 0.0, 2.0, 0, 0.48825340607534073 },
    { "sf_fresnel_s_c", 0.0, 2.0, 0, -0.34341567836369824 },
    { "sf_fresnel_c_c", 2.5, -2.5, 11017313.787134988, -11017313.787134988 },
    { "sf_fresnel_s_c", 2.5, -2.5, -11017313.287134988, -11017313.287134988 },
    { "sf_sinc_c", 0.6, 0.7, 1.0158745192160936, -0.14171002905512356 },
    { "sf_sinc_c", 0.75, 0.7, 0.97882371400005397, -0.17350579873993832 },
    { "sf_sinc_c", 3.0, -2.0, -0.4298743016327693, 0.91027199557341376 },
    { "sf_sinc_c", -40.0, 5.0, 1.208829319518125, 1.3883273007178918 },
    { "sf_sinc_c", 0.001, 2.0, 1.8134297843998028, -0.00097438263982562237 },
};
static const RRef R_REFS[] = {
    { "sf_erf_r", 0.49, 0.51166826118852304 },
    { "sf_erfc_r", 0.49, 0.4883317388114769 },
    { "sf_erf_r", 0.51, 0.52924361984117041 },
    { "sf_erfc_r", 0.51, 0.47075638015882959 },
    { "sf_erf_r", -2.7, -0.99986566726005943 },
    { "sf_erfc_r", -2.7, 1.9998656672600594 },
    { "sf_erf_r", 5.9, 0.99999999999999989 },
    { "sf_erfc_r", 5.9, 7.1904097835504783e-17 },
    { "sf_erf_r", 6.1, 1 },
    { "sf_erfc_r", 6.1, 6.3146021501937184e-18 },
    { "sf_erf_r", 20.0, 1 },
    { "sf_erfc_r", 20.0, 5.3958656116079012e-176 },
    { "sf_erf_r", 26.6, 1 },
    { "sf_erfc_r", 26.6, 1.0885125885442269e-309 },
    { "sf_erfi_r", 0.49, 0.60053514253298157 },
    { "sf_erfi_r", 0.51, 0.6295139420651148 },
    { "sf_erfi_r", 3.99, 1200593.4751780643 },
    { "sf_erfi_r", 4.01, 1401352.1882798546 },
    { "sf_erfi_r", 9.99, 1.2493843794925596e+42 },
    { "sf_erfi_r", 10.01, 1.8601003045443072e+42 },
    { "sf_erfi_r", -26.6, -4.1328960530520519e+305 },
    { "sf_erfi_r", 26.7, 8.4998672612689847e+307 },
    { "sf_fresnel_c_r", 0.99, 0.7797368571073986 },
    { "sf_fresnel_s_r", 0.99, 0.42826077993226197 },
    { "sf_fresnel_c_r", 1.01, 0.7797358102198938 },
    { "sf_fresnel_s_r", 1.01, 0.44825749017714234 },
    { "sf_fresnel_c_r", 3.99, 0.48845228197350049 },
    { "sf_fresnel_s_r", 3.99, 0.42114272525552493 },
    { "sf_fresnel_c_r", 4.01, 0.50839968557988013 },
    { "sf_fresnel_s_r", 4.01, 0.42114376749673127 },
    { "sf_fresnel_c_r", 7.99, 0.48990702618650817 },
    { "sf_fresnel_s_r", 7.99, 0.46146373697455201 },
    { "sf_fresnel_c_r", 8.01, 0.50969713855789456 },
    { "sf_fresnel_s_r", 8.01, 0.46146476440257278 },
    { "sf_fresnel_c_r", -2.5, -0.45741300964177706 },
    { "sf_fresnel_s_r", -2.5, -0.61918175581959289 },
    { "sf_fresnel_c_r", 100000.0, 0.49999999999999989 },
    { "sf_fresnel_s_r", 100000.0, 0.49999681690113817 },
    { "sf_fresnel_c_r", 123456.789, 0.49999923985222638 },
    { "sf_fresnel_s_r", 123456.789, 0.50000246370825019 },
    { "sf_inverse_erf_r", 0.3, 0.27246271472675432 },
    { "sf_inverse_erf_r", -0.874, -1.0819211672460041 },
    { "sf_inverse_erf_r", 0.876, 1.0876708437181992 },
    { "sf_inverse_erf_r", 0.99, 1.8213863677184494 },
    { "sf_inverse_erf_r", 0.9999999, 3.7665625816384707 },
    { "sf_inverse_erf_r", -0.999999999999, -5.042031898572696 },
    { "sf_inverse_erfc_r", 1e-300, 26.209469960516124 },
    { "sf_inverse_erfc_r", 1e-20, 6.6015806223551428 },
    { "sf_inverse_erfc_r", 0.001, 2.3267537655135246 },
    { "sf_inverse_erfc_r", 0.124, 1.0876708437181992 },
    { "sf_inverse_erfc_r", 0.126, 1.0819211672460041 },
    { "sf_inverse_erfc_r", 0.5, 0.47693627620446988 },
    { "sf_inverse_erfc_r", 1.0, 0 },
    { "sf_inverse_erfc_r", 1.5, -0.47693627620446988 },
    { "sf_inverse_erfc_r", 1.876, -1.0876708437181988 },
    { "sf_inverse_erfc_r", 1.999, -2.3267537655135464 },
    { "sf_sinc_r", 3.0, 0.047040002686622409 },
    { "sf_sinc_r", 1e-08, 1 },
    { "sf_sinc_r", 1000.5, 0.00099477656882080313 },
};

static double ulp_of(double a) {
    a = fabs(a);
    if (a < 2.2250738585072014e-308) return 4.9406564584124654e-324;
    int e;
    (void)frexp(a, &e);                 /* a = m 2^e, m in [1/2, 1) */
    return ldexp(1.0, e - 53);
}

typedef bool (*CK)(double, double, double*, double*);
typedef bool (*RK)(double, double*);

static CK ckernel(const char* k) {
    if (!strcmp(k, "sf_erf_c")) return sf_erf_c;
    if (!strcmp(k, "sf_erfc_c")) return sf_erfc_c;
    if (!strcmp(k, "sf_erfi_c")) return sf_erfi_c;
    if (!strcmp(k, "sf_fresnel_c_c")) return sf_fresnel_c_c;
    if (!strcmp(k, "sf_fresnel_s_c")) return sf_fresnel_s_c;
    if (!strcmp(k, "sf_sinc_c")) return sf_sinc_c;
    return NULL;
}
static RK rkernel(const char* k) {
    if (!strcmp(k, "sf_erf_r")) return sf_erf_r;
    if (!strcmp(k, "sf_erfc_r")) return sf_erfc_r;
    if (!strcmp(k, "sf_erfi_r")) return sf_erfi_r;
    if (!strcmp(k, "sf_fresnel_c_r")) return sf_fresnel_c_r;
    if (!strcmp(k, "sf_fresnel_s_r")) return sf_fresnel_s_r;
    if (!strcmp(k, "sf_inverse_erf_r")) return sf_inverse_erf_r;
    if (!strcmp(k, "sf_inverse_erfc_r")) return sf_inverse_erfc_r;
    if (!strcmp(k, "sf_sinc_r")) return sf_sinc_r;
    return NULL;
}

void test_kernel_spot_values_complex() {
    for (size_t i = 0; i < sizeof C_REFS / sizeof C_REFS[0]; i++) {
        const CRef* c = &C_REFS[i];
        double r, m;
        ASSERT_MSG(ckernel(c->k)(c->x, c->y, &r, &m), "%s(%g, %g) declined", c->k, c->x, c->y);
        double err = hypot(r - c->re, m - c->im);
        double mag = hypot(c->re, c->im);
        ASSERT_MSG(err <= 4.0 * ulp_of(mag), "%s(%g, %g): got %.17g + %.17g i, want %.17g + %.17g i (%.2f ulp)",
                   c->k, c->x, c->y, r, m, c->re, c->im, err / ulp_of(mag));
    }
}

void test_kernel_spot_values_real() {
    for (size_t i = 0; i < sizeof R_REFS / sizeof R_REFS[0]; i++) {
        const RRef* c = &R_REFS[i];
        double v;
        ASSERT_MSG(rkernel(c->k)(c->x, &v), "%s(%g) declined", c->k, c->x);
        double err = fabs(v - c->v);
        ASSERT_MSG(err <= 4.0 * ulp_of(c->v), "%s(%.17g): got %.17g, want %.17g (%.2f ulp)",
                   c->k, c->x, v, c->v, err / ulp_of(c->v));
    }
}

/* Near the axes the SMALL component must be right relative to itself, not
 * only relative to |f| (references at 400 digits: at 40, mpmath's own
 * Re erf(1e-300 + 2i) is wrong by 240 orders of magnitude). */
static const CRef AXIS_REFS[] = {
    { "sf_erf_c", 1e-300, 2.0, 6.1607415059355131e-299, 18.564802414575553 },
    { "sf_erf_c", 1.5, 1e-10, 0.96610514647531076, 1.1893028922362937e-11 },
    { "sf_erfc_c", 2.5, -1e-09, 0.00040695201744495892, 2.1782842303527099e-12 },
    { "sf_erfi_c", 3e-07, 1.2, 8.0203304101057222e-08, 0.9103139782296642 },
    { "sf_fresnel_c_c", 3.1, 1e-08, 0.56159390251135666, -8.1814971742502516e-09 },
    { "sf_fresnel_s_c", -2e-09, 2.2, 1.937166322257263e-09, -0.45570461212465707 },
    { "sf_erf_c", 0.05, 3.0, 449.97036886788084, 1562.0170564668124 },
};

void test_axis_components() {
    for (size_t i = 0; i < sizeof AXIS_REFS / sizeof AXIS_REFS[0]; i++) {
        const CRef* c = &AXIS_REFS[i];
        double r, m;
        ASSERT_MSG(ckernel(c->k)(c->x, c->y, &r, &m), "%s(%g, %g) declined", c->k, c->x, c->y);
        ASSERT_MSG(fabs(r - c->re) <= 4.0 * ulp_of(c->re) && fabs(m - c->im) <= 4.0 * ulp_of(c->im),
                   "%s(%g, %g): got %.17g + %.17g i, want %.17g + %.17g i (componentwise)",
                   c->k, c->x, c->y, r, m, c->re, c->im);
    }
    assert_eval_startswith("Re[Erf[1.*^-300 + 2. I]]", "6.16074");
}

/* Exact special points the kernels must answer exactly. */
void test_kernel_exact_points() {
    double v, r, m;
    ASSERT(sf_erf_r(0.0, &v) && v == 0.0);
    ASSERT(sf_erfc_r(0.0, &v) && v == 1.0);
    ASSERT(sf_erfc_r(30.0, &v) && v == 0.0);            /* underflow, not a decline */
    ASSERT(sf_erfc_r(-30.0, &v) && v == 2.0);
    ASSERT(sf_erfi_r(0.0, &v) && v == 0.0);
    ASSERT(sf_fresnel_c_r(0.0, &v) && v == 0.0);
    ASSERT(sf_fresnel_s_r(1e300, &v) && v == 0.5);        /* oscillation below 1/2 ulp */
    ASSERT(sf_sinc_r(0.0, &v) && v == 1.0);
    ASSERT(sf_inverse_erf_r(0.0, &v) && v == 0.0);
    ASSERT(sf_inverse_erfc_r(1.0, &v) && v == 0.0);
    ASSERT(sf_sinc_c(0.0, 0.0, &r, &m) && r == 1.0 && m == 0.0);
    ASSERT(sf_erf_c(0.0, 0.0, &r, &m) && r == 0.0 && m == 0.0);
    /* odd / conjugate symmetry hold exactly */
    double a, b;
    ASSERT(sf_erf_c(1.3, -0.4, &r, &m) && sf_erf_c(-1.3, 0.4, &a, &b) && a == -r && b == -m);
    ASSERT(sf_fresnel_c_c(2.1, 0.7, &r, &m) && sf_fresnel_c_c(2.1, -0.7, &a, &b) && a == r && b == -m);
}

/* Declines: the kernel answers false; the builtin then still answers (MPFR). */
void test_kernel_declines() {
    double v, r, m;
    ASSERT(!sf_inverse_erf_r(1.0, &v));
    ASSERT(!sf_inverse_erf_r(-1.0, &v));
    ASSERT(!sf_inverse_erf_r(1.5, &v));
    ASSERT(!sf_inverse_erfc_r(0.0, &v));
    ASSERT(!sf_inverse_erfc_r(2.0, &v));
    ASSERT(!sf_inverse_erfc_r(1e-310, &v));                /* subnormal */
    ASSERT(!sf_erfi_r(27.0, &v));                          /* overflows */
    ASSERT(!sf_erf_c(0.3, 30.0, &r, &m));                  /* e^{y^2} overflows */
    ASSERT(!sf_erfc_c(30.0, 1.0, &r, &m));
    ASSERT(!sf_sinc_c(1.0, 800.0, &r, &m));
    ASSERT(!sf_erf_r(NAN, &v));
    /* ... and the builtins fall through to arbitrary precision. */
    assert_eval_startswith("InverseErfc[1.*^-310]", "26.64");
    assert_eval_eq("InverseErf[1.]", "Infinity", 0);
    assert_eval_eq("InverseErfc[2.]", "-Infinity", 0);
    assert_eval_eq("InverseErf[1.5]", "InverseErf[1.5]", 0);
    assert_eval_startswith("Erfi[27.]", "8.33752");           /* MPFR, beyond double range */
    /* |z|^2 beyond any MPFR guard precision: erf z = 1 - erfc z with erfc from
     * arb (never an MPFR abort). erfc z ~ e^{-z^2}/(z Sqrt[Pi]) is ~4e-301 here,
     * so erf z is 1 plus a tiny COMPLEX part -- a Real 1.0 dropped it. */
    assert_eval_eq("Head[Erf[1.*^300 + 1.*^300 I]]", "Complex", 0);
    assert_eval_eq("Erf[1.*^300 + 1.*^300 I] == 1 - Erfc[1.*^300 + 1.*^300 I]", "True", 0);
    assert_eval_eq("Erf[-1.*^300 - 1.*^300 I] == -Erf[1.*^300 + 1.*^300 I]", "True", 0);
    assert_eval_eq("0 < Abs[Im[Erf[1.*^300 + 1.*^300 I]]] < 1.*^-300", "True", 0);
    /* huge |z| next to the real axis stays on the machine axis series */
    assert_eval_startswith("Im[FresnelC[1.*^10 + 1.*^-10 I]]", "3.67608");
    assert_eval_startswith("Re[FresnelC[1.*^-8 + 3.*^7 I]]", "1.15476");
}

/* scalar == unpacked List == packed List == NDArray == Compile (scalar and
 * rank-1), bit for bit, real and complex. */
static void check_surfaces(const char* head, const char* data, const char* type) {
    char buf[4096];
    snprintf(buf, sizeof buf,
        "Module[{xs = %s, sc, cf, cfa}, sc = %s /@ xs;"
        " cf = Compile[{{x, %s}}, %s[x]]; cfa = Compile[{{v, %s, 1}}, %s[v]];"
        " {Max[Abs[sc - %s[FromPackedArray[xs]]]], Max[Abs[sc - %s[ToPackedArray[xs]]]],"
        "  Max[Abs[sc - FromNDArray[%s[ToNDArray[xs]]]]], Max[Abs[sc - (cf /@ xs)]],"
        "  Max[Abs[sc - cfa[xs]]]}]",
        data, head, type, head, type, head, head, head, head);
    assert_eval_eq(buf, "{0.0, 0.0, 0.0, 0.0, 0.0}", 0);
}

void test_five_surfaces() {
    const char* re = "(SeedRandom[7]; RandomReal[{-3, 3}, 300])";
    const char* cx = "(SeedRandom[8]; RandomComplex[{-2.5 - 2.5 I, 2.5 + 2.5 I}, 300])";
    const char* heads[] = { "Erf", "Erfc", "Erfi", "FresnelC", "FresnelS", "Sinc" };
    for (size_t i = 0; i < 6; i++) {
        check_surfaces(heads[i], re, "_Real");
        check_surfaces(heads[i], cx, "_Complex");
    }
    check_surfaces("InverseErf", "(SeedRandom[9]; RandomReal[{-0.9999, 0.9999}, 300])", "_Real");
    check_surfaces("InverseErfc", "(SeedRandom[9]; RandomReal[{0.0001, 1.9999}, 300])", "_Real");
}

void test_compile_lowering() {
    const char* heads[] = { "Erf", "Erfc", "Erfi", "FresnelC", "FresnelS", "Sinc" };
    char buf[512];
    for (size_t i = 0; i < 6; i++) {
        snprintf(buf, sizeof buf,
            "Lookup[#, \"Compiled\"] & /@ {CompileDiagnostics[{{x, _Real}}, %s[x]],"
            " CompileDiagnostics[{{v, _Real, 1}}, %s[v]], CompileDiagnostics[{{z, _Complex}}, %s[z]],"
            " CompileDiagnostics[{{w, _Complex, 1}}, %s[w]]}", heads[i], heads[i], heads[i], heads[i]);
        assert_eval_eq(buf, "{True, True, True, True}", 0);
    }
    assert_eval_eq("Lookup[#, \"Compiled\"] & /@ {CompileDiagnostics[{{x, _Real}}, InverseErf[x]],"
                   " CompileDiagnostics[{{v, _Real, 1}}, InverseErfc[v]]}", "{True, True}", 0);
}

/* A machine-complex argument gives a machine-complex result. */
void test_scalar_shapes() {
    assert_eval_eq("Head /@ {Erf[0.5 + 0.7 I], Erfc[1.2 - 0.3 I], Erfi[0.3 + 2. I],"
                   " FresnelC[1.2 + 0.8 I], FresnelS[-0.4 + 1.9 I], Sinc[1. + 2. I]}",
                   "{Complex, Complex, Complex, Complex, Complex, Complex}", 0);
    assert_eval_eq("Precision /@ {Erf[0.5 + 0.7 I], FresnelS[2.5], InverseErfc[0.001]}",
                   "{MachinePrecision, MachinePrecision, MachinePrecision}", 0);
}

int main() {
    symtab_init();
    core_init();

    TEST(test_kernel_spot_values_complex);
    TEST(test_kernel_spot_values_real);
    TEST(test_axis_components);
    TEST(test_kernel_exact_points);
    TEST(test_kernel_declines);
    TEST(test_five_surfaces);
    TEST(test_compile_lowering);
    TEST(test_scalar_shapes);

    printf("All sf_machine_error tests passed.\n");
    return 0;
}
