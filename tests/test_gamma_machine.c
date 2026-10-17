/* Tests for the Gamma-family machine kernels (src/special_functions/sf_machine_gamma.c)
 * and their wiring into Gamma, LogGamma, PolyGamma, Beta and Pochhammer.
 *
 * Covers:
 *   - spot values of every kernel against mpmath (40 digits), as a relative error
 *     in units of 2^-53 with a per-case budget: 8 ulp in the tiled / near-real
 *     regimes, more where the function itself is ill-conditioned (|z psi(z)| large,
 *     e.g. Gamma[30 + 40 I]), where the budget is the condition number's;
 *   - the regime boundaries (digamma at 1, 2, 10; the tile edges Re z = 1/2, 5/2,
 *     10 and Im z = 19/4; the root disk) by continuity across each seam;
 *   - the declines (poles, n >= 1 in the left half-plane, out-of-range order),
 *     and that the builtin then still answers through its older path;
 *   - the scalar builtin results and their heads (Real vs Complex), including the
 *     one-argument PolyGamma, which now computes directly;
 *   - scalar == List == packed List == NDArray == Compile, bit for bit, real and
 *     complex, for all five heads;
 *   - LogGamma's branch: continuity across the reflection region and the
 *     conjugate symmetry off the negative real axis.
 *
 * Reference values: mpmath 1.3 at 40 digits (tools/sf_accuracy/gamma.py holds the
 * full >= 2000-point sweeps). */

#include "eval.h"
#include "parse.h"
#include "expr.h"
#include "symtab.h"
#include "attr.h"
#include "core.h"
#include "print.h"
#include "test_utils.h"
#include "sf_machine_gamma.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define ULP (1.0 / 9007199254740992.0)   /* 2^-53 */

/* ---- reference tables (mpmath, 40 digits): inputs, value, budget in ulp ---- */

static const double DIG_R[][3] = {
    { 0.37, -2.79530141089056405e+00, 8 },
    { 1.0, -5.77215664901532866e-01, 8 },
    { 1.25, -2.27453533376265421e-01, 8 },
    { 2.0, 4.22784335098467134e-01, 8 },
    { 2.5, 7.03156640645243192e-01, 8 },
    { 7.9, 2.00223848756357103e+00, 8 },
    { 9.999, 2.25164741720573547e+00, 8 },
    { 10.0, 2.25175258906672093e+00, 8 },
    { 12.0, 2.44266167997581185e+00, 8 },
    { 40.0, 3.67632737403484322e+00, 8 },
    { -2.5, 1.10315664064524310e+00, 8 },
    { -10.7, 1.33744812447064432e-01, 8 },
    { 0.001, -1.00057557193181026e+03, 8 },
};
static const double PG_R[][4] = {
    { 1, 0.37, 8.36047382779909753e+00, 8 },
    { 1, 30.5, 3.33302481128016082e-02, 8 },
    { 2, 0.37, -4.05303269975773901e+01, 8 },
    { 2, 2.5, -2.36204051641727392e-01, 8 },
    { 2, 7.9, -1.81790048383226044e-02, 8 },
    { 2, 12.0, -7.54720536899891169e-03, 8 },
    { 3, 0.9, 9.73901409377672778e+00, 8 },
    { 3, 27.7, 9.93184710741578178e-05, 8 },
    { 5, 1.5, 1.11135486024354968e+01, 8 },
    { 12, 3.25, -1.09612755215217831e+02, 8 },
};
static const double DIG_C[][5] = {
    { 0.6, 0.3, -1.17645876215802669e+00, 9.19347064506384593e-01, 8 },
    { 1.2, -2.7, 1.02103402271089383e+00, -1.31440857046995752e+00, 8 },
    { 3.9, 2.9, 1.49746181948375479e+00, 7.04144857704147031e-01, 8 },
    { -2.5, 0.7, 1.12900820397203550e+00, 2.83790499786108175e+00, 8 },
    { 1.46, 0.01, -1.53616048313680699e-03, 9.69093665299915544e-03, 8 },
    { 20.0, 15.0, 3.20283847355883822e+00, 6.55629097313798703e-01, 8 },
    { 0.5, 4.7, 1.54566096711156287e+00, 1.57079632679442671e+00, 8 },
    { 9.99, 4.74, 2.36182228179074150e+00, 4.62930540052214567e-01, 8 },
    { 0.75, 0.0001, -1.08586085327330650e+00, 2.54187961473275302e-04, 8 },
    { -40.3, 2.0, 3.70992778847961757e+00, 3.09260792324568401e+00, 8 },
};
static const double PG_C[][6] = {
    { 1, 0.7, 0.6, 8.37755540335082149e-01, -1.46651711104625160e+00, 8 },
    { 2, 0.71, 0.64, 1.03676940435073250e+00, 2.23300644727506725e+00, 8 },
    { 2, 22.0, -20.0, -8.40611346907010231e-05, -1.15668147465857103e-03, 8 },
    { 2, 3.9, 2.9, -8.48523908530176771e-03, 4.92373186520455933e-02, 8 },
    { 3, 0.2, 5.0, 2.94419435204701938e-03, 1.59712103180844646e-02, 8 },
    { 7, 1.5, -0.5, -1.08705141646149755e+02, 7.24010249745260666e+01, 8 },
};
static const double GAMMA_C[][5] = {
    { 0.6, 0.3, 1.16913001881741807e+00, -5.25917107255656657e-01, 8 },
    { 1.2, -2.7, 3.88396185823712062e-02, -6.13124366625638420e-02, 8 },
    { 3.9, 2.9, -1.32299370609796707e+00, -1.13060346015856794e+00, 8 },
    { -4.5, 3.0, -6.32915482239325550e-06, 2.17712585788771438e-05, 32 },
    { 12.0, 1.0, -2.92906378716227710e+07, 2.45550232059890442e+07, 16 },
    { 0.5, 10.0, 3.37872437623423601e-07, 1.68936983903891904e-07, 64 },
    { 2.0001, 0.0001, 1.00004227843334670e+00, 4.22866704796091036e-05, 8 },
    { 0.75, 4.74, -2.14598008174196131e-03, 2.39357509049717003e-04, 8 },
    { 9.99, -4.76, -1.08626705339641430e+04, 1.12124860997633165e+05, 64 },
    { 30.0, 40.0, 1.87419976730378030e+21, -1.51084450333286780e+21, 256 },
};
static const double LGAMMA_C[][5] = {
    { 0.6, 0.3, 2.48400002894635208e-01, -4.22717774331417662e-01, 8 },
    { 1.2, -2.7, -2.62307779289511522e+00, -1.00612905614266235e+00, 8 },
    { 3.9, 2.9, 5.54046358609709744e-01, 3.84873942279282932e+00, 8 },
    { -4.5, 3.0, -1.06943542765744599e+01, -1.07126607034147359e+01, 8 },
    { -4.5, -3.0, -1.06943542765744599e+01, 1.07126607034147359e+01, 8 },
    { 1.4988, -0.001, -1.20825820206447027e-01, -3.53677526132009767e-05, 8 },
    { 2.0385, -0.0002, 1.67513621003261663e-02, -8.94638837231431463e-05, 8 },
    { 1.0001, 0.0001, -5.77215656888838838e-05, -5.77051179508560700e-05, 8 },
    { 0.5, 30.0, -4.62049512706422263e+01, 7.20373104288057959e+01, 8 },
    { -25.5, 0.5, -5.93982648063846099e+01, -8.00522991100169605e+01, 8 },
    { 100.0, 200.0, 2.17773331290709280e+02, 9.92138792528162526e+02, 8 },
};
static const double BETA_R[][4] = {
    { 2.5, 1.7, 1.55722381342194183e-01, 8 },
    { -2.3, 1.1, -2.83801050746294359e-01, 8 },
    { 0.001, 4.5, 9.98036541235101936e+02, 8 },
    { 0.3, 0.3, 6.00962368373101530e+00, 8 },
    { 4.24, 0.864, 3.19873536726145102e-01, 8 },
    { 150.5, 10.25, 2.24891414496915626e-17, 8 },
};
static const double BETA_C[][7] = {
    { 1.37, 1.38, 2.81, -1.94, 2.76238328007806025e-02, -3.45934425319905994e-02, 8 },
    { 0.6, 0.3, 1.3, -0.4, 1.03906537719663228e+00, -3.57735060526675641e-01, 8 },
    { 2.5, 1.0, 1.5, 0.0, 1.53920094207259356e-01, -8.99803248902973085e-02, 8 },
    { -3.5, 2.0, 1.25, -0.5, -3.31240507631904851e-02, -2.09958321331376978e-02, 48 },
};
static const double POCH_R[][4] = {
    { 2.5, 3.0, 3.93750000000000000e+01, 8 },
    { 0.37, 3.0, 1.20135299999999989e+00, 8 },
    { -15.18, -18.0, 1.31165189011346989e-25, 8 },
    { 2.5, 2.5, 1.80540666735282009e+01, 8 },
    { -16.9, -11.6, 5.05123415611734479e-17, 8 },
    { 7.25, -4.0, 2.20642102995044151e-03, 8 },
    { 1.5, 40.0, 5.87716684727422530e+48, 8 },
};
static const double POCH_C[][7] = {
    { 1.2, -1.5, 3.0, 0.0, -6.40200000000000014e+00, -1.69049999999999976e+01, 8 },
    { 2.0, 1.0, 0.5, 0.2, 1.17499972120003870e+00, 5.28043395248772085e-01, 8 },
    { 0.6, 0.3, -2.0, 0.0, 9.17073170731707354e-01, 1.05365853658536568e+00, 8 },
    { 3.9, 2.9, 2.5, 0.0, 1.11281665249494335e+01, 6.88588920599773360e+01, 8 },
};

static void check_rel_r(const char* what, double got, double ref, double tol_ulp) {
    double e = fabs(got - ref) / fabs(ref) / ULP;
    ASSERT_MSG(e <= tol_ulp, "%s: got %.17g, want %.17g (%.2f ulp > %.0f)", what, got, ref, e, tol_ulp);
}
static void check_rel_c(const char* what, double gr, double gi, double rr, double ri, double tol_ulp) {
    double e = hypot(gr - rr, gi - ri) / hypot(rr, ri) / ULP;
    ASSERT_MSG(e <= tol_ulp, "%s: got %.17g%+.17gi, want %.17g%+.17gi (%.2f ulp > %.0f)",
               what, gr, gi, rr, ri, e, tol_ulp);
}
#define N_OF(a) (sizeof(a) / sizeof((a)[0]))

void test_kernel_spot_values(void) {
    char buf[96];
    for (size_t i = 0; i < N_OF(DIG_R); i++) {
        double v; snprintf(buf, sizeof buf, "digamma(%g)", DIG_R[i][0]);
        ASSERT_MSG(sfg_digamma(DIG_R[i][0], &v), "%s declined", buf);
        check_rel_r(buf, v, DIG_R[i][1], DIG_R[i][2]);
    }
    for (size_t i = 0; i < N_OF(PG_R); i++) {
        double v; snprintf(buf, sizeof buf, "polygamma(%g, %g)", PG_R[i][0], PG_R[i][1]);
        ASSERT_MSG(sfg_polygamma((int)PG_R[i][0], PG_R[i][1], &v), "%s declined", buf);
        check_rel_r(buf, v, PG_R[i][2], PG_R[i][3]);
    }
    for (size_t i = 0; i < N_OF(DIG_C); i++) {
        double r, s; snprintf(buf, sizeof buf, "digamma(%g%+gi)", DIG_C[i][0], DIG_C[i][1]);
        ASSERT_MSG(sfg_polygamma_c(0, DIG_C[i][0], DIG_C[i][1], &r, &s), "%s declined", buf);
        check_rel_c(buf, r, s, DIG_C[i][2], DIG_C[i][3], DIG_C[i][4]);
    }
    for (size_t i = 0; i < N_OF(PG_C); i++) {
        double r, s;
        snprintf(buf, sizeof buf, "polygamma(%g, %g%+gi)", PG_C[i][0], PG_C[i][1], PG_C[i][2]);
        ASSERT_MSG(sfg_polygamma_c((int)PG_C[i][0], PG_C[i][1], PG_C[i][2], &r, &s), "%s declined", buf);
        check_rel_c(buf, r, s, PG_C[i][3], PG_C[i][4], PG_C[i][5]);
    }
    for (size_t i = 0; i < N_OF(GAMMA_C); i++) {
        double r, s; snprintf(buf, sizeof buf, "gamma(%g%+gi)", GAMMA_C[i][0], GAMMA_C[i][1]);
        ASSERT_MSG(sfg_gamma_c(GAMMA_C[i][0], GAMMA_C[i][1], &r, &s), "%s declined", buf);
        check_rel_c(buf, r, s, GAMMA_C[i][2], GAMMA_C[i][3], GAMMA_C[i][4]);
    }
    for (size_t i = 0; i < N_OF(LGAMMA_C); i++) {
        double r, s; snprintf(buf, sizeof buf, "loggamma(%g%+gi)", LGAMMA_C[i][0], LGAMMA_C[i][1]);
        ASSERT_MSG(sfg_loggamma_c(LGAMMA_C[i][0], LGAMMA_C[i][1], &r, &s), "%s declined", buf);
        check_rel_c(buf, r, s, LGAMMA_C[i][2], LGAMMA_C[i][3], LGAMMA_C[i][4]);
    }
    for (size_t i = 0; i < N_OF(BETA_R); i++) {
        double v; snprintf(buf, sizeof buf, "beta(%g, %g)", BETA_R[i][0], BETA_R[i][1]);
        if (fabs(BETA_R[i][0]) > 170.0) {         /* |a| > 170: declines by design */
            ASSERT_MSG(!sfg_beta(BETA_R[i][0], BETA_R[i][1], &v), "%s should decline", buf);
            continue;
        }
        ASSERT_MSG(sfg_beta(BETA_R[i][0], BETA_R[i][1], &v), "%s declined", buf);
        check_rel_r(buf, v, BETA_R[i][2], BETA_R[i][3]);
    }
    for (size_t i = 0; i < N_OF(BETA_C); i++) {
        double r, s;
        snprintf(buf, sizeof buf, "beta(%g%+gi, %g%+gi)", BETA_C[i][0], BETA_C[i][1], BETA_C[i][2], BETA_C[i][3]);
        ASSERT_MSG(sfg_beta_c(BETA_C[i][0], BETA_C[i][1], BETA_C[i][2], BETA_C[i][3], &r, &s), "%s declined", buf);
        check_rel_c(buf, r, s, BETA_C[i][4], BETA_C[i][5], BETA_C[i][6]);
    }
    for (size_t i = 0; i < N_OF(POCH_R); i++) {
        double v; snprintf(buf, sizeof buf, "pochhammer(%g, %g)", POCH_R[i][0], POCH_R[i][1]);
        ASSERT_MSG(sfg_pochhammer(POCH_R[i][0], POCH_R[i][1], &v), "%s declined", buf);
        check_rel_r(buf, v, POCH_R[i][2], POCH_R[i][3]);
    }
    for (size_t i = 0; i < N_OF(POCH_C); i++) {
        double r, s;
        snprintf(buf, sizeof buf, "pochhammer(%g%+gi, %g%+gi)", POCH_C[i][0], POCH_C[i][1], POCH_C[i][2], POCH_C[i][3]);
        ASSERT_MSG(sfg_pochhammer_c(POCH_C[i][0], POCH_C[i][1], POCH_C[i][2], POCH_C[i][3], &r, &s), "%s declined", buf);
        check_rel_c(buf, r, s, POCH_C[i][4], POCH_C[i][5], POCH_C[i][6]);
    }
}

/* Continuity across every seam: the value just below a boundary and just above
 * it differ by no more than the derivative predicts plus a few ulp. Both sides
 * come from different code paths (fit vs recurrence, tile vs Stirling, root
 * series vs tile), so a coefficient or branch slip shows as a jump. */
static void seam_r(const char* what, double x) {
    double lo = nextafter(x, -INFINITY), hi = nextafter(x, INFINITY);
    double a, b, c, d;
    ASSERT(sfg_digamma(lo, &a) && sfg_digamma(x, &b) && sfg_digamma(hi, &c));
    ASSERT(sfg_polygamma(1, x, &d));                     /* the slope psi'(x) */
    double scale = fabs(b) > 1.0 ? fabs(b) : 1.0;
    double tol = 8 * ULP * scale;
    ASSERT_MSG(fabs(a - b) <= tol + d * (x - lo) && fabs(c - b) <= tol + d * (hi - x),
               "%s seam at %.17g: %.17g %.17g %.17g", what, x, a, b, c);
}
typedef bool (*cfun)(double, double, double*, double*);
static bool lg_c(double r, double i, double* a, double* b) { return sfg_loggamma_c(r, i, a, b); }
static bool g_c(double r, double i, double* a, double* b) { return sfg_gamma_c(r, i, a, b); }
static bool d_c(double r, double i, double* a, double* b) { return sfg_polygamma_c(0, r, i, a, b); }
static bool p2_c(double r, double i, double* a, double* b) { return sfg_polygamma_c(2, r, i, a, b); }
static void seam_c(const char* what, cfun f, double zr, double zi, int along_im) {
    double ar, ai, br, bi, pr, pi, mr, mi;
    double lr = along_im ? zr : nextafter(zr, -INFINITY), li = along_im ? nextafter(zi, -INFINITY) : zi;
    ASSERT(f(lr, li, &ar, &ai) && f(zr, zi, &br, &bi));
    /* the step itself moves f by |f'| |dz|: estimate f' by a central difference */
    const double h = 1e-5;
    ASSERT(f(along_im ? zr : zr + h, along_im ? zi + h : zi, &pr, &pi));
    ASSERT(f(along_im ? zr : zr - h, along_im ? zi - h : zi, &mr, &mi));
    double dfdz = hypot(pr - mr, pi - mi) / (2.0 * h);
    double dz = along_im ? zi - li : zr - lr;
    double jump = hypot(ar - br, ai - bi) - dfdz * dz;
    double e = jump / hypot(br, bi) / ULP;
    ASSERT_MSG(e <= 12.0, "%s seam at %g%+gi: jump of %.2f ulp beyond the slope", what, zr, zi, e);
}

void test_kernel_seams(void) {
    seam_r("digamma", 1.0);   seam_r("digamma", 2.0);   seam_r("digamma", 10.0);
    seam_r("digamma", 0.25);  seam_r("digamma", 3.0);
    const double ims[] = { 0.1, 0.9, 2.2, 4.0 };
    for (size_t k = 0; k < N_OF(ims); k++) {
        double y = ims[k];
        seam_c("loggamma Re=1/2", lg_c, 0.5, y, 0);  seam_c("loggamma Re=1", lg_c, 1.0, y, 0);
        seam_c("loggamma Re=5/2", lg_c, 2.5, y, 0);  seam_c("loggamma Re=10", lg_c, 10.0, y, 0);
        seam_c("gamma Re=1/2", g_c, 0.5, y, 0);      seam_c("gamma Re=10", g_c, 10.0, y, 0);
        seam_c("gamma Re=3/2", g_c, 1.5, y, 0);      seam_c("gamma Re=7", g_c, 7.0, y, 0);
        seam_c("digamma Re=1/2", d_c, 0.5, y, 0);    seam_c("digamma Re=10", d_c, 10.0, y, 0);
        seam_c("digamma Re=0", d_c, 0.0, y, 0);      seam_c("polygamma2 Re=11", p2_c, 11.0, y, 0);
    }
    /* the tile rows (Im = 1/4 + k/2) and the tile / Stirling edge Im = 19/4 */
    const double res[] = { 0.7, 1.3, 2.0, 3.1, 6.6 };
    for (size_t k = 0; k < N_OF(res); k++) {
        seam_c("loggamma Im=1/4", lg_c, res[k], 0.25, 1);  seam_c("loggamma Im=19/4", lg_c, res[k], 4.75, 1);
        seam_c("gamma Im=19/4", g_c, res[k], 4.75, 1);     seam_c("digamma Im=19/4", d_c, res[k], 4.75, 1);
        seam_c("gamma Im=3/4", g_c, res[k], 0.75, 1);      seam_c("digamma Im=1/4", d_c, res[k], 0.25, 1);
    }
    /* psi's root series (|z - x0| <= 0.35) against the tile next to it */
    seam_c("digamma root disk", d_c, 1.4616321449683622 + 0.35, 0.0001, 0);
}

void test_kernel_declines(void) {
    double v, r, s;
    ASSERT(!sfg_digamma(0.0, &v));
    ASSERT(!sfg_digamma(-3.0, &v));
    ASSERT(!sfg_digamma(NAN, &v));
    ASSERT(!sfg_polygamma(2, -1.5, &v));               /* n >= 1, negative axis */
    ASSERT(!sfg_polygamma(2, 0.0, &v));
    ASSERT(!sfg_polygamma(SFG_POLYGAMMA_MAX_N + 1, 2.0, &v));
    ASSERT(!sfg_polygamma_c(2, -0.5, 1.0, &r, &s));    /* n >= 1, left half-plane */
    ASSERT(!sfg_polygamma_c(0, -2.0, 0.0, &r, &s));     /* pole */
    ASSERT(!sfg_gamma_c(-3.0, 0.0, &r, &s));            /* pole */
    ASSERT(!sfg_gamma_c(200.0, 1.0, &r, &s));           /* overflow */
    ASSERT(!sfg_gamma_c(-40.5, 1.0, &r, &s));           /* left of the recurrence */
    ASSERT(!sfg_loggamma_c(-2.5, 0.0, &r, &s));         /* the cut: caller's real path */
    ASSERT(!sfg_loggamma_c(-40.5, 1.0, &r, &s));
    ASSERT(!sfg_beta(-1.0, 2.0, &v));                   /* the exact pole lattice */
    ASSERT(!sfg_beta(-0.5, -0.5, &v));                  /* 1/Gamma(a+b) = 0 */
    ASSERT(!sfg_pochhammer(-2.0, 0.5, &v));
    ASSERT(!sfg_pochhammer(2.0, -3.0, &v));             /* (1)(0)(-1) in the denominator */
    ASSERT(!sfg_pochhammer_c(1.0, 0.0, -2.0, 0.0, &r, &s));
}

/* ---- the builtins ----------------------------------------------------- */

static Expr* eval_str(const char* input) {
    Expr* e = parse_expression(input);
    ASSERT(e != NULL);
    Expr* r = evaluate(e);
    expr_free(e);
    return r;
}
static void eval_str_free(const char* input) { expr_free(eval_str(input)); }
static void assert_true(const char* input) {
    Expr* r = eval_str(input);
    ASSERT_MSG(r->type == EXPR_SYMBOL && strcmp(r->data.symbol.name, "True") == 0,
               "%s: expected True", input);
    expr_free(r);
}

void test_builtin_results(void) {
    /* machine real in, Real out; machine complex in, Complex out */
    assert_true("Head[PolyGamma[2.5]] === Real");
    assert_true("Head[PolyGamma[0, 2.5]] === Real");
    assert_true("PolyGamma[2.5] === PolyGamma[0, 2.5]");
    assert_true("Head[PolyGamma[2, 2.5]] === Real");
    assert_true("Head[PolyGamma[1.5 + 0.5 I]] === Complex");
    assert_true("Head[PolyGamma[2, 1.5 + 0.5 I]] === Complex");
    assert_true("Head[Beta[2.5, 1.5]] === Real");
    assert_true("Head[Beta[2.5 + I, 1.5]] === Complex");
    assert_true("Head[Pochhammer[2.5, 3]] === Real");
    assert_true("Head[Pochhammer[2.5 + I, 3]] === Complex");
    assert_true("Head[Gamma[2.5 + 1. I]] === Complex");
    assert_true("Head[LogGamma[2.5 + 1. I]] === Complex");
    /* exact values exactly as before */
    assert_eval_eq("Pochhammer[2.5, 0]", "1", 0);
    assert_true("Pochhammer[x, 3] === x (1 + x) (2 + x)");
    assert_eval_eq("Beta[2, 3]", "1/12", 0);
    assert_eval_eq("PolyGamma[x]", "PolyGamma[0, x]", 0);
    assert_true("PolyGamma[1, 1] === Pi^2/6");
    /* a pole still answers through the layered path */
    assert_eval_eq("PolyGamma[0, 0.]", "ComplexInfinity", 0);
    assert_eval_eq("Gamma[-2.]", "ComplexInfinity", 0);
    /* what the kernel declines still evaluates */
    assert_true("NumberQ[Gamma[-40.5 + 1. I]]");
    assert_true("NumberQ[LogGamma[-40.5 + 1. I]]");
    assert_true("NumberQ[PolyGamma[2, -1.5]]");
    assert_true("NumberQ[PolyGamma[2, -0.5 + 1. I]]");
    /* the scalar builtin returns the kernel's bits */
    {
        Expr* r = eval_str("PolyGamma[2, 2.5]");
        double v; ASSERT(sfg_polygamma(2, 2.5, &v));
        ASSERT(r->type == EXPR_REAL && r->data.real == v);
        expr_free(r);
    }
    {
        Expr* r = eval_str("Gamma[Complex[0.6, 0.3]]");
        double gr, gi; ASSERT(sfg_gamma_c(0.6, 0.3, &gr, &gi));
        ASSERT(r->type == EXPR_FUNCTION && r->data.function.arg_count == 2);
        ASSERT(r->data.function.args[0]->data.real == gr && r->data.function.args[1]->data.real == gi);
        expr_free(r);
    }
}

/* LogGamma's branch: the continuous sheet, conjugate-symmetric, agreeing with the
 * real-axis convention Im = -Pi Ceiling[-x] as Im z -> 0+. */
void test_loggamma_branch(void) {
    double r1, i1, r2, i2;
    ASSERT(sfg_loggamma_c(-4.5, 3.0, &r1, &i1) && sfg_loggamma_c(-4.5, -3.0, &r2, &i2));
    ASSERT(r1 == r2 && i1 == -i2);
    /* approaching the cut from above: Im -> -5 Pi at x = -4.5 */
    ASSERT(sfg_loggamma_c(-4.5, 1e-9, &r1, &i1));
    ASSERT_MSG(fabs(i1 + 5.0 * M_PI) < 1e-7, "Im LogGamma(-4.5 + 0i) = %.17g", i1);
    ASSERT(sfg_loggamma_c(-4.5 + 1e-7, 1e-9, &r2, &i2));
    ASSERT(fabs(i1 - i2) < 1e-6);
    /* a path from Re z = 9 to Re z = -9 at Im z = 0.7 stays continuous */
    double pr = 0, pi = 0;
    for (int k = 0; k <= 1800; k++) {
        double x = 9.0 - 0.01 * k, r, i;
        ASSERT(sfg_loggamma_c(x, 0.7, &r, &i));
        if (k) ASSERT_MSG(fabs(i - pi) < 0.05 && fabs(r - pr) < 0.1, "jump at x = %g", x);
        pr = r; pi = i;
    }
}

/* scalar == List == packed List == NDArray == Compile, bit for bit. */
void test_surfaces_bitmatch(void) {
    eval_str_free("gmXs = {0.37`, 1.4616321449683622`, 2.5`, 3.75`, 7.9`, 0.83`}");
    eval_str_free("gmZs = {0.6` + 0.3` I, 1.2` - 2.7` I, 3.9` + 2.9` I, 2.5` + 0.1` I, -1.3` + 0.4` I}");
    eval_str_free("gmChk[f_, a_, t_] := Module[{s = f /@ a}, "
                  "s === f[a] && s === f[ToPackedArray[a]] && s === Normal[f[NDArray[a]]] && "
                  "s === Compile[{{v, t, 1}}, f[v]][a]]");
    const char* heads[] = {
        "Gamma", "LogGamma", "PolyGamma", "PolyGamma[0, #]&", "PolyGamma[2, #]&",
        "Beta[#, 1.7`]&", "Beta[#, 1.3` - 0.4` I]&", "Pochhammer[#, 3]&",
        "Pochhammer[#, 2.5`]&", "Pochhammer[#, 0.5` + 0.2` I]&",
    };
    char buf[256];
    for (size_t i = 0; i < N_OF(heads); i++) {
        snprintf(buf, sizeof buf, "gmChk[%s, gmXs, _Real]", heads[i]);
        if (strstr(heads[i], " I]") == NULL) assert_true(buf);
        snprintf(buf, sizeof buf, "gmChk[%s, gmZs, _Complex]", heads[i]);
        if (strcmp(heads[i], "PolyGamma[2, #]&") != 0) assert_true(buf);  /* Re z < 0 declines for n = 2 */
    }
    /* Compile lowers every one of them, real and complex */
    assert_true("(\"Compiled\" /. CompileDiagnostics[{{x, _Real}}, PolyGamma[2, x]]) === True");
    assert_true("(\"Compiled\" /. CompileDiagnostics[{{v, _Real, 1}}, Beta[v, 2.5]]) === True");
    assert_true("(\"Compiled\" /. CompileDiagnostics[{{x, _Complex}}, Pochhammer[x, 3]]) === True");
    assert_true("(\"Compiled\" /. CompileDiagnostics[{{v, _Complex, 1}}, LogGamma[v]]) === True");
}

int main(void) {
    symtab_init();
    core_init();
    TEST(test_kernel_spot_values);
    TEST(test_kernel_seams);
    TEST(test_kernel_declines);
    TEST(test_builtin_results);
    TEST(test_loggamma_branch);
    TEST(test_surfaces_bitmatch);
    printf("All Gamma-family machine-kernel tests passed.\n");
    return 0;
}
