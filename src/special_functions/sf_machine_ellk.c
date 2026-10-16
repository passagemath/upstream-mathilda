/* Machine kernels for the complete elliptic integrals K(m) and E(m).
 *
 * ---- Real m < 1: polynomial-with-logarithm approximations ----------------
 *
 * With t = 1 - m in (0, 1], both integrals have the exact structure
 * (Abramowitz & Stegun 17.3.34 / 17.3.36, after Hastings)
 *
 *     K(m) = P_K(t) - log(t) Q_K(t)
 *     E(m) = P_E(t) - t log(t) Q_E(t)
 *
 * where the logarithm carries the whole singularity at m = 1 and the four
 * pieces are smooth on [0, 1]. This is the form of Cephes' ellpk / ellpe
 * (S. L. Moshier, Cephes Math Library, 1984-2000). The coefficients below are
 * NOT Cephes': they were re-derived by tools/sf_accuracy/fit_ellk.py as a
 * linear relative-error minimax fit (Lawson's IRLS, mpmath at 50 digits, points
 * down to t = 1e-300 so that Q(0) is pinned against |log t| ~ 690). Fitted
 * error 1.2e-18 (K) and 5.5e-19 (E) relative, i.e. far below one ulp, so the
 * double evaluation (two independent degree-10 Horner chains and one log) is
 * what sets the accuracy: about 1 ulp.
 *
 * This replaces Carlson R_F / R_F+R_D duplication for the complete integrals
 * (6-10 duplication steps, three square roots each), which the incomplete
 * integrals F, E(phi), Pi still use (elliptic.c).
 *
 * m < 0 by the imaginary-modulus transformation (DLMF 19.7.5):
 *     K(m) = K(m/(m-1)) / sqrt(1-m),   E(m) = sqrt(1-m) E(m/(m-1)),
 * and with t = 1 - m > 1 the transformed parameter has t' = 1 - m/(m-1) = 1/t,
 * so the same polynomials apply at t' = 1/t.
 *
 * ---- Complex m: the arithmetic-geometric mean -----------------------------
 *
 *     K(m) = pi / (2 M(1, sqrt(1 - m)))                       (DLMF 19.8.5)
 *     E(m) = K(m) (1 - sum_{n>=0} 2^(n-1) c_n^2),  c_0^2 = m,
 *            c_{n+1} = (a_n - b_n)/2                          (DLMF 19.8.6)
 *
 * in raw double re/im (no <complex.h>). b_0 is the principal square root, and
 * each b_{n+1} = sqrt(a_n b_n) takes the "right choice" of sign,
 * Re(a_{n+1} conj(b_{n+1})) >= 0 (Cox 1984, "The arithmetic-geometric mean of
 * Gauss", L'Enseignement Math. 30), which makes M the principal ("optimal")
 * AGM and K the principal branch, with its cut on m in [1, inf). Convergence
 * is quadratic: 4-6 square roots for |m| <= 1. Verified against mpmath
 * ellipk/ellipe over |Re m|, |Im m| <= 20 (tools/sf_accuracy/expint.py).
 */
#include "sf_machine_ellk.h"

#include <math.h>
#include <stddef.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* tools/sf_accuracy/fit_ellk.py --deg 10 --iters 31 --npts 300.
 * Highest degree first (Horner). */
static const double ELLK_KP[11] = {
    0.00012583821034900002,
    0.0021459306520971956,
    0.007770390820197493,
    0.009953864065215373,
    0.007065343010038764,
    0.00624743450972908,
    0.008796228067248696,
    0.014938198470189658,
    0.030885147587724145,
    0.09657359028221657,
    1.3862943611198906,
};
static const double ELLK_KQ[11] = {
    2.6674766742464312e-05,
    0.0008512207225535883,
    0.00567994648241944,
    0.015156381121135452,
    0.023758725527666336,
    0.030094412209904317,
    0.0373757773025507,
    0.048828000779783796,
    0.07031249952277752,
    0.12499999999975031,
    0.5,
};
static const double ELLK_EP[11] = {
    0.00045712599470521795,
    0.004094074465768602,
    0.00961742989065945,
    0.009301102859462696,
    0.006633637255346611,
    0.007355055345065001,
    0.011554018911979907,
    0.021831508307184586,
    0.05680519320454549,
    0.44314718056017904,
    1.0,
};
static const double ELLK_EQ[11] = {
    -8.158263740523923e-07,
    0.00011652341268041542,
    0.0020051860310761243,
    0.009075129893019965,
    0.019086899311931655,
    0.026889806978585146,
    0.033572637093773804,
    0.042722185395773396,
    0.05859372316040849,
    0.09374999992705718,
    0.24999999999997546,
};

/* Two degree-10 Horner chains, written out so the compiler interleaves the
 * independent recurrences (they share only t). */
static inline void ellk_pq(const double* P, const double* Q, double t,
                           double* p, double* q) {
    double a = P[0], b = Q[0];
    for (int i = 1; i < 11; i++) {
        a = a * t + P[i];
        b = b * t + Q[i];
    }
    *p = a; *q = b;
}

/* K at t = 1 - m in (0, 1]. */
static inline double ellk_k_t(double t) {
    double p, q;
    ellk_pq(ELLK_KP, ELLK_KQ, t, &p, &q);
    return p - log(t) * q;
}

/* E at t = 1 - m in [0, 1]; t = 0 is E(1) = 1. */
static inline double ellk_e_t(double t) {
    if (t == 0.0) return 1.0;
    double p, q;
    ellk_pq(ELLK_EP, ELLK_EQ, t, &p, &q);
    return p - t * log(t) * q;
}

bool sf_machine_ellk_k(double m, double* out) {
    if (!(m < 1.0)) return false;                 /* pole at 1, complex above, NaN */
    double t = 1.0 - m;
    double v;
    if (t <= 1.0) {
        v = ellk_k_t(t);
    } else {                                      /* m < 0: DLMF 19.7.5 */
        if (!isfinite(t)) return false;
        v = ellk_k_t(1.0 / t) / sqrt(t);
    }
    if (!isfinite(v)) return false;
    *out = v;
    return true;
}

bool sf_machine_ellk_e(double m, double* out) {
    if (!(m <= 1.0)) return false;                /* complex above 1, NaN */
    double t = 1.0 - m;
    double v;
    if (t <= 1.0) {
        v = ellk_e_t(t);
    } else {
        if (!isfinite(t)) return false;
        v = ellk_e_t(1.0 / t) * sqrt(t);
    }
    if (!isfinite(v)) return false;
    *out = v;
    return true;
}

/* ---- complex AGM ------------------------------------------------------- */

/* Principal square root of x + i y; a zero imaginary part of EITHER sign is
 * taken as +0, so the negative real axis maps to +i sqrt|x| (the value from
 * above, which for b_0 = sqrt(1 - m) is the m-from-below side of the cut). */
static inline void ellk_csqrt(double x, double y, double* sr, double* si) {
    if (y == 0.0) {
        if (x >= 0.0) { *sr = sqrt(x); *si = 0.0; }
        else          { *sr = 0.0;     *si = sqrt(-x); }
        return;
    }
    double r = hypot(x, y);
    if (x >= 0.0) {
        double t = sqrt(0.5 * (r + x));
        *sr = t; *si = y / (2.0 * t);
    } else {
        double t = sqrt(0.5 * (r - x));
        *sr = fabs(y) / (2.0 * t); *si = copysign(t, y);
    }
}

#define ELLK_AGM_MAXIT 40

/* M(1, sqrt(1-m)) -> (ar, ai); and, if want_sum, S = sum 2^(n-1) c_n^2. */
static bool ellk_agm(double mr, double mi, bool want_sum,
                     double* ar_out, double* ai_out, double* sr_out, double* si_out) {
    if (!isfinite(mr) || !isfinite(mi)) return false;
    /* Keep the products a_n b_n (|a b| ~ |m|) clear of overflow. */
    if (fabs(mr) > 1e300 || fabs(mi) > 1e300) return false;
    if (mr == 1.0 && mi == 0.0) return false;     /* K pole; E(1) = 1 is the caller's */
    double ar = 1.0, ai = 0.0, br, bi;
    ellk_csqrt(1.0 - mr, -mi, &br, &bi);
    double sr = 0.5 * mr, si = 0.5 * mi, p = 0.5;
    for (int n = 0; n < ELLK_AGM_MAXIT; n++) {
        double cr = 0.5 * (ar - br), ci = 0.5 * (ai - bi);
        double nar = 0.5 * (ar + br), nai = 0.5 * (ai + bi);
        double gr, gi;
        ellk_csqrt(ar * br - ai * bi, ar * bi + ai * br, &gr, &gi);
        if (nar * gr + nai * gi < 0.0) { gr = -gr; gi = -gi; }   /* right choice */
        if (want_sum) {
            p *= 2.0;
            sr += p * (cr * cr - ci * ci);
            si += p * (2.0 * cr * ci);
        }
        ar = nar; ai = nai; br = gr; bi = gi;
        /* c_{n+1} ~ c_n^2 / (4 a): once |c_n| < 1e-9 |a| the remaining terms
         * of both M and S are below 2^-60 relative. */
        if (fabs(cr) + fabs(ci) < 1e-9 * (fabs(ar) + fabs(ai))) {
            if (!isfinite(ar) || !isfinite(ai) || (ar == 0.0 && ai == 0.0)) return false;
            *ar_out = ar; *ai_out = ai;
            if (want_sum) { *sr_out = sr; *si_out = si; }
            return true;
        }
    }
    return false;                                 /* m at/near 1: no convergence */
}

/* K = (pi/2) / a by Smith's division (Smith 1962, CACM 5:435): one rounding
 * fewer per component than conj(a)/|a|^2, which measurably lowers the worst
 * case (E over the benchmark box: 5.0 -> 4.0 ulp) and cannot overflow |a|^2. */
static inline void ellk_half_pi_over(double ar, double ai, double* kr, double* ki) {
    const double h = 0.5 * M_PI;
    if (fabs(ar) >= fabs(ai)) {
        double r = ai / ar, den = ar + ai * r;
        *kr = h / den; *ki = -(h * r) / den;
    } else {
        double r = ar / ai, den = ai + ar * r;
        *kr = (h * r) / den; *ki = -h / den;
    }
}

/* A real m below the cut, held in a complex slot (a complex buffer element, or
 * Complex[m, 0.]), takes the REAL kernel: the value is real, and the real
 * kernel is both more accurate and the one a real scalar m uses -- so a complex
 * NDArray element 0.9 + 0. I and the scalar 0.9 (which is what Complex[0.9, 0.]
 * evaluates to) agree bit for bit. */
bool sf_machine_ellk_k_c(double mr, double mi, double* rr, double* ri) {
    if (mi == 0.0 && mr < 1.0) {
        if (!sf_machine_ellk_k(mr, rr)) return false;
        *ri = 0.0;
        return true;
    }
    double ar, ai, kr, ki;
    if (!ellk_agm(mr, mi, false, &ar, &ai, NULL, NULL)) return false;
    ellk_half_pi_over(ar, ai, &kr, &ki);
    if (!isfinite(kr) || !isfinite(ki)) return false;
    *rr = kr; *ri = ki;
    return true;
}

bool sf_machine_ellk_e_c(double mr, double mi, double* rr, double* ri) {
    if (mi == 0.0 && mr <= 1.0) {                /* real, as for K; E(1) = 1 */
        if (!sf_machine_ellk_e(mr, rr)) return false;
        *ri = 0.0;
        return true;
    }
    double ar, ai, sr, si;
    if (!ellk_agm(mr, mi, true, &ar, &ai, &sr, &si)) return false;
    double kr, ki;
    ellk_half_pi_over(ar, ai, &kr, &ki);
    double ur = 1.0 - sr, ui = -si;
    double er = kr * ur - ki * ui, ei = kr * ui + ki * ur;
    if (!isfinite(er) || !isfinite(ei)) return false;
    *rr = er; *ri = ei;
    return true;
}
