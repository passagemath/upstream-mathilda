/* Machine-precision kernels for the exponential-integral family.
 *
 * See sf_machine_expint.h for the contract. This file holds every double
 * kernel of the family -- real and complex -- so the scalar builtins, the
 * NDArray element-wise path and Compile[] all compute the same bits.
 *
 * ---------------------------------------------------------------------------
 * REAL ARGUMENTS
 *
 * Classical double-precision designs, with coefficients refitted here from
 * mpmath (tools/sf_accuracy/fit_expint.py -> sf_machine_expint_coef.h):
 *
 *   Si, Ci   (Cephes `sici`, S. L. Moshier 1989)
 *     0 < x <= 4   Si = x S(x^2),  Ci = gamma + log x + x^2 C(x^2)
 *                  (the Maclaurin series DLMF 6.6.5-6.6.6, economised)
 *     x > 4        Si = pi/2 - f cos x - g sin x,  Ci = f sin x - g cos x
 *                  with the auxiliary functions of DLMF 6.2.17-6.2.18 written
 *                  x f = 1 - w RF(w), x^2 g = 1 - w RG(w), w = 1/x^2, RF/RG
 *                  rational on [4, 8] and [8, inf).
 *   Shi, Chi (Cephes `shichi`)
 *     0 < x <= 4   Shi = x Sh(x^2), Chi = gamma + log x + x^2 Ch(x^2)
 *     x > 4        Shi = (Ei(x) + E1(x))/2, Chi = (Ei(x) - E1(x))/2
 *   Ei, E1   (Cody & Thacher, Math. Comp. 22 (1968) 641 and 23 (1969) 289)
 *     Ei, 0 < x <= 16    gamma + log x + x R(x)          (DLMF 6.6.2), R on
 *                        six pieces, each a positive-coefficient polynomial
 *     Ei, 16 < x <= 64   e^x/x (1 + u H(u)), u = 1/x, rational on 2 octaves
 *     Ei, x > 64         e^x/x sum k!/x^k (DLMF 6.12.2), truncated at 1e-17
 *     E1, 0 < x <= 0.6   -gamma - log x + x R(-x)
 *     E1, x > 0.6        e^-x/x (1 - u J(u)), rational on [0.6,2], [2,8], [8,inf)
 *   li(x) = Ei(log x), with e^{log x} replaced by x itself in the large-argument
 *   forms, so the rounding of log x is not amplified by an exp.
 *   E_n(x)  (Numerical Recipes 6.3 / Cephes `expn`): the power series for
 *   x <= 1 and the continued fraction (DLMF 8.19.17) for x > 1, the latter
 *   evaluated by the forward Wallis recurrences with rescaling -- one division
 *   at the end instead of two per level, which is what makes it fast.
 *
 * ---------------------------------------------------------------------------
 * COMPLEX ARGUMENTS
 *
 * Two methods, chosen per argument so neither is used where it loses digits:
 *
 *   power series   Ei: gamma + log z + sum z^k/(k k!);  Si: sum (-1)^k
 *                  z^(2k+1)/((2k+1)(2k+1)!);  Ci: gamma + log z + sum ...
 *                  evaluated by HORNER with tabulated coefficients (the term
 *                  recurrence t_k = t_{k-1} z/k accumulates O(k) roundings per
 *                  term and was measured at 5 ulp where Horner gives 2).
 *                  Used where the terms do not cancel: the loss is about
 *                  exp(|z| - Re z) for Ei and exp(|z| - |Im z|) for Si/Ci, so
 *                  the series serves |z| - Re z <= 2 (resp. |z| - |Im z| <= 2).
 *   continued fraction  e^w E1(w) = 1/(w + 1 - 1/(w + 3 - 4/(w + 5 - ...)))
 *                  (DLMF 6.9.1, even contraction), evaluated BACKWARD from a
 *                  depth fixed in advance. Forward (Lentz) evaluation multiplies
 *                  one rounding per level into the result and measured 55 ulp at
 *                  depth ~60; the backward recurrence damps them. The depth
 *                  n = 240/s + 8, s = |w| + Re w, is calibrated against mpmath
 *                  (the truncation error decays like exp(-4 Re sqrt(n w))).
 *                  Si/Ci then come from E1(+-iz) (DLMF 6.5.5-6.5.6), Ei from
 *                  E1(-z) (Ei(z) = -E1(-z) + i pi sgn Im z).
 *   |z| > 45       the asymptotic expansions (DLMF 6.12.2, 6.12.3), whose
 *                  smallest term is below 1e-19 there.
 * Shi(z) = -i Si(iz) and Chi(z) = Ci(iz) + (log z - log iz) are exact
 * rotations of the Si/Ci kernels; li(z) = Ei(log z).
 *
 * All complex arithmetic is spelled out on (re, im) doubles: no <complex.h>
 * (libgcc's __divdc3 is slow, and the C99 complex surface is a portability
 * hazard on MinGW).
 */
#include "sf_machine.h"
#include "sf_machine_expint.h"

#include <float.h>
#include <math.h>

#include "sf_machine_expint_coef.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define EX_EULER   0.57721566490153286061
#define EX_PI_2_HI 1.5707963267948966          /* double nearest pi/2 */
#define EX_PI_2_LO 6.123233995736766e-17       /* pi/2 - EX_PI_2_HI    */
#define EX_NELEM(a) ((int)(sizeof(a) / sizeof((a)[0])))

/* ------------------------------------------------------------------------- */
/*  Small helpers                                                            */
/* ------------------------------------------------------------------------- */

/* Horner: c[0] + c[1] x + ... + c[n-1] x^(n-1). */
static double ex_poly(const double* c, int n, double x) {
    double r = c[n - 1];
    for (int i = n - 2; i >= 0; i--) r = r * x + c[i];
    return r;
}

#define EX_RAT(P, Q, x) (ex_poly(P, EX_NELEM(P), x) / ex_poly(Q, EX_NELEM(Q), x))

/* A double result is usable: finite, and not an underflow to a subnormal/zero
 * that only the arbitrary-precision path can represent. */
static bool ex_ok(double v) {
    return isfinite(v) && (v == 0.0 || fabs(v) >= DBL_MIN);
}

/* ------------------------------------------------------------------------- */
/*  Si / Ci, real                                                             */
/* ------------------------------------------------------------------------- */

/* x f(x) and x^2 g(x) for x > 4 (the auxiliary functions, DLMF 6.2.17). */
static void sici_fg(double x, double* F, double* G) {
    double w = 1.0 / (x * x);
    if (x <= 8.0) {
        *F = 1.0 - w * EX_RAT(SI_RFA_P, SI_RFA_Q, w);
        *G = 1.0 - w * EX_RAT(SI_RGA_P, SI_RGA_Q, w);
    } else {
        *F = 1.0 - w * EX_RAT(SI_RFB_P, SI_RFB_Q, w);
        *G = 1.0 - w * EX_RAT(SI_RGB_P, SI_RGB_Q, w);
    }
}

bool sf_machine_si(double x, double* out) {
    if (isnan(x)) return false;
    double ax = fabs(x), v;
    if (ax <= 4.0) {
        double t = ax * ax * 0.125 - 1.0;
        v = ax * ex_poly(SI_S, EX_NELEM(SI_S), t);
    } else if (isinf(ax)) {
        v = EX_PI_2_HI;
    } else {
        double F, G;
        sici_fg(ax, &F, &G);
        double s = sin(ax), c = cos(ax);
        double corr = (F * c + G * s / ax) / ax;
        v = (EX_PI_2_HI - corr) + EX_PI_2_LO;
    }
    *out = (x < 0.0) ? -v : v;
    return true;
}

bool sf_machine_ci(double x, double* out) {
    if (!(x > 0.0) || isinf(x)) return false;     /* complex below, -Inf at 0 */
    double v;
    if (x <= 4.0) {
        double y = x * x, t = y * 0.125 - 1.0;
        v = (EX_EULER + log(x)) + y * ex_poly(SI_C, EX_NELEM(SI_C), t);
    } else {
        double F, G;
        sici_fg(x, &F, &G);
        double s = sin(x), c = cos(x);
        v = (F * s - G * c / x) / x;
    }
    *out = v;
    return ex_ok(v);
}

/* ------------------------------------------------------------------------- */
/*  Ei / E1, real                                                             */
/* ------------------------------------------------------------------------- */

/* H(x) = x e^-x Ei(x) for x > 16. */
typedef struct { double ua, ub; const double* P; int np; const double* Q; int nq; } ExRatT;
static const ExRatT EI_H_TAB[] = {
    { 0.03125,  0.0625,  EI_HD_P, EX_NELEM(EI_HD_P), EI_HD_Q, EX_NELEM(EI_HD_Q) },  /* [16, 32] */
    { 0.015625, 0.03125, EI_HE_P, EX_NELEM(EI_HE_P), EI_HE_Q, EX_NELEM(EI_HE_Q) },  /* [32, 64] */
};

static double ei_H(double x) {
    if (x > 64.0) {
        /* Asymptotic sum k!/x^k (DLMF 6.12.2): at x > 64 the terms fall below
         * 1e-17 long before they turn, so truncation is exact to double. */
        double u = 1.0 / x, term = 1.0, sum = 1.0;
        for (int k = 1; k < 200; k++) {
            term *= k * u;
            sum += term;
            if (term < 1e-17 * sum) break;
        }
        return sum;
    }
    const ExRatT* r = &EI_H_TAB[x <= 32.0 ? 0 : 1];
    double u = 1.0 / x;
    double t = (2.0 * u - r->ua - r->ub) / (r->ub - r->ua);
    double R = ex_poly(r->P, r->np, t) / ex_poly(r->Q, r->nq, t);
    return 1.0 + u * R;
}

/* Ei(x) = gamma + log x + x R(x) for 0 < x <= 16, R = sum x^(k-1)/(k k!)
 * absolutely monotonic, fitted piecewise in v = (x - a)/(b - a) so every
 * coefficient is positive and Horner never cancels. */
static double ei_R(double x) {
    if (x <= 2.0)  return ex_poly(EI_R0, EX_NELEM(EI_R0), x * 0.5);
    if (x <= 4.0)  return ex_poly(EI_R1, EX_NELEM(EI_R1), (x - 2.0) * 0.5);
    if (x <= 6.0)  return ex_poly(EI_R2, EX_NELEM(EI_R2), (x - 4.0) * 0.5);
    if (x <= 8.0)  return ex_poly(EI_R3, EX_NELEM(EI_R3), (x - 6.0) * 0.5);
    if (x <= 12.0) return ex_poly(EI_R4, EX_NELEM(EI_R4), (x - 8.0) * 0.25);
    return ex_poly(EI_R5, EX_NELEM(EI_R5), (x - 12.0) * 0.25);
}

/* Ei(x) for 0 < x <= 16, no exponential needed. */
static double ei_mid(double x) {
    return (EX_EULER + log(x)) + x * ei_R(x);
}

/* J(x) = x e^x E1(x) for x > 0.6. */
static double e1_J(double x) {
    double u = 1.0 / x;
    if (x <= 2.0)  return 1.0 - u * EX_RAT(E1_JA_P, E1_JA_Q, u);
    if (x <= 8.0)  return 1.0 - u * EX_RAT(E1_JB_P, E1_JB_Q, u);
    return 1.0 - u * EX_RAT(E1_JC_P, E1_JC_Q, u);
}

/* E1(x) for 0 < x <= 0.6: -gamma - log x + x R(-x). */
static double e1_small(double x) {
    double t = x * (2.0 / 0.6) - 1.0;
    return x * ex_poly(E1_R, EX_NELEM(E1_R), t) - (EX_EULER + log(x));
}

bool sf_machine_e1(double x, double* out) {
    if (!(x > 0.0) || isinf(x)) return false;     /* pole at 0, complex below */
    double v;
    if (x <= 0.6) v = e1_small(x);
    else {
        if (x > 745.0) return false;              /* underflows: MPFR's job */
        v = exp(-x) / x * e1_J(x);
    }
    *out = v;
    return ex_ok(v);
}

bool sf_machine_ei(double x, double* out) {
    if (isnan(x) || isinf(x) || x == 0.0) return false;   /* -Infinity at 0 */
    if (x < 0.0) {
        double e1;
        if (!sf_machine_e1(-x, &e1)) return false;
        *out = -e1;
        return true;
    }
    double v;
    if (x <= 16.0) v = ei_mid(x);
    else {
        double h = 0.5 * x;
        double eh = exp(h);                       /* e^x as e^(x/2)^2: no early overflow */
        v = (eh * ei_H(x) / x) * eh;
    }
    *out = v;
    return ex_ok(v);
}

/* ------------------------------------------------------------------------- */
/*  Shi / Chi, real                                                           */
/* ------------------------------------------------------------------------- */

/* Ei(x) and E1(x) for x > 4, for Shi = (Ei + E1)/2 and Chi = (Ei - E1)/2. */
static bool shichi_ei_e1(double x, double* ep, double* em) {
    if (x > 720.0) return false;                  /* Shi, Chi overflow */
    double eh = exp(0.5 * x);
    *ep = (x <= 16.0) ? ei_mid(x) : (eh * ei_H(x) / x) * eh;
    *em = 1.0 / (eh * eh) / x * e1_J(x);
    return true;
}

bool sf_machine_shi(double x, double* out) {
    if (isnan(x)) return false;
    double ax = fabs(x), v;
    if (ax <= 4.0) {
        v = ax * ex_poly(SH_S, EX_NELEM(SH_S), ax * ax * 0.0625);
    } else {
        double ep, em;
        if (!shichi_ei_e1(ax, &ep, &em)) return false;
        v = 0.5 * (ep + em);
    }
    *out = (x < 0.0) ? -v : v;
    return ex_ok(v);
}

bool sf_machine_chi(double x, double* out) {
    if (!(x > 0.0) || isinf(x)) return false;     /* complex below, -Inf at 0 */
    double v;
    if (x <= 4.0) {
        double y = x * x;
        v = (EX_EULER + log(x)) + y * ex_poly(SH_C, EX_NELEM(SH_C), y * 0.0625);
    } else {
        double ep, em;
        if (!shichi_ei_e1(x, &ep, &em)) return false;
        v = 0.5 * (ep - em);
    }
    *out = v;
    return ex_ok(v);
}

/* ------------------------------------------------------------------------- */
/*  li, real                                                                  */
/* ------------------------------------------------------------------------- */

/* li(x) = Ei(log x). In the large-argument forms Ei(y) = e^y/y H(y) and
 * E1(t) = e^-t/t J(t), the factor e^{+-y} IS x: using x itself keeps the
 * rounding of y = log x from being amplified by an exponential. */
bool sf_machine_li(double x, double* out) {
    /* li(0) = 0 by the limit, but the interpreter answers Indeterminate there,
     * so the machine path declines rather than be the only one to give a
     * number; li(1) is the pole. */
    if (!(x > 0.0) || x == 1.0 || isinf(x)) return false;
    double y = log(x), v;
    if (y > 0.0) {
        if (y <= 16.0) v = ei_mid(y);
        else           v = x / y * ei_H(y);
    } else {
        double t = -y;
        if (t <= 0.6) v = -e1_small(t);
        else          v = -(x / t * e1_J(t));
    }
    *out = v;
    return ex_ok(v);
}

/* ------------------------------------------------------------------------- */
/*  E_n(x), real                                                              */
/* ------------------------------------------------------------------------- */

bool sf_machine_expn(double nd, double x, double* out) {
    if (!(x > 0.0) || isinf(x)) return false;
    if (nd != floor(nd) || nd < 0.0 || nd > 1.0e6) return false;
    long n = (long)nd;
    double v;
    if (n == 0) {
        if (x > 745.0) return false;
        v = exp(-x) / x;
    } else if (n == 1) {
        return sf_machine_e1(x, out);
    } else if (x <= 1.0) {
        /* Power series (DLMF 8.19.8 / A&S 5.1.12):
         *   E_n(x) = (-x)^(n-1)/(n-1)! (psi(n) - log x) - sum_{k != n-1} (-x)^k/(k! (k-n+1)). */
        long nm1 = n - 1;
        double ans = 1.0 / (double)nm1;           /* the k = 0 term */
        double fact = 1.0;                        /* (-x)^k / k! */
        for (long k = 1; k < 400; k++) {
            fact *= -x / (double)k;
            double del;
            if (k != nm1) {
                del = -fact / (double)(k - nm1);
            } else {
                double psi = -EX_EULER;
                for (long i = 1; i <= nm1; i++) psi += 1.0 / (double)i;
                del = fact * (psi - log(x));
            }
            ans += del;
            if (fabs(del) < fabs(ans) * 1e-17 && k >= nm1) break;
        }
        v = ans;
    } else {
        /* Continued fraction (DLMF 8.19.17, even form):
         *   E_n(x) = e^-x / (x + n - 1 n/(x + n + 2 - 2(n+1)/(x + n + 4 - ...)))
         * The convergents A_k/B_k obey A_k = b_k A_{k-1} + a_k A_{k-2} (Wallis),
         * which is multiply-add only; the pair is rescaled every 8 levels so it
         * cannot overflow, and the convergence test (one division) is made at
         * the same cadence. */
        double Am2 = 1.0, Bm2 = 0.0;              /* A_{-1}, B_{-1} */
        double b = x + (double)n;
        double Am1 = 0.0, Bm1 = 1.0;              /* A_0 = 0 ... the CF is 1/(b0 + ...) */
        /* level 1: A_1 = 1, B_1 = b0 */
        double A = 1.0, B = b;
        Am2 = Am1; Bm2 = Bm1; Am1 = A; Bm1 = B;
        double prev = A / B;
        bool done = false;
        for (long i = 1; i < 5000; i++) {
            double a = -(double)i * (double)(n - 1 + i);
            b += 2.0;
            A = b * Am1 + a * Am2;
            B = b * Bm1 + a * Bm2;
            Am2 = Am1; Bm2 = Bm1; Am1 = A; Bm1 = B;
            if ((i & 7) == 0) {
                double s = 1.0 / B;
                Am2 *= s; Bm2 *= s; Am1 *= s; Bm1 = 1.0;
                double cur = Am1;
                if (fabs(cur - prev) <= 1e-17 * fabs(cur)) { done = true; break; }
                prev = cur;
            }
        }
        if (!done) return false;
        if (x > 745.0) return false;
        v = Am1 / Bm1 * exp(-x);
    }
    *out = v;
    return ex_ok(v);
}

/* ========================================================================= */
/*  COMPLEX                                                                  */
/* ========================================================================= */

typedef struct { double re, im; } ExC;

static inline ExC exc(double re, double im) { ExC z; z.re = re; z.im = im; return z; }
static inline ExC ex_add(ExC a, ExC b) { return exc(a.re + b.re, a.im + b.im); }
static inline ExC ex_sub(ExC a, ExC b) { return exc(a.re - b.re, a.im - b.im); }
static inline ExC ex_mul(ExC a, ExC b) {
    return exc(a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re);
}
static inline ExC ex_scale(ExC a, double s) { return exc(a.re * s, a.im * s); }
/* a / b by Smith's algorithm (no spurious overflow in |b|^2). */
static ExC ex_div(ExC a, ExC b) {
    if (fabs(b.re) >= fabs(b.im)) {
        double r = b.im / b.re, d = b.re + b.im * r;
        return exc((a.re + a.im * r) / d, (a.im - a.re * r) / d);
    } else {
        double r = b.re / b.im, d = b.re * r + b.im;
        return exc((a.re * r + a.im) / d, (a.im * r - a.re) / d);
    }
}
static inline ExC ex_inv(ExC b) { return ex_div(exc(1.0, 0.0), b); }
static ExC ex_exp(ExC z) {
    double m = exp(z.re);
    return exc(m * cos(z.im), m * sin(z.im));
}
static ExC ex_log(ExC z) {
    return exc(log(hypot(z.re, z.im)), atan2(z.im, z.re));
}
static inline double ex_abs(ExC z) { return hypot(z.re, z.im); }
static inline bool ex_finite(ExC z) { return isfinite(z.re) && isfinite(z.im); }

/* Complex Horner: c[0] + c[1] w + ... + c[n-1] w^(n-1). */
static ExC ex_cpoly(const double* c, int n, ExC w) {
    ExC r = exc(c[n - 1], 0.0);
    for (int i = n - 2; i >= 0; i--) {
        r = ex_mul(r, w);
        r.re += c[i];
    }
    return r;
}

/* How many terms of sum c_k r^k (c from `tab`, tab[k-1] for power k, a series
 * with decreasing c) are needed before a term drops below 1e-18 * max(1, r):
 * the series sums only run where cancellation is bounded by e^2, so this
 * relative cut is safely below an ulp of the result. */
static int ex_nterms(const double* tab, int nmax, double r, int step) {
    double rp = r, rstep = 1.0;
    for (int i = 0; i < step; i++) rstep *= r;
    double lim = 1e-18 * (r > 1.0 ? r : 1.0);
    for (int k = 0; k < nmax; k++) {
        if (tab[k] * rp < lim && k >= 2) return k + 1;
        rp *= rstep;
    }
    return nmax;
}

/* g(w) = e^w E1(w) by the backward continued fraction (see the header
 * comment). Valid for |w| + Re w > ~2, i.e. away from the negative real axis. */
static ExC ex_ge1_cf(ExC w) {
    double s = ex_abs(w) + w.re;
    int n = (int)(240.0 / s) + 8;
    if (n > 400) n = 400;
    ExC t = exc(w.re + (2.0 * n + 1.0), w.im);
    for (int k = n; k >= 1; k--) {
        /* t = w + (2k-1) - k^2/t */
        double kk = (double)k * (double)k;
        ExC q = ex_div(exc(kk, 0.0), t);
        t = exc(w.re + (2.0 * k - 1.0) - q.re, w.im - q.im);
    }
    return ex_inv(t);
}

/* Ei(z) power series: gamma + log z + sum z^k/(k k!). */
static ExC ex_ei_series(ExC z) {
    double r = ex_abs(z);
    int n = ex_nterms(EIN_C, EIN_C_N, r, 1);
    ExC s = ex_mul(ex_cpoly(EIN_C, n, z), z);
    ExC l = ex_log(z);
    return exc(EX_EULER + l.re + s.re, l.im + s.im);
}

/* Ei(z) asymptotic for |z| > 45 near the positive real axis:
 * e^z/z sum k!/z^k (+ i pi sgn Im z, exponentially small relative). */
static ExC ex_ei_asym(ExC z) {
    ExC u = ex_inv(z);
    ExC term = exc(1.0, 0.0), sum = exc(1.0, 0.0);
    double r = ex_abs(z);
    for (int k = 1; k < 200; k++) {
        term = ex_scale(ex_mul(term, u), (double)k);
        sum = ex_add(sum, term);
        if (ex_abs(term) < 1e-18 * ex_abs(sum) || (double)k > r) break;
    }
    ExC v = ex_mul(ex_mul(ex_exp(z), u), sum);
    if (z.im > 0.0) v.im += M_PI;
    else if (z.im < 0.0) v.im -= M_PI;
    return v;
}

bool sf_machine_ei_c(double re, double im, double* ore, double* oim) {
    if (!isfinite(re) || !isfinite(im)) return false;
    if (re == 0.0 && im == 0.0) return false;                 /* -Infinity */
    ExC z = exc(re, im), v;
    double r = ex_abs(z);
    if (r - re <= 2.0) {
        v = (r > 45.0) ? ex_ei_asym(z) : ex_ei_series(z);
    } else {
        /* Ei(z) = -E1(-z) + i pi sgn(Im z); on the cut (Im z = 0, Re z < 0)
         * the value is real, -E1(|x|). */
        ExC w = exc(-re, -im);
        ExC e1 = ex_mul(ex_ge1_cf(w), ex_exp(z));   /* E1(w) = e^-w g(w) = e^z g */
        v = exc(-e1.re, -e1.im);
        if (im > 0.0) v.im += M_PI;
        else if (im < 0.0) v.im -= M_PI;
    }
    if (!ex_finite(v)) return false;
    *ore = v.re; *oim = v.im;
    return true;
}

bool sf_machine_e1_c(double re, double im, double* ore, double* oim) {
    if (!isfinite(re) || !isfinite(im)) return false;
    if (re == 0.0 && im == 0.0) return false;
    /* E1(w) = -Ei(-w) + i pi sgn(Im w)... on the principal branch, i.e.
     * Ei(z) = -E1(-z) + i pi sgn(Im z) read backwards; for Im w = 0, w < 0 the
     * upper side (Im w = +0): E1(-x + i0) = -Ei(x) - i pi. */
    ExC w = exc(re, im), v;
    double r = ex_abs(w);
    if (r + re > 2.0) {
        v = ex_mul(ex_ge1_cf(w), ex_exp(exc(-re, -im)));
    } else {
        double er, ei;
        if (!sf_machine_ei_c(-re, -im, &er, &ei)) return false;
        v = exc(-er, -ei);
        /* -Ei(-w) + i pi sgn(Im w), where sgn(Im(-w)) was used inside Ei:
         * net: E1(w) = -Ei(-w) - i pi sgn(Im(-w)) ... with Im w = 0 meaning the
         * upper side of the cut for w < 0. */
        if (im > 0.0)      v.im -= M_PI;
        else if (im < 0.0) v.im += M_PI;
        else if (re < 0.0) v.im -= M_PI;
    }
    if (!ex_finite(v)) return false;
    *ore = v.re; *oim = v.im;
    return true;
}

/* ---- Si / Ci, complex ---------------------------------------------------- */

/* Si(z) = z sum SIN_C[k] (s z^2)^k with s = -1 (Si) or +1 (Shi). */
static ExC ex_sin_series(ExC z, double sgn) {
    double r = ex_abs(z);
    int n = ex_nterms(SIN_C, SIN_C_N, r, 2);
    ExC w = ex_mul(z, z);
    w = exc(sgn * w.re, sgn * w.im);
    return ex_mul(ex_cpoly(SIN_C, n, w), z);
}

/* sum_{k>=1} COS_C[k-1] (s z^2)^k: Ci(z) = gamma + log z + this with s = -1. */
static ExC ex_cos_series(ExC z, double sgn) {
    double r = ex_abs(z);
    int n = ex_nterms(COS_C, COS_C_N, r, 2);
    ExC w = ex_mul(z, z);
    w = exc(sgn * w.re, sgn * w.im);
    return ex_mul(ex_cpoly(COS_C, n, w), w);
}

/* f(z), g(z) asymptotic for |z| > 45, Re z >= 0 (DLMF 6.12.3-6.12.4):
 *   f ~ (1/z) sum (-1)^k (2k)!/z^(2k),  g ~ (1/z^2) sum (-1)^k (2k+1)!/z^(2k). */
static void ex_fg_asym(ExC z, ExC* f, ExC* g) {
    ExC u = ex_inv(z);
    ExC u2 = ex_mul(u, u);
    ExC tf = exc(1.0, 0.0), sf = tf, tg = tf, sg = tf;
    double r = ex_abs(z);
    for (int k = 1; k < 200; k++) {
        tf = ex_scale(ex_mul(tf, u2), -(double)(2 * k - 1) * (double)(2 * k));
        tg = ex_scale(ex_mul(tg, u2), -(double)(2 * k) * (double)(2 * k + 1));
        sf = ex_add(sf, tf);
        sg = ex_add(sg, tg);
        if ((ex_abs(tf) < 1e-18 * ex_abs(sf) && ex_abs(tg) < 1e-18 * ex_abs(sg))
            || (double)(2 * k + 1) > r) break;
    }
    *f = ex_mul(sf, u);
    *g = ex_mul(sg, u2);
}

/* Si and Ci together for Re z > 0 (or = 0 with z != 0) outside the series
 * region. Ci is returned WITHOUT any branch adjustment (principal, Re z > 0).
 * want: 1 = Si only, 2 = Ci only, 3 = both. */
static void ex_sici_far(ExC z, ExC* si, ExC* ci) {
    double r = ex_abs(z);
    ExC f, g;
    ExC sz = exc(sin(z.re) * cosh(z.im), cos(z.re) * sinh(z.im));
    ExC cz = exc(cos(z.re) * cosh(z.im), -sin(z.re) * sinh(z.im));
    if (r > 45.0) {
        ex_fg_asym(z, &f, &g);
        /* Si = pi/2 - f cos z - g sin z,  Ci = f sin z - g cos z */
        ExC a = ex_add(ex_mul(f, cz), ex_mul(g, sz));
        *si = exc((EX_PI_2_HI - a.re) + EX_PI_2_LO, -a.im);
        *ci = ex_sub(ex_mul(f, sz), ex_mul(g, cz));
        return;
    }
    /* E1(iz) = e^{-iz} g(iz),  E1(-iz) = e^{iz} g(-iz) */
    ExC iz = exc(-z.im, z.re), miz = exc(z.im, -z.re);
    ExC g1 = ex_ge1_cf(iz), g2 = ex_ge1_cf(miz);
    ExC e1p = ex_mul(g1, ex_exp(miz));     /* E1(iz)  */
    ExC e1m = ex_mul(g2, ex_exp(iz));      /* E1(-iz) */
    /* Si = pi/2 + (E1(iz) - E1(-iz))/(2i),  Ci = -(E1(iz) + E1(-iz))/2 */
    ExC d = ex_sub(e1p, e1m);
    *si = exc((EX_PI_2_HI + 0.5 * d.im) + EX_PI_2_LO, -0.5 * d.re);
    ExC s = ex_add(e1p, e1m);
    *ci = exc(-0.5 * s.re, -0.5 * s.im);
}

/* Which method for Si/Ci at z: the series where |z| - |Im z| <= 2 (terms do
 * not cancel) or |z| <= 2.5, and |z| <= 45; otherwise the far forms. */
static bool ex_sici_use_series(ExC z) {
    double r = ex_abs(z);
    if (r > 45.0) return false;
    return r <= 2.5 || r - fabs(z.im) <= 2.0;
}

bool sf_machine_si_c(double re, double im, double* ore, double* oim) {
    if (!isfinite(re) || !isfinite(im)) return false;
    ExC z = exc(re, im), v;
    if (ex_sici_use_series(z)) {
        v = ex_sin_series(z, -1.0);
    } else {
        /* Si is odd: fold onto Re z >= 0. */
        bool neg = re < 0.0;
        ExC zz = neg ? exc(-re, -im) : z, ci;
        ex_sici_far(zz, &v, &ci);
        if (neg) v = exc(-v.re, -v.im);
    }
    if (!ex_finite(v)) return false;
    *ore = v.re; *oim = v.im;
    return true;
}

/* Ci on the principal branch, cut (-inf, 0], upper side on the cut. */
static bool ex_ci(ExC z, ExC* out) {
    if (z.re == 0.0 && z.im == 0.0) return false;
    if (ex_sici_use_series(z)) {
        ExC s = ex_cos_series(z, -1.0), l = ex_log(z);
        *out = exc((EX_EULER + l.re) + s.re, l.im + s.im);
    } else if (z.re >= 0.0) {
        ExC si;
        ex_sici_far(z, &si, out);
    } else {
        /* Ci(z) = Ci(-z) + log z - log(-z) = Ci(-z) +- i pi (DLMF 6.4.4). */
        ExC si, c;
        ex_sici_far(exc(-z.re, -z.im), &si, &c);
        c.im += (z.im >= 0.0) ? M_PI : -M_PI;
        *out = c;
    }
    return ex_finite(*out);
}

bool sf_machine_ci_c(double re, double im, double* ore, double* oim) {
    if (!isfinite(re) || !isfinite(im)) return false;
    ExC v;
    if (!ex_ci(exc(re, im), &v)) return false;
    *ore = v.re; *oim = v.im;
    return true;
}

/* Shi(z) = -i Si(iz)  (DLMF 6.2.15 with 6.2.5). */
bool sf_machine_shi_c(double re, double im, double* ore, double* oim) {
    if (!isfinite(re) || !isfinite(im)) return false;
    ExC z = exc(re, im), v;
    ExC iz = exc(-im, re);
    if (ex_sici_use_series(iz)) {
        v = ex_sin_series(z, 1.0);          /* the same sum with +z^2, no rotation */
    } else {
        double sr, si;
        if (!sf_machine_si_c(iz.re, iz.im, &sr, &si)) return false;
        v = exc(si, -sr);                   /* -i (sr + i si) */
    }
    if (!ex_finite(v)) return false;
    *ore = v.re; *oim = v.im;
    return true;
}

/* Chi(z) = gamma + log z + sum z^(2k)/(2k (2k)!) = Ci(iz) + (log z - log iz).
 * log z - log(iz) is the constant -i pi/2 when arg z <= pi/2 and +3 i pi/2
 * when arg z > pi/2 (where arg(iz) wraps); taken exactly, not by subtraction. */
bool sf_machine_chi_c(double re, double im, double* ore, double* oim) {
    if (!isfinite(re) || !isfinite(im)) return false;
    if (re == 0.0 && im == 0.0) return false;
    ExC z = exc(re, im), v;
    ExC iz = exc(-im, re);
    if (ex_sici_use_series(iz)) {
        ExC s = ex_cos_series(z, 1.0), l = ex_log(z);
        v = exc((EX_EULER + l.re) + s.re, l.im + s.im);
    } else {
        if (!ex_ci(iz, &v)) return false;
        /* arg z > pi/2  <=>  Re z < 0 and Im z >= 0 (upper side on the cut) */
        bool wrap = (re < 0.0 && im >= 0.0);
        v.im += wrap ? 1.5 * M_PI : -0.5 * M_PI;
    }
    if (!ex_finite(v)) return false;
    *ore = v.re; *oim = v.im;
    return true;
}

/* li(z) = Ei(log z), principal log (the identity the interpreter uses). */
bool sf_machine_li_c(double re, double im, double* ore, double* oim) {
    if (!isfinite(re) || !isfinite(im)) return false;
    if (re == 0.0 && im == 0.0) return false;              /* Indeterminate */
    if (re == 1.0 && im == 0.0) return false;              /* pole */
    ExC l = ex_log(exc(re, im));
    return sf_machine_ei_c(l.re, l.im, ore, oim);
}
