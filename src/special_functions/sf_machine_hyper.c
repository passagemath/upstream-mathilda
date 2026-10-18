/* Machine-precision kernels for the hypergeometric family.
 *
 * See sf_machine_hyper.h for the contract. Everything here is plain C99 on
 * `double`, with complex values carried as explicit (re, im) pairs -- no
 * <complex.h> (libgcc's __divdc3/__muldc3 are slow and their NaN/Inf recovery is
 * not what a kernel that declines on non-finite wants).
 *
 * The one idea that runs through the file: a SERIES MEASURES ITS OWN LOSS.
 * Every power series is summed once in double while tracking its largest term;
 * peak/|sum| is exactly how much of the double's precision the summation threw
 * away. Up to HY_LOSS_OK the double sum is kept; beyond it the same series is
 * re-summed in double-double (Dekker 1971; Hida, Li & Bailey 2001), which has
 * ~106 bits to lose before the answer is touched; beyond HY_LOSS_MAX the kernel
 * declines. The transformations below exist to keep the series argument small
 * and the terms one-signed, so the double-double pass is the exception.
 *
 * References
 *   [AS]  Abramowitz & Stegun, Handbook of Mathematical Functions (1964), ch. 13
 *         and 15 -- 15.3.3-15.3.12 (2F1 transformations and the degenerate
 *         logarithmic cases), 13.1.27 (Kummer), 13.1.3 (U connection).
 *   [DLMF] NIST DLMF 13.2.9 (U, integer b), 13.7.2 & 13.7.3 (1F1 / U large-|z|
 *         asymptotics), 10.40.1 (I_nu Hankel asymptotic), 10.23.15-type Neumann
 *         normalisation for the Miller algorithm.
 *   [Forrey] R. C. Forrey, "Computing the hypergeometric function", J. Comput.
 *         Phys. 137 (1997) 79-100 -- choice of transformation by |w|, and the
 *         degenerate-case treatment.
 *   [MS]  N. Michel & M. V. Stoitsov, "Fast computation of the Gauss
 *         hypergeometric function with all its parameters complex", Comput.
 *         Phys. Commun. 178 (2008) 535-551 -- transformation selection.
 *   [Miller] J. C. P. Miller (1952), backward recurrence; Gautschi, SIAM Rev. 9
 *         (1967) 24-82, "Computational aspects of three-term recurrence
 *         relations", for the normalisation by a Neumann series.
 */
#include "sf_machine_hyper.h"
#include "sf_machine.h"

#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ------------------------------------------------------------------ */
/*  Tuning constants                                                   */
/* ------------------------------------------------------------------ */

#define HY_MAXTERMS   20000        /* double series term cap                   */
#define HY_DD_MAXTERMS 6000        /* double-double series term cap            */
#define HY_TOL        5.5511151231257827e-17   /* 2^-54: stop when a term is below */
#define HY_DD_TOL     8.6736173798840355e-19   /* 2^-60 for the dd re-sum          */
#define HY_LOSS_OK    2.0          /* peak/|sum| a double sum may lose (1 bit)  */
#define HY_LOSS_MAX   3.5184372088832e13       /* 2^45: beyond this even dd declines */
#define HY_KPEAK_DD   6            /* peak term deeper than this: re-sum in dd   */
#define HY_KAPPA_MAX  2.0          /* condition of a connection-formula combine */
#define HY_KAPPA_LAST 8.0          /* ... accepted only when nothing else converges */

/* ------------------------------------------------------------------ */
/*  Small real helpers                                                 */
/* ------------------------------------------------------------------ */

static bool hy_is_nonpos_int(double x) { return x <= 0.0 && x == floor(x); }

/* 1/Gamma(x), finite everywhere it is representable: 0 at the poles of Gamma.
 * Non-finite (and the caller declines) only when 1/Gamma itself overflows. */
static double hy_rgamma(double x) {
    if (hy_is_nonpos_int(x)) return 0.0;
    if (x > 171.0) return 0.0 * x + exp(-lgamma(x));      /* underflows gracefully */
    double g = tgamma(x);
    return 1.0 / g;
}

/* cos(pi x) with the argument reduced exactly (x mod 2 is exact in binary), so
 * an integer x gives exactly +-1 and a large x keeps its accuracy. */
static double hy_cospi(double x) {
    double r = fmod(fabs(x), 2.0);              /* exact */
    if (r == 0.0) return 1.0;
    if (r == 1.0) return -1.0;
    if (r == 0.5 || r == 1.5) return 0.0;
    if (r > 1.0) r = 2.0 - r;                   /* cos(pi r) symmetric about 1 */
    if (r <= 0.25) return cos(M_PI * r);
    if (r <= 0.75) return sin(M_PI * (0.5 - r));
    return -cos(M_PI * (1.0 - r));
}

/* ------------------------------------------------------------------ */
/*  Double-double arithmetic (Dekker 1971; Hida-Li-Bailey QD library)  */
/* ------------------------------------------------------------------ */

typedef struct { double hi, lo; } hdd;

static inline hdd dd_two_sum(double a, double b) {
    hdd r; r.hi = a + b;
    double bb = r.hi - a;
    r.lo = (a - (r.hi - bb)) + (b - bb);
    return r;
}
static inline hdd dd_quick(double a, double b) {          /* |a| >= |b| */
    hdd r; r.hi = a + b; r.lo = b - (r.hi - a); return r;
}
/* Dekker's exact product by Veltkamp splitting -- no FMA needed, so the result
 * is the same on every C99 target. Valid while |a|,|b| < 2^995. */
static inline void dd_split(double a, double* hi, double* lo) {
    double t = 134217729.0 * a;                 /* 2^27 + 1 */
    *hi = t - (t - a);
    *lo = a - *hi;
}
static inline hdd dd_two_prod(double a, double b) {
    hdd r; r.hi = a * b;
    double ah, al, bh, bl;
    dd_split(a, &ah, &al); dd_split(b, &bh, &bl);
    r.lo = ((ah * bh - r.hi) + ah * bl + al * bh) + al * bl;
    return r;
}
static inline hdd dd_from(double a) { hdd r; r.hi = a; r.lo = 0.0; return r; }
static inline hdd dd_neg(hdd a) { a.hi = -a.hi; a.lo = -a.lo; return a; }
static inline hdd dd_add(hdd a, hdd b) {
    hdd s = dd_two_sum(a.hi, b.hi);
    hdd t = dd_two_sum(a.lo, b.lo);
    s.lo += t.hi;
    s = dd_quick(s.hi, s.lo);
    s.lo += t.lo;
    return dd_quick(s.hi, s.lo);
}
static inline hdd dd_sub(hdd a, hdd b) { return dd_add(a, dd_neg(b)); }
static inline hdd dd_mul(hdd a, hdd b) {
    hdd p = dd_two_prod(a.hi, b.hi);
    p.lo += a.hi * b.lo + a.lo * b.hi;
    return dd_quick(p.hi, p.lo);
}
static inline hdd dd_mul_d(hdd a, double b) {
    hdd p = dd_two_prod(a.hi, b);
    p.lo += a.lo * b;
    return dd_quick(p.hi, p.lo);
}
static inline hdd dd_div(hdd a, hdd b) {
    double q1 = a.hi / b.hi;
    hdd r = dd_sub(a, dd_mul_d(b, q1));
    double q2 = r.hi / b.hi;
    r = dd_sub(r, dd_mul_d(b, q2));
    double q3 = r.hi / b.hi;
    hdd q = dd_quick(q1, q2);
    return dd_add(q, dd_from(q3));
}
/* (hi + k) + lo as an exact-ish double-double: the parameter's own low part
 * and the integer shift, with no rounding beyond the dd's own. */
static inline hdd dd_param(double hi, double lo, double k) {
    hdd s = dd_two_sum(hi, k);
    return dd_add(s, dd_from(lo));
}

/* Complex double-double. */
typedef struct { hdd re, im; } cdd;

static inline cdd cdd_mul(cdd x, cdd y) {
    cdd r;
    r.re = dd_sub(dd_mul(x.re, y.re), dd_mul(x.im, y.im));
    r.im = dd_add(dd_mul(x.re, y.im), dd_mul(x.im, y.re));
    return r;
}
static inline cdd cdd_div(cdd x, cdd y) {
    hdd den = dd_add(dd_mul(y.re, y.re), dd_mul(y.im, y.im));
    cdd n;
    n.re = dd_add(dd_mul(x.re, y.re), dd_mul(x.im, y.im));
    n.im = dd_sub(dd_mul(x.im, y.re), dd_mul(x.re, y.im));
    n.re = dd_div(n.re, den);
    n.im = dd_div(n.im, den);
    return n;
}
static inline cdd cdd_add(cdd x, cdd y) {
    cdd r; r.re = dd_add(x.re, y.re); r.im = dd_add(x.im, y.im); return r;
}

/* ------------------------------------------------------------------ */
/*  Raw double complex helpers                                         */
/* ------------------------------------------------------------------ */

static inline void c_mul(double ar, double ai, double br, double bi,
                         double* rr, double* ri) {
    double r = ar * br - ai * bi;
    double i = ar * bi + ai * br;
    *rr = r; *ri = i;
}
/* Smith's algorithm: no overflow for representable quotients. */
static inline void c_div(double ar, double ai, double br, double bi,
                         double* rr, double* ri) {
    if (fabs(br) >= fabs(bi)) {
        double t = bi / br, d = br + bi * t;
        double r = (ar + ai * t) / d, i = (ai - ar * t) / d;
        *rr = r; *ri = i;
    } else {
        double t = br / bi, d = br * t + bi;
        double r = (ar * t + ai) / d, i = (ai * t - ar) / d;
        *rr = r; *ri = i;
    }
}
static inline double c_abs(double r, double i) { return hypot(r, i); }
/* Principal log. Near |z| = 1, log(hypot) would lose the small modulus to
 * cancellation, so 0.5*log1p((x-1)(x+1)+y^2) is used there. */
static void c_log(double zr, double zi, double* lr, double* li) {
    double ax = fabs(zr), ay = fabs(zi);
    double big = ax > ay ? ax : ay, sml = ax > ay ? ay : ax;
    double m;
    if (big > 0.5 && big < 2.0)
        m = 0.5 * log1p((big - 1.0) * (big + 1.0) + sml * sml);
    else
        m = log(hypot(zr, zi));
    *lr = m;
    *li = atan2(zi, zr);
}
static inline void c_exp(double zr, double zi, double* rr, double* ri) {
    double e = exp(zr);
    if (zi == 0.0) { *rr = e; *ri = 0.0 * zi; return; }
    *rr = e * cos(zi);
    *ri = e * sin(zi);
}
/* w^s for real s, w given by its modulus and principal argument. The modulus
 * goes through pow() rather than exp(s log|w|): the latter multiplies the
 * rounding error of log|w| by |s log|w||, which cost several ulp in the
 * transformation prefactors (1-z)^-a. On the real axis (arg = 0) the result is
 * pow's ~0.5 ulp. */
static double hy_cospi(double x);
static inline void c_pow_ma(double mod, double arg, double s, double* rr, double* ri) {
    double m = pow(mod, s);
    if (arg == 0.0) { *rr = m; *ri = 0.0 * arg; return; }
    if (arg == M_PI || arg == -M_PI) {
        /* The negative real axis: the phase is exactly +-pi s, so reduce it
         * exactly (an integer s gives a real result, not one with a 1e-16
         * imaginary residue from sin(fl(pi) s)). */
        double sn = hy_cospi(s - 0.5);                 /* sin(pi s) */
        *rr = m * hy_cospi(s);
        *ri = (arg > 0.0 ? m : -m) * sn;
        return;
    }
    double ph = s * arg;
    *rr = m * cos(ph);
    *ri = m * sin(ph);
}

/* ------------------------------------------------------------------ */
/*  The generic pFq series                                             */
/* ------------------------------------------------------------------ */

/* Parameters of one series: hi parts plus optional low parts (a parameter that
 * is itself a rounded combination, like c-b, is carried as hi+lo so the
 * factor (a+k) does not inherit a SYSTEMATIC rounding error), complex parts
 * optional. NULL arrays mean zero. */
typedef struct {
    const double *ar, *alo, *ai;
    size_t p;
    const double *br, *blo, *bi;
    size_t q;
} HySer;

typedef struct {
    bool terminating;    /* finite: some upper parameter is a non-positive integer */
    long kmin;           /* never stop before this index                            */
} HySerInfo;

static double hs_ar(const HySer* s, size_t i)  { return s->ar[i]; }
static double hs_alo(const HySer* s, size_t i) { return s->alo ? s->alo[i] : 0.0; }
static double hs_ai(const HySer* s, size_t i)  { return s->ai ? s->ai[i] : 0.0; }
static double hs_br(const HySer* s, size_t j)  { return s->br[j]; }
static double hs_blo(const HySer* s, size_t j) { return s->blo ? s->blo[j] : 0.0; }
static double hs_bi(const HySer* s, size_t j)  { return s->bi ? s->bi[j] : 0.0; }

/* Classify the series: termination, the radius-of-convergence gate, and kmin --
 * the index past every negative parameter, before which a term that happens to
 * be tiny (a+k near 0) must not end the sum. Returns false for a lower-parameter
 * pole inside the summation range or a divergent / out-of-disc series. */
static bool hy_ser_prep(const HySer* s, double zabs, HySerInfo* inf) {
    inf->terminating = false;
    double nterm = 0.0, neg = 0.0;
    for (size_t i = 0; i < s->p; i++) {
        double a = hs_ar(s, i) + hs_alo(s, i), ai = hs_ai(s, i);
        if (ai == 0.0 && hs_alo(s, i) == 0.0 && hy_is_nonpos_int(a)) {
            if (!inf->terminating || -a < nterm) nterm = -a;
            inf->terminating = true;
        } else if (a < 0.0 && -a > neg) {
            neg = -a;
        }
    }
    for (size_t j = 0; j < s->q; j++) {
        double b = hs_br(s, j) + hs_blo(s, j), bi = hs_bi(s, j);
        if (bi == 0.0 && hs_blo(s, j) == 0.0 && hy_is_nonpos_int(b)) {
            /* (b)_k hits zero at k = -b: a pole unless the series ended first. */
            if (!inf->terminating || -b < nterm) return false;
        } else if (b < 0.0 && -b > neg) {
            neg = -b;
        }
    }
    if (inf->terminating) {
        if (nterm > (double)HY_MAXTERMS) return false;
    } else {
        if (s->p > s->q + 1) return false;                 /* divergent */
        if (s->p == s->q + 1 && !(zabs < 1.0)) return false;
    }
    double km = ceil(neg) + 1.0;
    if (km > (double)HY_MAXTERMS) return false;
    inf->kmin = (long)km;
    return true;
}

/* Real double sum. *peak is the largest |term| (term 0 = 1 included). */
/* Summation is compensated (Knuth's TwoSum carried in `comp`): for a slowly
 * decaying series the recursive-summation error, not the terms' own rounding,
 * is what dominates -- 2F1(1,2;3;0.89) sums ~380 terms and the plain sum was
 * 13 ulp off, the compensated one 1.4. */
static bool hy_sum_r(const HySer* s, const HySerInfo* inf, double z,
                     double* out, double* peak_out, long* kpk_out) {
    double term = 1.0, sum = 1.0, comp = 0.0, peak = 1.0;
    long kpk = 0;
    double rho_z = (s->p == s->q + 1) ? fabs(z) : 0.0;
    for (long k = 0; k < HY_MAXTERMS; k++) {
        double kd = (double)k;
        double num = z, den = kd + 1.0;
        for (size_t i = 0; i < s->p; i++) num *= (s->ar[i] + kd) + hs_alo(s, i);
        if (num == 0.0) break;                        /* terminated */
        for (size_t j = 0; j < s->q; j++) den *= (s->br[j] + kd) + hs_blo(s, j);
        if (den == 0.0) return false;
        double r = num / den;
        term *= r;
        hdd q = dd_two_sum(sum, term);
        sum = q.hi; comp += q.lo;
        double at = fabs(term);
        if (at > peak) { peak = at; kpk = k + 1; }
        if (!isfinite(sum)) return false;
        if (!inf->terminating && k + 1 >= inf->kmin) {
            double rho = fabs(r);
            if (rho < rho_z) rho = rho_z;
            if (rho < 1.0 && at <= fabs(sum) * HY_TOL * (1.0 - rho)) {
                *out = sum + comp; *peak_out = peak; *kpk_out = kpk; return true;
            }
        }
        if (k == HY_MAXTERMS - 1) return false;
    }
    *out = sum + comp; *peak_out = peak; *kpk_out = kpk;
    return true;
}

/* The same series in double-double. */
static bool hy_sum_r_dd(const HySer* s, const HySerInfo* inf, double z, double* out) {
    hdd term = dd_from(1.0), sum = dd_from(1.0);
    double rho_z = (s->p == s->q + 1) ? fabs(z) : 0.0;
    for (long k = 0; k < HY_DD_MAXTERMS; k++) {
        double kd = (double)k;
        hdd num = dd_from(z), den = dd_from(kd + 1.0);
        for (size_t i = 0; i < s->p; i++)
            num = dd_mul(num, dd_param(s->ar[i], hs_alo(s, i), kd));
        if (num.hi == 0.0) break;
        for (size_t j = 0; j < s->q; j++)
            den = dd_mul(den, dd_param(s->br[j], hs_blo(s, j), kd));
        if (den.hi == 0.0) return false;
        hdd r = dd_div(num, den);
        term = dd_mul(term, r);
        sum = dd_add(sum, term);
        double at = fabs(term.hi);
        if (!isfinite(sum.hi) || at > 1e290) return false;
        if (!inf->terminating && k + 1 >= inf->kmin) {
            double rho = fabs(r.hi);
            if (rho < rho_z) rho = rho_z;
            if (rho < 1.0 && at <= fabs(sum.hi) * HY_DD_TOL * (1.0 - rho)) {
                *out = sum.hi + sum.lo; return isfinite(*out);
            }
        }
        if (k == HY_DD_MAXTERMS - 1) return false;
    }
    *out = sum.hi + sum.lo;
    return isfinite(*out);
}

/* Real series with the loss-measured double-double fallback. *loss (optional)
 * receives peak/|sum| of the double pass. */
static bool hy_series_r(const HySer* s, double z, double* out) {
    HySerInfo inf;
    if (!hy_ser_prep(s, fabs(z), &inf)) return false;
    double v, peak;
    long kpk;
    if (!hy_sum_r(s, &inf, z, &v, &peak, &kpk)) return false;
    /* Two ways the double sum can fall short: cancellation (peak/|sum|), and
     * a peak deep in the series, where each term carries the rounding of the
     * ~kpk ratio products before it (1F1(1.4;1.7;29) summed in double was 10
     * ulp off with no cancellation at all). Either sends it to the dd re-sum. */
    if (peak <= HY_LOSS_OK * fabs(v) && kpk <= HY_KPEAK_DD) { *out = v; return true; }
    if (!(peak <= HY_LOSS_MAX * fabs(v))) return false;
    return hy_sum_r_dd(s, &inf, z, out);
}

/* Complex double sum (parameters may be complex). Magnitudes are max-norms. */
static inline double c_mag(double r, double i) {
    double a = fabs(r), b = fabs(i);
    return a > b ? a : b;
}

static bool hy_sum_c(const HySer* s, const HySerInfo* inf, double zr, double zi,
                     bool real_params, double* orr, double* ori, double* peak_out,
                     long* kpk_out) {
    double tr = 1.0, ti = 0.0, sr = 1.0, si = 0.0, cr = 0.0, ci = 0.0, peak = 1.0;
    long kpk = 0;
    double rho_z = (s->p == s->q + 1) ? c_abs(zr, zi) : 0.0;
    for (long k = 0; k < HY_MAXTERMS; k++) {
        double kd = (double)k, rr, ri;
        if (real_params) {
            double num = 1.0, den = kd + 1.0;
            for (size_t i = 0; i < s->p; i++) num *= (s->ar[i] + kd) + hs_alo(s, i);
            if (num == 0.0) break;
            for (size_t j = 0; j < s->q; j++) den *= (s->br[j] + kd) + hs_blo(s, j);
            if (den == 0.0) return false;
            double f = num / den;
            rr = zr * f; ri = zi * f;
        } else {
            double nr = zr, ni = zi, dr = kd + 1.0, di = 0.0;
            for (size_t i = 0; i < s->p; i++)
                c_mul(nr, ni, (s->ar[i] + kd) + hs_alo(s, i), hs_ai(s, i), &nr, &ni);
            if (nr == 0.0 && ni == 0.0) break;
            for (size_t j = 0; j < s->q; j++)
                c_mul(dr, di, (s->br[j] + kd) + hs_blo(s, j), hs_bi(s, j), &dr, &di);
            if (dr == 0.0 && di == 0.0) return false;
            c_div(nr, ni, dr, di, &rr, &ri);
        }
        if (rr == 0.0 && ri == 0.0) break;            /* z == 0 */
        c_mul(tr, ti, rr, ri, &tr, &ti);
        hdd qr = dd_two_sum(sr, tr), qi = dd_two_sum(si, ti);
        sr = qr.hi; cr += qr.lo; si = qi.hi; ci += qi.lo;
        double at = c_mag(tr, ti);
        if (at > peak) { peak = at; kpk = k + 1; }
        if (!isfinite(sr) || !isfinite(si)) return false;
        if (!inf->terminating && k + 1 >= inf->kmin) {
            double rho = c_abs(rr, ri);
            if (rho < rho_z) rho = rho_z;
            if (rho < 1.0 && at <= c_mag(sr, si) * HY_TOL * (1.0 - rho)) {
                *orr = sr + cr; *ori = si + ci; *peak_out = peak; *kpk_out = kpk;
                return true;
            }
        }
        if (k == HY_MAXTERMS - 1) return false;
    }
    *orr = sr + cr; *ori = si + ci; *peak_out = peak; *kpk_out = kpk;
    return true;
}

/* Real parameters, complex z: the ratio (prod (a_i+k)) / (prod (b_j+k) (k+1)) is
 * a real double-double, and the term update is one complex-dd by real-dd
 * product and one complex-dd by complex-double product -- a fraction of the
 * general complex-parameter path below. Same series, same stop rule. */
static bool hy_sum_c_dd_rp(const HySer* s, const HySerInfo* inf, double zr, double zi,
                           double* orr, double* ori) {
    hdd tr = dd_from(1.0), ti = dd_from(0.0), sr = dd_from(1.0), si = dd_from(0.0);
    double rho_z = (s->p == s->q + 1) ? c_abs(zr, zi) : 0.0;
    double az = c_abs(zr, zi);
    for (long k = 0; k < HY_DD_MAXTERMS; k++) {
        double kd = (double)k;
        hdd num = dd_from(1.0), den = dd_from(kd + 1.0);
        for (size_t i = 0; i < s->p; i++)
            num = dd_mul(num, dd_param(s->ar[i], hs_alo(s, i), kd));
        if (num.hi == 0.0) break;
        for (size_t j = 0; j < s->q; j++)
            den = dd_mul(den, dd_param(s->br[j], hs_blo(s, j), kd));
        if (den.hi == 0.0) return false;
        hdd f = dd_div(num, den);
        /* (tr + i ti) * f * (zr + i zi) */
        hdd ur = dd_mul(tr, f), ui = dd_mul(ti, f);
        tr = dd_sub(dd_mul_d(ur, zr), dd_mul_d(ui, zi));
        ti = dd_add(dd_mul_d(ur, zi), dd_mul_d(ui, zr));
        sr = dd_add(sr, tr); si = dd_add(si, ti);
        double at = c_mag(tr.hi, ti.hi);
        if (!isfinite(sr.hi) || !isfinite(si.hi) || at > 1e290) return false;
        if (!inf->terminating && k + 1 >= inf->kmin) {
            double rho = fabs(f.hi) * az;
            if (rho < rho_z) rho = rho_z;
            if (rho < 1.0 && at <= c_mag(sr.hi, si.hi) * HY_DD_TOL * (1.0 - rho)) {
                *orr = sr.hi + sr.lo; *ori = si.hi + si.lo;
                return isfinite(*orr) && isfinite(*ori);
            }
        }
        if (k == HY_DD_MAXTERMS - 1) return false;
    }
    *orr = sr.hi + sr.lo; *ori = si.hi + si.lo;
    return isfinite(*orr) && isfinite(*ori);
}

static bool hy_sum_c_dd(const HySer* s, const HySerInfo* inf, double zr, double zi,
                        double* orr, double* ori) {
    cdd term, sum, zz;
    term.re = dd_from(1.0); term.im = dd_from(0.0);
    sum = term;
    zz.re = dd_from(zr); zz.im = dd_from(zi);
    double rho_z = (s->p == s->q + 1) ? c_abs(zr, zi) : 0.0;
    for (long k = 0; k < HY_DD_MAXTERMS; k++) {
        double kd = (double)k;
        cdd num = zz, den;
        den.re = dd_from(kd + 1.0); den.im = dd_from(0.0);
        for (size_t i = 0; i < s->p; i++) {
            cdd f; f.re = dd_param(s->ar[i], hs_alo(s, i), kd); f.im = dd_from(hs_ai(s, i));
            num = cdd_mul(num, f);
        }
        if (num.re.hi == 0.0 && num.im.hi == 0.0) break;
        for (size_t j = 0; j < s->q; j++) {
            cdd f; f.re = dd_param(s->br[j], hs_blo(s, j), kd); f.im = dd_from(hs_bi(s, j));
            den = cdd_mul(den, f);
        }
        if (den.re.hi == 0.0 && den.im.hi == 0.0) return false;
        cdd r = cdd_div(num, den);
        term = cdd_mul(term, r);
        sum = cdd_add(sum, term);
        double at = c_mag(term.re.hi, term.im.hi);
        if (!isfinite(sum.re.hi) || !isfinite(sum.im.hi) || at > 1e290) return false;
        if (!inf->terminating && k + 1 >= inf->kmin) {
            double rho = c_abs(r.re.hi, r.im.hi);
            if (rho < rho_z) rho = rho_z;
            if (rho < 1.0 && at <= c_mag(sum.re.hi, sum.im.hi) * HY_DD_TOL * (1.0 - rho)) {
                *orr = sum.re.hi + sum.re.lo; *ori = sum.im.hi + sum.im.lo;
                return isfinite(*orr) && isfinite(*ori);
            }
        }
        if (k == HY_DD_MAXTERMS - 1) return false;
    }
    *orr = sum.re.hi + sum.re.lo; *ori = sum.im.hi + sum.im.lo;
    return isfinite(*orr) && isfinite(*ori);
}

static bool hy_ser_real_params(const HySer* s) {
    for (size_t i = 0; i < s->p; i++) if (hs_ai(s, i) != 0.0) return false;
    for (size_t j = 0; j < s->q; j++) if (hs_bi(s, j) != 0.0) return false;
    return true;
}

static bool hy_series_c(const HySer* s, double zr, double zi, double* orr, double* ori) {
    HySerInfo inf;
    if (!hy_ser_prep(s, c_abs(zr, zi), &inf)) return false;
    double vr, vi, peak;
    long kpk;
    if (!hy_sum_c(s, &inf, zr, zi, hy_ser_real_params(s), &vr, &vi, &peak, &kpk))
        return false;
    double m = c_mag(vr, vi);
    if (peak <= HY_LOSS_OK * m && kpk <= HY_KPEAK_DD) { *orr = vr; *ori = vi; return true; }
    if (!(peak <= HY_LOSS_MAX * m)) return false;
    if (hy_ser_real_params(s)) return hy_sum_c_dd_rp(s, &inf, zr, zi, orr, ori);
    return hy_sum_c_dd(s, &inf, zr, zi, orr, ori);
}

/* Convenience builders. */
static bool hy_ser1_r(size_t p, const double* a, const double* alo,
                      size_t q, const double* b, const double* blo,
                      double z, double* out) {
    HySer s = { a, alo, NULL, p, b, blo, NULL, q };
    return hy_series_r(&s, z, out);
}
static bool hy_ser1_c(size_t p, const double* a, const double* alo, const double* ai,
                      size_t q, const double* b, const double* blo, const double* bi,
                      double zr, double zi, double* rr, double* ri) {
    HySer s = { a, alo, ai, p, b, blo, bi, q };
    return hy_series_c(&s, zr, zi, rr, ri);
}

/* ------------------------------------------------------------------ */
/*  0F1(; b; z)                                                        */
/* ------------------------------------------------------------------ */

/* Miller backward recurrence for z = -y < 0, b > 0, in a form that never sees
 * the rounded t = 2 sqrt(y). With nu = b - 1 and u_n = J_{nu+n}(t)/(t/2)^n
 * (t = 2 sqrt y), the Bessel recurrence becomes
 *     u_{n-1} = (nu + n) u_n - y u_{n+1},
 * u_n is its minimal solution, and the Neumann series
 *     (t/2)^nu = sum_k (nu+2k) Gamma(nu+k)/k! J_{nu+2k}(t)        [Gautschi 1967]
 * normalises it. Dividing through by Gamma(nu+1) = Gamma(b):
 *     0F1(; b; -y) = u_0 / sum_k w_k y^k u_{2k},
 *     w_0 = 1,  w_{k+1}/w_k = (nu+2k+2)/(nu+2k) * (nu+k)/(k+1)  (k >= 1),
 *     w_1 = nu + 2,
 * every weight positive for nu > -1 (b > 0). No Gamma function is needed and the
 * recurrence is backward-stable, so the result is accurate relative to the
 * envelope of the oscillation -- unlike the alternating power series, which loses
 * ~2 sqrt(y) / ln 2 bits. The normalising sum is accumulated Horner-wise on the
 * way down, with rescaling to keep the dominant solution in range. */
static bool hy0f1_miller(double b, double y, double* out) {
    double nu = b - 1.0;
    double t = 2.0 * sqrt(y);
    double nstart = t + 30.0 + 8.0 * cbrt(t);
    if (nstart > 20000.0) return false;
    long N = (long)nstart;
    if (N & 1) N++;
    double up1 = 0.0, u = 1.0e-200;               /* u_{N+1}, u_N */
    /* S' accumulates sum_{k>=K} w_k y^k u_{2k} / (w_K y^K) at even n = 2K. */
    double S = u;                                  /* at n = N (K = N/2)       */
    for (long n = N; n >= 1; n--) {
        double um1 = (nu + (double)n) * u - y * up1;   /* u_{n-1} */
        up1 = u; u = um1;
        long m = n - 1;
        if ((m & 1) == 0) {                         /* even index m = 2K       */
            long K = m / 2;
            double ratio;                           /* w_{K+1}/w_K             */
            if (K == 0) ratio = nu + 2.0;
            else ratio = ((nu + 2.0 * K + 2.0) * (nu + (double)K)) /
                         ((nu + 2.0 * K) * ((double)K + 1.0));
            S = u + S * (ratio * y);
        }
        double au = fabs(u) > fabs(S) ? fabs(u) : fabs(S);
        if (au > 1.0e200) {                         /* rescale everything live */
            u *= 1.0e-200; up1 *= 1.0e-200; S *= 1.0e-200;
        }
    }
    double v = u / S;
    if (!isfinite(v)) return false;
    *out = v;
    return true;
}

/* Large positive z: 0F1(; b; y) = Gamma(b) y^((1-b)/2) I_nu(t), t = 2 sqrt(y),
 * nu = b - 1, with the Hankel expansion [DLMF 10.40.1]
 *     I_nu(t) ~ e^t / sqrt(2 pi t) * sum_k (-1)^k a_k(nu) / t^k.
 * t is formed in double-double (sqrt plus its exact residual) and e^t as
 * exp(t_hi)(1 + t_lo): e^t amplifies an error in t by t, so a rounded t would
 * cost ~t/2 ulp. Accepted only when the optimally truncated series has a least
 * term below 2^-54 of the sum. */
static bool hy0f1_asym_pos(double b, double y, double* out) {
    double nu = b - 1.0;
    double s = sqrt(y);
    hdd ss = dd_two_prod(s, s);                         /* s*s exactly */
    double corr = ((y - ss.hi) - ss.lo) / (2.0 * s);    /* sqrt residual */
    double th = 2.0 * s, tl = 2.0 * corr;
    if (th > 700.0) return false;
    double mu = 4.0 * nu * nu;
    double term = 1.0, sum = 1.0, last = 1.0;
    bool ok = false;
    for (int k = 1; k < 200; k++) {
        double j2 = (2.0 * k - 1.0);
        term *= -(mu - j2 * j2) / (8.0 * k * th);
        double at = fabs(term);
        if (term == 0.0) { ok = true; break; }
        if (at > last) { ok = (last <= fabs(sum) * HY_TOL); break; }
        sum += term;
        last = at;
        if (at <= fabs(sum) * HY_TOL) { ok = true; break; }
    }
    if (!ok) return false;
    double g = tgamma(b);
    if (!isfinite(g)) return false;
    double et = exp(th) * (1.0 + tl);
    double v = g * pow(y, -0.5 * nu) * et / sqrt(2.0 * M_PI * th) * sum;
    if (!isfinite(v)) return false;
    *out = v;
    return true;
}

bool sf_machine_hyper0f1(double b, double z, double* out) {
    if (!isfinite(b) || !isfinite(z)) return false;
    if (hy_is_nonpos_int(b)) return false;                 /* (b)_k pole */
    if (z == 0.0) { *out = 1.0; return true; }
    if (z > 0.0) {
        if (b > 0.0 && z > 900.0 && hy0f1_asym_pos(b, z, out)) return true;
        return hy_ser1_r(0, NULL, NULL, 1, &b, NULL, z, out);
    }
    /* z < 0: the series alternates. Keep the double sum where it loses at most
     * a bit; otherwise the dd re-sum (<= 1 ulp up to 2^45 of cancellation, i.e.
     * |z| up to ~250); beyond that Miller's recurrence for b > 0, which is
     * accurate to a few ulp of the oscillation's envelope. */
    HySer s = { NULL, NULL, NULL, 0, &b, NULL, NULL, 1 };
    HySerInfo inf;
    if (!hy_ser_prep(&s, -z, &inf)) return false;
    double v = 0.0, peak = INFINITY;
    long kpk = 0;
    if (hy_sum_r(&s, &inf, z, &v, &peak, &kpk) && peak <= HY_LOSS_OK * fabs(v) &&
        kpk <= HY_KPEAK_DD) {
        *out = v; return true;
    }
    if (peak <= HY_LOSS_MAX * fabs(v) && hy_sum_r_dd(&s, &inf, z, out)) return true;
    if (b > 0.0) return hy0f1_miller(b, -z, out);
    return false;
}

bool sf_machine_hyper0f1_c(double br, double bi, double zr, double zi,
                           double* rr, double* ri) {
    if (!isfinite(br) || !isfinite(bi) || !isfinite(zr) || !isfinite(zi)) return false;
    if (bi == 0.0 && hy_is_nonpos_int(br)) return false;
    if (zr == 0.0 && zi == 0.0) { *rr = 1.0; *ri = 0.0; return true; }
    return hy_ser1_c(0, NULL, NULL, NULL, 1, &br, NULL, &bi, zr, zi, rr, ri);
}

/* ------------------------------------------------------------------ */
/*  1F1(a; b; z)                                                       */
/* ------------------------------------------------------------------ */

/* Asymptotic series sum_s (al)_s (be)_s / s! t^s summed to its least term.
 * *err receives the magnitude of the last term kept when the series had to be
 * truncated (0 when it terminated or converged below 2^-56). */
static void hy_asym_sum(double al, double be, double t, double* S, double* err) {
    double term = 1.0, sum = 1.0, last = 1.0;
    *err = 0.0;
    for (int s = 0; s < 400; s++) {
        double next = term * (al + s) * (be + s) / (s + 1.0) * t;
        if (next == 0.0) { *S = sum; return; }
        double an = fabs(next);
        if (!isfinite(next)) { *S = sum; *err = INFINITY; return; }
        if (an > last) { *S = sum; *err = last; return; }
        sum += next;
        term = next;
        last = an;
        if (an <= fabs(sum) * 1.3877787807814457e-17) { *S = sum; return; }
    }
    *S = sum; *err = last;
}

/* DLMF 13.7.2 for real x, |x| large, both exponential sectors kept (the
 * subdominant one with the Stokes-line average cos factor), which makes the
 * formula exact in the terminating cases and correct to double beyond them:
 *   x = +y > 0:  M = G(b)/G(a)   e^y y^(a-b) S(1-a, b-a; 1/y)
 *                  + G(b)/G(b-a) cos(pi a) y^(-a) S(a, a-b+1; -1/y)
 *   x = -y < 0:  M = G(b)/G(b-a) y^(-a) S(a, a-b+1; 1/y)
 *                  + G(b)/G(a)   cos(pi(b-a)) e^(-y) y^(a-b) S(1-a, b-a; -1/y)
 * y^(a-b) is formed as y^a * y^-b so a rounded a-b cannot be amplified by ln y. */
static bool hy1f1_asym(double a, double b, double x, double* out) {
    double y = fabs(x);
    if (fabs(a) > 160.0 || fabs(b) > 160.0 || fabs(b - a) > 160.0) return false;
    double gb = tgamma(b);
    double rga = hy_rgamma(a), rgba = hy_rgamma(b - a);
    if (!isfinite(gb) || !isfinite(rga) || !isfinite(rgba)) return false;
    double S1, e1, S2, e2;
    double D = 0.0, E = 0.0, errD = 0.0, errE = 0.0;
    if (x > 0.0) {
        if (rga != 0.0) {
            if (y > 700.0) return false;
            hy_asym_sum(1.0 - a, b - a, 1.0 / y, &S2, &e2);
            double pre = gb * rga * exp(y) * pow(y, a) * pow(y, -b);
            D = pre * S2; errD = fabs(pre) * e2;
        }
        if (rgba != 0.0) {
            hy_asym_sum(a, a - b + 1.0, -1.0 / y, &S1, &e1);
            double pre = gb * rgba * hy_cospi(a) * pow(y, -a);
            E = pre * S1; errE = fabs(pre) * e1;
        }
    } else {
        if (rgba != 0.0) {
            hy_asym_sum(a, a - b + 1.0, 1.0 / y, &S1, &e1);
            double pre = gb * rgba * pow(y, -a);
            D = pre * S1; errD = fabs(pre) * e1;
        }
        if (rga != 0.0 && y < 745.0) {
            hy_asym_sum(1.0 - a, b - a, -1.0 / y, &S2, &e2);
            double pre = gb * rga * hy_cospi(b - a) * exp(-y) * pow(y, a) * pow(y, -b);
            E = pre * S2; errE = fabs(pre) * e2;
        }
    }
    double v = D + E;
    if (!isfinite(v) || v == 0.0) return false;
    if (errD + errE > fabs(v) * HY_TOL) return false;          /* not converged */
    if (fabs(D) + fabs(E) > 2.0 * fabs(v)) return false;       /* cancellation */
    *out = v;
    return true;
}

bool sf_machine_hyper1f1(double a, double b, double x, double* out) {
    if (!isfinite(a) || !isfinite(b) || !isfinite(x)) return false;
    bool aterm = hy_is_nonpos_int(a);
    if (hy_is_nonpos_int(b) && !(aterm && a >= b)) return false;   /* pole */
    if (x == 0.0 || a == 0.0) { *out = 1.0; return true; }
    if (a == b) { *out = exp(x); return isfinite(*out); }
    if (a == 1.0 && b == 2.0) {                     /* (e^x - 1)/x, elementary */
        *out = expm1(x) / x;
        return isfinite(*out);
    }
    if (aterm)                                      /* a polynomial in x */
        return hy_ser1_r(1, &a, NULL, 1, &b, NULL, x, out);
    if (fabs(x) > 20.0 && hy1f1_asym(a, b, x, out)) return true;
    if (x > 0.0)
        return hy_ser1_r(1, &a, NULL, 1, &b, NULL, x, out);
    /* x < 0: Kummer's transformation [AS 13.1.27] M(a,b,x) = e^x M(b-a,b,-x),
     * whose series is one-signed when b > a, b > 0. b - a is carried as hi+lo. */
    hdd bma = dd_two_sum(b, -a);
    double v;
    if (!hy_ser1_r(1, &bma.hi, &bma.lo, 1, &b, NULL, -x, &v)) return false;
    double r = exp(x) * v;
    if (!isfinite(r)) return false;
    *out = r;
    return true;
}

bool sf_machine_hyper1f1_c(double ar, double ai, double br, double bi,
                           double zr, double zi, double* rr, double* ri) {
    if (!isfinite(ar) || !isfinite(ai) || !isfinite(br) || !isfinite(bi) ||
        !isfinite(zr) || !isfinite(zi)) return false;
    bool real_ab = (ai == 0.0 && bi == 0.0);
    bool aterm = real_ab && hy_is_nonpos_int(ar);
    if (bi == 0.0 && hy_is_nonpos_int(br) && !(aterm && ar >= br)) return false;
    if ((zr == 0.0 && zi == 0.0) || (ar == 0.0 && ai == 0.0)) {
        *rr = 1.0; *ri = 0.0; return true;
    }
    if (ar == br && ai == bi) {
        c_exp(zr, zi, rr, ri);
        return isfinite(*rr) && isfinite(*ri);
    }
    if (real_ab && ar == 1.0 && br == 2.0) {
        /* (e^z - 1)/z with e^z - 1 = expm1(x) cos y - 2 sin^2(y/2) + i e^x sin y,
         * free of cancellation at small |z|. */
        double sh = sin(0.5 * zi);
        double em = expm1(zr);
        double nr = em * cos(zi) - 2.0 * sh * sh;
        double ni = exp(zr) * sin(zi);
        c_div(nr, ni, zr, zi, rr, ri);
        return isfinite(*rr) && isfinite(*ri);
    }
    if (zr >= 0.0 || aterm)
        return hy_ser1_c(1, &ar, NULL, &ai, 1, &br, NULL, &bi, zr, zi, rr, ri);
    /* Re z < 0: Kummer, M(a,b,z) = e^z M(b-a, b, -z). */
    hdd bma = dd_two_sum(br, -ar);
    double bmai = bi - ai, vr, vi, er, ei;
    if (!hy_ser1_c(1, &bma.hi, &bma.lo, &bmai, 1, &br, NULL, &bi, -zr, -zi, &vr, &vi))
        return false;
    c_exp(zr, zi, &er, &ei);
    c_mul(er, ei, vr, vi, rr, ri);
    return isfinite(*rr) && isfinite(*ri);
}

/* ------------------------------------------------------------------ */
/*  2F1(a, b; c; z)                                                    */
/* ------------------------------------------------------------------ */

/* Gamma(x) with an explicit decline value: NaN when it is a pole or overflows. */
static double hy_gamma(double x) {
    if (hy_is_nonpos_int(x)) return NAN;
    return tgamma(x);
}

/* psi(x) - psi(y). When x - y is a small integer j the difference is the finite
 * sum sum_{i<j} 1/(y+i) (no rounding of two large digammas to cancel); otherwise
 * the two digammas, and *err accumulates the scale (|psi(x)| + |psi(y)|) their
 * rounding is relative to. */
static bool hy_psi_diff(double x, double y, double* d, double* err) {
    double j = x - y;
    if (j == nearbyint(j) && fabs(j) <= 64.0) {
        double s = 0.0;
        if (j >= 0.0) for (double i = 0.0; i < j; i += 1.0) s += 1.0 / (y + i);
        else          for (double i = 0.0; i < -j; i += 1.0) s -= 1.0 / (x + i);
        *d = s;
        return isfinite(s);
    }
    double px, py;
    if (!sf_machine_digamma(x, &px) || !sf_machine_digamma(y, &py)) return false;
    *d = px - py;
    *err += fabs(px) + fabs(py);
    return true;
}

/* Connection of F(a,b;c;.) at 1-w, i.e. F(a, b; c; 1 - w) expressed through
 * series in w, for complex w (|w| < 1) given with its principal log (lwr, lwi),
 * and REAL parameters. [AS 15.3.6] in general, [AS 15.3.10-15.3.12] when
 * s = c - a - b is an integer m. Writes the value and a condition estimate
 * kappa = (sum of |pieces|)/|value| -- how much the combination amplifies the
 * pieces' own rounding. Callers must already have peeled off terminating a, b. */
static bool hy2f1_conn(double a, double b, double c,
                       double wr, double wi, double lwr, double lwi,
                       double* rr, double* ri, double* kappa) {
    double s = c - a - b;
    double m = nearbyint(s);
    double scale = fabs(a) + fabs(b) + fabs(c) + 1.0;
    bool degenerate = fabs(s - m) <= 4.0 * DBL_EPSILON * scale;
    double gc = hy_gamma(c);
    if (!isfinite(gc)) return false;

    if (!degenerate) {
        /* F = A1 F(a,b;1-s;w) + A2 w^s F(c-a,c-b;1+s;w)
         * A1 = G(c)G(s)/(G(c-a)G(c-b)),  A2 = G(c)G(-s)/(G(a)G(b)). */
        double A1 = gc * hy_gamma(s) * hy_rgamma(c - a) * hy_rgamma(c - b);
        double A2 = gc * hy_gamma(-s) * hy_rgamma(a) * hy_rgamma(b);
        if (!isfinite(A1) || !isfinite(A2)) return false;
        double f1r = 0.0, f1i = 0.0, f2r = 0.0, f2i = 0.0;
        if (A1 != 0.0) {
            hdd p1 = dd_two_sum(1.0, -s);                 /* 1 - s */
            double up[2] = { a, b };
            if (!hy_ser1_c(2, up, NULL, NULL, 1, &p1.hi, &p1.lo, NULL, wr, wi, &f1r, &f1i))
                return false;
        }
        if (A2 != 0.0) {
            hdd ca = dd_two_sum(c, -a), cb = dd_two_sum(c, -b), p2 = dd_two_sum(1.0, s);
            double up[2] = { ca.hi, cb.hi }, uplo[2] = { ca.lo, cb.lo };
            if (!hy_ser1_c(2, up, uplo, NULL, 1, &p2.hi, &p2.lo, NULL, wr, wi, &f2r, &f2i))
                return false;
            double pr, pi;
            c_pow_ma(c_abs(wr, wi), lwi, s, &pr, &pi);
            c_mul(f2r, f2i, pr, pi, &f2r, &f2i);
        }
        double t1r = A1 * f1r, t1i = A1 * f1i, t2r = A2 * f2r, t2i = A2 * f2i;
        double vr = t1r + t2r, vi = t1i + t2i;
        double mag = c_abs(vr, vi);
        if (!isfinite(vr) || !isfinite(vi) || mag == 0.0) return false;
        *kappa = (c_abs(t1r, t1i) + c_abs(t2r, t2i)) / mag;
        *rr = vr; *ri = vi;
        return true;
    }

    /* Degenerate: s = m integer. Unified form with M = |m|:
     *   m >= 0: (al, be) = (a+M, b+M),  G2 = G(c)/(G(a)G(b)),       factor (-1)^M w^M
     *   m <  0: (al, be) = (a, b),      G2 = G(c)/(G(a-M)G(b-M)),   factor (-1)^M
     *   log part = -factor * G2 * sum_n e_n w^n [ln w + h_n],
     *     e_n = (al)_n (be)_n / (n! (n+M)!),
     *     h_n = -psi(n+1) - psi(n+M+1) + psi(al+n) + psi(be+n);
     *   finite part (M >= 1):
     *     m > 0: G(M)G(c)/(G(a+M)G(b+M)) sum_{n<M} (a)_n (b)_n / (n! (1-M)_n) w^n
     *     m < 0: G(M)G(c)/(G(a)G(b)) w^-M sum_{n<M} (a-M)_n (b-M)_n / (n! (1-M)_n) w^n */
    long M = (long)fabs(m);
    if (M > 2000) return false;
    double Md = (double)M;
    bool pos = (m >= 0.0);
    double al = pos ? a + Md : a, be = pos ? b + Md : b;
    double G2 = gc * (pos ? hy_rgamma(a) * hy_rgamma(b)
                          : hy_rgamma(a - Md) * hy_rgamma(b - Md));
    double finr = 0.0, fini = 0.0, finmag = 0.0;
    if (M >= 1) {
        double G1 = tgamma(Md) * gc * (pos ? hy_rgamma(a + Md) * hy_rgamma(b + Md)
                                           : hy_rgamma(a) * hy_rgamma(b));
        if (!isfinite(G1)) return false;
        double pa = pos ? a : a - Md, pb = pos ? b : b - Md;
        double tr = 1.0, ti = 0.0, sr = 1.0, si = 0.0;
        for (long n = 0; n < M - 1; n++) {
            double nd = (double)n;
            double f = (pa + nd) * (pb + nd) / ((nd + 1.0) * (1.0 - Md + nd));
            c_mul(tr, ti, wr * f, wi * f, &tr, &ti);
            sr += tr; si += ti;
        }
        if (!pos) {                                     /* times w^-M */
            double pr, pi;
            c_pow_ma(c_abs(wr, wi), lwi, -Md, &pr, &pi);
            c_mul(sr, si, pr, pi, &sr, &si);
        }
        finr = G1 * sr; fini = G1 * si;
        finmag = c_abs(finr, fini);
    }
    double logr = 0.0, logi = 0.0, logmag = 0.0;
    if (G2 != 0.0) {
        if (!isfinite(G2)) return false;
        if (hy_is_nonpos_int(al) || hy_is_nonpos_int(be)) return false;
        /* h_n = [psi(al+n) - psi(n+1)] + [psi(be+n) - psi(n+M+1)], each bracket
         * carried by its own recurrence
         *   d1_{n+1} = d1_n + (1 - al)/((al+n)(n+1)),
         *   d2_{n+1} = d2_n + (M+1 - be)/((be+n)(n+M+1)),
         * so that a bracket that is identically zero (al = 1, be = M+1) stays
         * exactly zero and the others carry only their own rounding -- four
         * separately accumulated digammas would leave ~u sum|psi| of noise in
         * every h_n. err0 is the rounding of the two starting differences. */
        double d1, d2, err0 = 0.0;
        if (!hy_psi_diff(al, 1.0, &d1, &err0) || !hy_psi_diff(be, Md + 1.0, &d2, &err0))
            return false;
        /* e_0 = 1/M! */
        double e0 = 1.0 / tgamma(Md + 1.0);
        double tr = e0, ti = 0.0;                        /* e_n w^n */
        double s0r = 0.0, s0i = 0.0, s1r = 0.0, s1i = 0.0, peak1 = 0.0;
        double c0r = 0.0, c0i = 0.0, c1r = 0.0, c1i = 0.0;   /* compensation */
        double hscale = 0.0;
        double neg = 0.0;
        if (al < 0.0) neg = -al;
        if (be < 0.0 && -be > neg) neg = -be;
        long kmin = (long)ceil(neg) + 1;
        double rho_w = c_abs(wr, wi);
        bool done = false;
        for (long n = 0; n < HY_MAXTERMS; n++) {
            double nd = (double)n;
            double h = d1 + d2;
            hdd q;
            q = dd_two_sum(s0r, tr); s0r = q.hi; c0r += q.lo;
            q = dd_two_sum(s0i, ti); s0i = q.hi; c0i += q.lo;
            double hr = tr * h, hi = ti * h;
            q = dd_two_sum(s1r, hr); s1r = q.hi; c1r += q.lo;
            q = dd_two_sum(s1i, hi); s1i = q.hi; c1i += q.lo;
            double mh = c_mag(hr, hi);
            if (mh > peak1) peak1 = mh;
            hscale += c_mag(tr, ti) * (fabs(d1) + fabs(d2) + err0);
            /* advance e_n w^n and the two digamma differences */
            double f = (al + nd) * (be + nd) / ((nd + 1.0) * (nd + Md + 1.0));
            if (f == 0.0) { done = true; break; }
            c_mul(tr, ti, wr * f, wi * f, &tr, &ti);
            d1 += (1.0 - al) / ((al + nd) * (nd + 1.0));
            d2 += (Md + 1.0 - be) / ((be + nd) * (nd + Md + 1.0));
            if (!isfinite(tr) || !isfinite(s1r) || !isfinite(s1i)) return false;
            if (n + 1 >= kmin) {
                double rho = fabs(f) * rho_w;
                if (rho < rho_w) rho = rho_w;
                double at = c_mag(tr, ti);
                double hn = fabs(d1 + d2);
                if (rho < 1.0 &&
                    at * (1.0 + hn) <= HY_TOL * (1.0 - rho) *
                        (c_mag(s0r, s0i) * (fabs(lwr) + fabs(lwi)) + c_mag(s1r, s1i))) {
                    done = true; break;
                }
            }
        }
        if (!done) return false;
        s0r += c0r; s0i += c0i; s1r += c1r; s1i += c1i;
        /* inner = ln w * S0 + S1, formed in double-double: the two pieces
         * cancel (that is what kappa measures), so the product's own rounding
         * must not be what the cancellation amplifies. */
        hdd ir_ = dd_add(dd_sub(dd_two_prod(lwr, s0r), dd_two_prod(lwi, s0i)), dd_from(s1r));
        hdd ii_ = dd_add(dd_add(dd_two_prod(lwr, s0i), dd_two_prod(lwi, s0r)), dd_from(s1i));
        double inr = ir_.hi + ir_.lo, ini = ii_.hi + ii_.lo;
        /* S1's own terms change sign (h_n does); its sum is compensated, so
         * only a gross internal cancellation matters. */
        double lws0 = c_abs(lwr, lwi) * c_abs(s0r, s0i), as1 = c_abs(s1r, s1i);
        if (peak1 > 4.0 * (lws0 > as1 ? lws0 : as1)) return false;
        double innerpeak = lws0 + (as1 > hscale ? as1 : hscale);
        double fr = -G2, fi = 0.0;
        if (M & 1) fr = -fr;                            /* (-1)^M */
        if (pos && M > 0) {
            double pr, pi;
            c_pow_ma(c_abs(wr, wi), lwi, Md, &pr, &pi);
            c_mul(fr, fi, pr, pi, &fr, &fi);
        }
        c_mul(fr, fi, inr, ini, &logr, &logi);
        logmag = c_abs(fr, fi) * innerpeak;
    }
    double vr = finr + logr, vi = fini + logi;
    double mag = c_abs(vr, vi);
    if (!isfinite(vr) || !isfinite(vi) || mag == 0.0) return false;
    *kappa = (finmag + logmag) / mag;
    *rr = vr; *ri = vi;
    return true;
}

/* Choose the cheaper of the two Pfaff forms for real parameters:
 *   F(a,b;c;z) = (1-z)^-a F(a, c-b; c; z/(z-1)) = (1-z)^-b F(c-a, b; c; z/(z-1)).
 * Prefer one whose inner series terminates, then one whose inner parameters are
 * both positive (one-signed terms for 0 < w). Returns 0 for the a-form, 1 for b. */
static int hy_pfaff_pick(double a, double b, double c) {
    if (hy_is_nonpos_int(c - b)) return 0;
    if (hy_is_nonpos_int(c - a)) return 1;
    if (a > 0.0 && c - b > 0.0) return 0;
    if (b > 0.0 && c - a > 0.0) return 1;
    return 0;
}

/* Region selection for complex z, real parameters [Forrey 1997; MS 2008]:
 *   0: direct            w = z
 *   1: Pfaff             w = z/(z-1)
 *   2: 1-z connection    w = 1-z
 *   3: Pfaff + 1-w conn  w = 1/(1-z)
 * The direct and Pfaff series are cancellation-free compared with the
 * connection formulas, so they are preferred while their |w| is moderate. */
#define HY_2F1_SERIES_W 0.75     /* a plain series is preferred up to this |w| */
#define HY_2F1_SLOW_W   0.9      /* ... and, before continuation, up to this    */
#define HY_2F1_CONN_W   0.6      /* a connection formula up to this |w|        */
#define HY_2F1_MAX_W    0.97     /* nothing converges usefully beyond          */

/* Pfaff's form, (1-z)^-e F(e, c-f; c; z/(z-1)), summed directly. */
static bool hy2f1_pfaff_c(double a, double b, double c, double zr, double zi,
                          double m2, double li, double* rr, double* ri) {
    int pick = hy_pfaff_pick(a, b, c);
    double wr, wi;
    c_div(zr, zi, zr - 1.0, zi, &wr, &wi);
    double e = pick == 0 ? a : b;
    hdd cx = dd_two_sum(c, pick == 0 ? -b : -a);
    double ip[2], iplo[2];
    if (pick == 0) { ip[0] = a; iplo[0] = 0.0; ip[1] = cx.hi; iplo[1] = cx.lo; }
    else           { ip[0] = cx.hi; iplo[0] = cx.lo; ip[1] = b; iplo[1] = 0.0; }
    double fr, fi, pr, pi;
    if (!hy_ser1_c(2, ip, iplo, NULL, 1, &c, NULL, NULL, wr, wi, &fr, &fi)) return false;
    c_pow_ma(m2, li, -e, &pr, &pi);
    c_mul(fr, fi, pr, pi, rr, ri);
    return isfinite(*rr) && isfinite(*ri);
}

/* One Taylor step of the hypergeometric ODE
 *   z(1-z) w'' + [c - (a+b+1) z] w' - ab w = 0
 * from a centre z0 where (w, w') are known to z0 + h: with e_n = d_n h^n the
 * Taylor coefficients obey
 *   e_{n+2} = [ (n+a)(n+b) h^2 e_n - ((1-2z0) n + c - (a+b+1) z0)(n+1) h e_{n+1} ]
 *             / ( z0 (1-z0) (n+1)(n+2) ),
 * so w(z0+h) = sum e_n and h w'(z0+h) = sum n e_n. Both sums compensated; the
 * caller keeps |h| <= 0.6 min(|z0|, |1-z0|), so the terms fall at least as
 * 0.6^n. Returns false on a non-finite value or a sum that lost more than a bit. */
static bool hy2f1_tstep(double a, double b, double c, double z0r, double z0i,
                        double hr, double hi, double rho,
                        double* wr, double* wi, double* dr, double* di) {
    double pr, pi, Dr, Di;
    c_mul(z0r, z0i, 1.0 - z0r, -z0i, &pr, &pi);
    c_div(1.0, 0.0, pr, pi, &Dr, &Di);
    double alr = 1.0 - 2.0 * z0r, ali = -2.0 * z0i;
    double s1 = a + b + 1.0;
    double ber = c - s1 * z0r, bei = -s1 * z0i;
    double h2r, h2i;
    c_mul(hr, hi, hr, hi, &h2r, &h2i);
    double e0r = *wr, e0i = *wi, e1r, e1i;
    c_mul(*dr, *di, hr, hi, &e1r, &e1i);
    double sr = e0r, si = e0i, cr = 0.0, ci = 0.0;    /* sum e_n   */
    double tr = e1r, ti = e1i, ur = 0.0, ui = 0.0;    /* sum n e_n */
    hdd q = dd_two_sum(sr, e1r); sr = q.hi; cr += q.lo;
    q = dd_two_sum(si, e1i); si = q.hi; ci += q.lo;
    double peak = c_mag(e0r, e0i) > c_mag(e1r, e1i) ? c_mag(e0r, e0i) : c_mag(e1r, e1i);
    double nmin = fabs(a) + fabs(b) + 2.0;
    bool done = false;
    for (long n = 0; n < HY_MAXTERMS; n++) {
        double nd = (double)n;
        double f = (nd + a) * (nd + b);
        double t1r, t1i;
        c_mul(h2r, h2i, e0r * f, e0i * f, &t1r, &t1i);
        double qr = (alr * nd + ber) * (nd + 1.0), qi = (ali * nd + bei) * (nd + 1.0);
        double t2r, t2i, vr, vi;
        c_mul(qr, qi, hr, hi, &vr, &vi);
        c_mul(vr, vi, e1r, e1i, &t2r, &t2i);
        double den = (nd + 1.0) * (nd + 2.0);
        double e2r, e2i;
        c_mul(Dr / den, Di / den, t1r - t2r, t1i - t2i, &e2r, &e2i);
        q = dd_two_sum(sr, e2r); sr = q.hi; cr += q.lo;
        q = dd_two_sum(si, e2i); si = q.hi; ci += q.lo;
        double m2 = nd + 2.0;
        q = dd_two_sum(tr, m2 * e2r); tr = q.hi; ur += q.lo;
        q = dd_two_sum(ti, m2 * e2i); ti = q.hi; ui += q.lo;
        double m = c_mag(e2r, e2i);
        if (m > peak) peak = m;
        if (!isfinite(sr) || !isfinite(si) || !isfinite(tr) || !isfinite(ti)) return false;
        e0r = e1r; e0i = e1i; e1r = e2r; e1i = e2i;
        double em = c_mag(e0r, e0i) + m;
        if (nd > nmin && em * (nd + 3.0) <= HY_TOL * (1.0 - rho) *
                                             (c_mag(sr, si) + c_mag(tr, ti))) {
            done = true; break;
        }
    }
    if (!done) return false;
    sr += cr; si += ci; tr += ur; ti += ui;
    if (peak > HY_LOSS_OK * c_mag(sr, si)) return false;
    *wr = sr; *wi = si;
    c_div(tr, ti, hr, hi, dr, di);                     /* w' = (sum n e_n) / h */
    return true;
}

/* Analytic continuation of 2F1 along a path by successive Taylor steps -- the
 * one method that works wherever the transformations do not: the band around
 * z = e^(+-i pi/3) (where |z|, |1-z|, |z/(z-1)| and |1/(1-z)| are all near 1),
 * parameters whose connection formulas are ill-conditioned, and the stretch of
 * the cut z > 1 around z = 2. [Buhring 1987; MS 2008 sec. 4.4; the standard ODE
 * continuation.] Start at z0 = 1/2 +- 0.45 i (|z0| = |1-z0|), where w and
 * w' = (ab/c) F(a+1,b+1;c+1;z0) come from the direct series, on the side of the
 * real axis the target lies on -- the lower side for a cut point, whose value is
 * the limit from below. Each step goes straight toward the target by at most
 * 0.5 of the distance to the nearer singularity (0 or 1). */
static bool hy2f1_continue(double a, double b, double c, double zr, double zi,
                           double* rr, double* ri) {
    bool lower = (zi < 0.0) || (zi == 0.0 && signbit(zi));
    double z0r = 0.5, z0i = lower ? -0.45 : 0.45;
    double up[2] = { a, b };
    double wr, wi, dr, di;
    if (!hy_ser1_c(2, up, NULL, NULL, 1, &c, NULL, NULL, z0r, z0i, &wr, &wi)) return false;
    hdd a1 = dd_two_sum(a, 1.0), b1 = dd_two_sum(b, 1.0), c1 = dd_two_sum(c, 1.0);
    double up1[2] = { a1.hi, b1.hi }, up1lo[2] = { a1.lo, b1.lo };
    if (!hy_ser1_c(2, up1, up1lo, NULL, 1, &c1.hi, &c1.lo, NULL, z0r, z0i, &dr, &di))
        return false;
    double k = a * b / c;
    dr *= k; di *= k;
    for (int step = 0; step < 48; step++) {
        double R = c_abs(z0r, z0i), R1 = c_abs(1.0 - z0r, -z0i);
        if (R1 < R) R = R1;
        double hr = zr - z0r, hi = zi - z0i, ah = c_abs(hr, hi);
        if (ah <= 0.6 * R) {                           /* final step: evaluate */
            if (!hy2f1_tstep(a, b, c, z0r, z0i, hr, hi, ah / R, &wr, &wi, &dr, &di))
                return false;
            *rr = wr; *ri = wi;
            return true;
        }
        double sc = 0.5 * R / ah;
        hr *= sc; hi *= sc;
        if (!hy2f1_tstep(a, b, c, z0r, z0i, hr, hi, 0.5, &wr, &wi, &dr, &di))
            return false;
        z0r += hr; z0i += hi;
    }
    return false;
}

/* The connection formula on the smaller of |1-z| (m2) and |1/(1-z)| (m3). */
static bool hy2f1_conn_pick(double a, double b, double c, double zr, double zi,
                            double m2, double m3, double lr, double li,
                            double* vr, double* vi, double* kappa) {
    (void)zr; (void)zi;
    double omr = 1.0 - zr, omi = -zi;
    if (m2 <= m3) {
        if (!hy2f1_conn(a, b, c, omr, omi, lr, li, vr, vi, kappa)) return false;
        return isfinite(*vr) && isfinite(*vi);
    }
    /* (1-z)^-e F(e, c-f; c; w'), w' = z/(z-1), evaluated at 1-w' = 1/(1-z)
     * whose log is -log(1-z). */
    int pick = hy_pfaff_pick(a, b, c);
    double e = pick == 0 ? a : b, f = pick == 0 ? b : a;
    double w3r, w3i, ir, ii, pr, pi;
    c_div(1.0, 0.0, omr, omi, &w3r, &w3i);
    if (!hy2f1_conn(e, c - f, c, w3r, w3i, -lr, -li, &ir, &ii, kappa)) return false;
    c_pow_ma(m2, li, -e, &pr, &pi);
    c_mul(ir, ii, pr, pi, vr, vi);
    return isfinite(*vr) && isfinite(*vi);
}

static bool hy2f1_core(double a, double b, double c, double zr, double zi,
                       double* rr, double* ri) {
    /* 1 - z and its principal log, shared by every transformation (one branch
     * convention: on the cut z > 1 the caller's zi = -0 gives arg(1-z) = +pi,
     * i.e. the limit from below). */
    double omr = 1.0 - zr, omi = -zi;
    double lr, li;
    double m0 = c_abs(zr, zi);
    double m2 = c_abs(omr, omi);
    double m1 = m0 / m2, m3 = 1.0 / m2;
    double up[2] = { a, b };
    c_log(omr, omi, &lr, &li);

    /* Plain series first while cheap and cancellation-free: the direct one when
     * it is the smaller |w|, else Pfaff's. */
    if (m0 <= HY_2F1_SERIES_W && m0 <= m1)
        return hy_ser1_c(2, up, NULL, NULL, 1, &c, NULL, NULL, zr, zi, rr, ri);
    if (m1 <= HY_2F1_SERIES_W)
        return hy2f1_pfaff_c(a, b, c, zr, zi, m2, li, rr, ri);

    /* Connection formulas, accepted outright while the combination is well
     * conditioned (kappa <= HY_KAPPA_MAX) -- first only where they converge
     * fast; a worse-conditioned value is kept as the last resort below. */
    double kappa = INFINITY, vr = 0.0, vi = 0.0;
    bool have = false;
    double mc = m2 < m3 ? m2 : m3;
    if (mc <= HY_2F1_CONN_W) {
        have = hy2f1_conn_pick(a, b, c, zr, zi, m2, m3, lr, li, &vr, &vi, &kappa);
        if (have && kappa <= HY_KAPPA_MAX) { *rr = vr; *ri = vi; return true; }
    }
    /* A plain series that is still reasonably short is cancellation-free, so it
     * beats both a long continuation and an ill-conditioned connection. */
    if (m0 <= HY_2F1_SLOW_W && m0 <= m1)
        return hy_ser1_c(2, up, NULL, NULL, 1, &c, NULL, NULL, zr, zi, rr, ri);
    if (m1 <= HY_2F1_SLOW_W)
        return hy2f1_pfaff_c(a, b, c, zr, zi, m2, li, rr, ri);
    /* The e^(+-i pi/3) band, an ill-conditioned connection, the cut near 2:
     * continue the ODE solution from z0 = 1/2 +- 0.45 i. */
    if (hy2f1_continue(a, b, c, zr, zi, rr, ri)) return true;
    /* Then a slower connection series, or a plain one (cancellation-free). */
    if (!have && mc <= HY_2F1_MAX_W) {
        have = hy2f1_conn_pick(a, b, c, zr, zi, m2, m3, lr, li, &vr, &vi, &kappa);
        if (have && kappa <= HY_KAPPA_MAX) { *rr = vr; *ri = vi; return true; }
    }
    if (m0 <= HY_2F1_MAX_W && m0 <= m1)
        return hy_ser1_c(2, up, NULL, NULL, 1, &c, NULL, NULL, zr, zi, rr, ri);
    if (m1 <= HY_2F1_MAX_W)
        return hy2f1_pfaff_c(a, b, c, zr, zi, m2, li, rr, ri);
    /* Finally a connection value whose combination loses up to ~3 bits. */
    if (have && kappa <= HY_KAPPA_LAST) { *rr = vr; *ri = vi; return true; }
    return false;
}

/* Gauss's sum at z = 1: G(c)G(c-a-b)/(G(c-a)G(c-b)), c-a-b > 0. */
static bool hy2f1_at_one(double a, double b, double c, double* out) {
    double s = c - a - b;
    if (!(s > 0.0)) return false;
    double v = hy_gamma(c) * hy_gamma(s) * hy_rgamma(c - a) * hy_rgamma(c - b);
    if (!isfinite(v)) return false;
    *out = v;
    return true;
}

/* ---- Real 2F1: a "plan" so a buffer pays the parameter work once ----------
 *
 * The two series the real kernel sums on (-9, 0.9] -- the direct one and
 * Pfaff's -- are written as sum_k c_k w^k with the coefficient recurrence
 *     c_0 = 1,  c_{k+1} = c_k * g_k,  g_k = (p0+k)(p1+k) / ((q0+k)(k+1)),
 * which depends on the parameters only. A buffer of 10^6 points with fixed
 * (a, b, c) tabulates c_k and g_k once (Hy2Tab.c/.g) and then pays two
 * multiplies and a compensated add per term, with no division; a scalar call
 * computes the same c_k and g_k on the fly with the same expression, so the two
 * are bit-identical (and the Compile[] VM, which calls the scalar kernel, too). */
#define HY2_TAB_MAX 4096

typedef struct {
    double p0, p0lo, p1, p1lo, q0;   /* series parameters (hi + lo)            */
    double* c;                       /* c_k, k < n (NULL: on the fly)          */
    double* g;                       /* g_k                                    */
    int n;
    bool pos;                        /* every g_k > 0 (all k): one-signed terms */
    bool lean;                       /* pos, and every g_k <= 1: terms never grow */
} Hy2Tab;

static inline double hy2_g(const Hy2Tab* t, double kd) {
    double num = ((t->p0 + kd) + t->p0lo) * ((t->p1 + kd) + t->p1lo);
    double den = (t->q0 + kd) * (kd + 1.0);
    return num / den;
}

static void hy2_tab_fill(Hy2Tab* t) {
    t->c = malloc(sizeof(double) * 2 * HY2_TAB_MAX);
    if (!t->c) { t->g = NULL; t->n = 0; return; }
    t->g = t->c + HY2_TAB_MAX;
    double ck = 1.0;
    for (int k = 0; k < HY2_TAB_MAX; k++) {
        double g = hy2_g(t, (double)k);
        t->c[k] = ck; t->g[k] = g;
        ck = ck * g;
    }
    t->n = HY2_TAB_MAX;
}

/* The lean case (t->lean: every g_k in (0, 1]) and 0 < w < 1: all terms
 * positive and non-increasing, so there is no cancellation or deep peak to
 * measure, every partial sum dominates every later term (Fast2Sum is exact:
 * 3 flops instead of TwoSum's 6), the ratio bound is w itself, and the stop
 * test can run once per four terms. Even and odd terms go to two compensated
 * sums seeded with t_0 = 1 and t_1 (two independent dependency chains).
 *
 * hy2_lean_x4 below runs this very loop on four buffer elements in lock-step;
 * the arithmetic and the stop rule are written once (HY_LEAN_* macros) so the
 * scalar, buffer and Compile[] surfaces cannot drift apart. */
#define HY_LEAN_INIT(S0, C0, S1, C1, PE, PO, W2, THR, W, C1COEF)              \
    double S0 = 1.0, C0 = 0.0, S1 = (C1COEF) * (W), C1 = 0.0;                 \
    double W2 = (W) * (W), PE = W2, PO = (W) * W2;                            \
    double THR = HY_TOL * (1.0 - (W))
/* terms k..k+3 (k even, k >= 2) with coefficients ca..cd */
#define HY_LEAN_STEP(S0, C0, S1, C1, PE, PO, W2, ca, cb, cc, cd, T4, T5) do {  \
    double t2_ = (ca) * PE, t3_ = (cb) * PO; PE *= W2; PO *= W2;               \
    T4 = (cc) * PE; T5 = (cd) * PO; PE *= W2; PO *= W2;                        \
    double h_;                                                                 \
    h_ = S0 + t2_; C0 += t2_ - (h_ - S0); S0 = h_;                             \
    h_ = S1 + t3_; C1 += t3_ - (h_ - S1); S1 = h_;                             \
    h_ = S0 + T4;  C0 += T4  - (h_ - S0); S0 = h_;                             \
    h_ = S1 + T5;  C1 += T5  - (h_ - S1); S1 = h_;                             \
} while (0)
#define HY_LEAN_DONE(S0, S1, THR, T4, T5) ((T4) + (T5) <= ((S0) + (S1)) * (THR))

static bool hy2_tab_sum_pos(const Hy2Tab* t, long kmin, double w, double* out) {
    /* coefficients c_1.. on the fly (or from the table: the same numbers) */
    double c1 = (t->n > 1) ? t->c[1] : hy2_g(t, 0.0);
    HY_LEAN_INIT(s0, e0, s1, e1, pe, po, w2, thr, w, c1);
    double ck = c1;                                 /* c_{k-1} for k = 2 */
    for (long k = 2; k + 4 < HY_MAXTERMS; k += 4) {
        double ca, cb, cc, cd;
        if (k + 3 < t->n) {
            ca = t->c[k]; cb = t->c[k + 1]; cc = t->c[k + 2]; cd = t->c[k + 3];
        } else {
            ca = ck * hy2_g(t, (double)(k - 1));
            cb = ca * hy2_g(t, (double)k);
            cc = cb * hy2_g(t, (double)(k + 1));
            cd = cc * hy2_g(t, (double)(k + 2));
        }
        ck = cd;
        double t4, t5;
        HY_LEAN_STEP(s0, e0, s1, e1, pe, po, w2, ca, cb, cc, cd, t4, t5);
        if (k + 4 >= kmin && HY_LEAN_DONE(s0, s1, thr, t4, t5)) break;
    }
    *out = (s0 + s1) + (e0 + e1);
    return isfinite(*out);
}

/* sum_k c_k w^k for 0 <= |w| < 1, compensated; *peak = max |term|.
 * Two interleaved chains (even and odd k, each stepping by w^2) so neither the
 * power nor the running sum is a one-term-per-latency dependency; the pair is
 * tested together with the generic series' stop rule (ratio |g_k w|). The
 * coefficients come from the table while it lasts and from the same
 * recurrence beyond it, so table and on-the-fly sums are the same numbers. */
static bool hy2_tab_sum(const Hy2Tab* t, long kmin, double w, double* out, double* peak,
                        long* kpk_out) {
    double s0 = 1.0, c0 = 0.0, s1 = 0.0, c1 = 0.0, pk = 1.0;
    long kpk = 0;
    double aw = fabs(w), w2 = w * w;
    double p0 = w2, p1 = w;                 /* w^(k+2), w^(k+1) for k = 0 */
    double ck = 1.0;                        /* c_k, k = 0                 */
    for (long k = 0; k + 2 < HY_MAXTERMS; k += 2) {
        double g0 = (k < t->n) ? t->g[k] : hy2_g(t, (double)k);
        double ck1 = (k + 1 < t->n) ? t->c[k + 1] : ck * g0;
        double g1 = (k + 1 < t->n) ? t->g[k + 1] : hy2_g(t, (double)(k + 1));
        double ck2 = (k + 2 < t->n) ? t->c[k + 2] : ck1 * g1;
        double t1 = ck1 * p1, t2 = ck2 * p0;
        hdd q = dd_two_sum(s1, t1); s1 = q.hi; c1 += q.lo;
        q = dd_two_sum(s0, t2); s0 = q.hi; c0 += q.lo;
        double a1 = fabs(t1), a2 = fabs(t2);
        if (a1 > pk) { pk = a1; kpk = k + 1; }
        if (a2 > pk) { pk = a2; kpk = k + 2; }
        if (!isfinite(s0) || !isfinite(s1)) return false;
        if (ck2 == 0.0) break;                     /* terminated */
        if (k + 2 >= kmin) {
            double rho = fabs(g1) * aw;
            if (rho < aw) rho = aw;
            if (rho < 1.0 && a1 + a2 <= fabs(s0 + s1) * HY_TOL * (1.0 - rho)) break;
        }
        p0 *= w2; p1 *= w2;
        ck = ck2;
    }
    *out = (s0 + s1) + (c0 + c1);
    *peak = pk;
    *kpk_out = kpk;
    return isfinite(*out);
}

typedef struct {
    double a, b, c;
    int kind;            /* 0 general, 1 a or b terminates (polynomial), 2 a==c, 3 b==c */
    Hy2Tab dir;          /* (a, b; c)            */
    Hy2Tab pf;           /* Pfaff (e, c-f; c)    */
    double e;            /* Pfaff exponent        */
    long kmin_dir, kmin_pf;
    bool dir_ok, pf_ok;  /* series classification passed */
} Hy2Plan;

static long hy2_kmin(double p0, double p1, double q0) {
    double neg = 0.0;
    if (p0 < 0.0) neg = -p0;
    if (p1 < 0.0 && -p1 > neg) neg = -p1;
    if (q0 < 0.0 && -q0 > neg) neg = -q0;
    return (long)ceil(neg) + 1;
}

/* Returns false when the parameters alone decide a decline (a (c)_k pole). */
static bool hy2_plan_init(Hy2Plan* P, double a, double b, double c, bool tabulate) {
    P->a = a; P->b = b; P->c = c;
    bool ta = hy_is_nonpos_int(a), tb = hy_is_nonpos_int(b);
    if (hy_is_nonpos_int(c) && !ta && !tb) return false;
    P->kind = (ta || tb) ? 1 : (a == c ? 2 : (b == c ? 3 : 0));
    P->dir.p0 = a; P->dir.p0lo = 0.0; P->dir.p1 = b; P->dir.p1lo = 0.0; P->dir.q0 = c;
    P->dir.c = P->dir.g = NULL; P->dir.n = 0;
    int pick = hy_pfaff_pick(a, b, c);
    P->e = pick == 0 ? a : b;
    hdd cx = dd_two_sum(c, pick == 0 ? -b : -a);
    if (pick == 0) { P->pf.p0 = a; P->pf.p0lo = 0.0; P->pf.p1 = cx.hi; P->pf.p1lo = cx.lo; }
    else           { P->pf.p0 = cx.hi; P->pf.p0lo = cx.lo; P->pf.p1 = b; P->pf.p1lo = 0.0; }
    P->pf.q0 = c; P->pf.c = P->pf.g = NULL; P->pf.n = 0;
    /* g_k > 0 for every k exactly when all three parameters are positive. */
    P->dir.pos = (a > 0.0 && b > 0.0 && c > 0.0);
    P->pf.pos = (P->pf.p0 + P->pf.p0lo > 0.0 && P->pf.p1 + P->pf.p1lo > 0.0 && c > 0.0);
    /* g_k = (p0+k)(p1+k)/((q0+k)(k+1)) <= 1 for all k >= 0 exactly when
     * p0 p1 + (p0+p1) k <= q0 + (q0+1) k for all k, i.e. (sufficient, and
     * necessary as k -> 0 and k -> infinity) p0 p1 <= q0 and p0 + p1 <= q0 + 1. */
    P->dir.lean = P->dir.pos && a * b <= c && a + b <= c + 1.0;
    {
        double q0 = P->pf.p0 + P->pf.p0lo, q1 = P->pf.p1 + P->pf.p1lo;
        P->pf.lean = P->pf.pos && q0 * q1 <= c && q0 + q1 <= c + 1.0;
    }
    P->kmin_dir = hy2_kmin(a, b, c);
    P->kmin_pf = hy2_kmin(P->pf.p0 + P->pf.p0lo, P->pf.p1 + P->pf.p1lo, c);
    /* A non-terminating series is the generic one; a terminating Pfaff inner
     * series (c-b or c-a a non-positive integer) is fine for the table too, as
     * g_k = 0 ends it. Neither may cross a (c)_k pole: c is not one here. */
    P->dir_ok = P->kmin_dir <= HY_MAXTERMS;
    P->pf_ok = P->kmin_pf <= HY_MAXTERMS;
    if (tabulate && P->kind == 0) {
        if (P->dir_ok) hy2_tab_fill(&P->dir);
        if (P->pf_ok) hy2_tab_fill(&P->pf);
    }
    return true;
}

static void hy2_plan_free(Hy2Plan* P) {
    free(P->dir.c); free(P->pf.c);
    P->dir.c = P->pf.c = NULL;
}

/* One of the two tabulated series at w, with the double-double re-sum when
 * the double one loses more than a bit (one-signed terms make that rare). */
static bool hy2_series(const Hy2Tab* t, long kmin, double w, double* out) {
    /* The lean loop needs one-signed terms that never grow (peak at k = 0), so
     * there is neither cancellation nor deep-peak rounding to measure: all
     * g_k > 0 and g_k <= 1 for every k, decided once per plan (t->lean). */
    if (t->lean && w > 0.0 && w < 1.0)
        return hy2_tab_sum_pos(t, kmin, w, out);
    double v, peak;
    long kpk;
    if (!hy2_tab_sum(t, kmin, w, &v, &peak, &kpk)) return false;
    if (peak <= HY_LOSS_OK * fabs(v) && kpk <= HY_KPEAK_DD) { *out = v; return true; }
    if (!(peak <= HY_LOSS_MAX * fabs(v))) return false;
    double up[2] = { t->p0, t->p1 }, uplo[2] = { t->p0lo, t->p1lo };
    HySer s = { up, uplo, NULL, 2, &t->q0, NULL, NULL, 1 };
    HySerInfo inf;
    if (!hy_ser_prep(&s, fabs(w), &inf)) return false;
    return hy_sum_r_dd(&s, &inf, w, out);
}

/* (1-x)^-e. e = 1 and 2 (the common Pfaff exponents) by division, which is
 * as accurate as pow (<= 1 ulp) and several times cheaper. */
static inline double hy_pow_1mx(double x, double e) {
    double u = 1.0 - x;
    if (e == 1.0) return 1.0 / u;
    if (e == 2.0) return 1.0 / (u * u);
    return pow(u, -e);
}

static bool hy2_plan_eval(const Hy2Plan* P, double x, double* out) {
    double a = P->a, b = P->b, c = P->c;
    if (!isfinite(x)) return false;
    if (x == 0.0 || a == 0.0 || b == 0.0) { *out = 1.0; return true; }
    if (P->kind == 1) {                            /* polynomial, any x */
        double up[2] = { a, b };
        return hy_ser1_r(2, up, NULL, 1, &c, NULL, x, out);
    }
    if (x > 1.0) return false;                     /* complex: the cut */
    if (x == 1.0) return hy2f1_at_one(a, b, c, out);
    if (P->kind == 2) { *out = pow(1.0 - x, -b); return isfinite(*out); }
    if (P->kind == 3) { *out = pow(1.0 - x, -a); return isfinite(*out); }
    if (x > 0.0 && x <= 0.9 && P->dir_ok)          /* one-signed terms for a,b,c > 0 */
        return hy2_series(&P->dir, P->kmin_dir, x, out);
    if (x < 0.0 && x >= -9.0 && P->pf_ok) {
        /* Pfaff [AS 15.3.4]: w = x/(x-1) in (0, 0.9] */
        double w = x / (x - 1.0), v;
        if (!hy2_series(&P->pf, P->kmin_pf, w, &v)) return false;
        double r = hy_pow_1mx(x, P->e) * v;
        if (!isfinite(r)) return false;
        *out = r;
        return true;
    }
    /* x in (0.9, 1) or x < -9: the connection formulas, in the shared
     * (complex-arithmetic) core; the result is real. */
    double rr, ri;
    if (!hy2f1_core(a, b, c, x, 0.0, &rr, &ri)) return false;
    *out = rr;
    return isfinite(rr);
}

bool sf_machine_hyper2f1(double a, double b, double c, double x, double* out) {
    if (!isfinite(a) || !isfinite(b) || !isfinite(c)) return false;
    Hy2Plan P;
    if (!hy2_plan_init(&P, a, b, c, false)) return false;
    return hy2_plan_eval(&P, x, out);
}

/* The buffer driver: out[i] = pFq(a; b; x[i]) for a constant parameter set,
 * returning the number of leading elements computed (n on full success; a
 * smaller value is the index of the first element the kernel declined). Every
 * element is computed exactly as the scalar entry point computes it -- the 2F1
 * plan's tables hold the very coefficients the scalar path forms on the fly,
 * and the 1F1(1;2;x) closed form is the same expression -- so the buffer and
 * scalar surfaces agree bit-for-bit. */
/* Four buffer elements' lean sums (hy2_tab_sum_pos, table part) in lock-step.
 * Each lane performs exactly the scalar loop's operations in the scalar order
 * and stops where the scalar loop would, so its result is bit-identical; the
 * point is only that four independent dependency chains keep the FP units busy
 * where one chain is latency-bound. A lane that runs off the end of the table
 * (never, for the |w| <= 0.9 the lean path sees) is left to the scalar loop
 * (*ok = false for it). */
static void hy2_lean_x4(const Hy2Tab* t, long kmin, const double* w, int nl,
                        double* out, bool* ok) {
    /* Four lanes as plain scalars so every chain stays in a register; a
     * finished lane keeps computing harmlessly (its result is already taken). */
    double wA = w[0], wB = nl > 1 ? w[1] : wA, wC = nl > 2 ? w[2] : wA,
           wD = nl > 3 ? w[3] : wA;
    for (int l = 0; l < 4; l++) { ok[l] = false; out[l] = 0.0; }
    const double* CC = t->c;
    double c1 = CC[1];
    HY_LEAN_INIT(s0A, e0A, s1A, e1A, peA, poA, w2A, thrA, wA, c1);
    HY_LEAN_INIT(s0B, e0B, s1B, e1B, peB, poB, w2B, thrB, wB, c1);
    HY_LEAN_INIT(s0C, e0C, s1C, e1C, peC, poC, w2C, thrC, wC, c1);
    HY_LEAN_INIT(s0D, e0D, s1D, e1D, peD, poD, w2D, thrD, wD, c1);
    bool aA = true, aB = nl > 1, aC = nl > 2, aD = nl > 3;
    for (long k = 2; k + 3 < (long)t->n && (aA || aB || aC || aD); k += 4) {
        double ca = CC[k], cb = CC[k + 1], cc = CC[k + 2], cd = CC[k + 3];
        bool test = (k + 4 >= kmin);
        double t4, t5;
#define HY_LANE(L, IDX) do {                                                   \
        HY_LEAN_STEP(s0##L, e0##L, s1##L, e1##L, pe##L, po##L, w2##L,          \
                     ca, cb, cc, cd, t4, t5);                                  \
        if (a##L && test && HY_LEAN_DONE(s0##L, s1##L, thr##L, t4, t5)) {      \
            out[IDX] = (s0##L + s1##L) + (e0##L + e1##L);                      \
            ok[IDX] = isfinite(out[IDX]);                                      \
            a##L = false;                                                      \
        }                                                                      \
    } while (0)
        HY_LANE(A, 0); HY_LANE(B, 1); HY_LANE(C, 2); HY_LANE(D, 3);
#undef HY_LANE
    }
    /* Lanes still running at the end of the table: the scalar loop decides. */
    if (aA) ok[0] = false;
    if (nl > 1 && aB) ok[1] = false;
    if (nl > 2 && aC) ok[2] = false;
    if (nl > 3 && aD) ok[3] = false;
}

/* The lean-path test of hy2_series, for the batch driver. */
static bool hy2_is_lean(const Hy2Tab* t, double w) {
    return t->c && t->lean && w > 0.0 && w < 1.0;
}

size_t sf_machine_pfq_vec(const double* a, size_t p, const double* b, size_t q,
                          const double* x, double* out, size_t n) {
    if (p == 2 && q == 1) {
        if (!isfinite(a[0]) || !isfinite(a[1]) || !isfinite(b[0])) return 0;
        Hy2Plan P;
        if (!hy2_plan_init(&P, a[0], a[1], b[0], n >= 64)) return 0;
        size_t fail = n;
        if (P.kind == 0 && P.a != 0.0 && P.b != 0.0) {
            /* Elements on the lean direct (0 < x <= 0.9) or lean Pfaff
             * (-9 <= x < 0) path are queued per table AND per |w| bucket (so
             * the four lanes of a batch need similar numbers of terms) and
             * summed four at a time; everything else goes through
             * hy2_plan_eval as is. */
            enum { NB = 32 };
            size_t qi[2][NB][4]; double qw[2][NB][4]; int qn[2][NB];
            memset(qn, 0, sizeof qn);
            const Hy2Tab* tab[2] = { &P.dir, &P.pf };
            long km[2] = { P.kmin_dir, P.kmin_pf };
            for (size_t i = 0; i <= n; i++) {
                int fw = -1, fb = -1;          /* the (table, bucket) to flush */
                if (i < n) {
                    double xi = x[i];
                    int which = -1; double w = 0.0;
                    if (isfinite(xi) && xi > 0.0 && xi <= 0.9 && P.dir_ok &&
                        hy2_is_lean(&P.dir, xi)) { which = 0; w = xi; }
                    else if (isfinite(xi) && xi < 0.0 && xi >= -9.0 && P.pf_ok) {
                        w = xi / (xi - 1.0);
                        if (hy2_is_lean(&P.pf, w)) which = 1;
                    }
                    if (which < 0) {
                        if (!hy2_plan_eval(&P, xi, &out[i]) && i < fail) fail = i;
                        continue;
                    }
                    int bk = (int)(w * NB);
                    if (bk >= NB) bk = NB - 1;
                    int m = qn[which][bk];
                    qi[which][bk][m] = i; qw[which][bk][m] = w; qn[which][bk] = m + 1;
                    if (m + 1 < 4) continue;
                    fw = which; fb = bk;
                }
                for (int wq = 0; wq < 2; wq++) {
                    for (int bk = 0; bk < NB; bk++) {
                        int m = qn[wq][bk];
                        if (m == 0) continue;
                        if (i < n && !(wq == fw && bk == fb)) continue;
                        double v[4]; bool okl[4];
                        hy2_lean_x4(tab[wq], km[wq], qw[wq][bk], m, v, okl);
                        for (int l = 0; l < m; l++) {
                            size_t j = qi[wq][bk][l];
                            bool good;
                            if (!okl[l]) {                  /* scalar loop decides */
                                good = hy2_plan_eval(&P, x[j], &out[j]);
                            } else if (wq == 0) {
                                out[j] = v[l]; good = true;
                            } else {                        /* as hy2_plan_eval */
                                double r = hy_pow_1mx(x[j], P.e) * v[l];
                                good = isfinite(r);
                                if (good) out[j] = r;
                            }
                            if (!good && j < fail) fail = j;
                        }
                        qn[wq][bk] = 0;
                    }
                }
            }
        } else {
            for (size_t i = 0; i < n; i++)
                if (!hy2_plan_eval(&P, x[i], &out[i])) { fail = i; break; }
        }
        hy2_plan_free(&P);
        return fail;
    }
    if (p == 1 && q == 1 && a[0] == 1.0 && b[0] == 2.0) {
        size_t i = 0;
        for (; i < n; i++) {
            double xi = x[i];
            if (xi == 0.0) { out[i] = 1.0; continue; }
            double v = expm1(xi) / xi;             /* = sf_machine_hyper1f1's form */
            if (!isfinite(v)) break;
            out[i] = v;
        }
        return i;
    }
    size_t i = 0;
    for (; i < n; i++) if (!sf_machine_pfq(a, p, b, q, x[i], &out[i])) break;
    return i;
}

bool sf_machine_hyper2f1_c(double ar, double ai, double br, double bi,
                           double cr, double ci, double zr, double zi,
                           double* rr, double* ri) {
    if (!isfinite(ar) || !isfinite(ai) || !isfinite(br) || !isfinite(bi) ||
        !isfinite(cr) || !isfinite(ci) || !isfinite(zr) || !isfinite(zi)) return false;
    bool real_p = (ai == 0.0 && bi == 0.0 && ci == 0.0);
    if ((zr == 0.0 && zi == 0.0) || (ar == 0.0 && ai == 0.0) || (br == 0.0 && bi == 0.0)) {
        *rr = 1.0; *ri = 0.0; return true;
    }
    double up[2] = { ar, br }, upi[2] = { ai, bi };
    bool ta = (ai == 0.0 && hy_is_nonpos_int(ar)), tb = (bi == 0.0 && hy_is_nonpos_int(br));
    if (ci == 0.0 && hy_is_nonpos_int(cr) && !ta && !tb) return false;
    if (ta || tb)
        return hy_ser1_c(2, up, NULL, upi, 1, &cr, NULL, &ci, zr, zi, rr, ri);
    if (!real_p) {
        if (!(c_abs(zr, zi) <= 0.9)) return false;
        return hy_ser1_c(2, up, NULL, upi, 1, &cr, NULL, &ci, zr, zi, rr, ri);
    }
    if (zr == 1.0 && zi == 0.0) {
        double v;
        if (!hy2f1_at_one(ar, br, cr, &v)) return false;
        *rr = v; *ri = 0.0; return true;
    }
    if (ar == cr || br == cr) {                       /* (1-z)^-e */
        double e = (ar == cr) ? br : ar, lr, li;
        c_log(1.0 - zr, -zi, &lr, &li);
        c_pow_ma(c_abs(1.0 - zr, -zi), li, -e, rr, ri);
        return isfinite(*rr) && isfinite(*ri);
    }
    return hy2f1_core(ar, br, cr, zr, zi, rr, ri);
}

/* ------------------------------------------------------------------ */
/*  pFq                                                                */
/* ------------------------------------------------------------------ */

bool sf_machine_pfq(const double* a, size_t p, const double* b, size_t q,
                    double z, double* out) {
    if (!isfinite(z)) return false;
    for (size_t i = 0; i < p; i++) if (!isfinite(a[i])) return false;
    for (size_t j = 0; j < q; j++) if (!isfinite(b[j])) return false;
    /* The specialised kernels own their shapes, so every spelling of the same
     * function -- Hypergeometric1F1[a,b,z] or HypergeometricPFQ[{a},{b},z] --
     * reaches the same code and the same bits. */
    if (p == 0 && q == 1) return sf_machine_hyper0f1(b[0], z, out);
    if (p == 1 && q == 1) return sf_machine_hyper1f1(a[0], b[0], z, out);
    if (p == 2 && q == 1) return sf_machine_hyper2f1(a[0], a[1], b[0], z, out);
    if (z == 0.0) { *out = 1.0; return true; }
    if (p == 0 && q == 0) { *out = exp(z); return isfinite(*out); }
    if (p == 1 && q == 0) {
        if (!(z < 1.0)) return false;
        *out = pow(1.0 - z, -a[0]);
        return isfinite(*out);
    }
    HySer s = { a, NULL, NULL, p, b, NULL, NULL, q };
    return hy_series_r(&s, z, out);
}

bool sf_machine_pfq_c(const double* ar, const double* ai, size_t p,
                      const double* br, const double* bi, size_t q,
                      double zr, double zi, double* rr, double* ri) {
    if (!isfinite(zr) || !isfinite(zi)) return false;
    for (size_t i = 0; i < p; i++) if (!isfinite(ar[i]) || !isfinite(ai[i])) return false;
    for (size_t j = 0; j < q; j++) if (!isfinite(br[j]) || !isfinite(bi[j])) return false;
    if (p == 0 && q == 1) return sf_machine_hyper0f1_c(br[0], bi[0], zr, zi, rr, ri);
    if (p == 1 && q == 1)
        return sf_machine_hyper1f1_c(ar[0], ai[0], br[0], bi[0], zr, zi, rr, ri);
    if (p == 2 && q == 1)
        return sf_machine_hyper2f1_c(ar[0], ai[0], ar[1], ai[1], br[0], bi[0],
                                     zr, zi, rr, ri);
    if (zr == 0.0 && zi == 0.0) { *rr = 1.0; *ri = 0.0; return true; }
    HySer s = { ar, NULL, ai, p, br, NULL, bi, q };
    return hy_series_c(&s, zr, zi, rr, ri);
}

/* ------------------------------------------------------------------ */
/*  HypergeometricU and the Whittaker functions, complex z             */
/* ------------------------------------------------------------------ */

/* U(a,b,z) ~ z^-a sum_k (a)_k (a-b+1)_k / k! (-1/z)^k   [DLMF 13.7.3],
 * |arg z| < 3 pi/2; exact when it terminates. Accepted only when the least term
 * is below 2^-54 of the sum. */
static bool hyu_asym_c(double a, double b, double zr, double zi, double li,
                       double* rr, double* ri) {
    double a2 = a - b + 1.0;
    double tr, ti;                                   /* -1/z */
    c_div(-1.0, 0.0, zr, zi, &tr, &ti);
    double cr = 1.0, ci = 0.0, sr = 1.0, si = 0.0, last = 1.0;
    bool ok = false;
    /* a or a-b+1 = -m: the series is a FINITE sum (exactly U), so it must not be
     * cut at a growing term as if divergent -- sum all m+1 terms and decline only
     * on measured cancellation (peak term over result, sf_series_usable). */
    bool term = hy_is_nonpos_int(a) || hy_is_nonpos_int(a2);
    double peak = 1.0;
    int kmax = 400;
    if (term) {
        double m = hy_is_nonpos_int(a) ? -a : -a2;
        if (hy_is_nonpos_int(a) && hy_is_nonpos_int(a2) && -a2 < m) m = -a2;
        if (m > 100000.0) return false;
        kmax = (int)m + 1;
    }
    for (int k = 0; k < kmax; k++) {
        double f = (a + k) * (a2 + k) / (k + 1.0);
        if (f == 0.0) { ok = true; break; }
        double nr, ni;
        c_mul(cr, ci, tr * f, ti * f, &nr, &ni);
        double an = c_abs(nr, ni);
        if (!isfinite(an)) return false;
        if (term) {                                  /* finite sum: no truncation */
            if (an > peak) peak = an;
            sr += nr; si += ni; cr = nr; ci = ni;
            continue;
        }
        if (an > last) { ok = (last <= c_abs(sr, si) * HY_TOL); break; }
        sr += nr; si += ni;
        cr = nr; ci = ni; last = an;
        if (an <= c_abs(sr, si) * HY_TOL) { ok = true; break; }
    }
    if (!ok) return false;
    if (term && !sf_series_usable(peak, c_abs(sr, si))) return false;
    double pr, pi;
    c_pow_ma(c_abs(zr, zi), li, -a, &pr, &pi);
    c_mul(sr, si, pr, pi, rr, ri);
    return isfinite(*rr) && isfinite(*ri);
}

/* Non-integer b: [AS 13.1.3]
 *   U = G(1-b)/G(a-b+1) M(a,b,z) + G(b-1)/G(a) z^(1-b) M(a-b+1, 2-b, z). */
static bool hyu_conn_c(double a, double b, double zr, double zi, double li,
                       double* rr, double* ri) {
    double g1 = hy_gamma(1.0 - b) * hy_rgamma(a - b + 1.0);
    double g2 = hy_gamma(b - 1.0) * hy_rgamma(a);
    if (!isfinite(g1) || !isfinite(g2)) return false;
    double m1r = 0.0, m1i = 0.0, m2r = 0.0, m2i = 0.0;
    if (g1 != 0.0 && !sf_machine_hyper1f1_c(a, 0.0, b, 0.0, zr, zi, &m1r, &m1i)) return false;
    if (g2 != 0.0) {
        if (!sf_machine_hyper1f1_c(a - b + 1.0, 0.0, 2.0 - b, 0.0, zr, zi, &m2r, &m2i))
            return false;
        double pr, pi;
        c_pow_ma(c_abs(zr, zi), li, 1.0 - b, &pr, &pi);
        c_mul(m2r, m2i, pr, pi, &m2r, &m2i);
    }
    double t1r = g1 * m1r, t1i = g1 * m1i, t2r = g2 * m2r, t2i = g2 * m2i;
    double vr = t1r + t2r, vi = t1i + t2i, mag = c_abs(vr, vi);
    if (!isfinite(vr) || !isfinite(vi) || mag == 0.0) return false;
    if (c_abs(t1r, t1i) + c_abs(t2r, t2i) > HY_KAPPA_MAX * mag) return false;
    *rr = vr; *ri = vi;
    return true;
}

/* Integer b = n+1 >= 1: the confluent limit [DLMF 13.2.9]
 *   U(a,n+1,z) = (-1)^(n+1)/(n! G(a-n)) sum_k (a)_k/((n+1)_k k!) z^k
 *                   [ln z + psi(a+k) - psi(k+1) - psi(n+k+1)]
 *              + (n-1)!/G(a) z^-n sum_{k<n} (a-n)_k/((1-n)_k k!) z^k. */
static bool hyu_intb_c(double a, int n, double zr, double zi, double lr, double li,
                       double* rr, double* ri) {
    /* a = 0, -1, ...: both 1/Gamma prefactors vanish and this formula degenerates
     * to 0 (DLMF 13.2.9 excludes it); the polynomial is hyu_asym_c's finite sum,
     * so reaching here means that declined on cancellation -- decline too. */
    if (hy_is_nonpos_int(a)) return false;
    double finr = 0.0, fini = 0.0, peak = 0.0;
    if (n >= 1) {
        double rga = hy_rgamma(a);
        if (rga != 0.0) {
            double tr = 1.0, ti = 0.0, sr = 1.0, si = 0.0;
            for (int k = 1; k <= n - 1; k++) {
                double f = (a - n + (k - 1)) / ((double)(k - n) * (double)k);
                c_mul(tr, ti, zr * f, zi * f, &tr, &ti);
                sr += tr; si += ti;
            }
            double pr, pi;
            c_pow_ma(c_abs(zr, zi), li, -(double)n, &pr, &pi);
            c_mul(sr, si, pr, pi, &sr, &si);
            double g = tgamma((double)n) * rga;
            finr = g * sr; fini = g * si;
            peak = c_abs(finr, fini);
        }
    }
    double logr = 0.0, logi = 0.0;
    double rgan = hy_rgamma(a - (double)n);
    if (rgan != 0.0) {
        double pref = ((n & 1) ? 1.0 : -1.0) / tgamma((double)n + 1.0) * rgan;
        double pa, pk, pnk;
        if (!sf_machine_digamma(a, &pa) || !sf_machine_digamma(1.0, &pk) ||
            !sf_machine_digamma((double)(n + 1), &pnk)) return false;
        double tr = 1.0, ti = 0.0, sr = 0.0, si = 0.0, scale = 0.0;
        bool done = false;
        for (int k = 0; k < HY_MAXTERMS; k++) {
            double h = pa - pk - pnk;
            double ar_, ai_;                        /* term * (ln z + h) */
            c_mul(tr, ti, lr + h, li, &ar_, &ai_);
            sr += ar_; si += ai_;
            /* condition scale: the bracket's own rounding is ~u(|ln z|+|psi|s) */
            scale += c_abs(tr, ti) * (c_abs(lr, li) + fabs(pa) + fabs(pk) + fabs(pnk));
            double f = (a + k) / (double)(n + 1 + k) / (double)(k + 1);
            c_mul(tr, ti, zr * f, zi * f, &tr, &ti);
            if (!isfinite(tr) || !isfinite(ti)) return false;
            pa += 1.0 / (a + (double)k);
            pk += 1.0 / (double)(k + 1);
            pnk += 1.0 / (double)(n + k + 1);
            if (tr == 0.0 && ti == 0.0) { done = true; break; }
            if (k > 2 && (double)k > fabs(a) &&
                c_abs(tr, ti) * (c_abs(lr, li) + fabs(pa - pk - pnk)) <=
                    c_abs(sr, si) * HY_TOL) { done = true; break; }
        }
        if (!done) return false;
        logr = pref * sr; logi = pref * si;
        if (fabs(pref) * scale > peak) peak = fabs(pref) * scale;
    }
    double vr = finr + logr, vi = fini + logi, mag = c_abs(vr, vi);
    if (!isfinite(vr) || !isfinite(vi) || mag == 0.0) return false;
    if (peak > HY_KAPPA_MAX * mag) return false;
    *rr = vr; *ri = vi;
    return true;
}

bool sf_machine_hyperu_c(double a, double b, double zr, double zi,
                         double* rr, double* ri) {
    if (!isfinite(a) || !isfinite(b) || !isfinite(zr) || !isfinite(zi)) return false;
    if (zr == 0.0 && zi == 0.0) return false;
    double lr, li;
    c_log(zr, zi, &lr, &li);
    if (hyu_asym_c(a, b, zr, zi, li, rr, ri)) return true;
    double rb = nearbyint(b);
    if (b == rb) {
        if (rb >= 1.0) return hyu_intb_c(a, (int)rb - 1, zr, zi, lr, li, rr, ri);
        /* b <= 0: U(a,b,z) = z^(1-b) U(a-b+1, 2-b, z) */
        double vr, vi, pr, pi;
        if (!sf_machine_hyperu_c(a - b + 1.0, 2.0 - b, zr, zi, &vr, &vi)) return false;
        c_pow_ma(c_abs(zr, zi), li, 1.0 - b, &pr, &pi);
        c_mul(vr, vi, pr, pi, rr, ri);
        return isfinite(*rr) && isfinite(*ri);
    }
    return hyu_conn_c(a, b, zr, zi, li, rr, ri);
}

/* Real U(a, b, z), z > 0, a > 0, by its integral representation [DLMF 13.4.4]
 *   U(a,b,z) = z^-a / Gamma(a) * int_0^inf e^-s s^(a-1) (1 + s/z)^(b-a-1) ds
 * with the double-exponential substitution s = exp(u - e^-u) [Ooura & Mori,
 * J. Comput. Appl. Math. 112 (1999); the exp-sinh rule for an algebraic
 * endpoint and exponential decay] and the trapezoidal rule, step 1/8. This is
 * the moderate-z regime where the asymptotic series has not converged and the
 * connection formulas cancel -- Zhang & Jin's CHGU integrates there too. Every
 * node is positive, so the compensated sum has no cancellation; the one error
 * source is the rounding of the exponent's terms, ~u times their magnitude at
 * the peak, which is measured and must stay below HYQ_EMAG (else decline). */
#define HYQ_H     0.125
#define HYQ_EMAG  12.0

bool sf_machine_hyperu_quad(double a, double b, double z, double* out) {
    if (!(z > 0.0) || !(a > 0.0) || !isfinite(a) || !isfinite(b) || !isfinite(z))
        return false;
    if (a > 160.0) return false;
    double c = b - a - 1.0;
    double sum = 0.0, comp = 0.0, fmax = 0.0, emag = 0.0;
    int nodes = 0;
    for (int dir = 0; dir < 2; dir++) {
        int step = dir == 0 ? 1 : -1;
        for (int k = dir == 0 ? 0 : -1; ; k += step) {
            double u = (double)k * HYQ_H, em = exp(-u), s = exp(u - em);
            double t1 = a * (u - em), t2 = c * log1p(s / z);
            double f = exp(t1 - s + t2) * (1.0 + em);
            if (!isfinite(f)) return false;
            hdd q = dd_two_sum(sum, f); sum = q.hi; comp += q.lo;
            if (f > fmax) { fmax = f; emag = s + fabs(t1) + fabs(t2); }
            if (++nodes > 4000) return false;
            /* stop once a node is negligible and the integrand is past its
             * peak on this side (double-exponential decay from there on) */
            if (f <= 1.0e-18 * sum && (dir == 1 || s > a + fabs(c) + 1.0)) break;
        }
    }
    if (emag > HYQ_EMAG) return false;
    double ga = tgamma(a);
    double v = pow(z, -a) / ga * ((sum + comp) * HYQ_H);
    if (!isfinite(v) || v == 0.0) return false;
    *out = v;
    return true;
}

/* e^(-z/2) z^(m+1/2) as one complex factor (principal power). */
static void hy_whit_pref(double m, double zr, double zi, double* pr, double* pi) {
    double e = m + 0.5;
    double ar = atan2(zi, zr);
    if (zi == 0.0) {                       /* real axis: exact phase via c_pow_ma */
        double qr, qi;
        c_pow_ma(fabs(zr), ar, e, &qr, &qi);
        double ex = exp(-0.5 * zr);
        *pr = qr * ex; *pi = qi * ex;
        return;
    }
    double mod = pow(c_abs(zr, zi), e) * exp(-0.5 * zr);
    double ph = e * ar - 0.5 * zi;
    *pr = mod * cos(ph);
    *pi = mod * sin(ph);
}

bool sf_machine_whittaker_m_c(double k, double m, double zr, double zi,
                              double* rr, double* ri) {
    if (!isfinite(k) || !isfinite(m) || !isfinite(zr) || !isfinite(zi)) return false;
    if (zr == 0.0 && zi == 0.0) return false;
    double hr, hi, pr, pi;
    if (!sf_machine_hyper1f1_c(m - k + 0.5, 0.0, 2.0 * m + 1.0, 0.0, zr, zi, &hr, &hi))
        return false;
    hy_whit_pref(m, zr, zi, &pr, &pi);
    c_mul(hr, hi, pr, pi, rr, ri);
    return isfinite(*rr) && isfinite(*ri);
}

bool sf_machine_whittaker_w_c(double k, double m, double zr, double zi,
                              double* rr, double* ri) {
    if (!isfinite(k) || !isfinite(m) || !isfinite(zr) || !isfinite(zi)) return false;
    if (zr == 0.0 && zi == 0.0) return false;
    double am = fabs(m);                             /* W is even in m */
    double ur, ui, pr, pi;
    if (!sf_machine_hyperu_c(am - k + 0.5, 2.0 * am + 1.0, zr, zi, &ur, &ui)) return false;
    hy_whit_pref(am, zr, zi, &pr, &pi);
    c_mul(ur, ui, pr, pi, rr, ri);
    return isfinite(*rr) && isfinite(*ri);
}
