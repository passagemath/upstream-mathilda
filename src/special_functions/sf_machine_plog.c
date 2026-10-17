/* Machine-precision kernels: ProductLog (Lambert W0) real and complex, and the
 * integer-degree LegendreP / LegendreQ over a complex argument.  See
 * sf_machine_plog.h for the contracts.
 *
 * ProductLog
 * ----------
 * Corless, Gonnet, Hare, Jeffrey & Knuth, "On the Lambert W function",
 * Adv. Comput. Math. 5 (1996) 329-359:
 *   - Halley's iteration (their eq. 5.9)
 *         w <- w - f / (e^w (w + 1) - (w + 2) f / (2w + 2)),   f = w e^w - z;
 *   - the branch-point series (their eq. 4.22) in p = sqrt(2 (e z + 1)),
 *         W0 = -1 + p - p^2/3 + 11/72 p^3 - 43/540 p^4 + 769/17280 p^5 - ...,
 *     principal square root for the principal branch;
 *   - the asymptotic seed L1 - L2 + L2/L1 (their eq. 4.19), L1 = log z,
 *     L2 = log L1, away from the origin and the branch point;
 * plus the (3 + 6z + z^2)/(3 + 9z + 5z^2) Pade seed near the origin.
 *
 * NEAR THE BRANCH POINT z = -1/e the textbook iteration loses everything: W is
 * ill-conditioned there and, worse, f = w e^w - z is formed as the difference of
 * two numbers that are both ~ -1/e, so the residual it steers by is pure
 * rounding.  Two changes make the region exact to a few ulps:
 *   1. d = z + 1/e is formed with 1/e split into a double head and tail
 *      (INVE_HI + INVE_LO): z + INVE_HI is exact by Sterbenz's lemma for z
 *      within a factor of two of -1/e, so d carries full relative accuracy.
 *      This is also what decides the real/complex boundary correctly: the
 *      double nearest -1/e lies BELOW -1/e, so its W0 is complex
 *      (-1 + 8.2e-9 I), which the MPFR path has always said and a test of the
 *      rounded constant cannot.
 *   2. The iteration runs in u = w + 1 on
 *         H(u) = g(u) - e d,   g(u) = (u - 1) e^u + 1 = Sum_{k>=2} (k-1) u^k / k!,
 *      i.e. e (w e^w - z) with the two -1/e's cancelled ANALYTICALLY; g is summed
 *      from its Maclaurin series for |u| < 1, where the closed form cancels.
 *
 * Legendre
 * --------
 * P_n and Q_n by the three-term recurrence (DLMF 14.10.3), which is the form
 * the real kernels in sf_machine.c already use; Q_0 is the inverse hyperbolic
 * tangent (DLMF 14.5.22 with the principal logs of 14.5.21 off the cut),
 * evaluated with the cancellation-free formulas of Kahan's "Branch cuts for
 * complex elementary functions" (1987):
 *     Re Q_0 = log1p(4x / ((1 - x)^2 + y^2)) / 4    (|x| < 1/2),
 *            = (log|1 + z| - log|1 - z|) / 2       (otherwise),
 *     Im Q_0 = (arg(1 + z) - arg(1 - z)) / 2.
 */
#include "sf_machine_plog.h"
#include "sf_machine.h"

#include <math.h>

#ifndef M_E
#define M_E 2.71828182845904523536
#endif
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* 1/e = INVE_HI + INVE_LO to ~2^-110 (mpmath, 40 digits). */
#define INVE_HI  0.36787944117144233
#define INVE_LO  (-1.2428753672788363168e-17)

/* ------------------------------------------------------------------ */
/* Complex helpers on (re, im) pairs                                   */
/* ------------------------------------------------------------------ */

static void c_mul(double ar, double ai, double br, double bi, double* r, double* i) {
    double rr = ar * br - ai * bi;
    double ii = ar * bi + ai * br;
    *r = rr; *i = ii;
}

/* Smith's algorithm: a / b without overflow in the intermediate |b|^2. */
static void c_div(double ar, double ai, double br, double bi, double* r, double* i) {
    if (fabs(br) >= fabs(bi)) {
        double t = bi / br, den = br + bi * t;
        double rr = (ar + ai * t) / den;
        double ii = (ai - ar * t) / den;
        *r = rr; *i = ii;
    } else {
        double t = br / bi, den = br * t + bi;
        double rr = (ar * t + ai) / den;
        double ii = (ai * t - ar) / den;
        *r = rr; *i = ii;
    }
}

static void c_exp(double ar, double ai, double* r, double* i) {
    double m = exp(ar);
    *r = m * cos(ai); *i = m * sin(ai);
}

static void c_log(double ar, double ai, double* r, double* i) {
    *r = log(hypot(ar, ai)); *i = atan2(ai, ar);
}

/* Principal square root; a +0 imaginary part on the negative real axis gives
 * +i sqrt|x|, the upper side of the cut. */
static void c_sqrt(double x, double y, double* r, double* i) {
    if (x == 0.0 && y == 0.0) { *r = 0.0; *i = y; return; }
    double t = sqrt((fabs(x) + hypot(x, y)) * 0.5);
    if (x >= 0.0) { *r = t; *i = y / (2.0 * t); }
    else          { *r = fabs(y) / (2.0 * t); *i = copysign(t, y); }
}

static double c_abs(double r, double i) { return hypot(r, i); }

/* ------------------------------------------------------------------ */
/* ProductLog: branch-point region                                     */
/* ------------------------------------------------------------------ */

/* g(u) = (u - 1) e^u + 1 = Sum_{k>=2} (k-1)/k! u^k: c[k-2] = (k-1)/k!.
 * 24 terms reach 2^-60 relative at |u| = 1. */
static const double PLOG_G[] = {
    5.0e-1,                         /* k=2  1/2!   */
    3.33333333333333333e-1,     /* k=3  2/3!   */
    1.25e-1,                    /* k=4  3/4!   */
    3.33333333333333333e-2,     /* k=5  4/5!   */
    6.94444444444444444e-3,     /* k=6  5/6!   */
    1.19047619047619048e-3,     /* k=7  6/7!   */
    1.73611111111111111e-4,     /* k=8  7/8!   */
    2.20458553791887125e-5,     /* k=9  8/9!   */
    2.48015873015873016e-6,     /* k=10 9/10!  */
    2.50521083854417188e-7,     /* k=11 10/11! */
    2.29644326866549089e-8,     /* k=12 11/12! */
    1.92708526041859375e-9,     /* k=13 12/13! */
    1.49119692770486421e-10,     /* k=14 13/14! */
    1.07060292245477431e-11,     /* k=15 14/15! */
    7.16921599858107795e-13,     /* k=16 15/16! */
    4.49833160695283322e-14,     /* k=17 16/17! */
    2.6552651846596585e-15,     /* k=18 17/18! */
    1.47971434439237935e-16,     /* k=19 18/19! */
    7.80960348429311323e-18,     /* k=20 19/20! */
    3.91458821267825225e-19,     /* k=21 20/21! */
    1.86832619241462039e-20,     /* k=22 21/22! */
    8.50997437538750488e-22,     /* k=23 22/23! */
    3.7069964135210722e-23,     /* k=24 23/24! */
};
#define PLOG_GN ((int)(sizeof PLOG_G / sizeof PLOG_G[0]))

/* g(u) for complex u, |u| < 1 by the series, otherwise the closed form. */
static void plog_g(double ur, double ui, double* gr, double* gi) {
    if (c_abs(ur, ui) < 1.0) {
        double sr = PLOG_G[PLOG_GN - 1], si = 0.0;
        for (int k = PLOG_GN - 2; k >= 0; k--) {
            double tr, ti;
            c_mul(sr, si, ur, ui, &tr, &ti);
            sr = tr + PLOG_G[k]; si = ti;
        }
        double u2r, u2i;
        c_mul(ur, ui, ur, ui, &u2r, &u2i);
        c_mul(sr, si, u2r, u2i, gr, gi);
    } else {
        double er, ei, tr, ti;
        c_exp(ur, ui, &er, &ei);
        c_mul(ur - 1.0, ui, er, ei, &tr, &ti);
        *gr = tr + 1.0; *gi = ti;
    }
}

/* Halley in u = w + 1 on H(u) = g(u) - ed, ed = e (z + 1/e).  Seed: the
 * branch-point series in p = sqrt(2 ed) (principal root -> W0). */
static bool plog_branch(double edr, double edi, double* wr, double* wi) {
    double pr, pi_;
    c_sqrt(2.0 * edr, 2.0 * edi, &pr, &pi_);
    /* u0 = p - p^2/3 + 11/72 p^3 - 43/540 p^4 + 769/17280 p^5 (Horner) */
    double sr = 769.0 / 17280.0, si = 0.0, tr, ti;
    static const double cf[4] = { 1.0, -1.0 / 3.0, 11.0 / 72.0, -43.0 / 540.0 };
    for (int k = 3; k >= 0; k--) {
        c_mul(sr, si, pr, pi_, &tr, &ti);
        sr = tr + cf[k]; si = ti;
    }
    double ur, ui;
    c_mul(sr, si, pr, pi_, &ur, &ui);
    if (edr == 0.0 && edi == 0.0) { *wr = -1.0; *wi = 0.0; return true; }
    for (int it = 0; it < 30; it++) {
        double gr, gi, er, ei;
        plog_g(ur, ui, &gr, &gi);
        double hr = gr - edr, hi = gi - edi;
        c_exp(ur, ui, &er, &ei);
        double h1r, h1i, h2r, h2i;
        c_mul(ur, ui, er, ei, &h1r, &h1i);              /* H'  = u e^u     */
        c_mul(ur + 1.0, ui, er, ei, &h2r, &h2i);        /* H'' = (u+1) e^u */
        if (h1r == 0.0 && h1i == 0.0) break;
        double qr, qi, dr, di;
        c_mul(hr, hi, h2r, h2i, &qr, &qi);
        c_div(qr, qi, 2.0 * h1r, 2.0 * h1i, &qr, &qi);  /* H H'' / (2 H')  */
        c_div(hr, hi, h1r - qr, h1i - qi, &dr, &di);
        if (!isfinite(dr) || !isfinite(di)) return false;
        ur -= dr; ui -= di;
        if (c_abs(dr, di) <= 1e-15 * c_abs(ur, ui)) break;
    }
    if (!isfinite(ur) || !isfinite(ui)) return false;
    /* A converged u solves H = 0; check it rather than trust the count. */
    double gr, gi;
    plog_g(ur, ui, &gr, &gi);
    double scale = c_abs(edr, edi) + 1e-300;
    if (c_abs(gr - edr, gi - edi) > 1e-12 * scale) return false;
    *wr = ur - 1.0; *wi = ui;
    return true;
}

/* ------------------------------------------------------------------ */
/* ProductLog: real                                                    */
/* ------------------------------------------------------------------ */

#define PLOG_EPS     1.0e-16
#define PLOG_BRANCH  0.5          /* |e (z + 1/e)| below this: the u-iteration */

/* Lambert W, principal branch, real x >= -1/e.
 *
 * THE STARTING GUESS IS THE WHOLE ALGORITHM, and the asymptotic form used to be
 * applied outside its domain: `l1 - l2 + l2/l1` with l1 = log x needs l1
 * comfortably above 1, and at x = 1 it is log(0)/0 = NaN.  ProductLog[1.01]
 * once came back -338.392 from this kernel (Compile[] and the packed path).
 * log(1 + x) is the guess for the band the asymptotic form cannot serve; it is
 * an upper bound on W for x >= 0 and tight, so Halley converges in 3-4 steps.
 * Near -1/e the branch-point u-iteration above takes over (see the file
 * comment); x below -1/e (decided from the SPLIT constant) is complex and
 * declines. */
bool sf_machine_productlog(double x, double* out) {
    if (!isfinite(x)) return false;
    if (x == 0.0) { *out = 0.0; return true; }
    if (x < 0.0) {
        double d = (x + INVE_HI) + INVE_LO;               /* x + 1/e, exact-ish */
        if (d < 0.0) return false;                        /* complex below -1/e */
        double ed = M_E * d;
        if (ed < PLOG_BRANCH) {
            double wr, wi;
            if (!plog_branch(ed, 0.0, &wr, &wi)) return false;
            *out = wr;
            return true;
        }
    }

    double w;
    if (x < 1.0) {
        w = x * (1.0 - x + 1.5 * x * x);
    } else {
        double l1 = log(x);
        if (l1 > 1.0) {                                   /* x > e: asymptotic */
            double l2 = log(l1);
            w = l1 - l2 + l2 / l1;
        } else {                                          /* 1 <= x <= e */
            w = log(1.0 + x);
        }
    }
    for (int i = 0; i < 60; i++) {
        double e = exp(w), we = w * e, f = we - x;
        double denom = e * (w + 1.0) - (w + 2.0) * f / (2.0 * w + 2.0);
        if (denom == 0.0) break;
        double dw = f / denom;
        w -= dw;
        if (fabs(dw) <= PLOG_EPS * (fabs(w) + 1.0)) break;
    }
    if (!isfinite(w)) return false;

    /* VERIFY, then answer: an iteration that stalls or lands on the wrong
     * branch returns a finite number that is not W(x).  w + log w == log x for
     * w > 0 (w e^w overflows for large x); the direct form for w < 0, |w| < 1. */
    double resid, scale;
    if (w > 0.0) {
        resid = fabs((w + log(w)) - log(x));
        scale = fabs(log(x)) + 1.0;
    } else {
        resid = fabs(w * exp(w) - x);
        scale = fabs(x) + 1e-300;
    }
    if (resid > 1e-9 * scale) return false;
    *out = w;
    return true;
}

/* ------------------------------------------------------------------ */
/* ProductLog: complex                                                 */
/* ------------------------------------------------------------------ */

/* Is w in the range of W0?  The principal branch's image is bounded by
 * |Im w| < Pi and the curve u = -v cot v (Corless et al. Fig. 4); on the real
 * axis it is w >= -1.  A small slack keeps a value on the boundary curve (z on
 * the cut) from being rejected for its rounding. */
static bool plog_in_w0(double wr, double wi) {
    double v = fabs(wi);
    if (!(v < M_PI)) return false;
    double slack = 1e-9 * (1.0 + fabs(wr));
    if (v == 0.0) return wr >= -1.0 - slack;
    return wr > -v / tan(v) - slack;
}

bool sf_machine_productlog_c(double zr, double zi, double* wr, double* wi) {
    zi += 0.0;                                            /* -0 -> +0: the cut's upper side */
    if (!isfinite(zr) || !isfinite(zi)) return false;
    if (zi == 0.0) {
        double v;
        if (sf_machine_productlog(zr, &v)) { *wr = v; *wi = 0.0; return true; }
        if (zr >= 0.0) return false;                      /* real kernel declined */
    }
    if (zr == 0.0 && zi == 0.0) { *wr = 0.0; *wi = 0.0; return true; }

    /* Branch-point region: d = z + 1/e with the split constant. */
    double dr = (zr + INVE_HI) + INVE_LO, di = zi;
    double edr = M_E * dr, edi = M_E * di;
    if (c_abs(edr, edi) < PLOG_BRANCH)
        return plog_branch(edr, edi, wr, wi);

    /* Seed (as in Corless et al. / SciPy's lambertw): Pade near the origin,
     * else the asymptotic L1 - L2 + L2/L1. */
    double w_r, w_i;
    if (zr > -1.0 && zr < 1.5 && fabs(zi) < 1.0 && -2.5 * fabs(zi) - 0.2 < zr) {
        double z2r, z2i, nr, ni, de_r, de_i;
        c_mul(zr, zi, zr, zi, &z2r, &z2i);
        c_mul(zr, zi, 3.0 + 6.0 * zr + z2r, 6.0 * zi + z2i, &nr, &ni);
        de_r = 3.0 + 9.0 * zr + 5.0 * z2r; de_i = 9.0 * zi + 5.0 * z2i;
        c_div(nr, ni, de_r, de_i, &w_r, &w_i);
    } else {
        double l1r, l1i, l2r, l2i, qr, qi;
        c_log(zr, zi, &l1r, &l1i);
        c_log(l1r, l1i, &l2r, &l2i);
        c_div(l2r, l2i, l1r, l1i, &qr, &qi);
        w_r = l1r - l2r + qr; w_i = l1i - l2i + qi;
        if (!isfinite(w_r) || !isfinite(w_i)) { w_r = l1r; w_i = l1i; }
    }

    bool conv = false;
    for (int it = 0; it < 60; it++) {
        double er, ei, fr, fi;
        c_exp(w_r, w_i, &er, &ei);
        c_mul(w_r, w_i, er, ei, &fr, &fi);
        fr -= zr; fi -= zi;                               /* f = w e^w - z */
        double ar, ai, qr, qi, sr, si;
        c_mul(er, ei, w_r + 1.0, w_i, &ar, &ai);          /* e^w (w + 1) */
        c_mul(w_r + 2.0, w_i, fr, fi, &qr, &qi);
        c_div(qr, qi, 2.0 * w_r + 2.0, 2.0 * w_i, &qr, &qi);
        c_div(fr, fi, ar - qr, ai - qi, &sr, &si);
        if (!isfinite(sr) || !isfinite(si)) return false;
        w_r -= sr; w_i -= si;
        if (c_abs(sr, si) <= 1e-15 * c_abs(w_r, w_i)) { conv = true; break; }
    }
    if (!conv || !isfinite(w_r) || !isfinite(w_i)) return false;
    if (!plog_in_w0(w_r, w_i)) return false;

    /* Verify w + log w = log z (mod 2 Pi i): bounded where w e^w would overflow. */
    double lwr, lwi, lzr, lzi;
    c_log(w_r, w_i, &lwr, &lwi);
    c_log(zr, zi, &lzr, &lzi);
    double rr = w_r + lwr - lzr, ri = w_i + lwi - lzi;
    ri -= 2.0 * M_PI * nearbyint(ri / (2.0 * M_PI));
    if (c_abs(rr, ri) > 1e-10 * (c_abs(lzr, lzi) + 1.0)) return false;
    *wr = w_r; *wi = w_i;
    return true;
}

/* ------------------------------------------------------------------ */
/* Legendre P_n, Q_n over a complex argument                           */
/* ------------------------------------------------------------------ */

#define LEG_MAXDEG 1000.0

bool sf_machine_legendre_p_c(double n, double xr, double xi, double* rr, double* ri) {
    if (!isfinite(n) || n != floor(n)) return false;
    if (n < 0.0) n = -n - 1.0;                            /* P_{-n-1} = P_n */
    if (n > LEG_MAXDEG) return false;
    if (xi == 0.0) {
        double v;
        if (!sf_machine_legendre_p(n, xr, &v)) return false;
        *rr = v; *ri = 0.0;
        return true;
    }
    if (!isfinite(xr) || !isfinite(xi)) return false;
    int N = (int)n;
    double p0r = 1.0, p0i = 0.0, p1r = xr, p1i = xi;
    if (N == 0) { *rr = 1.0; *ri = 0.0; return true; }
    for (int k = 2; k <= N; k++) {
        double tr, ti;
        c_mul(xr, xi, p1r, p1i, &tr, &ti);
        double a = 2.0 * k - 1.0, b = k - 1.0;
        double pr = (a * tr - b * p0r) / (double)k;
        double pi_ = (a * ti - b * p0i) / (double)k;
        p0r = p1r; p0i = p1i; p1r = pr; p1i = pi_;
    }
    if (!isfinite(p1r) || !isfinite(p1i)) return false;
    *rr = p1r; *ri = p1i;
    return true;
}

/* Q_n = P_n Q_0 - W_{n-1} (Christoffel's form, W_{n-1} a polynomial), with P_n and
 * W_{n-1} both run up the three-term recurrence (W obeys it too, from
 * W_{-1} = 0, W_0 = 1).  Running the recurrence on Q itself instead seeds it
 * with Q_1 = x Q_0 - 1, which cancels badly wherever x Q_0 ~ 1 (|x| ~ 0.8) and
 * carries that error into every higher degree: 9 ulp at LegendreQ[2, -0.77 -
 * 0.05 I] where this form gives 1-2.  The one subtraction left is the
 * function's own cancellation near its zeros. */
bool sf_machine_legendre_q(double n, double x, double* out) {
    if (n != floor(n) || n < 0.0 || n > LEG_MAXDEG) return false;
    if (!(x > -1.0 && x < 1.0)) return false;     /* also rejects NaN */
    double q0 = atanh(x);
    int N = (int)n;
    if (N == 0) { *out = q0; return isfinite(q0); }
    double p0 = 1.0, p1 = x, w0 = 0.0, w1 = 1.0;
    for (int k = 1; k < N; k++) {
        double a = 2.0 * k + 1.0, b = (double)k, c = k + 1.0;
        double p2 = (a * x * p1 - b * p0) / c;
        double w2 = (a * x * w1 - b * w0) / c;
        p0 = p1; p1 = p2; w0 = w1; w1 = w2;
    }
    double v = p1 * q0 - w1;
    *out = v;
    return isfinite(v);
}

bool sf_machine_legendre_q_c(double n, double xr, double xi, double* rr, double* ri) {
    if (!isfinite(n) || n != floor(n) || n < 0.0 || n > LEG_MAXDEG) return false;
    if (!isfinite(xr) || !isfinite(xi)) return false;
    xi += 0.0;                                            /* -0 -> +0 */
    if (xi == 0.0 && xr > -1.0 && xr < 1.0) {
        double v;
        if (!sf_machine_legendre_q(n, xr, &v)) return false;
        *rr = v; *ri = 0.0;
        return true;
    }
    if (xi == 0.0 && (xr == 1.0 || xr == -1.0)) return false;   /* log singularity */

    /* Q_0 = atanh-type, Kahan's cancellation-free real part; the imaginary part
     * as (arg(1+z) - arg(1-z))/2 with +0 imaginary parts on the real axis, so a
     * real x > 1 gets Log[1 - x] = log(x-1) + i Pi exactly as the closed form. */
    double omx = 1.0 - xr;
    double q0r;
    if (fabs(xr) < 0.5) {
        /* |1+z| ~ |1-z|: the ratio form, no cancellation between two logs. */
        q0r = 0.25 * log1p(4.0 * xr / (omx * omx + xi * xi));
    } else {
        /* Near |x| = 1 the ratio 4x/den -> -1 and log1p would magnify its
         * rounding ~1/|1+x|-fold; 1 +- x is exact here (Sterbenz), so the
         * difference of the two logs is accurate instead. */
        q0r = 0.5 * (log(hypot(1.0 + xr, xi)) - log(hypot(omx, xi)));
    }
    double q0i = 0.5 * (atan2(xi + 0.0, 1.0 + xr) - atan2(0.0 - xi, omx));
    int N = (int)n;
    if (N == 0) {
        if (!isfinite(q0r) || !isfinite(q0i)) return false;
        *rr = q0r; *ri = q0i; return true;
    }
    double p0r = 1.0, p0i = 0.0, p1r = xr, p1i = xi;
    double w0r = 0.0, w0i = 0.0, w1r = 1.0, w1i = 0.0;
    for (int k = 1; k < N; k++) {
        double a = 2.0 * k + 1.0, b = (double)k, c = k + 1.0, tr, ti, ur, ui;
        c_mul(xr, xi, p1r, p1i, &tr, &ti);
        c_mul(xr, xi, w1r, w1i, &ur, &ui);
        double p2r = (a * tr - b * p0r) / c, p2i = (a * ti - b * p0i) / c;
        double w2r = (a * ur - b * w0r) / c, w2i = (a * ui - b * w0i) / c;
        p0r = p1r; p0i = p1i; p1r = p2r; p1i = p2i;
        w0r = w1r; w0i = w1i; w1r = w2r; w1i = w2i;
    }
    double qkr, qki;
    c_mul(p1r, p1i, q0r, q0i, &qkr, &qki);
    qkr -= w1r; qki -= w1i;
    if (!isfinite(qkr) || !isfinite(qki)) return false;
    *rr = qkr; *ri = qki;
    return true;
}
