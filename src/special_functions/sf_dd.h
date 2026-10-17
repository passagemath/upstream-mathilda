/* sf_dd.h -- double-double (~106-bit) real and complex arithmetic, header-only.
 *
 * A value is hi + lo with |lo| <= ulp(hi)/2. The kernels are the classical
 * error-free transformations (Knuth's two-sum, Dekker's two-product by
 * Veltkamp splitting -- no fused multiply-add, which the default x86-64 target
 * lacks, so the results are identical on every C99 target) and the QD-library
 * algorithms for exp / log / sin / cos (Hida, Li & Bailey, "Library for
 * double-double and quad-double arithmetic", 2007):
 *
 *   exp:  e^a = 2^k (1 + s),  s = expm1(r 2^-10) by Taylor, then squared up ten
 *         times as s <- s (s + 2), which keeps s's RELATIVE accuracy (squaring
 *         1 + s directly would not).
 *   log:  one Newton step x0 + (a - e^x0) / e^x0 from x0 = log(a.hi): the
 *         double start has 53 bits, the step doubles them.
 *   sin/cos: reduction by pi/2 (pi/2 as a double-double), Taylor to degree 31
 *         on |r| <= pi/4.
 *
 * Used where a double sum measurably loses bits to cancellation and the result
 * must still come back correct to double precision (sf_machine_meijerg.c's
 * log-case residue series). Valid while every magnitude stays inside
 * [2^-960, 2^995] (Veltkamp splitting needs the headroom); callers check.
 */
#ifndef MATHILDA_SF_DD_H
#define MATHILDA_SF_DD_H

#include <math.h>

typedef struct { double hi, lo; } sdd;
typedef struct { sdd re, im; } scdd;

/* ---- error-free transformations ---------------------------------------- */

static inline sdd sdd_two_sum(double a, double b) {
    sdd r; r.hi = a + b;
    double bb = r.hi - a;
    r.lo = (a - (r.hi - bb)) + (b - bb);
    return r;
}
static inline sdd sdd_quick(double a, double b) {        /* requires |a| >= |b| */
    sdd r; r.hi = a + b; r.lo = b - (r.hi - a); return r;
}
static inline sdd sdd_two_prod(double a, double b) {
    const double SPLIT = 134217729.0;                      /* 2^27 + 1 */
    double t = SPLIT * a, ah = t - (t - a), al = a - ah;
    t = SPLIT * b;
    double bh = t - (t - b), bl = b - bh;
    sdd r; r.hi = a * b;
    r.lo = ((ah * bh - r.hi) + ah * bl + al * bh) + al * bl;
    return r;
}

/* ---- real double-double ------------------------------------------------- */

static inline sdd sdd_from(double a) { sdd r; r.hi = a; r.lo = 0.0; return r; }
static inline sdd sdd_neg(sdd a) { a.hi = -a.hi; a.lo = -a.lo; return a; }
static inline sdd sdd_add(sdd a, sdd b) {
    sdd s = sdd_two_sum(a.hi, b.hi);
    sdd t = sdd_two_sum(a.lo, b.lo);
    s.lo += t.hi;
    s = sdd_quick(s.hi, s.lo);
    s.lo += t.lo;
    return sdd_quick(s.hi, s.lo);
}
static inline sdd sdd_sub(sdd a, sdd b) { return sdd_add(a, sdd_neg(b)); }
static inline sdd sdd_add_d(sdd a, double b) {
    sdd s = sdd_two_sum(a.hi, b);
    s.lo += a.lo;
    return sdd_quick(s.hi, s.lo);
}
static inline sdd sdd_mul(sdd a, sdd b) {
    sdd p = sdd_two_prod(a.hi, b.hi);
    p.lo += a.hi * b.lo + a.lo * b.hi;
    return sdd_quick(p.hi, p.lo);
}
static inline sdd sdd_mul_d(sdd a, double b) {
    sdd p = sdd_two_prod(a.hi, b);
    p.lo += a.lo * b;
    return sdd_quick(p.hi, p.lo);
}
static inline sdd sdd_ldexp(sdd a, int e) { a.hi = ldexp(a.hi, e); a.lo = ldexp(a.lo, e); return a; }
static inline sdd sdd_div(sdd a, sdd b) {
    double q1 = a.hi / b.hi;
    sdd r = sdd_sub(a, sdd_mul_d(b, q1));
    double q2 = r.hi / b.hi;
    r = sdd_sub(r, sdd_mul_d(b, q2));
    double q3 = r.hi / b.hi;
    return sdd_add_d(sdd_quick(q1, q2), q3);
}
static inline sdd sdd_div_d(sdd a, double b) { return sdd_div(a, sdd_from(b)); }

/* ln 2 and pi/2, pi as double-doubles. */
#define SDD_LN2_HI  6.931471805599452862e-01
#define SDD_LN2_LO  2.319046813846299558e-17
#define SDD_PIO2_HI 1.570796326794896558e+00
#define SDD_PIO2_LO 6.123233995736766036e-17
#define SDD_PI_HI   3.141592653589793116e+00
#define SDD_PI_LO   1.224646799147353207e-16

/* e^a. Overflow gives +inf in hi (callers test isfinite), deep underflow 0. */
static inline sdd sdd_exp(sdd a) {
    if (a.hi > 709.0)  { sdd r = { INFINITY, 0.0 }; return r; }
    if (a.hi < -708.0) { sdd r = { 0.0, 0.0 }; return r; }
    double k = nearbyint(a.hi / SDD_LN2_HI);
    sdd ln2 = { SDD_LN2_HI, SDD_LN2_LO };
    sdd r = sdd_sub(a, sdd_mul_d(ln2, k));
    r = sdd_ldexp(r, -10);                             /* |r| <= 3.4e-4 */
    /* expm1(r) = r (1 + r/2 (1 + r/3 (... (1 + r/9)))) */
    sdd t = sdd_from(1.0);
    for (int n = 9; n >= 2; n--) t = sdd_add_d(sdd_div_d(sdd_mul(r, t), (double)n), 1.0);
    sdd s = sdd_mul(r, t);
    for (int i = 0; i < 10; i++) s = sdd_mul(s, sdd_add_d(s, 2.0));
    return sdd_ldexp(sdd_add_d(s, 1.0), (int)k);
}

/* ln a for a > 0 (a.hi > 0). */
static inline sdd sdd_log(sdd a) {
    double x0 = log(a.hi);
    sdd e = sdd_exp(sdd_from(x0));
    return sdd_add_d(sdd_div(sdd_sub(a, e), e), x0);
}

/* sin and cos of a double-double angle of moderate size (|a| < 2^20). */
static inline void sdd_sincos(sdd a, sdd* s, sdd* c) {
    double k = nearbyint(a.hi / SDD_PIO2_HI);
    sdd pio2 = { SDD_PIO2_HI, SDD_PIO2_LO };
    sdd r = sdd_sub(a, sdd_mul_d(pio2, k));            /* |r| <= pi/4 (+ulp) */
    sdd r2 = sdd_mul(r, r);
    /* sin r = r (1 - r^2/(2*3) (1 - r^2/(4*5) (...))), to degree 31;
     * cos r = 1 - r^2/(1*2) (1 - r^2/(3*4) (...)), to degree 30. */
    sdd ts = sdd_from(1.0), tc = sdd_from(1.0);
    for (int n = 15; n >= 1; n--) {
        ts = sdd_add_d(sdd_neg(sdd_div_d(sdd_mul(r2, ts), (double)(2 * n) * (double)(2 * n + 1))), 1.0);
        tc = sdd_add_d(sdd_neg(sdd_div_d(sdd_mul(r2, tc), (double)(2 * n - 1) * (double)(2 * n))), 1.0);
    }
    sdd sr = sdd_mul(r, ts), cr = tc;
    long q = (long)fmod(k, 4.0);
    if (q < 0) q += 4;
    switch (q) {
        case 0:  *s = sr;          *c = cr;          break;
        case 1:  *s = cr;          *c = sdd_neg(sr); break;
        case 2:  *s = sdd_neg(sr); *c = sdd_neg(cr); break;
        default: *s = sdd_neg(cr); *c = sr;          break;
    }
}

/* ---- complex double-double ---------------------------------------------- */

static inline scdd scdd_from(double re, double im) {
    scdd r; r.re = sdd_from(re); r.im = sdd_from(im); return r;
}
static inline scdd scdd_add(scdd a, scdd b) {
    scdd r; r.re = sdd_add(a.re, b.re); r.im = sdd_add(a.im, b.im); return r;
}
static inline scdd scdd_mul(scdd a, scdd b) {
    scdd r;
    r.re = sdd_sub(sdd_mul(a.re, b.re), sdd_mul(a.im, b.im));
    r.im = sdd_add(sdd_mul(a.re, b.im), sdd_mul(a.im, b.re));
    return r;
}
static inline scdd scdd_mul_d(scdd a, double b) {
    scdd r; r.re = sdd_mul_d(a.re, b); r.im = sdd_mul_d(a.im, b); return r;
}
/* 1/a = conj(a) / |a|^2. */
static inline scdd scdd_inv(scdd a) {
    sdd d = sdd_add(sdd_mul(a.re, a.re), sdd_mul(a.im, a.im));
    scdd r; r.re = sdd_div(a.re, d); r.im = sdd_neg(sdd_div(a.im, d)); return r;
}

/* Principal log of the double z = zr + i zi (z != 0), on the branch of
 * sf_machine_meijerg.c's c_log: a zero imaginary part is read as +0, so the
 * negative real axis has arg +pi. The argument is atan2 refined by one Newton
 * step against the double-double sin/cos. */
static inline scdd scdd_log_d(double zr, double zi) {
    scdd L;
    if (zi == 0.0) {
        L.re = sdd_log(sdd_from(fabs(zr)));
        L.im = (zr > 0.0) ? sdd_from(0.0) : (sdd){ SDD_PI_HI, SDD_PI_LO };
        return L;
    }
    sdd m2 = sdd_add(sdd_two_prod(zr, zr), sdd_two_prod(zi, zi));
    L.re = sdd_mul_d(sdd_log(m2), 0.5);
    double t0 = atan2(zi, zr);
    sdd s, c;
    sdd_sincos(sdd_from(t0), &s, &c);
    /* z e^{-i t0} = u + i v, and arg z = t0 + atan(v/u) with v/u ~ 1e-16. */
    sdd u = sdd_add(sdd_mul_d(c, zr), sdd_mul_d(s, zi));
    sdd v = sdd_sub(sdd_mul_d(c, zi), sdd_mul_d(s, zr));
    L.im = sdd_add_d(sdd_div(v, u), t0);
    return L;
}

static inline scdd scdd_exp(scdd a) {
    sdd m = sdd_exp(a.re), s, c;
    sdd_sincos(a.im, &s, &c);
    scdd r; r.re = sdd_mul(m, c); r.im = sdd_mul(m, s); return r;
}

#endif /* MATHILDA_SF_DD_H */
