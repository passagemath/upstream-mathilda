/* Mathilda -- Nielsen's generalized polylogarithm S_{n,p}(z) = PolyLog[n, p, z]
 * at machine precision.
 *
 * The double-precision twin of flint_num_nielsen_polylog (src/flint_num_bridge.c),
 * and the kernel the scalar builtin, the NDArray / packed map and Compile[] all
 * call, so the three surfaces agree bit for bit. It evaluates exactly the calls
 * the interpreter evaluates (polylog_nielsen in polylog.c) and declines the rest:
 *
 *   integer n <= 0    any z: (z d/dz)^(-n) (-Log(1 - z))^p / Gamma(p + 1), as a
 *                     polynomial in L = -Log(1 - z) and w = 1/(1 - z)
 *                     (z dL/dz = w - 1, z dw/dz = w^2 - w)
 *   |z| <= 1, z != 1  any n
 *   z == 1            Re n > 0
 *   |z| > 1           integer n >= 1 and a positive integer p
 *
 * with Re p > 0 throughout. p == 1 is Li_{n+1}(z), through the PolyLog kernels.
 * Otherwise: the Stirling series for an integer p and |z| <= 3/4, the Miller
 * series for any p at |z| <= 1/2, and tanh-sinh quadrature of
 *
 *   S_{n,p}(z) = 1/(Gamma(n) Gamma(p+1)) Int_0^1 (-Log t)^(n-1) (-Log(1 - z t))^p dt/t
 *
 * (Re n > 0), split at t = 1/z for a real z > 1, on the principal branch. A
 * quadrature that does not converge to 2e-14 relative, or a series that would
 * need more than 400 terms, DECLINES: the caller's arb path answers instead.
 */
#include "sf_machine_nielsen.h"
#include "sf_machine.h"        /* sf_machine_polylog */
#include "sf_machine_zeta.h"   /* sf_polylog_c */

#include <complex.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

typedef double _Complex cplx;

/* A complex from its parts, keeping a +0 imaginary part +0 (x + 0.0 * I can
 * pick up a -0 through the arithmetic, which moves a value across a cut). */
static cplx CMPLX_(double re, double im) {
    cplx z;
    ((double*)&z)[0] = re;
    ((double*)&z)[1] = im;
    return z;
}

static bool nl_finite(cplx v) { return isfinite(creal(v)) && isfinite(cimag(v)); }

/* 1 - z, with a +0 imaginary part when z is real. `1.0 - z` on a cplx makes it
 * -0.0, and clog then takes the LOWER side of the cut at z > 1, where the
 * principal value (and the reference) take + I Pi. */
static cplx nl_one_minus(cplx z) {
    double zi = cimag(z);
    return CMPLX_(1.0 - creal(z), zi == 0.0 ? 0.0 : -zi);
}

/* Log(1 + u) without the cancellation of forming 1 + u for small u. */
static cplx nl_clog1p(cplx u) {
    double ur = creal(u), ui = cimag(u);
    double re = 0.5 * log1p(2.0 * ur + ur * ur + ui * ui);
    return re + atan2(ui, 1.0 + ur) * I;
}

/* cos(pi x), sin(pi x) exactly zero at the half-integers and integers. */
static void nl_sincospi(double x, double* sn, double* cs) {
    double r = x - 2.0 * floor(x / 2.0);          /* [0, 2) */
    if (r == 0.0)      { *sn = 0.0;  *cs = 1.0;  return; }
    if (r == 0.5)      { *sn = 1.0;  *cs = 0.0;  return; }
    if (r == 1.0)      { *sn = 0.0;  *cs = -1.0; return; }
    if (r == 1.5)      { *sn = -1.0; *cs = 0.0;  return; }
    *sn = sin(M_PI * r); *cs = cos(M_PI * r);
}

/* b^e on the principal branch; a negative real base with a real exponent takes
 * the phase e^(I Pi e) exactly, so a half-integer e leaves no rounding dust in
 * the part that vanishes. */
static cplx nl_pow(cplx b, cplx e) {
    if (cimag(b) == 0.0 && creal(b) < 0.0) {
        /* arg = +Pi whatever the sign of the zero imaginary part: -Log of a
         * positive number carries -0.0 there, and cpow would take -Pi. */
        if (cimag(e) == 0.0) {
            double sn, cs, m = pow(-creal(b), creal(e));
            nl_sincospi(creal(e), &sn, &cs);
            return CMPLX_(m * cs, m * sn);
        }
        return cexp(e * CMPLX_(log(-creal(b)), M_PI));
    }
    if (cimag(b) == 0.0 && creal(b) > 0.0 && cimag(e) == 0.0)
        return pow(creal(b), creal(e));
    /* A complex base with a small integer exponent: binary powering. cpow goes
     * through exp(e Log b), which leaves rounding dust in a part that is zero --
     * (-I Pi)^2 came back as -Pi^2 - 6e-16 I, so the real PolyLog[0, 2, 2.] =
     * Log(-1)^2 / 2 = -Pi^2/2 printed as a complex number. */
    if (cimag(e) == 0.0 && creal(e) == floor(creal(e)) && fabs(creal(e)) <= 64.0) {
        int k = (int)creal(e);
        unsigned m = (unsigned)(k < 0 ? -k : k);
        cplx r = 1.0, sq = b;
        while (m) {
            if (m & 1u) r *= sq;
            m >>= 1;
            if (m) sq *= sq;
        }
        return k < 0 ? 1.0 / r : r;
    }
    return cpow(b, e);
}

/* 1 / Gamma(x) for complex x, via the real lgamma when x is real. */
static cplx nl_rgamma(cplx x) {
    if (cimag(x) == 0.0) {
        double r = creal(x);
        if (r <= 0.0 && r == floor(r)) return 0.0;
        double g = tgamma(r);
        return isfinite(g) ? 1.0 / g : 0.0;
    }
    /* Lanczos (g = 7, n = 9) on the right half plane, reflection on the left. */
    static const double c[9] = {
        0.99999999999980993, 676.5203681218851, -1259.1392167224028,
        771.32342877765313, -176.61502916214059, 12.507343278686905,
        -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7 };
    if (creal(x) < 0.5)
        return csin(M_PI * x) / M_PI * (1.0 / nl_rgamma(1.0 - x));
    cplx z = x - 1.0, a = c[0];
    for (int i = 1; i < 9; i++) a += c[i] / (z + (double)i);
    cplx t = z + 7.5;
    cplx g = sqrt(2.0 * M_PI) * nl_pow(t, z + 0.5) * cexp(-t) * a;
    return 1.0 / g;
}

/* ---- integer n <= 0: the closed form ---------------------------------- */

#define NL_MAXM 20

static bool nl_nonpositive(int m, cplx p, cplx z, cplx* out) {
    /* c[j][k] multiplies L^(p - j) w^k; start (1/Gamma(p+1)) L^p. */
    static cplx cc[2][NL_MAXM + 2][NL_MAXM + 2];
    memset(cc, 0, sizeof cc);
    cc[0][0][0] = nl_rgamma(p + 1.0);
    int cur = 0;
    for (int it = 0; it < m; it++) {
        int nx = 1 - cur;
        memset(cc[nx], 0, sizeof cc[nx]);
        for (int j = 0; j <= it; j++)
            for (int k = 0; k <= it; k++) {
                cplx a = cc[cur][j][k];
                if (a == 0.0) continue;
                cplx q = p - (double)j;
                /* q L^(q-1) (w^(k+1) - w^k)  +  k L^q (w^(k+1) - w^k) */
                cc[nx][j + 1][k + 1] += a * q;
                cc[nx][j + 1][k]     -= a * q;
                cc[nx][j][k + 1]     += a * (double)k;
                cc[nx][j][k]         -= a * (double)k;
            }
        cur = nx;
    }
    cplx omz = nl_one_minus(z);
    cplx L = -clog(omz), w = 1.0 / omz, s = 0.0;
    for (int j = 0; j <= m; j++) {
        cplx Lq = nl_pow(L, p - (double)j), wk = 1.0;
        for (int k = 0; k <= m; k++) {
            s += cc[cur][j][k] * Lq * wk;
            wk *= w;
        }
    }
    *out = s;
    return nl_finite(s);
}

/* ---- the series ------------------------------------------------------- */

/* Stirling: S_{n,p}(z) = Sum_{k>=p} (c(k, p)/k!) z^k / k^n, integer p. */
static bool nl_stirling(cplx n, int p, cplx z, cplx* out) {
    double az = cabs(z);
    int kmax = az == 0.0 ? p + 1 : (int)(55.0 / -log2(az)) + 4 * p + 40;
    if (kmax > 4000 || p > 60) return false;
    double e[62];
    memset(e, 0, sizeof e);
    e[0] = 1.0;
    cplx zk = 1.0, s = 0.0;
    for (int k = 0; k < kmax; k++) {
        for (int j = p; j >= 1; j--) e[j] = ((double)k * e[j] + e[j - 1]) / (double)(k + 1);
        e[0] = (double)k * e[0] / (double)(k + 1);
        zk *= z;
        int kp = k + 1;
        if (kp < p) continue;
        cplx t = zk * e[p] / nl_pow((double)kp, n);
        s += t;
        if (kp > 4 * p + 20 && cabs(t) < 1e-18 * cabs(s)) break;
    }
    *out = s;
    return nl_finite(s);
}

/* Miller: (1/Gamma(p+1)) Sum_k c_k z^(p+k) / (p+k)^n, c the coefficients of
 * (-Log(1-w)/w)^p by c_k = (1/k) Sum_j ((p+1) j - k) c_{k-j} / (j + 1). */
static bool nl_miller(cplx n, cplx p, cplx z, cplx* out) {
    double az = cabs(z);
    int K = az == 0.0 ? 1 : (int)(55.0 / -log2(az)) + 10;
    if (K > 400) return false;
    cplx c[401];
    c[0] = 1.0;
    for (int k = 1; k <= K; k++) {
        cplx s = 0.0;
        for (int j = 1; j <= k; j++) s += ((p + 1.0) * (double)j - (double)k) * c[k - j] / (double)(j + 1);
        c[k] = s / (double)k;
    }
    cplx zk = 1.0, s = 0.0;
    for (int k = 0; k <= K; k++) {
        s += c[k] * zk / nl_pow(p + (double)k, n);
        zk *= z;
    }
    *out = s * nl_pow(z, p) * nl_rgamma(p + 1.0);
    return nl_finite(*out);
}

/* ---- tanh-sinh quadrature --------------------------------------------- */

typedef struct {
    cplx n1, p, z;
    double c, omc;   /* split point 1/z and 1 - 1/z (real z > 1) */
    int piece;       /* 0 whole, 1 left [0, 1/z], 2 right [1/z, 1] */
    bool p_int; int p_i;
} nl_q;

static cplx nl_integrand(const nl_q* q, double tau, double omtau) {
    double t, omt;
    cplx w, lw;
    if (q->piece == 0) {
        t = tau; omt = omtau;
        cplx u = q->z * t;
        lw = cabs(u) < 0.5 ? -nl_clog1p(-u) : -clog(omtau - (q->z - 1.0) * tau);
    } else if (q->piece == 1) {
        t = q->c * tau; omt = 1.0 - t;
        lw = tau < 0.5 ? -log1p(-tau) : -log(omtau);
    } else {
        t = q->c + q->omc * tau; omt = q->omc * omtau;
        w = CMPLX_((1.0 - creal(q->z)) * tau, 0.0);   /* negative real: + I Pi */
        lw = -clog(w);
    }
    double lt = log1p(omt / t);            /* -Log t, accurate as t -> 1 */
    cplx pw;
    if (q->p_int) {                        /* exact repeated product */
        pw = 1.0;
        for (int i = 0; i < q->p_i; i++) pw *= lw;
    } else {
        pw = nl_pow(lw, q->p);
    }
    cplx v = nl_pow(lt, q->n1) * pw / t;
    if (q->piece == 1) v *= q->c;
    if (q->piece == 2) v *= q->omc;
    return v;
}

static bool nl_tanh_sinh(const nl_q* q, cplx* out) {
    const double khmax = 4.0;
    cplx prev = 0.0, S = 0.0;
    double h = 0.5;
    /* level 0: all nodes at step h; later levels add the odd multiples. */
    cplx sum = M_PI / 2.0 * nl_integrand(q, 0.5, 0.5) * 0.5;
    for (int level = 0; level <= 8; level++) {
        int kstep = level == 0 ? 1 : 2, k0 = level == 0 ? 1 : 1;
        for (int k = k0; (double)k * h <= khmax; k += kstep) {
            double kh = (double)k * h;
            double s2 = M_PI * sinh(kh), e = exp(-s2);
            double t = 1.0 / (1.0 + e), omt = e * t;
            double wgt = M_PI * cosh(kh) * omt * t;   /* (pi/2) cosh / cosh^2(s) / 2 * 2 */
            if (omt == 0.0) break;
            cplx a = nl_integrand(q, t, omt), b = nl_integrand(q, omt, t);
            if (nl_finite(a)) sum += wgt * a;
            if (nl_finite(b)) sum += wgt * b;
        }
        S = h * sum;
        if (level >= 3) {
            double d = cabs(S - prev);
            if (d <= 2e-14 * cabs(S)) { *out = S; return nl_finite(S); }
        }
        prev = S;
        h *= 0.5;
    }
    return false;
}

static bool nl_quad(cplx n, cplx p, cplx z, cplx* out) {
    nl_q q;
    memset(&q, 0, sizeof q);
    q.n1 = n - 1.0; q.p = p; q.z = z;
    q.p_int = cimag(p) == 0.0 && creal(p) == floor(creal(p)) && fabs(creal(p)) < 64;
    q.p_i = q.p_int ? (int)creal(p) : 0;
    cplx S;
    if (cimag(z) == 0.0 && creal(z) > 1.0) {
        q.c = 1.0 / creal(z); q.omc = 1.0 - q.c;
        cplx a, b;
        q.piece = 1; if (!nl_tanh_sinh(&q, &a)) return false;
        q.piece = 2; if (!nl_tanh_sinh(&q, &b)) return false;
        S = a + b;
    } else {
        q.piece = 0; if (!nl_tanh_sinh(&q, &S)) return false;
    }
    *out = S * nl_rgamma(n) * nl_rgamma(p + 1.0);
    return nl_finite(*out);
}

/* ---- the kernel -------------------------------------------------------- */

static bool nl_int_value(double re, double im, int* out) {
    if (im != 0.0 || re != floor(re) || fabs(re) > 1e6) return false;
    *out = (int)re;
    return true;
}

bool sf_nielsen_c(double nr, double ni, double pr, double pi,
                  double zr, double zi, double* rr, double* ri) {
    if (!isfinite(nr) || !isfinite(ni) || !isfinite(pr) || !isfinite(pi) ||
        !isfinite(zr) || !isfinite(zi)) return false;
    int pv = 0, nv = 0;
    bool p_int = nl_int_value(pr, pi, &pv), n_int = nl_int_value(nr, ni, &nv);
    if (!(pr > 0.0) || (p_int && pv < 1)) return false;
    cplx n = CMPLX_(nr, ni), p = CMPLX_(pr, pi), z = CMPLX_(zr, zi), v;
    bool ok;
    if (n_int && nv <= 0) {
        if (-nv > NL_MAXM) return false;
        ok = nl_nonpositive(-nv, p, z, &v);
    } else {
        double az = cabs(z);
        bool z_one = zr == 1.0 && zi == 0.0;
        if (!((az <= 1.0 && !z_one) || (z_one && nr > 0.0) || (p_int && n_int && !z_one)))
            return false;
        if (p_int && pv == 1) {
            /* Li_{n+1}(z), through the same kernels as PolyLog[s, z] */
            double s = nr + 1.0, a, b;
            if (ni != 0.0 || zi != 0.0 || !(fabs(zr) <= 1.0)) {
                if (!sf_polylog_c(s, ni, zr, zi, &a, &b)) return false;
            } else {
                if (!sf_machine_polylog(s, zr, &a)) return false;
                b = 0.0;
            }
            *rr = a; *ri = b;
            return true;
        }
        if (p_int && az <= 0.75)                 ok = nl_stirling(n, pv, z, &v);
        else if (az < 1.0 && (az <= 0.5 || !(nr > 0.0))) ok = nl_miller(n, p, z, &v);
        else if (nr > 0.0)                       ok = nl_quad(n, p, z, &v);
        else                                     ok = false;
    }
    if (!ok) return false;
    double a = creal(v), b = cimag(v);
    /* A real call whose value is real -- an integer p at z <= 1, or any real p
     * at 0 <= z <= 1 -- carries no imaginary part (the quadrature and the
     * principal powers can leave rounding dust there). */
    if (ni == 0.0 && pi == 0.0 && zi == 0.0 &&
        ((p_int && zr <= 1.0) || (zr >= 0.0 && zr <= 1.0)))
        b = 0.0;
    *rr = a; *ri = b;
    return isfinite(a) && isfinite(b);
}
