/* sf_machine_dirichlet.c -- machine-precision Dirichlet L-functions.
 * See sf_machine_dirichlet.h for the method and the decline contract.
 * No <complex.h>: complex values travel as (re, im) double pairs.
 */
#include "sf_machine_dirichlet.h"
#include "sf_machine_zeta.h"   /* sf_zeta_r / sf_zeta_c (conductor 1) */
#include "sf_machine_gamma.h"  /* sfg_loggamma_c                      */
#include <math.h>
#include <stdlib.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#ifndef M_LN2
#define M_LN2 0.693147180559945309417
#endif

/* Euler-Maclaurin order: J Bernoulli terms, B_{2j}/(2j)!. */
#define DL_EM_J 10
static const double DL_BC[DL_EM_J] = {
    1.0 / 12.0,
    -1.0 / 720.0,
    1.0 / 30240.0,
    -1.0 / 1209600.0,
    1.0 / 47900160.0,
    -691.0 / 1307674368000.0,
    1.0 / 74724249600.0,
    -3617.0 / 10670622842880000.0,
    43867.0 / 5109094217170944000.0,
    -174611.0 / 802857662698291200000.0,
};

/* Ceiling on direct-sum terms before declining to arbitrary precision. */
#define DL_MAX_TERMS 20000000.0

/* Tail start M (in units of f). The first omitted Euler-Maclaurin term is
 * bounded by
 *     2 |(s)_{2J+1}| / (2 Pi)^(2J+2) * M^-(2J+1+sigma),
 * with |(s)_{2J+1}| <= Gamma(|s| + 2J + 1)/Gamma(|s|); M is the smallest
 * integer (>= 2) putting that below 1e-17. sigma >= 1/2 here (the functional
 * equation covers the left half-plane). */
static double dl_tail_start(double sabs, double sigma) {
    double a = sabs < 0.5 ? 0.5 : sabs;
    double lpoch = lgamma(a + 2.0 * DL_EM_J + 1.0) - lgamma(a);
    double num = M_LN2 + lpoch - (2.0 * DL_EM_J + 2.0) * log(2.0 * M_PI) + 39.1439465808987777;
    double M = ceil(exp(num / (2.0 * DL_EM_J + 1.0 + sigma)));
    return M < 2.0 ? 2.0 : M;
}

/* (e^u - 1)/u for complex u, accurate for small |u|. */
static void dl_phi(double ur, double ui, double* pr, double* pi) {
    double m = fabs(ur) + fabs(ui);
    if (m < 1e-5) {
        /* 1 + u/2 + u^2/6 + u^3/24 */
        double u2r = ur * ur - ui * ui, u2i = 2.0 * ur * ui;
        double u3r = u2r * ur - u2i * ui, u3i = u2r * ui + u2i * ur;
        *pr = 1.0 + ur / 2.0 + u2r / 6.0 + u3r / 24.0;
        *pi = ui / 2.0 + u2i / 6.0 + u3i / 24.0;
        return;
    }
    /* expm1(u) = (expm1(ur) cos ui - 2 sin^2(ui/2)) + I e^ur sin ui */
    double sh = sin(0.5 * ui);
    double er = expm1(ur) * cos(ui) - 2.0 * sh * sh;
    double ei = exp(ur) * sin(ui);
    double d = ur * ur + ui * ui;
    *pr = (er * ur + ei * ui) / d;
    *pi = (ei * ur - er * ui) / d;
}

static double dl_phi_r(double u) {
    if (fabs(u) < 1e-5) return 1.0 + u / 2.0 + u * u / 6.0 + u * u * u / 24.0;
    return expm1(u) / u;
}

/* L(s, chi*) for a primitive NON-principal chi* (f >= 2), Re s >= 1/2, by
 * Euler-Maclaurin. `conj` uses conj(chi*). Complex arithmetic throughout. */
static bool dl_em_c(const DchTable* T, bool conj, double sr, double si,
                    double* vr, double* vi) {
    uint64_t f = T->f;
    const double* cr = T->cre;
    const double* ci = T->cim;
    double sg = conj ? -1.0 : 1.0;
    double M = dl_tail_start(hypot(sr, si), sr);
    if (M * (double)f > DL_MAX_TERMS) return false;
    uint64_t Mi = (uint64_t)M;

    /* Direct sum over n = m f + a, 1 <= n <= M f. */
    double Dr = 0.0, Di = 0.0;
    for (uint64_t a = 1; a < f; a++) {
        double xr = cr[a], xi = sg * ci[a];
        if (xr == 0.0 && xi == 0.0) continue;
        double Sr = 0.0, Si = 0.0;
        for (uint64_t m = 0; m < Mi; m++) {
            double ln = log((double)(m * f + a));
            double mag = exp(-sr * ln), ang = -si * ln;
            Sr += mag * cos(ang); Si += mag * sin(ang);
        }
        Dr += xr * Sr - xi * Si;
        Di += xr * Si + xi * Sr;
    }

    /* Pochhammer coefficients C_j = c_j (s)_{2j-1}. */
    double Cr[DL_EM_J], Ci[DL_EM_J];
    {
        double pr = sr, pi = si;                       /* (s)_1 */
        for (int j = 0; j < DL_EM_J; j++) {
            Cr[j] = DL_BC[j] * pr; Ci[j] = DL_BC[j] * pi;
            /* (s)_{2j+1} = (s)_{2j-1} (s + 2j - 1)(s + 2j) */
            double a1r = sr + 2.0 * j + 1.0, a2r = sr + 2.0 * j + 2.0;
            double tr = a1r * a2r - si * si, ti = si * (a1r + a2r);
            double nr = pr * tr - pi * ti, ni = pr * ti + pi * tr;
            pr = nr; pi = ni;
        }
    }

    /* Tail: sum_a chi(a) [ -M^(1-s) l_a phi(u_a)  (via M^(1-s) below)
     *                     + x^-s (1/2 + (1/x) sum_j C_j x^-(2j-2)) ]. */
    double lM = log(M);
    double Ir = 0.0, Ii = 0.0;      /* integral part, without the M^(1-s) factor */
    double Br = 0.0, Bi = 0.0;      /* boundary + Bernoulli part                 */
    double fM = (double)f * M;
    for (uint64_t a = 1; a < f; a++) {
        double xr = cr[a], xi = sg * ci[a];
        if (xr == 0.0 && xi == 0.0) continue;
        double la = log1p((double)a / fM);
        double ur = (1.0 - sr) * la, ui = -si * la;
        double phr, phi_;
        dl_phi(ur, ui, &phr, &phi_);
        double tr = -la * phr, ti = -la * phi_;
        Ir += xr * tr - xi * ti;
        Ii += xr * ti + xi * tr;

        double x = M + (double)a / (double)f;
        double lx = lM + la;
        double mag = exp(-sr * lx), ang = -si * lx;
        double Pr = mag * cos(ang), Pi = mag * sin(ang);
        double y = 1.0 / (x * x);
        double hr = Cr[DL_EM_J - 1], hi = Ci[DL_EM_J - 1];
        for (int j = DL_EM_J - 2; j >= 0; j--) { hr = hr * y + Cr[j]; hi = hi * y + Ci[j]; }
        hr = 0.5 + hr / x; hi = hi / x;
        double qr = Pr * hr - Pi * hi, qi = Pr * hi + Pi * hr;
        Br += xr * qr - xi * qi;
        Bi += xr * qi + xi * qr;
    }
    /* M^(1-s) */
    {
        double mag = exp((1.0 - sr) * lM), ang = -si * lM;
        double mr = mag * cos(ang), mi = mag * sin(ang);
        double tr = mr * Ir - mi * Ii, ti = mr * Ii + mi * Ir;
        Br += tr; Bi += ti;
    }
    /* times f^-s */
    {
        double lf = log((double)f);
        double mag = exp(-sr * lf), ang = -si * lf;
        double fr = mag * cos(ang), fi = mag * sin(ang);
        double tr = fr * Br - fi * Bi, ti = fr * Bi + fi * Br;
        Dr += tr; Di += ti;
    }
    *vr = Dr; *vi = Di;
    return isfinite(Dr) && isfinite(Di);
}

/* The same for a REAL character and real s >= 1/2, in real arithmetic. */
static bool dl_em_r(const DchTable* T, double s, double* v) {
    uint64_t f = T->f;
    const double* cr = T->cre;
    double M = dl_tail_start(fabs(s), s);
    if (M * (double)f > DL_MAX_TERMS) return false;
    uint64_t Mi = (uint64_t)M;
    double D = 0.0;
    for (uint64_t a = 1; a < f; a++) {
        double x = cr[a];
        if (x == 0.0) continue;
        double S = 0.0;
        for (uint64_t m = 0; m < Mi; m++) S += exp(-s * log((double)(m * f + a)));
        D += x * S;
    }
    double C[DL_EM_J];
    {
        double p = s;
        for (int j = 0; j < DL_EM_J; j++) {
            C[j] = DL_BC[j] * p;
            p *= (s + 2.0 * j + 1.0) * (s + 2.0 * j + 2.0);
        }
    }
    double lM = log(M), fM = (double)f * M;
    double I = 0.0, B = 0.0;
    for (uint64_t a = 1; a < f; a++) {
        double xc = cr[a];
        if (xc == 0.0) continue;
        double la = log1p((double)a / fM);
        I += xc * (-la * dl_phi_r((1.0 - s) * la));
        double x = M + (double)a / (double)f;
        double P = exp(-s * (lM + la));
        double y = 1.0 / (x * x);
        double h = C[DL_EM_J - 1];
        for (int j = DL_EM_J - 2; j >= 0; j--) h = h * y + C[j];
        B += xc * P * (0.5 + h / x);
    }
    B += exp((1.0 - s) * lM) * I;
    D += exp(-s * log((double)f)) * B;
    *v = D;
    return isfinite(D);
}

/* L(s, chi*) for primitive chi* (any conductor), complex s. */
static bool dl_prim_c(const DchTable* T, bool conj, double sr, double si,
                      double* vr, double* vi) {
    if (T->f == 1) return sf_zeta_c(sr, si, vr, vi);
    if (sr >= 0.5) return dl_em_c(T, conj, sr, si, vr, vi);
    /* Functional equation. Trivial zeros: (s + d)/2 a non-positive integer. */
    double d = (double)T->odd;
    double hr = 0.5 * (sr + d), hi = 0.5 * si;
    if (hi == 0.0 && hr <= 0.0 && hr == floor(hr)) { *vr = 0.0; *vi = 0.0; return true; }
    double Lr, Li;
    if (!dl_em_c(T, !conj, 1.0 - sr, -si, &Lr, &Li)) return false;
    double g1r, g1i, g2r, g2i;
    if (si == 0.0) {
        /* Real s: the Gamma ratio is real. The complex LogGamma declines on the
         * negative real axis, where (s + d)/2 lies, so take log|Gamma| and carry
         * the sign of Gamma(h) as a phase of Pi. */
        g1r = lgamma(0.5 * (1.0 - sr + d)); g1i = 0.0;   /* argument > 1/4 */
        g2r = lgamma(hr);
        g2i = (hr < 0.0 && fmod(ceil(-hr), 2.0) != 0.0) ? M_PI : 0.0;
    } else {
        if (!sfg_loggamma_c(0.5 * (1.0 - sr + d), -0.5 * si, &g1r, &g1i)) return false;
        if (!sfg_loggamma_c(hr, hi, &g2r, &g2i)) return false;
    }
    double lfp = log((double)T->f / M_PI);
    /* log of (f/Pi)^(1/2-s) Gamma(..)/Gamma(..) */
    double er = (0.5 - sr) * lfp + g1r - g2r;
    double ei = -si * lfp + g1i - g2i;
    if (er > 700.0) return false;
    double mag = exp(er);
    double wr = mag * cos(ei), wi = mag * sin(ei);
    /* Root number: eps(conj chi) = conj(eps(chi)), since tau(conj chi) =
     * chi(-1) conj(tau(chi)) and chi(-1) = (-1)^d. */
    double epr = T->eps_re, epi = conj ? -T->eps_im : T->eps_im;
    double tr = epr * wr - epi * wi, ti = epr * wi + epi * wr;
    *vr = tr * Lr - ti * Li;
    *vi = tr * Li + ti * Lr;
    return isfinite(*vr) && isfinite(*vi);
}

static bool dl_prim_r(const DchTable* T, double s, double* v) {
    if (T->f == 1) return sf_zeta_r(s, v);
    if (s >= 0.5) return dl_em_r(T, s, v);
    double d = (double)T->odd;
    double h = 0.5 * (s + d);
    if (h <= 0.0 && h == floor(h)) { *v = 0.0; return true; }
    double L;
    if (!dl_em_r(T, 1.0 - s, &L)) return false;
    double g1 = 0.5 * (1.0 - s + d);           /* > 1/4 */
    double lr = lgamma(g1) - lgamma(h) + (0.5 - s) * log((double)T->f / M_PI);
    if (lr > 700.0) return false;
    /* sign of Gamma(h): positive for h > 0; (-1)^ceil(-h) for h < 0. */
    double sgn = 1.0;
    if (h < 0.0) { double c = ceil(-h); if (fmod(c, 2.0) != 0.0) sgn = -1.0; }
    *v = sgn * exp(lr) * L;                      /* eps = 1 for a real character */
    return isfinite(*v);
}

/* chi*(p) as machine re/im from its exponent over the order. */
static void dl_unit(const DchTable* T, uint64_t e, double* cr, double* ci) {
    uint64_t O = T->order;
    e %= O;
    if ((4 * e) % O == 0) {
        static const double ur[4] = { 1.0, 0.0, -1.0, 0.0 };
        static const double ui[4] = { 0.0, 1.0, 0.0, -1.0 };
        int q = (int)(4 * e / O);
        *cr = ur[q]; *ci = ui[q];
        return;
    }
    double num = (double)e;
    if (2 * e > O) num -= (double)O;
    double th = 2.0 * M_PI * num / (double)O;
    *cr = cos(th); *ci = sin(th);
}

bool sf_dirichlet_l_tab_c(const DchTable* T, double sr, double si, double* vr, double* vi) {
    if (!T || !T->estar) return false;
    if (!isfinite(sr) || !isfinite(si)) return false;
    double Lr, Li;
    if (!dl_prim_c(T, false, sr, si, &Lr, &Li)) return false;
    for (int i = 0; i < T->neul; i++) {
        double cr, ci;
        dl_unit(T, T->eul_e[i], &cr, &ci);
        double lp = log((double)T->eul_p[i]);
        double mag = exp(-sr * lp), ang = -si * lp;
        double pr = mag * cos(ang), pi = mag * sin(ang);
        double fr = 1.0 - (cr * pr - ci * pi), fi = -(cr * pi + ci * pr);
        double nr = Lr * fr - Li * fi, ni = Lr * fi + Li * fr;
        Lr = nr; Li = ni;
    }
    *vr = Lr; *vi = Li;
    return isfinite(Lr) && isfinite(Li);
}

bool sf_dirichlet_l_tab_r(const DchTable* T, double s, double* v) {
    if (!T || !T->estar || !T->real) return false;
    if (!isfinite(s)) return false;
    double L;
    if (!dl_prim_r(T, s, &L)) return false;
    for (int i = 0; i < T->neul; i++) {
        double c = (T->eul_e[i] % T->order == 0) ? 1.0 : -1.0;
        L *= 1.0 - c * exp(-s * log((double)T->eul_p[i]));
    }
    *v = L;
    return isfinite(L);
}

/* (k, j) as doubles -> table, or NULL. */
static const DchTable* dl_table_d(double k, double j, DchTable** owned) {
    *owned = NULL;
    if (!(k >= 1.0 && j >= 1.0 && k < 9.0e15 && j < 9.0e15)) return NULL;
    if (k != floor(k) || j != floor(j)) return NULL;
    return dch_get((int64_t)k, (int64_t)j, true, owned, NULL);
}

bool sf_dirichlet_l_c(double k, double j, double sr, double si, double* vr, double* vi) {
    DchTable* owned;
    const DchTable* T = dl_table_d(k, j, &owned);
    if (!T) return false;
    bool ok;
    if (si == 0.0 && T->real) { ok = sf_dirichlet_l_tab_r(T, sr, vr); *vi = 0.0; }
    else ok = sf_dirichlet_l_tab_c(T, sr, si, vr, vi);
    dch_free(owned);
    return ok;
}

bool sf_dirichlet_l_r(double k, double j, double s, double* v) {
    DchTable* owned;
    const DchTable* T = dl_table_d(k, j, &owned);
    if (!T) return false;
    bool ok = sf_dirichlet_l_tab_r(T, s, v);
    dch_free(owned);
    return ok;
}
