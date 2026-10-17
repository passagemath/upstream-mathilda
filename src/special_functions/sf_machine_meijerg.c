/* sf_machine_meijerg.c -- machine-precision Meijer G-function.
 *
 * See sf_machine_meijerg.h for the method (Slater's theorem with regularized
 * pFq terms) and the decline contract. The per-term pFq values come from the
 * shared hypergeometric kernels (sf_machine_pfq / sf_machine_pfq_c), so every
 * Slater term is evaluated by exactly the code HypergeometricPFQ uses.
 */
#include "sf_machine_meijerg.h"
#include "sf_machine.h"        /* sf_machine_pfq (real)              */
#include "sf_machine_hyper.h"  /* sf_machine_pfq_c (complex)         */
#include "sf_machine_gamma.h"  /* sfg_gamma_c, sfg_loggamma_c        */
#include "flint_num_bridge.h"  /* flint_mg_eval_d: the arb fallback */
#include <math.h>
#include <string.h>
#include <stdlib.h>
#ifdef MATHILDA_THREADS
#include <pthread.h>
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#ifndef M_LN2
#define M_LN2 0.693147180559945309417
#endif

/* ------------------------------------------------------------------ */
/*  (re, im) arithmetic                                                */
/* ------------------------------------------------------------------ */

static inline void c_mul(double ar, double ai, double br, double bi,
                         double* rr, double* ri) {
    double r = ar * br - ai * bi;
    double i = ar * bi + ai * br;
    *rr = r; *ri = i;
}

static inline void c_div(double ar, double ai, double br, double bi,
                         double* rr, double* ri) {
    /* Smith's algorithm: no overflow for large |b|. */
    if (fabs(br) >= fabs(bi)) {
        double t = bi / br, d = br + bi * t;
        *rr = (ar + ai * t) / d;
        *ri = (ai - ar * t) / d;
    } else {
        double t = br / bi, d = br * t + bi;
        *rr = (ar * t + ai) / d;
        *ri = (ai * t - ar) / d;
    }
}

static inline double c_abs(double r, double i) { return hypot(r, i); }

static inline void c_exp(double zr, double zi, double* rr, double* ri) {
    double m = exp(zr);
    if (zi == 0.0) { *rr = m; *ri = 0.0; return; }
    *rr = m * cos(zi);
    *ri = m * sin(zi);
}

/* Principal log; a signed-zero imaginary part is read as +0, so the negative
 * real axis has arg +pi (Mathematica's branch, continuous from above). */
static inline void c_log(double zr, double zi, double* lr, double* li) {
    if (zi == 0.0) zi = 0.0;               /* -0.0 -> +0.0 */
    *lr = log(c_abs(zr, zi));
    *li = atan2(zi, zr);
}

/* sin(pi x) and cos(pi x) with the argument reduced EXACTLY to |f| <= 1/4 first,
 * so sin(pi k) is exactly 0 and the pi/sin(pi d) factors are as accurate as the
 * difference d itself. */
static void sincospi(double x, double* s, double* c) {
    double n = nearbyint(2.0 * x);
    double f = x - 0.5 * n;                /* exact: |f| <= 1/4 */
    double sf = sin(M_PI * f), cf = cos(M_PI * f);
    double k4 = fmod(n, 4.0);
    if (k4 < 0) k4 += 4.0;
    switch ((int)k4) {
        case 0:  *s = sf;  *c = cf;  break;
        case 1:  *s = cf;  *c = -sf; break;
        case 2:  *s = -sf; *c = -cf; break;
        default: *s = -cf; *c = sf;  break;
    }
}

/* sin(pi (x + i y)). */
static void c_sinpi(double x, double y, double* rr, double* ri) {
    double s, c;
    sincospi(x, &s, &c);
    if (y == 0.0) { *rr = s; *ri = 0.0; return; }
    *rr = s * cosh(M_PI * y);
    *ri = c * sinh(M_PI * y);
}

static inline bool is_nonpos_int(double r, double i) {
    return i == 0.0 && r <= 0.0 && r == nearbyint(r);
}

static inline bool is_int(double r, double i) {
    return i == 0.0 && r == nearbyint(r);
}

/* ------------------------------------------------------------------ */
/*  The coefficient of one Slater term                                 */
/* ------------------------------------------------------------------ */

/* A running product kept as a finite complex mantissa plus a power-of-two
 * exponent, so a coefficient that is the ratio of two huge Gamma values stays
 * representable until the very end. */
typedef struct { double r, i; long e; bool zero; } Prod;

static void prod_init(Prod* P) { P->r = 1.0; P->i = 0.0; P->e = 0; P->zero = false; }

static void prod_norm(Prod* P) {
    double a = c_abs(P->r, P->i);
    if (a == 0.0 || !isfinite(a)) return;
    int ex;
    (void)frexp(a, &ex);
    if (ex > 64 || ex < -64) {
        P->r = ldexp(P->r, -ex); P->i = ldexp(P->i, -ex); P->e += ex;
    }
}

static void prod_mul(Prod* P, double r, double i) {
    c_mul(P->r, P->i, r, i, &P->r, &P->i);
    prod_norm(P);
}

static void prod_div(Prod* P, double r, double i) {
    c_div(P->r, P->i, r, i, &P->r, &P->i);
    prod_norm(P);
}

/* Multiply by Gamma(z) (inv = false) or 1/Gamma(z) (inv = true). A pole of the
 * numerator Gamma is an inconsistent parameter set: returns false. A pole under
 * 1/Gamma zeroes the product. Large arguments, where Gamma itself overflows,
 * go through LogGamma. */
static bool prod_gamma(Prod* P, double zr, double zi, bool inv) {
    if (is_nonpos_int(zr, zi)) {
        if (inv) { P->zero = true; return true; }
        return false;
    }
    double gr, gi;
    if (c_abs(zr, zi) < 150.0 && sfg_gamma_c(zr, zi, &gr, &gi) &&
        isfinite(gr) && isfinite(gi) && (gr != 0.0 || gi != 0.0)) {
        if (inv) prod_div(P, gr, gi); else prod_mul(P, gr, gi);
        return true;
    }
    double lr, li;
    if (!sfg_loggamma_c(zr, zi, &lr, &li) || !isfinite(lr) || !isfinite(li)) return false;
    if (inv) { lr = -lr; li = -li; }
    /* exp(lr + i li) = 2^(lr / ln 2) * e^(i li): split off the binary exponent. */
    double l2 = lr / M_LN2;
    double ie = floor(l2);
    double mag = exp((l2 - ie) * M_LN2);
    double er = mag * cos(li), ei = mag * sin(li);
    c_mul(P->r, P->i, er, ei, &P->r, &P->i);
    P->e += (long)ie;
    prod_norm(P);
    return true;
}

static bool prod_value(const Prod* P, double* rr, double* ri) {
    if (P->zero) { *rr = 0.0; *ri = 0.0; return true; }
    if (P->e > 2000 || P->e < -2000) return false;
    *rr = ldexp(P->r, (int)P->e);
    *ri = ldexp(P->i, (int)P->e);
    return isfinite(*rr) && isfinite(*ri);
}

/* ------------------------------------------------------------------ */
/*  One Slater expansion                                               */
/* ------------------------------------------------------------------ */

/* Build S for G^{m,n}_{p,q}(A; B; .) expanded about 0 (the "z-series").
 * Returns false when the parameters are inconsistent (a numerator Gamma pole:
 * a_j - b_h a positive integer for j <= n, h <= m). */
static bool build_series(SfMgSeries* S,
                         const double* Ar, const double* Ai, size_t n, size_t p,
                         const double* Br, const double* Bi, size_t m, size_t q) {
    memset(S, 0, sizeof *S);
    S->nt = m;
    S->pu = p;
    S->ql = q ? q - 1 : 0;
    S->sgn = (((long)p - (long)m - (long)n) % 2 != 0) ? -1.0 : 1.0;
    S->real = true;
    for (size_t j = 0; j < p; j++) if (Ai[j] != 0.0) S->real = false;
    for (size_t j = 0; j < q; j++) if (Bi[j] != 0.0) S->real = false;

    /* The log case: two of b_1..b_m an integer apart. */
    for (size_t h = 0; h < m; h++)
        for (size_t j = h + 1; j < m; j++)
            if (is_int(Br[j] - Br[h], Bi[j] - Bi[h])) S->log_case = true;

    for (size_t h = 0; h < m; h++) {
        S->er[h] = Br[h]; S->ei[h] = Bi[h];
        /* Upper parameters 1 + b_h - a_k, lower 1 + b_h - b_j (j != h). */
        size_t nl = 0;
        for (size_t k = 0; k < p; k++) {
            S->ur[h][k] = 1.0 + Br[h] - Ar[k];
            S->ui[h][k] = Bi[h] - Ai[k];
        }
        for (size_t j = 0; j < q; j++) {
            if (j == h) continue;
            S->lr[h][nl] = 1.0 + Br[h] - Br[j];
            S->li[h][nl] = Bi[h] - Bi[j];
            nl++;
        }
        Prod C; prod_init(&C);
        /* Numerator Gammas Gamma(1 + b_h - a_j), j <= n. */
        for (size_t j = 0; j < n; j++)
            if (!prod_gamma(&C, 1.0 + Br[h] - Ar[j], Bi[h] - Ai[j], false)) return false;
        if (S->log_case) { S->cr[h] = 0.0; S->ci[h] = 0.0; continue; }
        /* 1/Gamma(a_j - b_h), j > n. */
        for (size_t j = n; j < p; j++)
            if (!prod_gamma(&C, Ar[j] - Br[h], Ai[j] - Bi[h], true)) return false;
        /* A lower parameter 1 + b_h - b_j (necessarily j > m: j <= m would be
         * the log case) that is a non-positive integer needs the regularized
         * form; the common case does not, and keeps the textbook coefficient. */
        bool shift = false;
        for (size_t l = 0; l < nl; l++)
            if (is_nonpos_int(S->lr[h][l], S->li[h][l])) shift = true;
        if (!shift) {
            /* Gamma(b_j - b_h), j <= m, j != h; 1/Gamma(1 + b_h - b_j), j > m:
             * the latter pairs with the term's own lower parameter. */
            for (size_t j = 0; j < q; j++) {
                if (j == h) continue;
                bool ok = (j < m) ? prod_gamma(&C, Br[j] - Br[h], Bi[j] - Bi[h], false)
                                  : prod_gamma(&C, 1.0 + Br[h] - Br[j], Bi[h] - Bi[j], true);
                if (!ok) return false;
            }
            S->kx[h] = 0;
            if (!prod_value(&C, &S->cr[h], &S->ci[h])) return false;
            if (S->ci[h] != 0.0) S->real = false;
            continue;
        }
        /* Regularized route. pi / sin(pi (b_j - b_h)), j <= m, j != h, is the
         * Gamma pair Gamma(b_j - b_h) Gamma(1 + b_h - b_j); the term is then
         * C * R(u; l; x) with R the REGULARIZED pFq, every lower parameter's
         * 1/Gamma applied below. */
        for (size_t j = 0; j < m; j++) {
            if (j == h) continue;
            double sr, si;
            c_sinpi(Br[j] - Br[h], Bi[j] - Bi[h], &sr, &si);
            prod_mul(&C, M_PI, 0.0);
            prod_div(&C, sr, si);
        }
        /* A non-positive integer lower parameter -N is shifted out exactly:
         *   R(u; -N, l'; x) = prod_k (u_k)_{N+1} x^{N+1} R(u+N+1; N+2, l'+N+1; x)
         * repeated until no such parameter remains. */
        int kx = 0;
        for (int guard = 0; guard < 64; guard++) {
            size_t hit = (size_t)-1;
            for (size_t l = 0; l < nl; l++)
                if (is_nonpos_int(S->lr[h][l], S->li[h][l])) { hit = l; break; }
            if (hit == (size_t)-1) break;
            long N = (long)(-S->lr[h][hit]);
            for (size_t k = 0; k < p; k++) {
                for (long t = 0; t <= N; t++) prod_mul(&C, S->ur[h][k] + (double)t, S->ui[h][k]);
                S->ur[h][k] += (double)(N + 1);
            }
            for (size_t l = 0; l < nl; l++) {
                if (l == hit) S->lr[h][l] = (double)(N + 2);
                else          S->lr[h][l] += (double)(N + 1);
            }
            kx += (int)(N + 1);
        }
        S->kx[h] = kx;
        /* 1/Gamma of every (now pole-free) lower parameter. */
        for (size_t l = 0; l < nl; l++)
            if (!prod_gamma(&C, S->lr[h][l], S->li[h][l], true)) return false;
        if (!prod_value(&C, &S->cr[h], &S->ci[h])) return false;
        if (S->ci[h] != 0.0) S->real = false;
    }
    S->valid = true;
    return true;
}

bool sf_mg_plan(const double* ar, const double* ai, size_t n, size_t p,
                const double* br, const double* bi, size_t m, size_t q,
                double r, SfMgPlan* P) {
    memset(P, 0, sizeof *P);
    if (p > SF_MG_MAXP || q > SF_MG_MAXP || n > p || m > q) return false;
    if (!isfinite(r) || r == 0.0) return false;
    P->n = n; P->p = p; P->m = m; P->q = q; P->r = r;
    P->real_params = true;
    for (size_t j = 0; j < p; j++) {
        P->ar[j] = ar[j]; P->ai[j] = ai ? ai[j] : 0.0;
        if (!isfinite(P->ar[j]) || !isfinite(P->ai[j])) return false;
        if (P->ai[j] != 0.0) P->real_params = false;
    }
    for (size_t j = 0; j < q; j++) {
        P->br[j] = br[j]; P->bi[j] = bi ? bi[j] : 0.0;
        if (!isfinite(P->br[j]) || !isfinite(P->bi[j])) return false;
        if (P->bi[j] != 0.0) P->real_params = false;
    }
    if (p <= q && !build_series(&P->zs, P->ar, P->ai, n, p, P->br, P->bi, m, q))
        P->zs.valid = false;
    if (p >= q) {
        /* Dual: G^{n,m}_{q,p}(1 - b; 1 - a; 1/w). */
        double Ar[SF_MG_MAXP], Ai[SF_MG_MAXP], Br[SF_MG_MAXP], Bi[SF_MG_MAXP];
        for (size_t j = 0; j < q; j++) { Ar[j] = 1.0 - P->br[j]; Ai[j] = -P->bi[j]; }
        for (size_t j = 0; j < p; j++) { Br[j] = 1.0 - P->ar[j]; Bi[j] = -P->ai[j]; }
        if (!build_series(&P->ds, Ar, Ai, m, q, Br, Bi, n, p)) P->ds.valid = false;
    }
    P->ok = true;
    return true;
}

/* ------------------------------------------------------------------ */
/*  Evaluation                                                         */
/* ------------------------------------------------------------------ */

/* pFq(u; l; x) for one term, real lane. */
static bool term_pfq_r(const SfMgSeries* S, size_t h, double x, double* F) {
    if (S->pu == 0 && S->ql == 0) { *F = exp(x); return isfinite(*F); }
    return sf_machine_pfq(S->ur[h], S->pu, S->lr[h], S->ql, x, F);
}

/* pFq(u; l; x) for one term, complex lane. */
static bool term_pfq_c(const SfMgSeries* S, size_t h, double xr, double xi,
                       double* Fr, double* Fi) {
    if (S->pu == 0 && S->ql == 0) { c_exp(xr, xi, Fr, Fi); return isfinite(*Fr) && isfinite(*Fi); }
    if (S->pu == 1 && S->ql == 0) {
        /* 1F0(u;;x) = (1 - x)^(-u), principal branch: the closed form is valid
         * on the whole plane bar the cut, so |x| = 1 needs no series. */
        if (xr == 1.0 && xi == 0.0) return false;
        double lr, li, er, ei;
        c_log(1.0 - xr, -xi, &lr, &li);
        c_mul(-S->ur[h][0], -S->ui[h][0], lr, li, &er, &ei);
        c_exp(er, ei, Fr, Fi);
        return isfinite(*Fr) && isfinite(*Fi);
    }
    return sf_machine_pfq_c(S->ur[h], S->ui[h], S->pu, S->lr[h], S->li[h], S->ql,
                            xr, xi, Fr, Fi);
}

/* Sum one expansion at log w = (lr, li). */
static bool eval_series(const SfMgSeries* S, double lr, double li, bool real_lane,
                        double* rr, double* ri) {
    if (!S->valid || S->log_case) return false;
    if (S->nt == 0) { *rr = 0.0; *ri = 0.0; return true; }
    double wr, wi;
    c_exp(lr, li, &wr, &wi);
    double xr = S->sgn * wr, xi = S->sgn * wi;
    double sr = 0.0, si = 0.0, maxabs = 0.0;
    for (size_t h = 0; h < S->nt; h++) {
        if (S->cr[h] == 0.0 && S->ci[h] == 0.0) continue;
        double tr, ti;
        if (real_lane) {
            double F;
            if (!term_pfq_r(S, h, xr, &F)) return false;
            /* w^{e} for real w > 0 and real e. */
            double pw = exp(S->er[h] * lr);
            tr = S->cr[h] * pw * F;
            for (int k = 0; k < S->kx[h]; k++) tr *= xr;
            ti = 0.0;
        } else {
            double Fr, Fi, er, ei, pr, pi;
            if (!term_pfq_c(S, h, xr, xi, &Fr, &Fi)) return false;
            c_mul(S->er[h], S->ei[h], lr, li, &er, &ei);
            c_exp(er, ei, &pr, &pi);
            c_mul(S->cr[h], S->ci[h], pr, pi, &tr, &ti);
            c_mul(tr, ti, Fr, Fi, &tr, &ti);
            for (int k = 0; k < S->kx[h]; k++) c_mul(tr, ti, xr, xi, &tr, &ti);
        }
        if (!isfinite(tr) || !isfinite(ti)) return false;
        double a = c_abs(tr, ti);
        if (a > maxabs) maxabs = a;
        sr += tr; si += ti;
    }
    double s = c_abs(sr, si);
    if (maxabs > 0.0) {
        if (s == 0.0) return false;                     /* total cancellation */
        if (log2(maxabs / s) > SF_MG_MAX_LOST_BITS) return false;
    }
    *rr = sr; *ri = si;
    return true;
}

bool sf_mg_eval_machine(const SfMgPlan* P, double zr, double zi,
                        double* rr, double* ri, bool* is_c) {
    if (!P->ok || !isfinite(zr) || !isfinite(zi)) return false;
    if (zr == 0.0 && zi == 0.0) return false;           /* the exact path owns z = 0 */
    double lr, li;
    c_log(zr, zi, &lr, &li);
    lr /= P->r; li /= P->r;
    const SfMgSeries* S;
    double sl = 1.0;                                    /* +1: w, -1: 1/w */
    if (P->p < P->q) S = &P->zs;
    else if (P->p > P->q) { S = &P->ds; sl = -1.0; }
    else {
        if (lr < 0.0) S = &P->zs;
        else if (lr > 0.0) { S = &P->ds; sl = -1.0; }
        else {
            /* |w| = 1 exactly: only the shapes whose term function continues
             * analytically past the unit circle (1F0, 2F1). */
            if (P->p > 2) return false;
            S = &P->zs;
        }
    }
    bool real_lane = P->real_params && S->real && li == 0.0 && zi == 0.0 && zr > 0.0;
    if (real_lane && eval_series(S, sl * lr, 0.0, true, rr, ri)) {
        *ri = 0.0; *is_c = false; return true;
    }
    if (!eval_series(S, sl * lr, sl * li, false, rr, ri)) return false;
    *is_c = true;
    if (real_lane) { *ri = 0.0; *is_c = false; }       /* real by construction */
    return true;
}

bool sf_mg_eval(const SfMgPlan* P, double zr, double zi,
                double* rr, double* ri, bool* is_c) {
    if (sf_mg_eval_machine(P, zr, zi, rr, ri, is_c)) return true;
    if (!P->ok || !isfinite(zr) || !isfinite(zi)) return false;
    if (zr == 0.0 && zi == 0.0) return false;
    return flint_mg_eval_d(P->ar, P->ai, P->n, P->p, P->br, P->bi, P->m, P->q,
                           P->r, zr, zi, rr, ri, is_c);
}

/* ------------------------------------------------------------------ */
/*  Flattened n-ary layout (ndkernels.c, Compile[])                    */
/* ------------------------------------------------------------------ */

double sf_mg_shape_code(size_t n, size_t p, size_t m) {
    return (double)n + 64.0 * (double)p + 4096.0 * (double)m;
}

/* A small, append-only cache of plans keyed by the parameter block, so a
 * compiled loop over z builds the plan once. Entries are never replaced or
 * freed while the process runs (plans are immutable once published), which is
 * what makes the lookup safe from the threaded fused-map loops without holding
 * the lock during evaluation; once the cache is full, further parameter sets
 * build their plan per call on the stack. */
#define MG_MEMO_SLOTS 16
#define MG_KEY_MAX (2 + 4 * SF_MG_MAXP)
typedef struct { size_t len; double key[MG_KEY_MAX]; SfMgPlan* plan; } MgMemo;
static MgMemo g_memo[MG_MEMO_SLOTS];
static size_t g_memo_n = 0;
#ifdef MATHILDA_THREADS
static pthread_mutex_t g_memo_lock = PTHREAD_MUTEX_INITIALIZER;
#define MG_LOCK()   pthread_mutex_lock(&g_memo_lock)
#define MG_UNLOCK() pthread_mutex_unlock(&g_memo_lock)
#else
#define MG_LOCK()   ((void)0)
#define MG_UNLOCK() ((void)0)
#endif

static const SfMgPlan* memo_find(const double* key, size_t len) {
    const SfMgPlan* hit = NULL;
    MG_LOCK();
    for (size_t s = 0; s < g_memo_n; s++)
        if (g_memo[s].len == len && memcmp(g_memo[s].key, key, len * sizeof(double)) == 0) {
            hit = g_memo[s].plan; break;
        }
    MG_UNLOCK();
    return hit;
}

static void memo_add(const double* key, size_t len, const SfMgPlan* P) {
    MG_LOCK();
    if (g_memo_n < MG_MEMO_SLOTS) {
        SfMgPlan* c = malloc(sizeof *c);
        if (c) {
            *c = *P;
            g_memo[g_memo_n].len = len;
            memcpy(g_memo[g_memo_n].key, key, len * sizeof(double));
            g_memo[g_memo_n].plan = c;
            g_memo_n++;
        }
    }
    MG_UNLOCK();
}

bool sf_mg_eval_flat(const double* re, const double* im, size_t count,
                     double* rr, double* ri) {
    if (count < 3 || im[0] != 0.0 || im[1] != 0.0) return false;
    double code = re[0];
    if (!(code >= 0.0) || code != nearbyint(code)) return false;
    long c = (long)code;
    size_t n = (size_t)(c % 64), p = (size_t)((c / 64) % 64), m = (size_t)(c / 4096);
    if (count < p + 3) return false;
    size_t q = count - p - 3;
    if (p > SF_MG_MAXP || q > SF_MG_MAXP) return false;
    /* Key: code, r, then (re, im) of every parameter. */
    double key[MG_KEY_MAX];
    size_t len = 0;
    key[len++] = code; key[len++] = re[1];
    for (size_t j = 0; j < p + q; j++) { key[len++] = re[2 + j]; key[len++] = im[2 + j]; }
    const SfMgPlan* P = memo_find(key, len);
    SfMgPlan local;
    if (!P) {
        if (!sf_mg_plan(re + 2, im + 2, n, p, re + 2 + p, im + 2 + p, m, q, re[1], &local))
            return false;
        memo_add(key, len, &local);
        P = &local;
    }
    bool is_c;
    return sf_mg_eval(P, re[count - 1], im[count - 1], rr, ri, &is_c);
}
