/*
 * fi_rbd_lattice.c — the Howgrave-Graham approximate-divisor detector of
 * lattice-augmented rational base descent ([LRBD] = lattice_rbd.tex,
 * Theorem 2 and Appendix A).
 *
 * Given N and an integer Q within X of a multiple v M of an unknown divisor
 * M >= N^beta, the polynomials
 *
 *     g_i(x) = N^(t-i) (Q + x)^i          0 <= i <= t,
 *     h_j(x) = x^j (Q + x)^t              1 <= j < d - t,
 *
 * all vanish modulo M^t at x0 = vM - Q.  Their coefficient vectors, column k
 * scaled by X^k, span a d-dimensional lattice of determinant
 * N^(t(t+1)/2) X^(d(d-1)/2) (the basis is triangular).  LLL finds a vector of
 * norm at most 2^((d-1)/4) det^(1/d); once that is below M^t / sqrt(d), its
 * polynomial vanishes at x0 over the integers (Howgrave-Graham's lemma), and
 * x0 is read off as an integer root.  fi_rbd_hg_params solves that
 * inequality for the tolerance or for the dimension.
 *
 * Integer roots are found p-adically, which needs nothing beyond GMP: take a
 * small prime p not dividing the leading coefficient with h square-free mod
 * p, find the roots mod p by evaluation, Hensel-lift each (Newton, doubling
 * the precision) past 2 * bound, and keep the symmetric residues that are
 * genuine integer roots.
 */
#include "fi_rbd.h"

#include <math.h>
#include <stdlib.h>

/* ------------------------------------------------------------------ */
/* Parameters                                                          */
/* ------------------------------------------------------------------ */

/* log2 of the largest tolerance the LLL guarantee covers at (t, d):
 *   (d-1)/4 + [t(t+1)/2 L + d(d-1)/2 lx] / d  <  t beta L - log2(d)/2. */
static double hg_max_log2X(int t, int d, double L, double beta) {
    double rhs = t * beta * L - 0.5 * log2((double)d) - (d - 1) / 4.0 -
                 (double)t * (t + 1) * L / (2.0 * d);
    return rhs * 2.0 / (d - 1);
}

int fi_rbd_hg_params(int* t, int* d, double* log2X_out, double log2N, double beta,
                     double log2X) {
    if (!(beta > 0) || beta > 1 || log2N < 8) return 0;
    if (log2X >= 0) {
        for (int dd = 2; dd <= 48; dd++) {
            int best_t = 0;
            double best = -1e300;
            for (int tt = 1; tt < dd; tt++) {
                double m = hg_max_log2X(tt, dd, log2N, beta);
                if (m > best) { best = m; best_t = tt; }
            }
            if (best >= log2X) {
                *t = best_t;
                *d = dd;
                if (log2X_out) *log2X_out = best;
                return 1;
            }
        }
        return 0;
    }
    int dd = *d > 1 ? *d : 12;
    int best_t = 0;
    double best = -1e300;
    for (int tt = 1; tt < dd; tt++) {
        double m = hg_max_log2X(tt, dd, log2N, beta);
        if (m > best) { best = m; best_t = tt; }
    }
    if (best < 1) return 0;
    *t = best_t;
    *d = dd;
    if (log2X_out) *log2X_out = best;
    return 1;
}

/* ------------------------------------------------------------------ */
/* Integer roots                                                       */
/* ------------------------------------------------------------------ */

/* Small primes for the p-adic root search, from 1009 up. */
static int next_prime_from(uint32_t* p) {
    for (uint32_t c = *p + 1; c < 65536; c++) {
        int ok = c > 1;
        for (uint32_t q = 2; q * q <= c && ok; q++) if (c % q == 0) ok = 0;
        if (ok) { *p = c; return 1; }
    }
    return 0;
}

/* Polynomials over F_p as uint32 arrays (degree tracked separately). */
static int fp_deg(const uint32_t* f, int n) {
    while (n >= 0 && f[n] == 0) n--;
    return n;
}

static uint32_t fp_inv(uint32_t a, uint32_t p) {
    /* a^(p-2) mod p */
    uint64_t r = 1, b = a % p;
    uint32_t e = p - 2;
    while (e) {
        if (e & 1) r = r * b % p;
        b = b * b % p;
        e >>= 1;
    }
    return (uint32_t)r;
}

/* f <- f mod g in place (deg g >= 0); returns the new degree of f. */
static int fp_rem(uint32_t* f, int df, const uint32_t* g, int dg, uint32_t p) {
    uint32_t inv = fp_inv(g[dg], p);
    while (df >= dg && df >= 0) {
        uint64_t c = (uint64_t)f[df] * inv % p;
        if (c) {
            int sh = df - dg;
            for (int i = 0; i <= dg; i++)
                f[sh + i] = (uint32_t)((f[sh + i] + (uint64_t)(p - c) * g[i]) % p);
        }
        df = fp_deg(f, df - 1);
    }
    return df;
}

/* Degree of gcd(f, f') over F_p (f of degree n, leading coefficient nonzero
 * mod p); 0 means f is square-free mod p. */
static int fp_gcd_deriv_deg(const uint32_t* f, int n, uint32_t p, uint32_t* wa, uint32_t* wb) {
    for (int i = 0; i <= n; i++) wa[i] = f[i];
    for (int i = 1; i <= n; i++) wb[i - 1] = (uint32_t)((uint64_t)f[i] * (uint32_t)i % p);
    int da = n, db = fp_deg(wb, n - 1);
    if (db < 0) return n;                   /* f' = 0 mod p */
    uint32_t* A = wa;
    uint32_t* Bp = wb;
    while (db >= 0) {
        da = fp_rem(A, da, Bp, db, p);
        uint32_t* tmp = A; A = Bp; Bp = tmp;
        int td = da; da = db; db = td;
    }
    return da;
}

/* h(x) mod m (m > 0), result in [0, m). */
static void eval_mod(mpz_t r, const mpz_t* h, int deg, const mpz_t x, const mpz_t m) {
    mpz_set(r, h[deg]);
    for (int i = deg - 1; i >= 0; i--) {
        mpz_mul(r, r, x);
        mpz_add(r, r, h[i]);
        mpz_mod(r, r, m);
    }
    mpz_mod(r, r, m);
}

/* h'(x) mod m. */
static void eval_deriv_mod(mpz_t r, const mpz_t* h, int deg, const mpz_t x, const mpz_t m,
                           mpz_t t) {
    mpz_mul_ui(r, h[deg], (unsigned long)deg);
    for (int i = deg - 1; i >= 1; i--) {
        mpz_mul(r, r, x);
        mpz_mul_ui(t, h[i], (unsigned long)i);
        mpz_add(r, r, t);
        mpz_mod(r, r, m);
    }
    mpz_mod(r, r, m);
}

/* Exact h(x) over Z, == 0 ? */
static int is_root(const mpz_t* h, int deg, const mpz_t x, mpz_t r) {
    mpz_set(r, h[deg]);
    for (int i = deg - 1; i >= 0; i--) {
        mpz_mul(r, r, x);
        mpz_add(r, r, h[i]);
    }
    return mpz_sgn(r) == 0;
}

/* ---- Square-free part over Z (primitive PRS) ----------------------- */

/* Degree of f[0..n] (-1 for the zero polynomial). */
static int zdeg(mpz_t* f, int n) {
    while (n >= 0 && mpz_sgn(f[n]) == 0) n--;
    return n;
}

/* f <- f / content(f), with a positive leading coefficient. */
static void zprim(mpz_t* f, int n, mpz_t c) {
    mpz_set_ui(c, 0);
    for (int i = 0; i <= n; i++) mpz_gcd(c, c, f[i]);
    if (mpz_sgn(c) == 0) return;
    if (mpz_sgn(f[n]) < 0) mpz_neg(c, c);
    if (mpz_cmp_ui(c, 1) != 0)
        for (int i = 0; i <= n; i++) mpz_divexact(f[i], f[i], c);
}

/* a <- prem(a, b) = lc(b)^(da-db+1) a mod b; returns the new degree. */
static int zprem(mpz_t* a, int da, mpz_t* b, int db, mpz_t t) {
    while (da >= db && da >= 0) {
        /* a <- lc(b) a - lc(a) x^(da-db) b */
        mpz_set(t, a[da]);
        for (int i = 0; i <= da; i++) mpz_mul(a[i], a[i], b[db]);
        for (int i = 0; i <= db; i++) mpz_submul(a[i + da - db], t, b[i]);
        da = zdeg(a, da - 1);
    }
    return da;
}

/* The square-free part s = h / gcd(h, h') of h (degree n >= 1), primitive,
 * into s[0..]; returns its degree.  The roots of s are those of h, each
 * simple, so the p-adic lifting below applies.  Small degrees only (the
 * detector's polynomials have d <= 48 coefficients). */
static int zsquarefree(mpz_t* s, const mpz_t* h, int n) {
    size_t sz = (size_t)(n + 1);
    mpz_t* a = malloc(3 * sz * sizeof(mpz_t));
    if (!a) { for (int i = 0; i <= n; i++) mpz_set(s[i], h[i]); return n; }
    mpz_t* b = a + sz;
    mpz_t* q = b + sz;
    for (size_t i = 0; i < 3 * sz; i++) mpz_init(a[i]);
    mpz_t t, c;
    mpz_inits(t, c, NULL);
    /* g = gcd(h, h') by the primitive PRS. */
    for (int i = 0; i <= n; i++) mpz_set(a[i], h[i]);
    for (int i = 1; i <= n; i++) mpz_mul_ui(b[i - 1], h[i], (unsigned long)i);
    int da = n, db = n - 1;
    zprim(a, da, c);
    zprim(b, db, c);
    while (db > 0) {
        da = zprem(a, da, b, db, t);
        if (da < 0) break;
        zprim(a, da, c);
        /* swap (a, da) <-> (b, db) */
        for (int i = 0; i <= (da > db ? da : db); i++) mpz_swap(a[i], b[i]);
        int td = da; da = db; db = td;
    }
    /* The gcd is b when the last remainder vanished (da < 0) and db > 0;
     * a non-zero constant remainder means gcd = 1. */
    int dg = (da < 0) ? db : 0;
    int ds = n;
    for (int i = 0; i <= n; i++) mpz_set(s[i], h[i]);
    if (dg > 0) {
        /* s = h / g exactly: g is primitive and divides h, so (Gauss) the
         * quotient is integral; long division with exact steps. */
        mpz_t* g = b;
        for (int i = 0; i <= n; i++) mpz_set(a[i], h[i]);
        zprim(a, n, c);
        ds = n - dg;
        for (int i = 0; i <= ds; i++) mpz_set_ui(q[i], 0);
        for (int k = n; k >= dg; k--) {
            if (mpz_sgn(a[k]) == 0) continue;
            mpz_divexact(t, a[k], g[dg]);
            mpz_set(q[k - dg], t);
            for (int i = 0; i <= dg; i++) mpz_submul(a[i + k - dg], t, g[i]);
        }
        for (int i = 0; i <= ds; i++) mpz_set(s[i], q[i]);
    }
    zprim(s, ds, c);
    mpz_clears(t, c, NULL);
    for (size_t i = 0; i < 3 * sz; i++) mpz_clear(a[i]);
    free(a);
    return ds;
}

int fi_rbd_int_roots(mpz_t* roots, const mpz_t* h_in, int deg, const mpz_t bound) {
    while (deg > 0 && mpz_sgn(h_in[deg]) == 0) deg--;
    if (deg < 1) return 0;
    int nroots = 0;
    /* Strip the factor x^v: 0 is then a root. */
    int v = 0;
    while (v < deg && mpz_sgn(h_in[v]) == 0) v++;
    if (v > 0) mpz_set_ui(roots[nroots++], 0);
    if (deg - v < 1) return nroots;
    /* Work on the square-free part, so every root is simple. */
    mpz_t* hs = malloc((size_t)(deg - v + 1) * sizeof(mpz_t));
    if (!hs) return nroots;
    for (int i = 0; i <= deg - v; i++) mpz_init(hs[i]);
    int dfull = deg - v;
    deg = zsquarefree(hs, h_in + v, dfull);
    const mpz_t* h = (const mpz_t*)hs;
    if (deg < 1) {
        for (int i = 0; i <= dfull; i++) mpz_clear(hs[i]);
        free(hs);
        return nroots;
    }

    uint32_t* f = malloc((size_t)(deg + 1) * 3 * sizeof(uint32_t));
    if (!f) {
        for (int i = 0; i <= dfull; i++) mpz_clear(hs[i]);
        free(hs);
        return nroots;
    }
    uint32_t* wa = f + deg + 1;
    uint32_t* wb = wa + deg + 1;
    mpz_t m, m2, r, dr, t, x, lim;
    mpz_inits(m, m2, r, dr, t, x, lim, NULL);
    mpz_mul_2exp(lim, bound, 1);
    mpz_add_ui(lim, lim, 1);                /* lift past 2 bound + 1 */

    uint32_t p = 1008;
    int tries = 0;
    while (tries < 48 && next_prime_from(&p)) {
        tries++;
        if (mpz_fdiv_ui(h[deg], p) == 0) continue;
        for (int i = 0; i <= deg; i++) f[i] = (uint32_t)mpz_fdiv_ui(h[i], p);
        if (fp_gcd_deriv_deg(f, deg, p, wa, wb) != 0) continue;
        /* Square-free mod p: every root mod p is simple and lifts uniquely. */
        for (uint32_t r0 = 0; r0 < p && nroots < deg + v; r0++) {
            uint64_t acc = f[deg];
            for (int i = deg - 1; i >= 0; i--) acc = (acc * r0 + f[i]) % p;
            if (acc) continue;
            mpz_set_ui(x, r0);
            mpz_set_ui(m, p);
            while (mpz_cmp(m, lim) <= 0) {
                mpz_mul(m2, m, m);
                eval_mod(r, h, deg, x, m2);
                eval_deriv_mod(dr, h, deg, x, m2, t);
                if (!mpz_invert(dr, dr, m2)) break;
                mpz_mul(r, r, dr);
                mpz_sub(x, x, r);
                mpz_mod(x, x, m2);
                mpz_swap(m, m2);
            }
            /* Symmetric residue. */
            mpz_tdiv_q_2exp(t, m, 1);
            if (mpz_cmp(x, t) > 0) mpz_sub(x, x, m);
            if (mpz_cmpabs(x, bound) > 0) continue;
            if (!is_root(h, deg, x, r)) continue;
            mpz_set(roots[nroots++], x);
        }
        break;
    }
    mpz_clears(m, m2, r, dr, t, x, lim, NULL);
    free(f);
    for (int i = 0; i <= dfull; i++) mpz_clear(hs[i]);
    free(hs);
    return nroots;
}

/* ------------------------------------------------------------------ */
/* Detector                                                            */
/* ------------------------------------------------------------------ */

static int hg_detect(mpz_t g, mpz_t x0, const mpz_t N, const mpz_t Q, const mpz_t X, int t,
                     int d, int native) {
    if (t < 1 || d < t + 1 || d > 64 || mpz_sgn(X) <= 0) return 0;
    size_t dd = (size_t)d;
    mpz_t* M = malloc(dd * dd * sizeof(mpz_t));
    mpz_t* Qp = malloc((size_t)(t + 1) * sizeof(mpz_t));
    mpz_t* Np = malloc((size_t)(t + 1) * sizeof(mpz_t));
    mpz_t* Xp = malloc(dd * sizeof(mpz_t));
    mpz_t* h = malloc(dd * sizeof(mpz_t));
    mpz_t* roots = malloc(dd * sizeof(mpz_t));
    if (!M || !Qp || !Np || !Xp || !h || !roots) {
        free(M); free(Qp); free(Np); free(Xp); free(h); free(roots);
        return 0;
    }
    for (size_t i = 0; i < dd * dd; i++) mpz_init(M[i]);
    for (int i = 0; i <= t; i++) { mpz_init(Qp[i]); mpz_init(Np[i]); }
    for (size_t i = 0; i < dd; i++) { mpz_init(Xp[i]); mpz_init(h[i]); mpz_init(roots[i]); }
    mpz_t c, y;
    mpz_inits(c, y, NULL);

    mpz_set_ui(Qp[0], 1);
    mpz_set_ui(Np[0], 1);
    for (int i = 1; i <= t; i++) { mpz_mul(Qp[i], Qp[i - 1], Q); mpz_mul(Np[i], Np[i - 1], N); }
    mpz_set_ui(Xp[0], 1);
    for (size_t i = 1; i < dd; i++) mpz_mul(Xp[i], Xp[i - 1], X);

    /* Rows i = 0..t: N^(t-i) (Q + x)^i; coefficient of x^k is
     * N^(t-i) C(i, k) Q^(i-k). */
    for (int i = 0; i <= t; i++)
        for (int k = 0; k <= i; k++) {
            mpz_bin_uiui(c, (unsigned long)i, (unsigned long)k);
            mpz_mul(c, c, Qp[i - k]);
            mpz_mul(c, c, Np[t - i]);
            mpz_mul(M[(size_t)i * dd + (size_t)k], c, Xp[k]);
        }
    /* Rows t + j, j = 1..d-t-1: x^j (Q + x)^t. */
    for (int j = 1; j < d - t; j++)
        for (int k = 0; k <= t; k++) {
            mpz_bin_uiui(c, (unsigned long)t, (unsigned long)k);
            mpz_mul(c, c, Qp[t - k]);
            mpz_mul(M[(size_t)(t + j) * dd + (size_t)(k + j)], c, Xp[k + j]);
        }

    int found = 0;
    if ((native ? fi_rbd_lll_native(M, d, d) : fi_rbd_lll(M, d, d)) == 0) {
        /* The short vectors come first; the target polynomial is almost
         * always the first, but a few more cost little. */
        int rows_to_try = d < 4 ? d : 4;
        for (int r = 0; r < rows_to_try && !found; r++) {
            int deg = -1;
            for (int k = 0; k < d; k++) {
                mpz_divexact(h[k], M[(size_t)r * dd + (size_t)k], Xp[k]);
                if (mpz_sgn(h[k])) deg = k;
            }
            if (deg < 1) continue;
            int nr = fi_rbd_int_roots(roots, (const mpz_t*)h, deg, X);
            for (int i = 0; i < nr && !found; i++) {
                mpz_add(y, Q, roots[i]);
                if (mpz_sgn(y) <= 0) continue;
                mpz_gcd(c, y, N);
                if (mpz_cmp_ui(c, 1) > 0 && mpz_cmp(c, N) < 0) {
                    mpz_set(g, c);
                    mpz_set(x0, roots[i]);
                    found = 1;
                }
            }
        }
    }

    mpz_clears(c, y, NULL);
    for (size_t i = 0; i < dd * dd; i++) mpz_clear(M[i]);
    for (int i = 0; i <= t; i++) { mpz_clear(Qp[i]); mpz_clear(Np[i]); }
    for (size_t i = 0; i < dd; i++) { mpz_clear(Xp[i]); mpz_clear(h[i]); mpz_clear(roots[i]); }
    free(M); free(Qp); free(Np); free(Xp); free(h); free(roots);
    return found;
}

int fi_rbd_hg_detect(mpz_t g, mpz_t x0, const mpz_t N, const mpz_t Q, const mpz_t X, int t,
                     int d) {
    return hg_detect(g, x0, N, Q, X, t, d, 0);
}

int fi_rbd_hg_detect_native(mpz_t g, mpz_t x0, const mpz_t N, const mpz_t Q, const mpz_t X,
                            int t, int d) {
    return hg_detect(g, x0, N, Q, X, t, d, 1);
}
