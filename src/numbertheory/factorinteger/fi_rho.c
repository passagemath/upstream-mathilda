/*
 * fi_rho.c — Pollard rho (Brent's variant) on bignums.
 *
 * Two-limb inputs (below 2^126) go to the Montgomery engine in fi_u128.c;
 * larger ones run here on GMP.  Brent's cycle finding with the products of
 * |x - y| accumulated 256 at a time between gcds, and Brent's backtrack when
 * a batch overshoots to gcd = n.  Each start uses f(y) = y^2 + c with a
 * fresh c; the budget is `iterations` steps per start (default 2^24, which
 * finds factors up to ~2^40 with high probability) over `starts` starts.
 */
#include "fi_internal.h"

static int rho_mpz(mpz_t d, const mpz_t n, unsigned long c, uint64_t max_iters) {
    mpz_t x, y, ys, q, g, t;
    mpz_inits(x, y, ys, q, g, t, NULL);
    mpz_set_ui(y, 2);
    mpz_set_ui(q, 1);
    mpz_set_ui(g, 1);
    const uint64_t m = 256;
    uint64_t done = 0;
#define STEP(v) do { mpz_mul(t, v, v); mpz_add_ui(t, t, c); mpz_tdiv_r(v, t, n); } while (0)
    for (uint64_t r = 1; mpz_cmp_ui(g, 1) == 0 && done < max_iters; r <<= 1) {
        mpz_set(x, y);
        for (uint64_t i = 0; i < r; i++) STEP(y);
        done += r;
        for (uint64_t k = 0; k < r && mpz_cmp_ui(g, 1) == 0; k += m) {
            mpz_set(ys, y);
            uint64_t lim = (r - k < m) ? r - k : m;
            for (uint64_t i = 0; i < lim; i++) {
                STEP(y);
                mpz_sub(t, x, y);
                mpz_mul(q, q, t);
                mpz_tdiv_r(q, q, n);
            }
            done += lim;
            mpz_gcd(g, q, n);
        }
    }
    if (mpz_cmp(g, n) == 0 || mpz_sgn(g) == 0) {
        mpz_set_ui(g, 1);
        for (uint64_t i = 0; i < 2 * m && mpz_cmp_ui(g, 1) == 0; i++) {
            STEP(ys);
            mpz_sub(t, x, ys);
            mpz_gcd(g, t, n);
        }
    }
#undef STEP
    int ok = mpz_cmp_ui(g, 1) > 0 && mpz_cmp(g, n) < 0;
    if (ok) mpz_set(d, g);
    mpz_clears(x, y, ys, q, g, t, NULL);
    return ok;
}

int fi_split_rho(mpz_t d, const mpz_t n, FiCtx* c, uint64_t max_iters, uint64_t starts) {
    (void)c;
    if (!max_iters) max_iters = UINT64_C(1) << 24;
    if (!starts) starts = 4;
    size_t bits = mpz_sizeinbase(n, 2);
    for (uint64_t s = 0; s < starts; s++) {
        unsigned long cc = (unsigned long)(1 + s);
        if (bits > 64 && bits <= 126) {
            if (fi_u128_rho(d, n, cc, max_iters)) return 1;
        } else if (rho_mpz(d, n, cc, max_iters)) {
            return 1;
        }
    }
    return 0;
}
