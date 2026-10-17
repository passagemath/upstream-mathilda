/*
 * fi_fermat.c — Fermat's method (residue-sieved) and Lehman's method on
 * bignums.  Word-sized inputs use the versions in fi_u64.c.
 *
 * Fermat looks for a with a^2 - n a perfect square, starting at ceil(sqrt n).
 * A square must be a quadratic residue modulo every small m, so for each of
 * nine coprime moduli a table marks the residues a mod m for which a^2 - n
 * can be a square; a step then costs nine word increments and lookups, and
 * only the ~1 in 3000 candidates passing every table pay for a bignum
 * square test.  Budget: `iterations` values of a (default 2^26).
 *
 * Lehman (1974) is Fermat applied to 4kn for k = 1, 2, ... after trial
 * division to n^(1/3); it is O(n^(1/3)) and complete when run to the end.
 * On bignums it runs with a budget, as a multiplier-Fermat that also finds
 * p, q with p/q close to a small rational.
 */
#include "fi_internal.h"

#include <math.h>
#include <stdlib.h>

static const unsigned FERMAT_MOD[] = { 64, 63, 65, 11, 17, 19, 23, 29, 31 };
#define NFM (sizeof FERMAT_MOD / sizeof FERMAT_MOD[0])

/* Fermat on kn, scanning a from ceil(sqrt(kn)) for at most `iters` steps.
 * On a hit, d = gcd(a - b, n). */
static int fermat_scan(mpz_t d, const mpz_t n, const mpz_t kn, uint64_t iters) {
    unsigned char* ok[NFM];
    unsigned cur[NFM];
    mpz_t a0, b2, a, r, g;
    mpz_inits(a0, b2, a, r, g, NULL);
    mpz_sqrtrem(a0, b2, kn);
    if (mpz_sgn(b2) != 0) mpz_add_ui(a0, a0, 1);
    for (size_t i = 0; i < NFM; i++) {
        unsigned m = FERMAT_MOD[i];
        unsigned char* sq = calloc(m, 1);
        for (unsigned x = 0; x < m; x++) sq[(x * x) % m] = 1;
        ok[i] = calloc(m, 1);
        unsigned nm = (unsigned)mpz_fdiv_ui(kn, m);
        for (unsigned x = 0; x < m; x++) ok[i][x] = sq[((x * x) % m + m - nm) % m];
        free(sq);
        cur[i] = (unsigned)mpz_fdiv_ui(a0, m);
    }
    int found = 0;
    for (uint64_t it = 0; it < iters && !found; it++) {
        size_t i = 0;
        while (i < NFM && ok[i][cur[i]]) i++;
        if (i == NFM) {
            mpz_add_ui(a, a0, (unsigned long)it);
            mpz_mul(b2, a, a);
            mpz_sub(b2, b2, kn);
            if (mpz_perfect_square_p(b2)) {
                mpz_sqrt(r, b2);
                mpz_sub(b2, a, r);
                mpz_gcd(g, b2, n);
                if (mpz_cmp_ui(g, 1) > 0 && mpz_cmp(g, n) < 0) { mpz_set(d, g); found = 1; }
                else {
                    mpz_add(b2, a, r);
                    mpz_gcd(g, b2, n);
                    if (mpz_cmp_ui(g, 1) > 0 && mpz_cmp(g, n) < 0) { mpz_set(d, g); found = 1; }
                }
            }
        }
        for (i = 0; i < NFM; i++) if (++cur[i] == FERMAT_MOD[i]) cur[i] = 0;
    }
    for (size_t i = 0; i < NFM; i++) free(ok[i]);
    mpz_clears(a0, b2, a, r, g, NULL);
    return found;
}

int fi_split_fermat(mpz_t d, const mpz_t n, FiCtx* c, uint64_t max_iters) {
    (void)c;
    if (!max_iters) max_iters = UINT64_C(1) << 26;
    return fermat_scan(d, n, n, max_iters);
}

int fi_split_lehman(mpz_t d, const mpz_t n, FiCtx* c, uint64_t max_iters) {
    (void)c;
    if (!max_iters) max_iters = UINT64_C(1) << 24;
    /* Trial division to min(n^(1/3), 2^20). */
    mpz_t cb, kn;
    mpz_inits(cb, kn, NULL);
    mpz_root(cb, n, 3);
    uint64_t tb = mpz_cmp_ui(cb, FI_SMALL_PRIME_LIMIT) < 0 ? mpz_get_ui(cb) : FI_SMALL_PRIME_LIMIT;
    size_t np;
    const uint32_t* P = fi_small_primes(&np);
    int found = 0;
    for (size_t i = 1; i < np && P[i] <= tb; i++) {
        if (mpz_divisible_ui_p(n, P[i])) { mpz_set_ui(d, P[i]); found = 1; break; }
    }
    /* Then Fermat on 4kn with a short window per k; the window for k is
     * n^(1/6) / (4 sqrt k) + 1 values of a, so the budget is spread over the
     * k that Lehman's bound needs. */
    if (!found) {
        mpz_t s6;
        mpz_init(s6);
        mpz_root(s6, n, 6);
        double sixth = mpz_get_d(s6);
        mpz_clear(s6);
        uint64_t spent = 0;
        for (uint64_t k = 1; spent < max_iters && !found; k++) {
            uint64_t win = (uint64_t)(sixth / (4.0 * sqrt((double)k))) + 2;
            mpz_mul_ui(kn, n, 4 * (unsigned long)k);
            found = fermat_scan(d, n, kn, win);
            spent += win + 16;
        }
    }
    mpz_clears(cb, kn, NULL);
    return found;
}
