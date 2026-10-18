/*
 * fi_fermat.c — Fermat's method (residue-sieved) and Lehman's method on
 * bignums.  Word-sized inputs use the versions in fi_u64.c.
 *
 * Fermat looks for a with a^2 - n a perfect square, starting at ceil(sqrt n).
 * A square must be a quadratic residue modulo every small m, so for each of
 * nine coprime moduli a table marks the residues a mod m for which a^2 - n
 * can be a square; a step then costs nine word increments and lookups, and
 * only the ~1 in 3000 candidates passing every table pay for a bignum
 * square test.  Budget: `iterations` values of a (default 2^26; 2^24 for
 * Lehman, whose budget also covers the per-multiplier set-up).
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

/* uint64_t into / times an mpz, portably: unsigned long is 32 bits on LLP64
 * Windows, so the value is assembled from two 32-bit halves. */
static void fm_set_u64(mpz_t z, uint64_t v) {
    mpz_set_ui(z, (unsigned long)(v >> 32));
    mpz_mul_2exp(z, z, 32);
    mpz_add_ui(z, z, (unsigned long)(v & UINT64_C(0xffffffff)));
}

static void fm_mul_u64(mpz_t z, const mpz_t x, uint64_t v) {
    if (v <= UINT64_C(0xffffffff)) { mpz_mul_ui(z, x, (unsigned long)v); return; }
    mpz_t t;
    mpz_init(t);
    fm_set_u64(t, v);
    mpz_mul(z, x, t);
    mpz_clear(t);
}

/* Fermat on kn: test a = ceil(sqrt(kn)) + start, ..., + start + iters - 1.
 * On a hit, d = gcd(a - b, n).  The residue tables live on the stack (the
 * largest modulus is 65), so a call does no malloc -- it is made once per
 * (multiplier, round) by Lehman below, and a TimeConstrained SIGPROF landing
 * here never unwinds out of the allocator. */
#define FM_MAXMOD 65
#define FERMAT_CALL_COST 32   /* set-up of one call, in units of one step */
/* x^2 mod m and the squares mod m, per modulus: they do not depend on kn, so
 * they are built once and a call only re-derives its kn-shifted filter. */
static unsigned char FM_XX[NFM][FM_MAXMOD];
static unsigned char FM_SQ[NFM][FM_MAXMOD];
static int fm_tables_ready = 0;

static void fermat_tables_init(void) {
    if (fm_tables_ready) return;
    for (size_t i = 0; i < NFM; i++) {
        unsigned m = FERMAT_MOD[i];
        for (unsigned x = 0; x < m; x++) {
            FM_XX[i][x] = (unsigned char)((x * x) % m);
            FM_SQ[i][(x * x) % m] = 1;
        }
    }
    fm_tables_ready = 1;
}

static int fermat_scan(mpz_t d, const mpz_t n, const mpz_t kn, uint64_t start,
                       uint64_t iters) {
    unsigned char ok[NFM][FM_MAXMOD];
    unsigned cur[NFM];
    mpz_t a0, b2, a, r, g;
    mpz_inits(a0, b2, a, r, g, NULL);
    mpz_sqrtrem(a0, b2, kn);
    if (mpz_sgn(b2) != 0) mpz_add_ui(a0, a0, 1);
    if (start) {
        fm_set_u64(r, start);
        mpz_add(a0, a0, r);
    }
    fermat_tables_init();
    for (size_t i = 0; i < NFM; i++) {
        unsigned m = FERMAT_MOD[i];
        unsigned nm = (unsigned)mpz_fdiv_ui(kn, m);
        for (unsigned x = 0; x < m; x++) {
            unsigned t = FM_XX[i][x] + m - nm;     /* x^2 - kn mod m, in [0, 2m) */
            ok[i][x] = FM_SQ[i][t >= m ? t - m : t];
        }
        cur[i] = (unsigned)mpz_fdiv_ui(a0, m);
    }
    int found = 0;
    for (uint64_t it = 0; it < iters && !found; it++) {
        size_t i = 0;
        while (i < NFM && ok[i][cur[i]]) i++;
        if (i == NFM) {
            fm_set_u64(a, it);
            mpz_add(a, a, a0);
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
    mpz_clears(a0, b2, a, r, g, NULL);
    return found;
}

int fi_split_fermat(mpz_t d, const mpz_t n, FiCtx* c, uint64_t max_iters) {
    (void)c;
    if (!max_iters) max_iters = UINT64_C(1) << 26;
    return fermat_scan(d, n, n, 0, max_iters);
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
    /* Then Fermat on 4kn, k = 1, 2, ...  Lehman's window for k is
     * W_k = n^(1/6) / (4 sqrt k) values of a past ceil(sqrt(4kn)).
     *
     * A factor pair p, q with v p - u q = delta (p/q close to u/v) is found
     * at multiplier k = uv, at offset ~ delta^2 / (4 sqrt(kn)) -- that is, at
     * the FRACTION delta^2 / n^(2/3) of W_k, independent of k.  Scanning each
     * k's whole window in turn (as Lehman's complete algorithm does) spends
     * the entire budget on k = 1 once n passes ~45 digits (W_1 > 2^24), so
     * the "small rational ratio" factors at k >= 2 were never reached.
     *
     * Instead the (k, fraction) plane is swept in rounds of a growing level
     * c: in round c, multiplier k is scanned to the fraction min(1, c/k) of
     * its window, for every k whose share is at least one value of a,
     *     k <= K(c) = min((c W_1)^(2/3), n^(1/3)).
     * The cost of a round is ~ 4 W_1 sqrt(c), so c *= 4 doubles it; each
     * round resumes every k where the previous round stopped (the previous
     * target is recomputed, so no per-k state is stored).  Level c covers
     * every ratio u/v with uv <= K(c) and delta^2 <= n^(2/3) c / uv: the
     * near-1 ratios of plain Fermat at small c * W_1, small-rational ratios
     * after a handful of rounds, and at c >= n^(1/3) the whole of Lehman's
     * complete search.  "Iterations" bounds the total, set-up included. */
    if (!found) {
        mpz_t s6;
        mpz_init(s6);
        mpz_root(s6, n, 6);
        double w1 = mpz_get_d(s6) / 4.0;
        mpz_clear(s6);
        if (w1 < 1.0) w1 = 1.0;
        double klim = mpz_get_d(cb);                 /* Lehman's k bound */
        if (klim < 1.0) klim = 1.0;
        uint64_t spent = 0;
        double c = 1.0 / w1, prev_c = 0.0;
        uint64_t prev_K = 0;
        for (;;) {
            double kd = pow(c * w1, 2.0 / 3.0);
            if (kd > klim) kd = klim;
            if (kd > 9.0e18) kd = 9.0e18;
            uint64_t K = kd < 1.0 ? 1 : (uint64_t)kd;
            for (uint64_t k = 1; k <= K && !found; k++) {
                double wk = w1 / sqrt((double)k) + 2.0;
                double f = c / (double)k;
                uint64_t want = (uint64_t)ceil((f < 1.0 ? f : 1.0) * wk);
                uint64_t have = 0;
                if (k <= prev_K) {
                    double pf = prev_c / (double)k;
                    have = (uint64_t)ceil((pf < 1.0 ? pf : 1.0) * wk);
                }
                if (want < 1) want = 1;
                if (want <= have) continue;
                if (spent + FERMAT_CALL_COST >= max_iters) break;
                uint64_t len = want - have;
                if (len > max_iters - spent - FERMAT_CALL_COST)
                    len = max_iters - spent - FERMAT_CALL_COST;
                mpz_mul_ui(kn, n, 4);
                fm_mul_u64(kn, kn, k);
                found = fermat_scan(d, n, kn, have, len);
                spent += len + FERMAT_CALL_COST;
            }
            if (found || spent + FERMAT_CALL_COST >= max_iters) break;
            if ((double)K >= klim && c >= klim) break;   /* Lehman complete */
            prev_c = c;
            prev_K = K;
            c *= 4.0;
        }
    }
    mpz_clears(cb, kn, NULL);
    return found;
}
