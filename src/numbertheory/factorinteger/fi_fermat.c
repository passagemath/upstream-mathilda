/*
 * fi_fermat.c — Fermat's method (residue-sieved) and Lehman's method on
 * bignums.  Word-sized inputs use the versions in fi_u64.c.
 *
 * Fermat looks for a with a^2 - n a perfect square, starting at ceil(sqrt n).
 * A square must be a quadratic residue modulo every small m, so for each of
 * nine coprime moduli the residues a mod m for which a^2 - n can be a square
 * are known; only the ~1 in 3000 candidates passing every modulus pay for a
 * bignum square test.  The filter runs 64 values of a at a time: per modulus
 * a table of m 64-bit masks gives, for each residue r of the block's first a,
 * which of the 64 values pass, so a block costs nine loads and ANDs instead
 * of 64 x nine counter steps (~0.5 ns per value of a, against ~10 ns for
 * the per-value counter walk of v0.464).  Budget: `iterations` values of a
 * (default 2^26; 2^28 for Lehman, whose budget also covers the
 * per-multiplier set-up).
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
#define FM_MAXMOD 65

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

/* The filter masks, for every modulus m and every residue nm = kn mod m:
 * FM_MASK[FM_BASE[i] + nm m + r] has bit j set when a = r + j (mod m) can
 * make a^2 - kn a square mod m.  They do not depend on anything else, so
 * they are built once (sum m^2 = 15392 words, 120 KB, ~40k word operations)
 * and a call only picks its nine rows: the per-call set-up is a bignum
 * square root and four word residues, not a table build. */
#define FM_NMASK (64 * 64 + 63 * 63 + 65 * 65 + 11 * 11 + 17 * 17 + 19 * 19 + \
                  23 * 23 + 29 * 29 + 31 * 31)
static uint64_t FM_MASK[FM_NMASK];
static unsigned FM_BASE[NFM];
static int fm_tables_ready = 0;

static void fermat_tables_init(void) {
    if (fm_tables_ready) return;
    unsigned base = 0;
    for (size_t i = 0; i < NFM; i++) {
        unsigned m = FERMAT_MOD[i];
        unsigned char sq[FM_MAXMOD] = { 0 }, ok[FM_MAXMOD];
        for (unsigned x = 0; x < m; x++) sq[(x * x) % m] = 1;
        FM_BASE[i] = base;
        for (unsigned nm = 0; nm < m; nm++) {
            for (unsigned x = 0; x < m; x++) ok[x] = sq[((x * x) % m + m - nm) % m];
            /* Row r = 0 directly, then each next residue shifts one value
             * out and the next one in; branch-free, ok[] being 0/1. */
            uint64_t* w = FM_MASK + base + nm * m;
            uint64_t w0 = 0;
            for (unsigned j = 0, x = 0; j < 64; j++) {
                w0 |= (uint64_t)ok[x] << j;
                x = x + 1 == m ? 0 : x + 1;
            }
            w[0] = w0;
            unsigned in = 64 % m;                  /* (r + 64) mod m, r = 0 */
            for (unsigned r = 1; r < m; r++) {
                w0 = (w0 >> 1) | ((uint64_t)ok[in] << 63);
                w[r] = w0;
                if (++in == m) in = 0;
            }
        }
        base += m * m;
    }
    fm_tables_ready = 1;
}

/* The nine moduli split into two products below 2^32 (unsigned long is 32
 * bits on LLP64 Windows), so a bignum is reduced by all nine with two
 * mpz_fdiv_ui calls. */
#define FM_PROD1 (64UL * 63UL * 65UL * 11UL * 17UL)   /* 49008960 */
#define FM_PROD2 (19UL * 23UL * 29UL * 31UL)          /* 392863 */
#define FM_SPLIT 5                                    /* moduli in FM_PROD1 */

static void fm_residues(unsigned out[NFM], const mpz_t z) {
    unsigned long r1 = mpz_fdiv_ui(z, FM_PROD1), r2 = mpz_fdiv_ui(z, FM_PROD2);
    for (size_t i = 0; i < NFM; i++)
        out[i] = (unsigned)((i < FM_SPLIT ? r1 : r2) % FERMAT_MOD[i]);
}

/* Set-up of one fermat_scan call (a bignum square root and four residues),
 * in units of one value of a: ~1 us against ~0.5 ns per value. */
#define FERMAT_CALL_COST 2048

/* Test the candidate a: on a^2 - kn = b^2, d = a nontrivial gcd(a -+ b, n). */
static int fermat_hit(mpz_t d, const mpz_t n, const mpz_t kn, const mpz_t a,
                      mpz_t b2, mpz_t r, mpz_t g) {
    mpz_mul(b2, a, a);
    mpz_sub(b2, b2, kn);
    if (!mpz_perfect_square_p(b2)) return 0;
    mpz_sqrt(r, b2);
    mpz_sub(b2, a, r);
    mpz_gcd(g, b2, n);
    if (mpz_cmp_ui(g, 1) > 0 && mpz_cmp(g, n) < 0) { mpz_set(d, g); return 1; }
    mpz_add(b2, a, r);
    mpz_gcd(g, b2, n);
    if (mpz_cmp_ui(g, 1) > 0 && mpz_cmp(g, n) < 0) { mpz_set(d, g); return 1; }
    return 0;
}

/* Fermat on kn: test a = ceil(sqrt(kn)) + start, ..., + start + iters - 1,
 * 64 values per block.  On a hit, d = gcd(a - b, n).  A call does no malloc
 * beyond GMP's own -- it is made once per (multiplier, round) by Lehman
 * below, and a TimeConstrained SIGPROF landing here never unwinds out of the
 * allocator. */
static int fermat_scan(mpz_t d, const mpz_t n, const mpz_t kn, uint64_t start,
                       uint64_t iters) {
    const uint64_t* row[NFM];
    unsigned cur[NFM], step[NFM], nm[NFM];
    mpz_t a0, b2, a, r, g;
    mpz_inits(a0, b2, a, r, g, NULL);
    mpz_sqrtrem(a0, b2, kn);
    if (mpz_sgn(b2) != 0) mpz_add_ui(a0, a0, 1);
    if (start) {
        fm_set_u64(r, start);
        mpz_add(a0, a0, r);
    }
    fermat_tables_init();
    fm_residues(nm, kn);
    fm_residues(cur, a0);
    for (size_t i = 0; i < NFM; i++) {
        unsigned m = FERMAT_MOD[i];
        row[i] = FM_MASK + FM_BASE[i] + nm[i] * m;
        step[i] = 64 % m;
    }
    int found = 0;
    for (uint64_t done = 0; done < iters && !found; done += 64) {
        uint64_t bits = ~UINT64_C(0);
        for (size_t i = 0; i < NFM; i++) {
            bits &= row[i][cur[i]];
            unsigned c = cur[i] + step[i];
            cur[i] = c >= FERMAT_MOD[i] ? c - FERMAT_MOD[i] : c;
        }
        if (iters - done < 64) bits &= (UINT64_C(1) << (iters - done)) - 1;
        while (bits && !found) {
            int j = fi_ctz64(bits);
            bits &= bits - 1;
            fm_set_u64(a, done + (uint64_t)j);
            mpz_add(a, a, a0);
            found = fermat_hit(d, n, kn, a, b2, r, g);
        }
    }
    mpz_clears(a0, b2, a, r, g, NULL);
    return found;
}

int fi_split_fermat(mpz_t d, const mpz_t n, FiCtx* c, uint64_t max_iters) {
    (void)c;
    if (!max_iters) max_iters = UINT64_C(1) << 26;
    return fermat_scan(d, n, n, 0, max_iters);
}

/* Lehman schedule constants (see fi_split_lehman).  LEHMAN_K1_WEIGHT: how far
 * ahead of the level c the k = 1 window runs; LEHMAN_GROWTH: level growth per
 * round; LEHMAN_MIN_SHARE: the smallest extension a multiplier is given, in
 * values of a, so that no call is dominated by its set-up. */
#define LEHMAN_K1_WEIGHT 8.0
#define LEHMAN_GROWTH 2.0
#define LEHMAN_MIN_SHARE ((double)FERMAT_CALL_COST)

int fi_split_lehman(mpz_t d, const mpz_t n, FiCtx* c, uint64_t max_iters) {
    (void)c;
    /* 2^28 values of a: with the 64-wide filter an exhausted run costs
     * ~0.1 s, less than 2^24 did with the per-value filter (v0.464). */
    if (!max_iters) max_iters = UINT64_C(1) << 28;
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
     * c: in round c, multiplier k >= 2 is scanned to the fraction min(1, c/k)
     * of its window, and k = 1 to min(1, w c) with w = LEHMAN_K1_WEIGHT.
     * Weighting k = 1 is the prior: p/q near 1 (plain Fermat) is the common
     * structured case, and with w = 8 the k = 1 prefix takes ~5/6 of every
     * round (sum_{k>=2} k^(-3/2) ~ 1.61), so a near-square factor costs
     * ~1.2x what plain Fermat pays for it, not the 2.6x of an unweighted
     * sweep.  A multiplier joins once its share reaches LEHMAN_MIN_SHARE
     * values (or its whole window): k <= K(c) = max((c W_1 / MIN)^(2/3), c),
     * capped at n^(1/3).  The level grows LEHMAN_GROWTH-fold per round; each
     * round resumes every k where the previous one stopped (the previous
     * target is recomputed, so no per-k state is stored).  The sweep order
     * does not depend on the budget, so a larger "Iterations" only ever
     * scans further along the same order (it never finds less).  Level c
     * covers every ratio u/v with uv <= K(c) and delta^2 <= n^(2/3) c / uv,
     * and at c >= n^(1/3) the whole of Lehman's complete search.
     * "Iterations" bounds the total, set-up included. */
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
        double c = LEHMAN_MIN_SHARE / (LEHMAN_K1_WEIGHT * w1), prev_c = 0.0;
        uint64_t prev_K = 0;
        for (;;) {
            double kd = pow(c * w1 / LEHMAN_MIN_SHARE, 2.0 / 3.0);
            if (kd < c) kd = c;
            if (kd > klim) kd = klim;
            if (kd > 9.0e18) kd = 9.0e18;
            uint64_t K = kd < 1.0 ? 1 : (uint64_t)kd;
            for (uint64_t k = 1; k <= K && !found; k++) {
                double wk = w1 / sqrt((double)k) + 2.0;
                double wt = k == 1 ? LEHMAN_K1_WEIGHT : 1.0 / (double)k;
                double f = c * wt;
                uint64_t want = (uint64_t)ceil((f < 1.0 ? f : 1.0) * wk);
                uint64_t have = 0;
                if (k <= prev_K) {
                    double pf = prev_c * wt;
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
            c *= LEHMAN_GROWTH;
        }
    }
    mpz_clears(cb, kn, NULL);
    return found;
}
