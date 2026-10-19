/* primecount.c -- pi(x), the prime-counting function, with several selectable
 * algorithms (see primecount.h).
 *
 * Shared infrastructure:
 *   - a small sieve table of the primes below 10^6 (also serves Prime[n]);
 *   - gen_primes(bound): the primes up to a modest bound (<= sqrt(x));
 *   - PiTable: a bit-packed sieve up to B with block prefix counts giving
 *     O(1) pi(v) lookups for v <= B;
 *   - phi(v, b): the partial sieve function #{n <= v : n coprime to the first b
 *     primes}, via a PhiTiny wheel base plus Lehmer's p_b^2 > v prune (which
 *     turns the otherwise-exponential Legendre recursion into a bounded one,
 *     using PiTable for the pruned leaves).
 *
 * Methods:
 *   Sieve     -- segmented sieve of Eratosthenes (ground truth, small x).
 *   Lucy      -- Lucy_Hedgehog DP, O(x^3/4) time, O(sqrt x) space.
 *   Legendre  -- pi(x) = phi(x, pi(sqrt x)) + pi(sqrt x) - 1.
 *   Meissel   -- pi(x) = phi(x, a) + a - 1 - P2,  a = pi(x^1/3).
 *   Lehmer    -- Meissel refined with a = pi(x^1/4) plus a P3 correction.
 *   LMO / DR  -- Lagarias-Miller-Odlyzko / Deleglise-Rivat (segmented).
 *
 * The methods deliberately overlap in range so they cross-validate each other;
 * tests check all of them against one another and against the known pi(10^k). */

#include "primecount.h"
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>

/* ------------------------------------------------------------------ */
/* integer roots                                                       */
/* ------------------------------------------------------------------ */

static int64_t isqrt64(int64_t n) {
    if (n < 0) return 0;
    int64_t r = (int64_t)sqrtl((long double)n);
    while (r > 0 && r * r > n) r--;
    while ((r + 1) * (r + 1) <= n) r++;
    return r;
}

static int64_t icbrt64(int64_t n) {
    if (n < 0) return 0;
    int64_t r = (int64_t)cbrtl((long double)n);
    while (r > 0 && r * r * r > n) r--;
    while ((r + 1) * (r + 1) * (r + 1) <= n) r++;
    return r;
}

static int64_t iroot4(int64_t n) {
    int64_t r = isqrt64(isqrt64(n));
    while (r > 0 && r * r * r * r > n) r--;
    while ((r + 1) * (r + 1) * (r + 1) * (r + 1) <= n) r++;
    return r;
}

/* ------------------------------------------------------------------ */
/* small-prime table (primes < 10^6)                                   */
/* ------------------------------------------------------------------ */

#define SMALL_PRIME_LIMIT 1000000

static uint32_t *small_primes = NULL;
static uint32_t  small_prime_count = 0;

void primecount_init(void) {
    if (small_primes) return;
    uint8_t *is_prime = calloc((size_t)SMALL_PRIME_LIMIT + 1, sizeof(uint8_t));
    if (!is_prime) return;
    for (uint32_t i = 2; i <= SMALL_PRIME_LIMIT; i++) is_prime[i] = 1;
    for (uint32_t p = 2; (uint64_t)p * p <= SMALL_PRIME_LIMIT; p++)
        if (is_prime[p])
            for (uint32_t i = p * p; i <= SMALL_PRIME_LIMIT; i += p) is_prime[i] = 0;
    for (uint32_t i = 2; i <= SMALL_PRIME_LIMIT; i++)
        if (is_prime[i]) small_prime_count++;
    small_primes = malloc((size_t)small_prime_count * sizeof(uint32_t));
    if (small_primes) {
        uint32_t idx = 0;
        for (uint32_t i = 2; i <= SMALL_PRIME_LIMIT; i++)
            if (is_prime[i]) small_primes[idx++] = i;
    } else {
        small_prime_count = 0;
    }
    free(is_prime);
}

int64_t  primecount_small_table_size(void) { primecount_init(); return small_prime_count; }
uint32_t primecount_small_prime(int64_t n) { return small_primes[n - 1]; }

/* pi(x) for x <= SMALL_PRIME_LIMIT via binary search over the table. */
static int64_t pi_small(int64_t x) {
    if (x < 2) return 0;
    int64_t lo = 0, hi = (int64_t)small_prime_count - 1;
    while (lo <= hi) {
        int64_t mid = lo + (hi - lo) / 2;
        if ((int64_t)small_primes[mid] <= x) lo = mid + 1;
        else hi = mid - 1;
    }
    return lo;
}

/* ------------------------------------------------------------------ */
/* generate primes up to `bound` (bound is small, <= sqrt(x))          */
/* ------------------------------------------------------------------ */

/* Returns a malloc'd ascending array of all primes <= bound; *count set.
 * Caller frees.  Returns NULL on allocation failure. */
static uint32_t *gen_primes(int64_t bound, int64_t *count) {
    *count = 0;
    if (bound < 2) { uint32_t *e = malloc(sizeof(uint32_t)); return e; }
    uint8_t *sieve = calloc((size_t)bound + 1, 1);
    if (!sieve) return NULL;
    for (int64_t i = 2; i <= bound; i++) sieve[i] = 1;
    for (int64_t p = 2; p * p <= bound; p++)
        if (sieve[p])
            for (int64_t i = p * p; i <= bound; i += p) sieve[i] = 0;
    int64_t c = 0;
    for (int64_t i = 2; i <= bound; i++) if (sieve[i]) c++;
    uint32_t *primes = malloc((size_t)(c > 0 ? c : 1) * sizeof(uint32_t));
    if (!primes) { free(sieve); return NULL; }
    int64_t idx = 0;
    for (int64_t i = 2; i <= bound; i++) if (sieve[i]) primes[idx++] = (uint32_t)i;
    free(sieve);
    *count = c;
    return primes;
}

/* ------------------------------------------------------------------ */
/* PiTable: bit-packed prime sieve up to B with O(1) pi(v) lookups      */
/* ------------------------------------------------------------------ */

typedef struct {
    uint64_t *words;   /* bit i set  <=>  i is prime, for i in [0, B]   */
    uint32_t *cnt;     /* cnt[w] = number of primes in words[0 .. w-1]  */
    int64_t   B;
    int64_t   nwords;
} PiTable;

static int pitable_build(PiTable *t, int64_t B) {
    t->words = NULL; t->cnt = NULL; t->B = B;
    t->nwords = (B >> 6) + 1;
    t->words = calloc((size_t)t->nwords, sizeof(uint64_t));
    t->cnt   = malloc((size_t)t->nwords * sizeof(uint32_t));
    if (!t->words || !t->cnt) { free(t->words); free(t->cnt); return 0; }

    /* mark every integer in [2, B] as prime, then sieve out composites */
    for (int64_t w = 0; w < t->nwords; w++) t->words[w] = ~0ULL;
    t->words[0] &= ~3ULL;                      /* 0 and 1 are not prime */
    int64_t hi = t->nwords * 64 - 1;
    for (int64_t i = B + 1; i <= hi; i++)      /* clear padding bits past B */
        t->words[i >> 6] &= ~(1ULL << (i & 63));

    int64_t r = isqrt64(B);
    for (int64_t p = 2; p <= r; p++) {
        if (t->words[p >> 6] & (1ULL << (p & 63))) {
            for (int64_t m = p * p; m <= B; m += p)
                t->words[m >> 6] &= ~(1ULL << (m & 63));
        }
    }

    uint32_t acc = 0;
    for (int64_t w = 0; w < t->nwords; w++) {
        t->cnt[w] = acc;
        acc += (uint32_t)__builtin_popcountll(t->words[w]);
    }
    return 1;
}

static void pitable_free(PiTable *t) { free(t->words); free(t->cnt); t->words = NULL; t->cnt = NULL; }

/* pi(v) = number of primes <= v, for 0 <= v <= B. */
static int64_t pitable_count(const PiTable *t, int64_t v) {
    if (v < 2) return 0;
    if (v > t->B) v = t->B;
    int64_t w = v >> 6, r = v & 63;
    uint64_t mask = (r == 63) ? ~0ULL : ((1ULL << (r + 1)) - 1);
    return (int64_t)t->cnt[w] + __builtin_popcountll(t->words[w] & mask);
}

/* ------------------------------------------------------------------ */
/* phi(v, b): PhiTiny wheel base + Lehmer p^2 > v prune                 */
/* ------------------------------------------------------------------ */

/* Wheel over the first PHI_C primes. */
#define PHI_C 7
static const int PHI_C_PRIMES[PHI_C] = {2, 3, 5, 7, 11, 13, 17};
static int64_t  PHI_P = 0;          /* product of the first PHI_C primes */
static int64_t  PHI_PHI_P = 0;      /* totient of PHI_P (= 92160, exceeds 16 bits) */
static uint32_t *phi_prefix = NULL; /* phi_prefix[r] = #{1..r coprime to PHI_P} */

static void phi_tiny_init(void) {
    if (phi_prefix) return;
    int64_t P = 1;
    for (int i = 0; i < PHI_C; i++) P *= PHI_C_PRIMES[i];
    PHI_P = P;
    phi_prefix = malloc((size_t)(P + 1) * sizeof(uint32_t));
    if (!phi_prefix) return;
    uint8_t *coprime = malloc((size_t)(P + 1));
    for (int64_t i = 0; i <= P; i++) coprime[i] = 1;
    for (int i = 0; i < PHI_C; i++) {
        int p = PHI_C_PRIMES[i];
        for (int64_t m = p; m <= P; m += p) coprime[m] = 0;
    }
    uint32_t acc = 0;
    phi_prefix[0] = 0;
    for (int64_t i = 1; i <= P; i++) { acc += coprime[i]; phi_prefix[i] = acc; }
    PHI_PHI_P = phi_prefix[P];
    free(coprime);
}

/* phi(v, PHI_C) via the wheel. */
static int64_t phi_wheel(int64_t v) {
    if (v <= 0) return 0;
    return (v / PHI_P) * PHI_PHI_P + phi_prefix[v % PHI_P];
}

/* phi(v, b) for small b (< PHI_C) by inclusion-exclusion over first b primes. */
static int64_t phi_small_b(int64_t v, int b, const uint32_t *primes) {
    /* iterate the 2^b squarefree subsets of the first b primes */
    int64_t total = 0;
    for (int64_t mask = 0; mask < (1LL << b); mask++) {
        int64_t prod = 1;
        int bits = 0;
        int ok = 1;
        for (int i = 0; i < b; i++) {
            if (mask & (1LL << i)) {
                if (prod > v / primes[i]) { ok = 0; break; }  /* prod*prime > v */
                prod *= primes[i];
                bits++;
            }
        }
        if (ok) total += (bits & 1) ? -(v / prod) : (v / prod);
    }
    return total;
}

/* phi(v, b) = #{n <= v : n not divisible by any of the first b primes}.
 * `primes` is 0-indexed (primes[0] = 2); `pt` covers all pruned leaves. */
static int64_t phi_rec(int64_t v, int64_t b, const uint32_t *primes, const PiTable *pt) {
    if (v == 0) return 0;
    if (b == 0) return v;
    if (b <= PHI_C) {
        if (b == PHI_C) return phi_wheel(v);
        return phi_small_b(v, (int)b, primes);
    }
    /* Lehmer prune: if p_b^2 > v then n <= v coprime to the first b primes is
     * either 1 or a single prime in (p_b, v], so phi(v, b) = pi(v) - b + 1. */
    int64_t pb = primes[b - 1];
    if (pb * pb > v) return pitable_count(pt, v) - b + 1;
    return phi_rec(v, b - 1, primes, pt) - phi_rec(v / pb, b - 1, primes, pt);
}

/* ------------------------------------------------------------------ */
/* Method: Lucy_Hedgehog (exact, O(x^3/4))                              */
/* ------------------------------------------------------------------ */

static int64_t count_lucy(int64_t n) {
    if (n < 2) return 0;
    int64_t r = isqrt64(n);
    int64_t *small = malloc((size_t)(r + 1) * sizeof(int64_t));
    int64_t *large = malloc((size_t)(r + 1) * sizeof(int64_t));
    if (!small || !large) { free(small); free(large); return -1; }

    for (int64_t i = 1; i <= r; i++) { small[i] = i - 1; large[i] = n / i - 1; }
    for (int64_t p = 2; p <= r; p++) {
        if (small[p] == small[p - 1]) continue;     /* p composite */
        int64_t sp = small[p - 1];
        int64_t p2 = p * p;
        int64_t imax = n / p2; if (imax > r) imax = r;
        for (int64_t i = 1; i <= imax; i++) {
            int64_t d = i * p;
            int64_t prev = (d <= r) ? large[d] : small[n / d];
            large[i] -= prev - sp;
        }
        for (int64_t v = r; v >= p2; v--) small[v] -= small[v / p] - sp;
    }
    int64_t result = large[1];
    free(small); free(large);
    return result;
}

/* ------------------------------------------------------------------ */
/* Method: segmented sieve of Eratosthenes                             */
/* ------------------------------------------------------------------ */

#define SIEVE_MAX 10000000000LL     /* 10^10: above this the sieve is too slow */

static int64_t count_sieve(int64_t x) {
    if (x < 2) return 0;
    if (x > SIEVE_MAX) return -1;
    int64_t sqrtx = isqrt64(x);
    int64_t bp_count;
    uint32_t *bp = gen_primes(sqrtx, &bp_count);
    if (!bp) return -1;

    const int64_t SEG = 1 << 18;
    uint8_t *seg = malloc((size_t)SEG);
    if (!seg) { free(bp); return -1; }

    int64_t count = 0;
    for (int64_t low = 2; low <= x; low += SEG) {
        int64_t high = low + SEG - 1; if (high > x) high = x;
        int64_t len = high - low + 1;
        memset(seg, 1, (size_t)len);
        for (int64_t i = 0; i < bp_count; i++) {
            int64_t p = bp[i];
            if (p * p > high) break;
            int64_t start = ((low + p - 1) / p) * p;
            if (start < p * p) start = p * p;
            for (int64_t j = start; j <= high; j += p) seg[j - low] = 0;
        }
        for (int64_t k = 0; k < len; k++) if (seg[k]) count++;
    }
    free(seg); free(bp);
    return count;
}

/* ------------------------------------------------------------------ */
/* Methods: Legendre, Meissel, Lehmer (PiTable + pruned phi)            */
/* ------------------------------------------------------------------ */

#define LEGENDRE_MAX 1000000000LL    /* 10^9:  PiTable spans [0, x]      */
#define MEISSEL_MAX  10000000000000LL/* 10^13: PiTable spans [0, x^2/3]  */
#define LEHMER_MAX   1000000000000LL /* 10^12: PiTable spans [0, x^3/4]  */

static int64_t count_legendre(int64_t x) {
    if (x < 2) return 0;
    if (x > LEGENDRE_MAX) return -1;
    phi_tiny_init();
    int64_t sqrtx = isqrt64(x);
    int64_t a;
    uint32_t *primes = gen_primes(sqrtx, &a);     /* a = pi(sqrt x) */
    if (!primes) return -1;
    PiTable pt;
    if (!pitable_build(&pt, x)) { free(primes); return -1; }   /* prune needs pi up to x */
    int64_t result = phi_rec(x, a, primes, &pt) + a - 1;
    pitable_free(&pt);
    free(primes);
    return result;
}

static int64_t count_meissel(int64_t x) {
    if (x < 2) return 0;
    if (x > MEISSEL_MAX) return -1;
    phi_tiny_init();
    int64_t sqrtx = isqrt64(x);
    int64_t cbrtx = icbrt64(x);
    int64_t nprimes;
    uint32_t *primes = gen_primes(sqrtx, &nprimes);  /* primes up to sqrt x */
    if (!primes) return -1;
    PiTable pt;
    int64_t B = x / cbrtx;                            /* ~ x^2/3 >= x/p for p > x^1/3 */
    if (!pitable_build(&pt, B)) { free(primes); return -1; }

    int64_t a = pitable_count(&pt, cbrtx);            /* pi(x^1/3) */
    int64_t b = pitable_count(&pt, sqrtx);            /* pi(sqrt x) */

    int64_t result = phi_rec(x, a, primes, &pt) + a - 1;
    /* P2 = sum_{i=a+1}^{b} (pi(x/p_i) - i + 1) */
    for (int64_t i = a + 1; i <= b; i++) {
        int64_t p = primes[i - 1];
        result -= pitable_count(&pt, x / p) - i + 1;
    }
    pitable_free(&pt);
    free(primes);
    return result;
}

static int64_t count_lehmer(int64_t x) {
    if (x < 2) return 0;
    if (x > LEHMER_MAX) return -1;
    phi_tiny_init();
    int64_t sqrtx  = isqrt64(x);
    int64_t cbrtx  = icbrt64(x);
    int64_t root4x = iroot4(x);
    int64_t nprimes;
    uint32_t *primes = gen_primes(sqrtx, &nprimes);
    if (!primes) return -1;
    PiTable pt;
    /* P2 needs pi(x/p_i) for p_i > x^1/4, i.e. arguments up to x^3/4. */
    int64_t B = x / root4x;
    if (!pitable_build(&pt, B)) { free(primes); return -1; }

    int64_t a = pitable_count(&pt, root4x);   /* pi(x^1/4) */
    int64_t b = pitable_count(&pt, sqrtx);    /* pi(sqrt x) */
    int64_t c = pitable_count(&pt, cbrtx);    /* pi(x^1/3) */

    int64_t result = phi_rec(x, a, primes, &pt)
                   + (b + a - 2) * (b - a + 1) / 2;
    for (int64_t i = a + 1; i <= b; i++) {
        int64_t xi = x / primes[i - 1];
        result -= pitable_count(&pt, xi);
        if (i <= c) {
            int64_t bi = pitable_count(&pt, isqrt64(xi));
            for (int64_t j = i; j <= bi; j++)
                result += (j - 1) - pitable_count(&pt, xi / primes[j - 1]);
        }
    }
    pitable_free(&pt);
    free(primes);
    return result;
}

/* ------------------------------------------------------------------ */
/* Method: Lagarias-Miller-Odlyzko (segmented special leaves)          */
/* ------------------------------------------------------------------ */

/* Fenwick (binary-indexed) tree over positions [0, n) holding 0/1 counts,
 * used per segment to query phi(v, b) = #{unsieved positions in [0, v]}. */
static void bit_build_ones(int32_t *bit, int64_t n) {
    /* bit[1..n], all underlying values 1; O(n) build */
    for (int64_t i = 1; i <= n; i++) bit[i] = 1;
    for (int64_t i = 1; i <= n; i++) {
        int64_t j = i + (i & -i);
        if (j <= n) bit[j] += bit[i];
    }
}
static void bit_update(int32_t *bit, int64_t n, int64_t pos, int32_t delta) {
    for (int64_t i = pos + 1; i <= n; i += i & -i) bit[i] += delta;
}
static int64_t bit_query(const int32_t *bit, int64_t pos) {   /* sum over [0, pos] */
    int64_t s = 0;
    for (int64_t i = pos + 1; i > 0; i -= i & -i) s += bit[i];
    return s;
}

/* Computes pi(x) by the LMO method:
 *   pi(x) = phi(x, a) + a - 1 - P2(x, a),  a = pi(x^1/3),
 * with phi(x, a) = S1 (ordinary leaves) + S2 (special leaves).  S2 is evaluated
 * over [1, z] (z = x/y) in segments carrying a Fenwick tree, so its memory is
 * O(sqrt z); P2 reuses a PiTable up to z. */
static int64_t count_lmo(int64_t x) {
    if (x < 2) return 0;
    if (x > PI_COUNT_MAX) return -1;
    phi_tiny_init();

    int64_t y     = icbrt64(x);
    int64_t z     = x / y;
    int64_t sqrtx = isqrt64(x);
    const int c   = PHI_C;          /* PhiTiny base: first 7 primes */

    int64_t nprimes;
    uint32_t *primes = gen_primes(sqrtx, &nprimes);
    if (!primes) return -1;
    PiTable pt;
    if (!pitable_build(&pt, z)) { free(primes); return -1; }

    int64_t a   = pitable_count(&pt, y);      /* pi(x^1/3) */
    int64_t bpi = pitable_count(&pt, sqrtx);  /* pi(sqrt x) */
    int64_t pc  = primes[c - 1];              /* p_c = 17 */

    /* mu(n) and least-prime-factor(n) for n in [1, y]. */
    int8_t  *mu  = malloc((size_t)(y + 1) * sizeof(int8_t));
    int64_t *lpf = malloc((size_t)(y + 1) * sizeof(int64_t));
    if (!mu || !lpf) { free(mu); free(lpf); pitable_free(&pt); free(primes); return -1; }
    for (int64_t i = 0; i <= y; i++) { mu[i] = 1; lpf[i] = INT64_MAX; }
    for (int64_t p = 2; p <= y; p++) {
        if (lpf[p] == INT64_MAX) {              /* p is prime */
            for (int64_t m = p; m <= y; m += p) { if (lpf[m] == INT64_MAX) lpf[m] = p; mu[m] = (int8_t)-mu[m]; }
            int64_t p2 = p * p;
            for (int64_t m = p2; m <= y; m += p2) mu[m] = 0;
        }
    }

    /* S1: ordinary leaves -- sum over n <= y with lpf(n) > p_c. */
    int64_t S1 = 0;
    for (int64_t n = 1; n <= y; n++)
        if (mu[n] != 0 && lpf[n] > pc)
            S1 += mu[n] * phi_wheel(x / n);

    /* S2: special leaves over [1, z], segmented with a Fenwick tree. */
    int64_t S2 = 0;
    const int64_t SEG = 1 << 18;
    int32_t *bit    = malloc((size_t)(SEG + 1) * sizeof(int32_t));
    uint8_t *sieve  = malloc((size_t)SEG);
    int64_t *phi    = calloc((size_t)(a + 1), sizeof(int64_t)); /* phi[k] = phi(low-1, k) */
    int64_t *leaf_m = malloc((size_t)(a + 1) * sizeof(int64_t));
    if (!bit || !sieve || !phi || !leaf_m) {
        free(bit); free(sieve); free(phi); free(leaf_m);
        free(mu); free(lpf); pitable_free(&pt); free(primes); return -1;
    }
    for (int64_t k = 0; k <= a; k++) leaf_m[k] = y;

    for (int64_t low = 1; low <= z; low += SEG) {
        int64_t high = low + SEG; if (high > z + 1) high = z + 1;
        int64_t len = high - low;
        memset(sieve, 1, (size_t)len);
        bit_build_ones(bit, len);

        /* base: remove the first c primes from this segment */
        for (int k = 1; k <= c; k++) {
            int64_t p = primes[k - 1];
            int64_t start = ((low + p - 1) / p) * p;
            for (int64_t mlt = start; mlt < high; mlt += p) {
                int64_t pos = mlt - low;
                if (sieve[pos]) { sieve[pos] = 0; bit_update(bit, len, pos, -1); }
            }
        }

        /* levels b = c+1 .. a-1 (1-indexed prime index) */
        for (int64_t b = c + 1; b <= a - 1; b++) {
            int64_t prime = primes[b - 1];      /* p_b */
            int64_t lim = y / prime;
            while (leaf_m[b] > lim) {
                int64_t m = leaf_m[b];
                if (mu[m] != 0 && prime < lpf[m]) {
                    int64_t v = x / (prime * m);
                    if (v >= high) break;       /* leaf belongs to a later segment */
                    if (v >= low)
                        S2 -= mu[m] * (phi[b - 1] + bit_query(bit, v - low));
                }
                leaf_m[b]--;
            }
            phi[b - 1] += bit_query(bit, len - 1);   /* carry phi(., b-1) to next segment */
            int64_t start = ((low + prime - 1) / prime) * prime;
            for (int64_t mlt = start; mlt < high; mlt += prime) {
                int64_t pos = mlt - low;
                if (sieve[pos]) { sieve[pos] = 0; bit_update(bit, len, pos, -1); }
            }
        }
    }

    int64_t phi_xa = S1 + S2;

    /* P2 = sum_{i=a+1}^{pi(sqrt x)} (pi(x/p_i) - i + 1) */
    int64_t P2 = 0;
    for (int64_t i = a + 1; i <= bpi; i++)
        P2 += pitable_count(&pt, x / primes[i - 1]) - i + 1;

    int64_t result = phi_xa + a - 1 - P2;

    free(bit); free(sieve); free(phi); free(leaf_m);
    free(mu); free(lpf); pitable_free(&pt); free(primes);
    return result;
}

/* ------------------------------------------------------------------ */
/* Method: Deleglise-Rivat                                             */
/* ------------------------------------------------------------------ */

/* Deleglise-Rivat, with the standard tuning that makes it sublinear in
 * practice rather than only in theory:
 *
 *   y = alpha * x^(1/3),  alpha = 0.0008 ln(x)^3 (clamped to [1, sqrt(x)/x^(1/3)])
 *
 * pi(x) = phi(x, a) + a - 1 - P2(x, a) with a = pi(y); P3 vanishes because
 * p_{a+1} > y >= x^(1/3). phi(x, a) = S1 + S2 over the ordinary leaves
 * (n <= y, squarefree, lpf(n) > p_c: mu(n) phi(x/n, c) off the PhiTiny wheel)
 * and the special leaves n = p_b m (m <= y < p_b m, squarefree m, lpf(m) >
 * p_b): -mu(m) phi(v, b-1) with v = floor(x / (p_b m)) < z = x / y.
 *
 * Leaves are split by how phi(v, b-1) can be had:
 *   trivial  v < p_b            phi = 1
 *   easy     v < p_{b-1}^2      phi = pi(v) - (b-1) + 1   (one PiTable lookup)
 *   hard     otherwise          the segmented sieve over [1, z]
 * Hard leaves need p_{b-1}^2 <= v < z, so only levels with p_{b-1} <= sqrt(z)
 * ever touch the sieve; every level above that is answered from the PiTable
 * alone, and there (p_b > sqrt(z) >= sqrt(y)) a valid m is necessarily a prime,
 * so the trivial ones are COUNTED (one pi difference per level) and only the
 * easy ones are visited. Within the sieve, levels with p_b > sqrt(y) likewise
 * walk the prime list rather than every m in (y/p_b, y].
 *
 * The sieve is a bitset per segment; a level's hard queries arrive in
 * increasing v (m decreases), so phi(v, b-1) is a running popcount sweep --
 * O(1) per removal, O(segment/64) per level and segment -- instead of LMO's
 * Fenwick tree with its O(log) on every crossed-off number. Quotients are
 * floor(floor(x/p_b)/m) from a precomputed 1/m and an exact one-step
 * correction, so no 64-bit division runs per leaf.
 *
 * Measured against the y = x^(1/3) version this replaces (one core, x86-64):
 * 10^11 0.54 s -> 0.04 s, 10^12 4.4 s -> 0.20 s, 10^13 40 s -> 1.06 s.
 * Cross-checked against Meissel and LMO on 400 random x in [10^6, 10^11] and
 * on cube/square edges for every alpha from 1 to 64, and against the known
 * pi(10^k) and Mathematica's pi(5 * 10^13). */

double primecount_dr_alpha = -1.0;   /* > 0 forces alpha (tests); < 0 automatic */

static double dr_alpha(int64_t x) {
    if (primecount_dr_alpha > 0) return primecount_dr_alpha;
    double l = log((double)x);
    double a = 0.0008 * l * l * l;
    return a < 1.0 ? 1.0 : a;
}

/* floor(xb / m) from a reciprocal estimate, corrected exactly (xb < 2^53). */
static inline int64_t dr_div(int64_t xb, int64_t m, const double* inv) {
    int64_t q = (int64_t)((double)xb * inv[m]);
    if (q * m > xb) q--;
    else if ((q + 1) * m <= xb) q++;
    return q;
}

/* phi(v, b-1) for a hard leaf: the carried count below the segment plus a
 * running popcount up to v (queries at one level arrive in increasing v). */
static inline int64_t dr_sweep(const uint64_t* bits, int64_t pos, int64_t* wpos,
                               int64_t* acc) {
    int64_t wi = pos >> 6;
    while (*wpos < wi) *acc += __builtin_popcountll(bits[(*wpos)++]);
    uint64_t mask = ((pos & 63) == 63) ? ~0ULL : ((1ULL << ((pos & 63) + 1)) - 1);
    return *acc + __builtin_popcountll(bits[wi] & mask);
}

static int64_t count_dr(int64_t x) {
    if (x < 2) return 0;
    if (x > PI_COUNT_MAX) return -1;
    primecount_init();
    if (x <= SMALL_PRIME_LIMIT) return pi_small(x);
    phi_tiny_init();

    int64_t cb = icbrt64(x), sqrtx = isqrt64(x);
    int64_t y = (int64_t)(dr_alpha(x) * (double)cb);
    if (y < cb) y = cb;
    if (y > sqrtx) y = sqrtx;
    int64_t z = x / y;
    const int c = PHI_C;

    int64_t nprimes;
    uint32_t *primes = gen_primes(sqrtx, &nprimes);     /* y <= sqrtx */
    if (!primes) return -1;
    PiTable pt;
    if (!pitable_build(&pt, z)) { free(primes); return -1; }
    int64_t a   = pitable_count(&pt, y);
    int64_t bpi = pitable_count(&pt, sqrtx);
    int64_t pc  = primes[c - 1];

    /* mu(n), least prime factor lpf(n) (lpf(1) = "infinity") and 1/n, n <= y. */
    int8_t   *mu  = malloc((size_t)(y + 1));
    uint32_t *lpf = malloc((size_t)(y + 1) * sizeof(uint32_t));
    double   *inv = malloc((size_t)(y + 1) * sizeof(double));
    if (!mu || !lpf || !inv) {
        free(mu); free(lpf); free(inv); pitable_free(&pt); free(primes); return -1;
    }
    for (int64_t i = 0; i <= y; i++) { mu[i] = 1; lpf[i] = 0; inv[i] = i ? 1.0 / (double)i : 0.0; }
    for (int64_t p = 2; p <= y; p++) {
        if (lpf[p] == 0) {
            for (int64_t m = p; m <= y; m += p) { if (!lpf[m]) lpf[m] = (uint32_t)p; mu[m] = (int8_t)-mu[m]; }
            if (p <= y / p) for (int64_t m = p * p; m <= y; m += p * p) mu[m] = 0;
        }
    }
    lpf[1] = UINT32_MAX;

    /* S1: ordinary leaves. */
    int64_t S1 = 0;
    for (int64_t n = 1; n <= y; n++)
        if (mu[n] != 0 && (int64_t)lpf[n] > pc) S1 += mu[n] * phi_wheel(x / n);

    /* bh: the last level that can hold a hard leaf (p_{b-1} <= sqrt(z)). */
    int64_t sqz = isqrt64(z);
    int64_t bh = c;
    while (bh + 1 <= a - 1 && (int64_t)primes[bh - 1] <= sqz) bh++;

    int64_t S2 = 0;

    /* Levels above bh: m is a prime q in (p_b, y], mu(q) = -1. */
    for (int64_t b = bh + 1; b <= a - 1; b++) {
        int64_t pb = primes[b - 1];
        int64_t qmax = x / (pb * pb);               /* q <= qmax  <=>  v >= p_b */
        if (qmax > y) qmax = y;
        int64_t i_nt = (qmax > pb) ? pitable_count(&pt, qmax) : b;
        S2 += a - i_nt;                             /* trivial: phi = 1 each */
        int64_t xb = x / pb;
        for (int64_t i = b + 1; i <= i_nt; i++)     /* easy */
            S2 += pitable_count(&pt, dr_div(xb, primes[i - 1], inv)) - b + 2;
    }

    /* Levels c+1 .. bh: one segmented pass over [1, z]. */
    if (bh >= c + 1) {
        const int64_t SEG = (int64_t)1 << 19;       /* bits per segment */
        uint64_t *bits = malloc((size_t)(SEG / 64) * sizeof(uint64_t));
        int64_t  *phic = calloc((size_t)(bh + 1), sizeof(int64_t)); /* phi(low-1, b) */
        int64_t  *lm   = malloc((size_t)(bh + 1) * sizeof(int64_t)); /* next m (or prime index) */
        if (!bits || !phic || !lm) {
            free(bits); free(phic); free(lm);
            free(mu); free(lpf); free(inv); pitable_free(&pt); free(primes); return -1;
        }
        for (int64_t b = 0; b <= bh; b++)
            lm[b] = (b >= 1 && (int64_t)primes[b - 1] > y / (int64_t)primes[b - 1]) ? a : y;

        for (int64_t low = 1; low <= z; low += SEG) {
            int64_t high = low + SEG; if (high > z + 1) high = z + 1;
            int64_t len = high - low, nw = (len + 63) / 64;
            for (int64_t w = 0; w < nw; w++) bits[w] = ~0ULL;
            if (len & 63) bits[nw - 1] = (1ULL << (len & 63)) - 1;
            for (int k = 0; k < c; k++) {               /* the wheel primes */
                int64_t p = primes[k];
                for (int64_t m = ((low + p - 1) / p) * p; m < high; m += p) {
                    int64_t pos = m - low; bits[pos >> 6] &= ~(1ULL << (pos & 63));
                }
            }
            for (int64_t b = c + 1; b <= bh; b++) {
                int64_t pb = primes[b - 1], pp = primes[b - 2], pp2 = pp * pp;
                int64_t xb = x / pb, lim = y / pb;
                int64_t wpos = 0, acc = 0;
                if (pb > lim) {
                    /* p_b > sqrt(y): walk primes q, descending (v ascending). */
                    while (lm[b] > b) {
                        int64_t v = dr_div(xb, primes[lm[b] - 1], inv);
                        if (v >= high) break;
                        S2 += (v < pb) ? 1
                            : (pp2 > v) ? pitable_count(&pt, v) - b + 2
                            : phic[b - 1] + dr_sweep(bits, v - low, &wpos, &acc);
                        lm[b]--;
                    }
                } else {
                    while (lm[b] > lim) {
                        int64_t m = lm[b];
                        if (mu[m] != 0 && pb < (int64_t)lpf[m]) {
                            int64_t v = dr_div(xb, m, inv);
                            if (v >= high) break;
                            int64_t phv = (v < pb) ? 1
                                : (pp2 > v) ? pitable_count(&pt, v) - b + 2
                                : phic[b - 1] + dr_sweep(bits, v - low, &wpos, &acc);
                            S2 -= mu[m] * phv;
                        }
                        lm[b]--;
                    }
                }
                while (wpos < nw) acc += __builtin_popcountll(bits[wpos++]);
                phic[b - 1] += acc;                       /* carry phi(., b-1) */
                for (int64_t m = ((low + pb - 1) / pb) * pb; m < high; m += pb) {
                    int64_t pos = m - low; bits[pos >> 6] &= ~(1ULL << (pos & 63));
                }
            }
        }
        free(bits); free(phic); free(lm);
    }

    /* P2 = sum_{i=a+1}^{pi(sqrt x)} (pi(x/p_i) - i + 1) */
    int64_t P2 = 0;
    for (int64_t i = a + 1; i <= bpi; i++)
        P2 += pitable_count(&pt, x / primes[i - 1]) - i + 1;

    int64_t result = S1 + S2 + a - 1 - P2;
    free(mu); free(lpf); free(inv); pitable_free(&pt); free(primes);
    return result;
}


/* ------------------------------------------------------------------ */
/* the n-th prime from an anchor (Prime[n])                            */
/* ------------------------------------------------------------------ */

/* Given c = pi(x), return the n-th prime by sieving the gap between x and it
 * -- forward when n > c, backward when n <= c -- in segments of 2^18 odd and
 * even numbers alike, with the base primes up to sqrt of the far end. The gap
 * is ~|n - c| ln x numbers; from a Riemann-R estimate that is O(sqrt x), so
 * this costs a few milliseconds where stepping mpz_nextprime cost a primality
 * test per prime. Returns -1 on allocation failure or overflow. */
int64_t primecount_nth_from(int64_t x, int64_t c, int64_t n) {
    if (n <= 0) return -1;
    const int64_t SEG = (int64_t)1 << 18;
    double lx = log((double)(x > 3 ? x : 3));
    int64_t span = (int64_t)(((double)(n > c ? n - c : c - n) + 64.0) * (lx + 2.0)) + SEG;
    int64_t far = (n > c) ? x + span : x;
    if (far < x) return -1;                              /* overflow */
    int64_t nbp;
    uint32_t *bp = gen_primes(isqrt64(far + 2 * span) + 1, &nbp);
    if (!bp) return -1;
    uint8_t *seg = malloc((size_t)SEG);
    if (!seg) { free(bp); return -1; }

    int64_t result = -1;
    int64_t k = c;
    if (n > c) {
        for (int64_t low = x + 1; result < 0; low += SEG) {
            int64_t high = low + SEG;                     /* [low, high) */
            if (high < low) break;
            memset(seg, 1, (size_t)SEG);
            for (int64_t i = 0; i < nbp; i++) {
                int64_t p = bp[i];
                if (p * p >= high) break;
                int64_t st = ((low + p - 1) / p) * p;
                if (st < p * p) st = p * p;
                for (int64_t m = st; m < high; m += p) seg[m - low] = 0;
            }
            for (int64_t j = 0; j < SEG; j++) {
                int64_t v = low + j;
                if (v >= 2 && seg[j] && ++k == n) { result = v; break; }
            }
        }
    } else {
        for (int64_t high = x + 1; result < 0 && high > 2; high -= SEG) {
            int64_t low = high - SEG; if (low < 2) low = 2;   /* [low, high) */
            int64_t len = high - low;
            memset(seg, 1, (size_t)len);
            for (int64_t i = 0; i < nbp; i++) {
                int64_t p = bp[i];
                if (p * p >= high) break;
                int64_t st = ((low + p - 1) / p) * p;
                if (st < p * p) st = p * p;
                for (int64_t m = st; m < high; m += p) seg[m - low] = 0;
            }
            for (int64_t j = len - 1; j >= 0; j--) {
                if (seg[j]) {                              /* the k-th prime */
                    if (k == n) { result = low + j; break; }
                    k--;
                }
            }
        }
    }
    free(seg); free(bp);
    return result;
}

/* ------------------------------------------------------------------ */
/* dispatch                                                            */
/* ------------------------------------------------------------------ */

int64_t prime_count(int64_t x, PrimeCountMethod method) {
    if (x < 2) return 0;

    if (method == PC_AUTOMATIC) {
        primecount_init();
        if (x <= SMALL_PRIME_LIMIT) return pi_small(x);
        if (x > PI_COUNT_MAX) return -1;
        /* The tuned Deleglise-Rivat wins from the small-table limit up
         * (measured, one core: 4 ms at 10^9 against Lucy's 30 ms, 0.20 s at
         * 10^12 against 2.5 s). */
        return count_dr(x);
    }

    /* Fast exact path for tiny x shared by all methods. */
    if (x <= SMALL_PRIME_LIMIT && method != PC_SIEVE) { primecount_init(); return pi_small(x); }

    switch (method) {
        case PC_LUCY:     return (x > PI_COUNT_MAX) ? -1 : count_lucy(x);
        case PC_SIEVE:    return count_sieve(x);
        case PC_LEGENDRE: return count_legendre(x);
        case PC_MEISSEL:  return count_meissel(x);
        case PC_LEHMER:   return count_lehmer(x);
        case PC_LMO:      return (x > PI_COUNT_MAX) ? -1 : count_lmo(x);
        case PC_DR:       return (x > PI_COUNT_MAX) ? -1 : count_dr(x);
        default:          return -1;
    }
}
