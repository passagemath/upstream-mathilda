/* nt_dlog.c -- discrete logarithms in a cyclic subgroup of (Z/nZ)^*.
 *
 * The engine behind MultiplicativeOrder[k, n, {r1, r2, ...}]: the smallest
 * m >= 1 with k^m == r_i (mod n) for some i.  Given d = ord(k) and its
 * factorisation, each target is solved by
 *
 *   Pohlig-Hellman  over the prime powers q^e || d, digit by digit, the
 *                   residues recombined by CRT;
 *   BSGS            for each prime q of up to DL_BSGS_BITS (u64) / DL_BSGS_BITS_MPZ
 *                   (mpz) bits: one baby-step table per q, shared by every digit
 *                   of every target -- the common case pays it once;
 *   Pollard rho     above that: Teske's r-adding walk with Brent cycle
 *                   detection, exponents tracked mod q.
 *
 * Arithmetic runs on one of two back-ends: a u64 one (Montgomery form for odd
 * n, a 128-bit remainder for even n) whenever n < 2^63 and the compiler has a
 * 128-bit integer, else GMP.  The Pohlig-Hellman driver itself is in GMP: it
 * does O(e log d) operations per prime, against the O(sqrt q) of the solver.
 *
 * Correctness does not rest on the solvers: every accepted exponent is checked
 * with one mpz_powm before it is returned.  Membership is decided as follows:
 * a BSGS miss certifies h not in <gamma>; a non-degenerate rho collision gives
 * the only candidate, so a failed check certifies non-membership too.  The one
 * uncertified case is a rho run whose collisions keep being degenerate, which
 * can only persist when h lies outside <gamma> -- possible only when (Z/nZ)^*
 * is not cyclic; there the solver gives up after DL_RHO_ATTEMPTS bounded walks
 * (failure probability far below 2^-100 at the sizes rho is used for). */

#include "numbertheory_internal.h"
#include "core.h"        /* tc_deadline_passed */
#include <gmp.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

typedef enum { DL_FOUND, DL_NONMEMBER, DL_ABORT } DlStatus;

/* BSGS while q has at most this many bits; Pollard rho above.  36 bits keeps
 * the table at 2^18 entries (3 MB) on the u64 back-end; the mpz back-end stores
 * every baby step as an mpz_t, so its threshold is lower. */
#define DL_BSGS_BITS      36
#define DL_BSGS_BITS_MPZ  32
/* Rho: number of r-adding multipliers, and the poll period for
 * TimeConstrained (a clean return, not a siglongjmp mid-frame). */
#define DL_RHO_R         32
#define DL_POLL_MASK     ((UINT64_C(1) << 14) - 1)
/* Non-cyclic groups only: bounded rho walks before declaring non-membership. */
#define DL_RHO_ATTEMPTS  4

/* =====================================================================
 * Shared helpers
 * ===================================================================*/

/* Deterministic xorshift for the walk seeds; the answer is unique, so the
 * seeds affect only the running time. */
static uint64_t dl_rand(uint64_t* s) {
    uint64_t x = *s;
    x ^= x << 13; x ^= x >> 7; x ^= x << 17;
    return *s = x;
}

/* A uniform-ish draw in [0, q) as an mpz. */
static void dl_rand_mpz(mpz_t out, const mpz_t q, uint64_t* s) {
    mpz_set_ui(out, (unsigned long)(dl_rand(s) & 0xffffffffu));
    mpz_mul_2exp(out, out, 32);
    mpz_add_ui(out, out, (unsigned long)(dl_rand(s) & 0xffffffffu));
    mpz_mul_2exp(out, out, 32);
    mpz_add_ui(out, out, (unsigned long)(dl_rand(s) & 0xffffffffu));
    mpz_mod(out, out, q);
}

/* Open-addressing table: key -> baby-step index.  val == DL_EMPTY marks a
 * free slot. */
#define DL_EMPTY UINT32_MAX
typedef struct {
    uint64_t* key;
    uint32_t* val;
    uint64_t  mask;
    int       shift;
} DlHash;

static bool dl_hash_init(DlHash* H, uint64_t m) {
    int bits = 4;
    while ((UINT64_C(1) << bits) < 2 * m) bits++;
    uint64_t size = UINT64_C(1) << bits;
    H->key = (uint64_t*)malloc((size_t)size * sizeof(uint64_t));
    H->val = (uint32_t*)malloc((size_t)size * sizeof(uint32_t));
    if (!H->key || !H->val) { free(H->key); free(H->val); return false; }
    memset(H->val, 0xff, (size_t)size * sizeof(uint32_t));
    H->mask = size - 1;
    H->shift = 64 - bits;
    return true;
}

static void dl_hash_free(DlHash* H) { free(H->key); free(H->val); }

static inline uint64_t dl_hash_slot(const DlHash* H, uint64_t key) {
    return (key * UINT64_C(0x9E3779B97F4A7C15)) >> H->shift;
}

static inline void dl_hash_put(DlHash* H, uint64_t key, uint32_t v) {
    uint64_t s = dl_hash_slot(H, key);
    while (H->val[s] != DL_EMPTY) s = (s + 1) & H->mask;
    H->key[s] = key;
    H->val[s] = v;
}

/* ceil(sqrt(q)) for q >= 1. */
static uint64_t dl_isqrt_ceil(const mpz_t q) {
    mpz_t r;
    mpz_init(r);
    mpz_sqrt(r, q);
    unsigned long m = mpz_get_ui(r);     /* q < 2^36 here: fits 32 bits */
    mpz_mul(r, r, r);
    if (mpz_cmp(r, q) < 0) m++;
    mpz_clear(r);
    return (uint64_t)m;
}

/* =====================================================================
 * u64 back-end: n < 2^63
 * ===================================================================*/

#if defined(__SIZEOF_INT128__) && GMP_NUMB_BITS == 64
#define DL_HAVE_U64 1
__extension__ typedef unsigned __int128 dl_u128;

/* Values live in Montgomery form (x R mod n, R = 2^64) when n is odd; for an
 * even n they are plain residues and products go through a 128-bit remainder.
 * Either way the representation is a bijection on [0, n), so equality tests
 * and hashing work directly on it. */
typedef struct {
    uint64_t n, ninv, r2;
    int odd;
} DlMont;

static inline uint64_t dl_mul(const DlMont* M, uint64_t a, uint64_t b) {
    dl_u128 t = (dl_u128)a * b;
    if (!M->odd) return (uint64_t)(t % M->n);
    uint64_t lo = (uint64_t)t, hi = (uint64_t)(t >> 64);
    uint64_t m = lo * M->ninv;
    uint64_t mh = (uint64_t)(((dl_u128)m * M->n) >> 64);
    /* lo + low(m n) == 0 mod 2^64, carrying exactly when lo != 0; n < 2^63
     * keeps the sum below 2^64. */
    uint64_t r = hi + mh + (lo != 0);
    return r >= M->n ? r - M->n : r;
}

static uint64_t dl_get_u64(const mpz_t x) {
    return mpz_size(x) ? (uint64_t)mpz_getlimbn(x, 0) : 0;
}

static void dl_set_u64(mpz_t x, uint64_t v) {
    mpz_import(x, 1, -1, sizeof v, 0, 0, &v);
}

static void dl_mont_init(DlMont* M, const mpz_t n) {
    M->n = dl_get_u64(n);
    M->odd = (int)(M->n & 1);
    M->ninv = 0;
    M->r2 = 0;
    if (!M->odd) return;
    uint64_t inv = M->n;                         /* n^-1 mod 2^64 by Newton */
    for (int i = 0; i < 5; i++) inv *= 2 - M->n * inv;
    M->ninv = (uint64_t)0 - inv;
    mpz_t r;
    mpz_init_set_ui(r, 1);
    mpz_mul_2exp(r, r, 128);
    mpz_mod(r, r, n);
    M->r2 = dl_get_u64(r);
    mpz_clear(r);
}

/* mpz residue in [0, n) -> back-end representation, and back. */
static uint64_t dl_to(const DlMont* M, const mpz_t x) {
    uint64_t v = dl_get_u64(x);
    return M->odd ? dl_mul(M, v, M->r2) : v;
}

static void dl_from(const DlMont* M, mpz_t out, uint64_t v) {
    dl_set_u64(out, M->odd ? dl_mul(M, v, 1) : v);
}

static uint64_t dl_one(const DlMont* M) {
    return M->odd ? dl_mul(M, 1, M->r2) : 1 % M->n;
}

static uint64_t dl_pow(const DlMont* M, uint64_t b, uint64_t e) {
    uint64_t r = dl_one(M);
    while (e) {
        if (e & 1) r = dl_mul(M, r, b);
        e >>= 1;
        if (e) b = dl_mul(M, b, b);
    }
    return r;
}

static uint64_t dl_addmod(uint64_t a, uint64_t b, uint64_t q) {
    uint64_t s = a + b;                          /* a, b < q < 2^63 */
    return s >= q ? s - q : s;
}
#endif

/* =====================================================================
 * Per-prime solver context: gamma = k^(d/q) of prime order q.
 * ===================================================================*/

typedef struct {
    mpz_t    q, gamma;
    int      u64;        /* arithmetic back-end in use */
    int      built;      /* BSGS table present */
    uint64_t m;          /* baby steps */
    DlHash   H;
#ifdef DL_HAVE_U64
    DlMont   M;
    uint64_t giant;      /* gamma^-m */
#endif
    mpz_t    giant_z;    /* gamma^-m (mpz back-end) */
    mpz_t*   baby;       /* gamma^i, i < m (mpz back-end) */
} DlPrime;

static void dl_prime_init(DlPrime* P, const mpz_t q, const mpz_t gamma,
                          const mpz_t n) {
    mpz_init_set(P->q, q);
    mpz_init_set(P->gamma, gamma);
    mpz_init(P->giant_z);
    P->built = 0;
    P->baby = NULL;
    P->m = 0;
    P->u64 = 0;
#ifdef DL_HAVE_U64
    if (mpz_sizeinbase(n, 2) <= 63) {
        P->u64 = 1;
        dl_mont_init(&P->M, n);
    }
#else
    (void)n;
#endif
}

static void dl_prime_clear(DlPrime* P) {
    if (P->built) {
        dl_hash_free(&P->H);
        if (P->baby) {
            for (uint64_t i = 0; i < P->m; i++) mpz_clear(P->baby[i]);
            free(P->baby);
        }
    }
    mpz_clears(P->q, P->gamma, P->giant_z, NULL);
}

static bool dl_use_bsgs(const DlPrime* P) {
    return mpz_sizeinbase(P->q, 2) <= (P->u64 ? DL_BSGS_BITS : DL_BSGS_BITS_MPZ);
}

/* Build the baby-step table on first use. */
static DlStatus dl_bsgs_build(DlPrime* P, const mpz_t n) {
    if (P->built) return DL_FOUND;
    P->m = dl_isqrt_ceil(P->q);
    if (!dl_hash_init(&P->H, P->m)) return DL_ABORT;
    mpz_t t;
    mpz_init(t);
#ifdef DL_HAVE_U64
    if (P->u64) {
        const DlMont* M = &P->M;
        uint64_t g = dl_to(M, P->gamma), y = dl_one(M);
        for (uint64_t i = 0; i < P->m; i++) {
            dl_hash_put(&P->H, y, (uint32_t)i);
            y = dl_mul(M, y, g);
        }
        /* y == gamma^m; the giant step is its inverse. */
        dl_from(M, t, y);
        mpz_invert(t, t, n);
        P->giant = dl_to(M, t);
        mpz_clear(t);
        P->built = 1;
        return DL_FOUND;
    }
#endif
    P->baby = (mpz_t*)malloc((size_t)P->m * sizeof(mpz_t));
    if (!P->baby) { dl_hash_free(&P->H); mpz_clear(t); return DL_ABORT; }
    mpz_set_ui(t, 1);
    for (uint64_t i = 0; i < P->m; i++) {
        mpz_init_set(P->baby[i], t);
        dl_hash_put(&P->H, mpz_size(t) ? (uint64_t)mpz_getlimbn(t, 0) : 0,
                    (uint32_t)i);
        mpz_mul(t, t, P->gamma);
        mpz_mod(t, t, n);
    }
    mpz_invert(P->giant_z, t, n);
    mpz_clear(t);
    P->built = 1;
    return DL_FOUND;
}

/* Solve gamma^x == h by baby-step giant-step; h is reduced mod n. */
static DlStatus dl_bsgs(DlPrime* P, const mpz_t h, const mpz_t n, mpz_t x) {
    DlStatus st = dl_bsgs_build(P, n);
    if (st != DL_FOUND) return st;
    /* m^2 >= q, so j < m giant steps cover every exponent in [0, q). */
#ifdef DL_HAVE_U64
    if (P->u64) {
        const DlMont* M = &P->M;
        const DlHash* H = &P->H;
        uint64_t y = dl_to(M, h);
        for (uint64_t j = 0; j < P->m; j++) {
            for (uint64_t s = dl_hash_slot(H, y); H->val[s] != DL_EMPTY;
                 s = (s + 1) & H->mask) {
                if (H->key[s] == y) {
                    dl_set_u64(x, j * P->m + H->val[s]);
                    return DL_FOUND;
                }
            }
            y = dl_mul(M, y, P->giant);
            if ((j & DL_POLL_MASK) == DL_POLL_MASK && tc_deadline_passed())
                return DL_ABORT;
        }
        return DL_NONMEMBER;
    }
#endif
    const DlHash* H = &P->H;
    mpz_t y;
    mpz_init_set(y, h);
    for (uint64_t j = 0; j < P->m; j++) {
        uint64_t key = mpz_size(y) ? (uint64_t)mpz_getlimbn(y, 0) : 0;
        for (uint64_t s = dl_hash_slot(H, key); H->val[s] != DL_EMPTY;
             s = (s + 1) & H->mask) {
            if (H->key[s] == key && mpz_cmp(P->baby[H->val[s]], y) == 0) {
                mpz_set_ui(x, (unsigned long)P->m);
                mpz_mul_ui(x, x, (unsigned long)j);
                mpz_add_ui(x, x, (unsigned long)H->val[s]);
                mpz_clear(y);
                return DL_FOUND;
            }
        }
        mpz_mul(y, y, P->giant_z);
        mpz_mod(y, y, n);
        if ((j & DL_POLL_MASK) == DL_POLL_MASK && tc_deadline_passed()) {
            mpz_clear(y);
            return DL_ABORT;
        }
    }
    mpz_clear(y);
    return DL_NONMEMBER;
}

/* From a rho collision gamma^a1 h^b1 == gamma^a2 h^b2: x = (a2-a1)/(b1-b2)
 * mod q.  Returns 0 for a degenerate collision (b1 == b2), else 1 with x
 * set to the only candidate. */
static int dl_rho_solve(mpz_t x, const mpz_t a1, const mpz_t b1,
                        const mpz_t a2, const mpz_t b2, const mpz_t q) {
    mpz_t db;
    mpz_init(db);
    mpz_sub(db, b1, b2);
    mpz_mod(db, db, q);
    int ok = mpz_sgn(db) != 0 && mpz_invert(db, db, q);
    if (ok) {
        mpz_sub(x, a2, a1);
        mpz_mul(x, x, db);
        mpz_mod(x, x, q);
    }
    mpz_clear(db);
    return ok;
}


/* Bound on one rho walk when (Z/nZ)^* is not cyclic: 16 sqrt(q) steps, far
 * past the ~1.25 sqrt(q) a member needs.  Unbounded (0) in the cyclic case. */
static uint64_t dl_rho_cap(const mpz_t q, int cyclic) {
    if (cyclic || mpz_sizeinbase(q, 2) > 120) return 0;
    mpz_t r;
    mpz_init(r);
    mpz_sqrt(r, q);
    uint64_t c = (uint64_t)mpz_get_d(r) * 16 + 1024;
    mpz_clear(r);
    return c;
}

#ifdef DL_HAVE_U64
/* Pollard rho on the u64 back-end (q < n < 2^63). */
static DlStatus dl_rho_u64(DlPrime* P, const mpz_t hz, const mpz_t n,
                           int cyclic, mpz_t x) {
    const DlMont* M = &P->M;
    uint64_t q = dl_get_u64(P->q);
    uint64_t g = dl_to(M, P->gamma), h = dl_to(M, hz);
    uint64_t seed = q ^ UINT64_C(0x2545F4914F6CDD1D);
    uint64_t cap = dl_rho_cap(P->q, cyclic);
    uint64_t ma[DL_RHO_R], mb[DL_RHO_R], mul[DL_RHO_R];
    mpz_t a1, b1, a2, b2, chk;
    mpz_inits(a1, b1, a2, b2, chk, NULL);
    DlStatus st = DL_NONMEMBER;
    for (int attempt = 0; cyclic || attempt < DL_RHO_ATTEMPTS; attempt++) {
        for (int s = 0; s < DL_RHO_R; s++) {
            ma[s] = dl_rand(&seed) % q;
            mb[s] = dl_rand(&seed) % q;
            mul[s] = dl_mul(M, dl_pow(M, g, ma[s]), dl_pow(M, h, mb[s]));
        }
        uint64_t ta = dl_rand(&seed) % q, tb = dl_rand(&seed) % q;
        uint64_t ty = dl_mul(M, dl_pow(M, g, ta), dl_pow(M, h, tb));
        uint64_t hy = ty, ha = ta, hb = tb;
        uint64_t power = 1, lam = 0, steps = 0;
        /* Brent: the tortoise parks at each power of two, the hare walks. */
        for (;;) {
            unsigned s = (unsigned)((hy * UINT64_C(0x9E3779B97F4A7C15)) >> 59);
            hy = dl_mul(M, hy, mul[s]);
            ha = dl_addmod(ha, ma[s], q);
            hb = dl_addmod(hb, mb[s], q);
            lam++;
            steps++;
            if (hy == ty) break;
            if (lam == power) { ty = hy; ta = ha; tb = hb; power <<= 1; lam = 0; }
            if ((steps & DL_POLL_MASK) == 0 && tc_deadline_passed()) {
                st = DL_ABORT;
                goto done;
            }
            if (cap && steps > cap) break;
        }
        if (hy != ty) continue;                  /* capped walk */
        dl_set_u64(a1, ta); dl_set_u64(b1, tb);
        dl_set_u64(a2, ha); dl_set_u64(b2, hb);
        if (!dl_rho_solve(x, a1, b1, a2, b2, P->q)) continue;   /* degenerate */
        mpz_powm(chk, P->gamma, x, n);
        st = mpz_cmp(chk, hz) == 0 ? DL_FOUND : DL_NONMEMBER;
        goto done;
    }
done:
    mpz_clears(a1, b1, a2, b2, chk, NULL);
    return st;
}
#endif

/* Pollard rho on the mpz back-end: the same walk with GMP state. */
static DlStatus dl_rho_mpz(DlPrime* P, const mpz_t h, const mpz_t n,
                           int cyclic, mpz_t x) {
    uint64_t seed = (uint64_t)mpz_getlimbn(P->q, 0) ^ UINT64_C(0x2545F4914F6CDD1D);
    if (!seed) seed = 1;
    uint64_t cap = dl_rho_cap(P->q, cyclic);
    mpz_t ma[DL_RHO_R], mb[DL_RHO_R], mul[DL_RHO_R];
    mpz_t ty, ta, tb, hy, ha, hb, t;
    for (int s = 0; s < DL_RHO_R; s++) mpz_inits(ma[s], mb[s], mul[s], NULL);
    mpz_inits(ty, ta, tb, hy, ha, hb, t, NULL);
    DlStatus st = DL_NONMEMBER;
    for (int attempt = 0; cyclic || attempt < DL_RHO_ATTEMPTS; attempt++) {
        for (int s = 0; s < DL_RHO_R; s++) {
            dl_rand_mpz(ma[s], P->q, &seed);
            dl_rand_mpz(mb[s], P->q, &seed);
            mpz_powm(mul[s], P->gamma, ma[s], n);
            mpz_powm(t, h, mb[s], n);
            mpz_mul(mul[s], mul[s], t);
            mpz_mod(mul[s], mul[s], n);
        }
        dl_rand_mpz(ta, P->q, &seed);
        dl_rand_mpz(tb, P->q, &seed);
        mpz_powm(ty, P->gamma, ta, n);
        mpz_powm(t, h, tb, n);
        mpz_mul(ty, ty, t);
        mpz_mod(ty, ty, n);
        mpz_set(hy, ty); mpz_set(ha, ta); mpz_set(hb, tb);
        uint64_t power = 1, lam = 0, steps = 0;
        int met = 0;
        for (;;) {
            uint64_t key = mpz_size(hy) ? (uint64_t)mpz_getlimbn(hy, 0) : 0;
            unsigned s = (unsigned)((key * UINT64_C(0x9E3779B97F4A7C15)) >> 59);
            mpz_mul(hy, hy, mul[s]);
            mpz_mod(hy, hy, n);
            mpz_add(ha, ha, ma[s]);
            if (mpz_cmp(ha, P->q) >= 0) mpz_sub(ha, ha, P->q);
            mpz_add(hb, hb, mb[s]);
            if (mpz_cmp(hb, P->q) >= 0) mpz_sub(hb, hb, P->q);
            lam++;
            steps++;
            if (mpz_cmp(hy, ty) == 0) { met = 1; break; }
            if (lam == power) {
                mpz_set(ty, hy); mpz_set(ta, ha); mpz_set(tb, hb);
                power <<= 1;
                lam = 0;
            }
            if ((steps & DL_POLL_MASK) == 0 && tc_deadline_passed()) {
                st = DL_ABORT;
                goto done;
            }
            if (cap && steps > cap) break;
        }
        if (!met) continue;
        if (!dl_rho_solve(x, ta, tb, ha, hb, P->q)) continue;
        mpz_powm(t, P->gamma, x, n);
        st = mpz_cmp(t, h) == 0 ? DL_FOUND : DL_NONMEMBER;
        goto done;
    }
done:
    for (int s = 0; s < DL_RHO_R; s++) mpz_clears(ma[s], mb[s], mul[s], NULL);
    mpz_clears(ty, ta, tb, hy, ha, hb, t, NULL);
    return st;
}

/* gamma^x == h for h of order dividing q: x in [0, q). */
static DlStatus dl_solve_prime(DlPrime* P, const mpz_t h, const mpz_t n,
                               int cyclic, mpz_t x) {
    if (mpz_cmp_ui(h, 1) == 0) { mpz_set_ui(x, 0); return DL_FOUND; }
    if (dl_use_bsgs(P)) return dl_bsgs(P, h, n, x);
#ifdef DL_HAVE_U64
    if (P->u64) return dl_rho_u64(P, h, n, cyclic, x);
#endif
    return dl_rho_mpz(P, h, n, cyclic, x);
}

/* =====================================================================
 * Pohlig-Hellman driver
 * ===================================================================*/

bool nt_dlog_min(const mpz_t k_in, const mpz_t n, const mpz_t d,
                 const mpz_t* q, const unsigned long* e, size_t nq,
                 int cyclic, const mpz_t* targets, size_t nt, mpz_t m_out) {
    if (mpz_cmp_ui(n, 1) == 0) {                 /* every residue is 0 == 1 */
        mpz_set_ui(m_out, 1);
        return nt > 0;
    }
    mpz_t k, kinv, t, X, Mod, xq, qpow, ex, hh, tmp, delta, best;
    mpz_inits(k, kinv, t, X, Mod, xq, qpow, ex, hh, tmp, delta, best, NULL);
    mpz_mod(k, k_in, n);
    mpz_invert(kinv, k, n);                      /* k is a unit (caller) */

    DlPrime* ctx = nq ? (DlPrime*)malloc(nq * sizeof(DlPrime)) : NULL;
    char* live = nq ? (char*)calloc(nq, 1) : NULL;
    bool found = false, aborted = (nq && (!ctx || !live));

    for (size_t ti = 0; ti < nt && !aborted; ti++) {
        mpz_mod(t, targets[ti], n);
        mpz_gcd(tmp, t, n);
        if (mpz_cmp_ui(tmp, 1) != 0) continue;   /* not a unit: never a power */
        mpz_powm(tmp, t, d, n);
        if (mpz_cmp_ui(tmp, 1) != 0) continue;   /* order does not divide d */

        mpz_set_ui(X, 0);
        mpz_set_ui(Mod, 1);
        bool member = true;
        for (size_t j = 0; j < nq && member; j++) {
            if (!live[j]) {                      /* gamma = k^(d/q), order q */
                mpz_divexact(ex, d, q[j]);
                mpz_powm(tmp, k, ex, n);
                dl_prime_init(&ctx[j], q[j], tmp, n);
                live[j] = 1;
            }
            /* Digit i of x mod q^e: (t k^-xq)^(d/q^(i+1)) == gamma^digit. */
            mpz_set_ui(xq, 0);
            mpz_set_ui(qpow, 1);
            mpz_divexact(ex, d, q[j]);
            for (unsigned long i = 0; i < e[j]; i++) {
                mpz_powm(tmp, kinv, xq, n);
                mpz_mul(tmp, tmp, t);
                mpz_mod(tmp, tmp, n);
                mpz_powm(hh, tmp, ex, n);
                if (!cyclic) {                   /* hh must have order | q */
                    mpz_powm(tmp, hh, q[j], n);
                    if (mpz_cmp_ui(tmp, 1) != 0) { member = false; break; }
                }
                DlStatus st = dl_solve_prime(&ctx[j], hh, n, cyclic, delta);
                if (st == DL_ABORT) { aborted = true; member = false; break; }
                if (st == DL_NONMEMBER) { member = false; break; }
                mpz_addmul(xq, delta, qpow);
                mpz_mul(qpow, qpow, q[j]);
                if (i + 1 < e[j]) mpz_divexact(ex, ex, q[j]);
            }
            if (!member) break;
            /* CRT: X == previous X (mod Mod), X == xq (mod qpow). */
            mpz_sub(tmp, xq, X);
            mpz_invert(hh, Mod, qpow);
            mpz_mul(tmp, tmp, hh);
            mpz_mod(tmp, tmp, qpow);
            mpz_addmul(X, tmp, Mod);
            mpz_mul(Mod, Mod, qpow);
        }
        if (!member) continue;
        mpz_mod(X, X, d);
        if (mpz_sgn(X) == 0) mpz_set(X, d);      /* smallest POSITIVE exponent */
        mpz_powm(tmp, k, X, n);
        if (mpz_cmp(tmp, t) != 0) continue;      /* certificate */
        if (!found || mpz_cmp(X, best) < 0) mpz_set(best, X);
        found = true;
    }

    for (size_t j = 0; j < nq; j++) if (live && live[j]) dl_prime_clear(&ctx[j]);
    free(ctx);
    free(live);
    if (found && !aborted) mpz_set(m_out, best);
    mpz_clears(k, kinv, t, X, Mod, xq, qpow, ex, hh, tmp, delta, best, NULL);
    return found && !aborted;
}
