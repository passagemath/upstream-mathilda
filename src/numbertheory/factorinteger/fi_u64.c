/*
 * fi_u64.c — the word-size factoring engine.
 *
 * Everything below 2^64 is factored here without touching GMP: Montgomery
 * multiplication on single words, a deterministic Miller–Rabin test, and the
 * splitters Pollard–Brent rho, SQUFOF, Lehman and Fermat, plus a micro-ECM
 * (Montgomery curves with Suyama's parametrisation, so 12 divides every group
 * order; prime-power ladder stage 1; baby-step/giant-step stage 2).  Micro-ECM
 * is the Automatic splitter for 40..64-bit composites: a 32-bit factor costs
 * rho ~2^16 modular multiplications, micro-ECM a handful of ~3000-multiply
 * curves.
 *
 * Montgomery convention: R = 2^64, n odd, ninv = n^-1 mod 2^64.  A residue x
 * is held as xR mod n.  mulredc(a, b) = a b R^-1 mod n for a, b < n, computed
 * as hi(ab) - hi(m n) with m = lo(ab) ninv: the low words cancel exactly, so
 * there is no carry to track and it is valid for every odd n < 2^64.
 */
#include "fi_internal.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

/* ------------------------------------------------------------------ */
/* 64x64 -> 128 multiply                                               */
/* ------------------------------------------------------------------ */

#if defined(__SIZEOF_INT128__)
__extension__ typedef unsigned __int128 fi_uint128;
void fi_mul64(uint64_t a, uint64_t b, uint64_t* hi, uint64_t* lo) {
    fi_uint128 p = (fi_uint128)a * b;
    *hi = (uint64_t)(p >> 64);
    *lo = (uint64_t)p;
}
static inline uint64_t mulhi(uint64_t a, uint64_t b) {
    return (uint64_t)(((fi_uint128)a * b) >> 64);
}
#else
void fi_mul64(uint64_t a, uint64_t b, uint64_t* hi, uint64_t* lo) {
    uint64_t al = a & 0xffffffffu, ah = a >> 32, bl = b & 0xffffffffu, bh = b >> 32;
    uint64_t ll = al * bl, lh = al * bh, hl = ah * bl, hh = ah * bh;
    uint64_t mid = (ll >> 32) + (lh & 0xffffffffu) + (hl & 0xffffffffu);
    *lo = (ll & 0xffffffffu) | (mid << 32);
    *hi = hh + (lh >> 32) + (hl >> 32) + (mid >> 32);
}
static inline uint64_t mulhi(uint64_t a, uint64_t b) {
    uint64_t h, l;
    fi_mul64(a, b, &h, &l);
    return h;
}
#endif

/* ------------------------------------------------------------------ */
/* Montgomery arithmetic                                               */
/* ------------------------------------------------------------------ */

typedef struct { uint64_t n, ninv, one, r2; } Mont;

static inline uint64_t inv64(uint64_t n) {
    uint64_t x = n;                     /* 3 correct bits: n*n == 1 (mod 8) */
    for (int i = 0; i < 5; i++) x *= 2 - n * x;
    return x;
}
static inline uint64_t addm(uint64_t a, uint64_t b, uint64_t n) {
    uint64_t s = a + b;
    if (s < a || s >= n) s -= n;
    return s;
}
static inline uint64_t subm(uint64_t a, uint64_t b, uint64_t n) {
    return a >= b ? a - b : a - b + n;
}
static inline uint64_t mulredc(uint64_t a, uint64_t b, const Mont* M) {
    uint64_t th, tl;
    fi_mul64(a, b, &th, &tl);
    uint64_t mh = mulhi(tl * M->ninv, M->n);
    uint64_t r = th - mh;
    if (th < mh) r += M->n;
    return r;
}
static void mont_setup(Mont* M, uint64_t n) {
    M->n = n;
    M->ninv = inv64(n);
    M->one = (0 - n) % n;               /* R mod n */
    uint64_t x = M->one;                /* R^2 mod n = R * 2^64: 64 doublings */
    for (int i = 0; i < 64; i++) x = addm(x, x, n);
    M->r2 = x;
}
static inline uint64_t to_mont(uint64_t a, const Mont* M) { return mulredc(a % M->n, M->r2, M); }
static inline uint64_t from_mont(uint64_t a, const Mont* M) { return mulredc(a, 1, M); }

static uint64_t mont_pow(uint64_t b, uint64_t e, const Mont* M) {
    uint64_t r = M->one;
    while (e) {
        if (e & 1) r = mulredc(r, b, M);
        b = mulredc(b, b, M);
        e >>= 1;
    }
    return r;
}

/* ------------------------------------------------------------------ */
/* Small helpers                                                       */
/* ------------------------------------------------------------------ */

static inline int ctz64(uint64_t x) {
#if defined(__GNUC__)
    return __builtin_ctzll(x);
#else
    int c = 0;
    while (!(x & 1)) { x >>= 1; c++; }
    return c;
#endif
}

static inline int clz64(uint64_t x) {
#if defined(__GNUC__)
    return __builtin_clzll(x);
#else
    int c = 0;
    while (!(x >> 63)) { x <<= 1; c++; }
    return c;
#endif
}

uint64_t fi_u64_gcd(uint64_t a, uint64_t b) {
    if (!a) return b;
    if (!b) return a;
    int s = ctz64(a | b);
    a >>= ctz64(a);
    do {
        b >>= ctz64(b);
        if (a > b) { uint64_t t = a; a = b; b = t; }
        b -= a;
    } while (b);
    return a << s;
}

uint64_t fi_u64_isqrt(uint64_t n) {
    uint64_t r = (uint64_t)sqrt((double)n);
    /* The double estimate can be off by one either way near 2^64. */
    while (r > 0xffffffffu || r * r > n) r--;
    while ((r + 1) <= 0xffffffffu && (r + 1) * (r + 1) <= n) r++;
    return r;
}

/* Perfect-square test with quadratic-residue filters mod 64, 63, 65, 11 (the
 * classic Cohen filter: rejects ~99.4% of non-squares without a sqrt). */
static unsigned char sq64[64], sq63[63], sq65[65], sq11[11];
static int sq_ready = 0;
static void sq_init(void) {
    for (int i = 0; i < 64; i++) sq64[(i * i) % 64] = 1;
    for (int i = 0; i < 63; i++) sq63[(i * i) % 63] = 1;
    for (int i = 0; i < 65; i++) sq65[(i * i) % 65] = 1;
    for (int i = 0; i < 11; i++) sq11[(i * i) % 11] = 1;
    sq_ready = 1;
}
static inline int is_square_u64(uint64_t x, uint64_t* root) {
    if (!sq64[x & 63]) return 0;
    uint64_t r = x % 45045;             /* 63 * 65 * 11 */
    if (!sq63[r % 63] || !sq65[r % 65] || !sq11[r % 11]) return 0;
    uint64_t s = fi_u64_isqrt(x);
    if (s * s != x) return 0;
    if (root) *root = s;
    return 1;
}

/* ------------------------------------------------------------------ */
/* Primality                                                           */
/* ------------------------------------------------------------------ */

static int mr_witness(uint64_t a, uint64_t d, int s, const Mont* M) {
    uint64_t am = to_mont(a, M);
    if (am == 0) return 1;              /* a == 0 mod n: not a witness */
    uint64_t x = mont_pow(am, d, M);
    uint64_t mone = M->n - M->one;
    if (x == M->one || x == mone) return 1;
    for (int i = 1; i < s; i++) {
        x = mulredc(x, x, M);
        if (x == mone) return 1;
        if (x == M->one) return 0;
    }
    return 0;
}

static const uint8_t SMALL_P[] = { 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53 };

int fi_u64_is_prime(uint64_t n) {
    if (n < 2) return 0;
    if (!(n & 1)) return n == 2;
    for (size_t i = 0; i < sizeof SMALL_P; i++) {
        if (n == SMALL_P[i]) return 1;
        if (n % SMALL_P[i] == 0) return 0;
    }
    if (n < 59 * 59) return 1;
    Mont M;
    mont_setup(&M, n);
    uint64_t d = n - 1;
    int s = ctz64(d);
    d >>= s;
    /* Deterministic: {2,7,61} below 2^32 (Jaeschke); the 7-base set of
     * Sinclair for the full 64-bit range. */
    if (n < UINT64_C(4294967296)) {
        static const uint64_t B[] = { 2, 7, 61 };
        for (int i = 0; i < 3; i++) if (!mr_witness(B[i], d, s, &M)) return 0;
        return 1;
    }
    static const uint64_t B[] = { 2, 325, 9375, 28178, 450775, 9780504, 1795265022 };
    for (int i = 0; i < 7; i++) if (!mr_witness(B[i], d, s, &M)) return 0;
    return 1;
}

/* ------------------------------------------------------------------ */
/* Pollard–Brent rho                                                   */
/* ------------------------------------------------------------------ */

/* f(y) = y^2 + c in the Montgomery domain.  Products of |x - y| are batched
 * 128 at a time between gcds; on a gcd of n the last batch is replayed one
 * step at a time from the saved ys (Brent's backtrack). */
uint64_t fi_u64_rho(uint64_t n, uint64_t c, uint64_t max_iters) {
    if (!(n & 1)) return 2;
    Mont M;
    mont_setup(&M, n);
    uint64_t cm = to_mont(c ? c : 1, &M);
    uint64_t y = to_mont(2, &M), x = y, ys = y, q = M.one, g = 1;
    const uint64_t m = 128;
    uint64_t done = 0;
    for (uint64_t r = 1; g == 1 && done < max_iters; r <<= 1) {
        x = y;
        for (uint64_t i = 0; i < r; i++) y = addm(mulredc(y, y, &M), cm, n);
        done += r;
        for (uint64_t k = 0; k < r && g == 1; k += m) {
            ys = y;
            uint64_t lim = (r - k < m) ? r - k : m;
            for (uint64_t i = 0; i < lim; i++) {
                y = addm(mulredc(y, y, &M), cm, n);
                q = mulredc(q, x > y ? x - y : y - x, &M);
            }
            done += lim;
            g = fi_u64_gcd(q, n);
        }
    }
    if (g == n || g == 0) {
        g = 1;
        for (uint64_t i = 0; i < 4 * m && g == 1; i++) {
            ys = addm(mulredc(ys, ys, &M), cm, n);
            g = fi_u64_gcd(x > ys ? x - ys : ys - x, n);
        }
    }
    return (g > 1 && g < n) ? g : 0;
}

/* ------------------------------------------------------------------ */
/* SQUFOF (Shanks), sequential multipliers                             */
/* ------------------------------------------------------------------ */

static uint64_t squfof_one(uint64_t n, uint64_t k, uint64_t max_iters) {
    uint64_t hi, kn;
    fi_mul64(n, k, &hi, &kn);
    if (hi || kn >> 62) return 0;       /* need kN < 2^62 for the form bounds */
    uint64_t s = fi_u64_isqrt(kn);
    if (s * s == kn) {
        uint64_t g = fi_u64_gcd(n, s);
        return (g > 1 && g < n) ? g : 0;
    }
    uint64_t P = s, Qp = 1, Q = kn - s * s;
    uint64_t lim = max_iters ? max_iters : 4 * (uint64_t)sqrt(2.0 * sqrt((double)kn)) + 64;
    uint64_t r = 0;
    uint64_t i;
    for (i = 1; i < lim; i++) {
        uint64_t b = (s + P) / Q;
        uint64_t Pn = b * Q - P;
        /* Q' = Qprev + b (P - P'); P - P' may be negative: signed. */
        int64_t dq = (int64_t)Qp + (int64_t)b * ((int64_t)P - (int64_t)Pn);
        Qp = Q; Q = (uint64_t)dq; P = Pn;
        /* After step i, Q holds Q_{i+1}: test the even-indexed forms. */
        if ((i & 1) && is_square_u64(Q, &r)) break;
    }
    if (i >= lim) return 0;
    /* Reverse cycle from the square root form. */
    uint64_t b = (s - P) / r;
    P = b * r + P;
    Qp = r;
    Q = (kn - P * P) / Qp;
    for (uint64_t j = 0; j < lim; j++) {
        uint64_t bb = (s + P) / Q;
        uint64_t Pn = bb * Q - P;
        if (Pn == P) break;
        int64_t dq = (int64_t)Qp + (int64_t)bb * ((int64_t)P - (int64_t)Pn);
        Qp = Q; Q = (uint64_t)dq; P = Pn;
    }
    uint64_t g = fi_u64_gcd(n, P);
    return (g > 1 && g < n) ? g : 0;
}

/* Multipliers of Gower & Wagstaff (squarefree products of 3, 5, 7, 11). */
static const uint16_t SQUFOF_MULT[] = {
    1, 3, 5, 7, 11, 15, 21, 33, 35, 55, 77, 105, 165, 231, 385, 1155
};

uint64_t fi_u64_squfof(uint64_t n, uint64_t max_iters) {
    if (!sq_ready) sq_init();
    if (!(n & 1)) return 2;
    for (size_t i = 0; i < sizeof SQUFOF_MULT / sizeof SQUFOF_MULT[0]; i++) {
        uint64_t g = squfof_one(n, SQUFOF_MULT[i], max_iters);
        if (g) return g;
    }
    return 0;
}

/* ------------------------------------------------------------------ */
/* Fermat and Lehman                                                   */
/* ------------------------------------------------------------------ */

uint64_t fi_u64_fermat(uint64_t n, uint64_t max_iters) {
    if (!sq_ready) sq_init();
    if (!(n & 1)) return 2;
    if (n >> 62) return 0;              /* a^2 must fit a word */
    uint64_t a = fi_u64_isqrt(n);
    if (a * a == n) return a;
    a++;
    uint64_t b2 = a * a - n, r;
    for (uint64_t i = 0; i < max_iters && a <= n / 2; i++) {
        if (is_square_u64(b2, &r)) {
            uint64_t d = a - r;
            return (d > 1 && d < n) ? d : 0;
        }
        b2 += 2 * a + 1;
        a++;
    }
    return 0;
}

/* Lehman (1974): O(n^(1/3)).  Trial division to n^(1/3), then for each
 * k <= n^(1/3) a short window of a with a^2 - 4kn square.  Requires
 * n < 2^42 so 4 k n stays below 2^63.  Trial division uses the primes with
 * precomputed inverses (n divisible by p <=> n p^-1 mod 2^64 <= floor(
 * (2^64-1)/p)), so it costs a multiply per prime instead of a division. */
#define LEH_PMAX 16384u
static uint32_t leh_p[2000];
static uint64_t leh_inv[2000], leh_lim[2000];
static size_t leh_n = 0;
static void leh_init(void) {
    size_t np;
    const uint32_t* P = fi_small_primes(&np);
    for (size_t i = 1; i < np && P[i] < LEH_PMAX; i++) {
        leh_p[leh_n] = P[i];
        leh_inv[leh_n] = inv64(P[i]);
        leh_lim[leh_n] = UINT64_MAX / P[i];
        leh_n++;
    }
}

uint64_t fi_u64_lehman(uint64_t n) {
    if (!sq_ready) sq_init();
    if (!leh_n) leh_init();
    if (!(n & 1)) return 2;
    if (n >> 42) return 0;
    uint64_t c = (uint64_t)cbrt((double)n) + 1;
    for (size_t i = 0; i < leh_n && leh_p[i] <= c; i++)
        if (n * leh_inv[i] <= leh_lim[i]) return leh_p[i] == n ? 0 : leh_p[i];
    double sixth = pow((double)n, 1.0 / 6.0) * 0.25;
    double dn4 = 4.0 * (double)n;
    for (uint64_t k = 1; k <= c; k++) {
        uint64_t fkn = 4 * k * n;
        double sk = sqrt(dn4 * (double)k);
        uint64_t a0 = (uint64_t)sk;
        if (a0 * a0 < fkn) a0++;
        else if ((a0 - 1) * (a0 - 1) >= fkn) a0--;
        uint64_t a1 = (uint64_t)(sk + sixth / sqrt((double)k)) + 1;
        /* a must have the parity making a^2 - 4kn a square: for odd k, a is
         * even unless ...; checking both parities is cheap with the filter. */
        uint64_t b2 = a0 * a0 - fkn;
        for (uint64_t a = a0; a <= a1; a++) {
            uint64_t b;
            if (is_square_u64(b2, &b)) {
                uint64_t g = fi_u64_gcd(a + b, n);
                if (g > 1 && g < n) return g;
            }
            b2 += 2 * a + 1;
        }
    }
    return 0;
}

/* ------------------------------------------------------------------ */
/* Micro-ECM                                                           */
/* ------------------------------------------------------------------ */

typedef struct { uint64_t X, Z; } Pt;

static inline Pt ec_dbl(Pt P, uint64_t a24, const Mont* M) {
    uint64_t n = M->n;
    uint64_t s = addm(P.X, P.Z, n), d = subm(P.X, P.Z, n);
    uint64_t s2 = mulredc(s, s, M), d2 = mulredc(d, d, M);
    uint64_t t = subm(s2, d2, n);
    Pt R;
    R.X = mulredc(s2, d2, M);
    R.Z = mulredc(t, addm(d2, mulredc(a24, t, M), n), M);
    return R;
}

/* P + Q given D = P - Q. */
static inline Pt ec_add(Pt P, Pt Q, Pt D, const Mont* M) {
    uint64_t n = M->n;
    uint64_t u = mulredc(subm(P.X, P.Z, n), addm(Q.X, Q.Z, n), M);
    uint64_t v = mulredc(addm(P.X, P.Z, n), subm(Q.X, Q.Z, n), M);
    uint64_t sp = addm(u, v, n), sm = subm(u, v, n);
    Pt R;
    R.X = mulredc(D.Z, mulredc(sp, sp, M), M);
    R.Z = mulredc(D.X, mulredc(sm, sm, M), M);
    return R;
}

/* Montgomery ladder: k P for k >= 1. */
static Pt ec_mul(Pt P, uint64_t k, uint64_t a24, const Mont* M) {
    if (k == 1) return P;
    Pt R0 = P, R1 = ec_dbl(P, a24, M);
    int top = 63 - clz64(k);
    for (int b = top - 1; b >= 0; b--) {
        if ((k >> b) & 1) { R0 = ec_add(R1, R0, P, M); R1 = ec_dbl(R1, a24, M); }
        else              { R1 = ec_add(R1, R0, P, M); R0 = ec_dbl(R0, a24, M); }
    }
    return R0;
}

/* Inverse of a mod n (plain integers) by extended Euclid; returns 0 and sets
 * *g to gcd(a, n) when it is not 1. */
static uint64_t inv_mod(uint64_t a, uint64_t n, uint64_t* g) {
    int64_t t0 = 0, t1 = 1;
    uint64_t r0 = n, r1 = a % n;
    while (r1) {
        uint64_t qq = r0 / r1, r2 = r0 - qq * r1;
        int64_t t2 = t0 - (int64_t)qq * t1;
        r0 = r1; r1 = r2; t0 = t1; t1 = t2;
    }
    *g = r0;
    if (r0 != 1) return 0;
    return t0 < 0 ? (uint64_t)(t0 + (int64_t)n) : (uint64_t)t0;
}

/* Prime powers <= B1 packed into word-sized multipliers, and a prime bitmap
 * for stage 2, built once. */
#define UECM_MAXB2 65536
static unsigned char uecm_isprime[UECM_MAXB2 + 1];
static int uecm_ready = 0;
static void uecm_init(void);

/* Prime bitmap up to 65536, shared with the two-limb micro-ECM. */
const unsigned char* fi_prime_bitmap(uint32_t* limit) {
    if (!uecm_ready) uecm_init();
    *limit = UECM_MAXB2;
    return uecm_isprime;
}

static void uecm_init(void) {
    memset(uecm_isprime, 1, sizeof uecm_isprime);
    uecm_isprime[0] = uecm_isprime[1] = 0;
    for (uint32_t i = 2; i * i <= UECM_MAXB2; i++)
        if (uecm_isprime[i])
            for (uint32_t j = i * i; j <= UECM_MAXB2; j += i) uecm_isprime[j] = 0;
    uecm_ready = 1;
}

#define UECM_D 210
static const uint8_t UECM_J[] = {   /* j < D/2 with gcd(j, D) = 1 */
    1, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73,
    79, 83, 89, 97, 101, 103
};
#define UECM_NJ (sizeof UECM_J / sizeof UECM_J[0])

/* One curve.  Returns a factor, 0 for none, or n when every prime of n was
 * caught at once (caller moves on to the next curve). */
static uint64_t uecm_curve(uint64_t sigma, uint32_t B1, uint32_t B2, const Mont* M) {
    uint64_t n = M->n, g;
    /* Suyama: u = s^2 - 5, v = 4 s, x0 = u^3, z0 = v^3,
     * (A + 2)/4 = (v - u)^3 (3u + v) / (16 u^3 v). */
    uint64_t s = to_mont(sigma, M);
    uint64_t u = subm(mulredc(s, s, M), to_mont(5, M), n);
    uint64_t v = to_mont((4 * sigma) % n, M);
    uint64_t u3 = mulredc(mulredc(u, u, M), u, M);
    uint64_t v3 = mulredc(mulredc(v, v, M), v, M);
    uint64_t vmu = subm(v, u, n);
    uint64_t num = mulredc(mulredc(mulredc(vmu, vmu, M), vmu, M),
                           addm(addm(addm(u, u, n), u, n), v, n), M);
    uint64_t den = mulredc(mulredc(to_mont(16, M), u3, M), v, M);
    uint64_t di = inv_mod(from_mont(den, M), n, &g);
    if (!di) return (g > 1 && g < n) ? g : 0;
    uint64_t a24 = mulredc(num, to_mont(di, M), M);

    Pt P = { u3, v3 };
    /* Stage 1: multiply by every prime power <= B1, packed into words. */
    uint64_t acc = 1;
    for (uint32_t p = 2; p <= B1; p++) {
        if (!uecm_isprime[p]) continue;
        uint64_t q = p;
        while (q * p <= B1) q *= p;
        uint64_t h, l;
        fi_mul64(acc, q, &h, &l);
        if (h) { P = ec_mul(P, acc, a24, M); acc = q; }
        else acc = l;
    }
    P = ec_mul(P, acc, a24, M);
    g = fi_u64_gcd(from_mont(P.Z, M), n);
    if (g != 1) return g;               /* g == n: caller retries */

    /* Stage 2: primes q in (B1, B2] as q = m D +- j. */
    Pt baby[UECM_NJ];
    Pt P2 = ec_dbl(P, a24, M);
    {
        Pt prev = P, cur = ec_add(P2, P, P, M);  /* 3P */
        size_t bi = 0;
        baby[bi++] = P;                          /* j = 1 */
        for (uint32_t j = 3; j < UECM_D / 2 && bi < UECM_NJ; j += 2) {
            if (j == UECM_J[bi]) baby[bi++] = cur;
            Pt nx = ec_add(cur, P2, prev, M);    /* (j+2)P */
            prev = cur; cur = nx;
        }
    }
    Pt G = ec_mul(P, UECM_D, a24, M);
    /* First giant step whose window [mD - D/2, mD + D/2] reaches past B1. */
    uint32_t m0 = (B1 + UECM_D / 2) / UECM_D;
    if (m0 < 1) m0 = 1;
    Pt Mc = ec_mul(G, m0, a24, M);
    Pt Mp = (m0 > 1) ? ec_mul(G, m0 - 1, a24, M) : G;   /* (m-1)G: the add difference */
    uint64_t prod = M->one;
    for (uint32_t m = m0; (uint64_t)m * UECM_D - UECM_D / 2 <= B2; m++) {
        for (size_t k = 0; k < UECM_NJ; k++) {
            uint32_t j = UECM_J[k];
            uint32_t q1 = m * UECM_D - j, q2 = m * UECM_D + j;
            int hit = (q1 > B1 && q1 <= B2 && uecm_isprime[q1]) ||
                      (q2 > B1 && q2 <= B2 && uecm_isprime[q2]);
            if (hit) {
                uint64_t t = subm(mulredc(Mc.X, baby[k].Z, M), mulredc(baby[k].X, Mc.Z, M), n);
                prod = mulredc(prod, t, M);
            }
        }
        Pt Mn;
        if (m == 1) Mn = ec_dbl(Mc, a24, M);       /* 2G from G */
        else Mn = ec_add(Mc, G, Mp, M);
        Mp = Mc; Mc = Mn;
    }
    g = fi_u64_gcd(from_mont(prod, M), n);
    return g;
}

uint64_t fi_u64_ecm(uint64_t n, uint64_t* rng, int max_curves) {
    if (!uecm_ready) uecm_init();
    if (!(n & 1)) return 2;
    int bits = 64 - clz64(n);
    uint32_t B1 = bits <= 44 ? 47 : bits <= 52 ? 85 : bits <= 58 ? 125 : 165;
    uint32_t B2 = 25 * B1;
    Mont M;
    mont_setup(&M, n);
    for (int c = 0; c < max_curves; c++) {
        /* Deterministic sigma sequence; xorshift keeps it reproducible per n. */
        uint64_t x = *rng;
        x ^= x << 13; x ^= x >> 7; x ^= x << 17;
        *rng = x;
        uint64_t sigma = 6 + (x % (UINT64_C(1) << 30));
        uint64_t g = uecm_curve(sigma, B1, B2, &M);
        if (g > 1 && g < n) return g;
    }
    return 0;
}

/* ------------------------------------------------------------------ */
/* Pollard p-1 and Williams p+1 on a word                              */
/* ------------------------------------------------------------------ */

/* Both walk the prime powers q <= B1 in batches of 32 with one gcd per
 * batch; when a batch overshoots (gcd = n: every prime caught at once) it is
 * replayed one prime power at a time, then one prime at a time, from the
 * state saved before the batch.  Stage 1 only: on a word B1 = 10^6 already
 * covers what a stage 2 would add for most inputs. */

/* x <- step(x, q): a^q for p-1, V_q(V) for p+1. */
static uint64_t pm1_step(uint64_t a, uint64_t q, const Mont* M, uint64_t aux) {
    (void)aux;
    return mont_pow(a, q, M);
}

/* Lucas V_k(v) by the binary ladder (V_{2k} = V_k^2 - 2,
 * V_{2k+1} = V_k V_{k+1} - V_1); `two` is 2 in Montgomery form. */
static uint64_t pp1_step(uint64_t v, uint64_t k, const Mont* M, uint64_t two) {
    uint64_t n = M->n;
    uint64_t a = v, b = subm(mulredc(v, v, M), two, n);    /* V_1, V_2 */
    int top = 63 - clz64(k);
    for (int i = top - 1; i >= 0; i--) {
        if ((k >> i) & 1) { a = subm(mulredc(a, b, M), v, n); b = subm(mulredc(b, b, M), two, n); }
        else              { b = subm(mulredc(a, b, M), v, n); a = subm(mulredc(a, a, M), two, n); }
    }
    return a;
}

typedef uint64_t (*StepFn)(uint64_t, uint64_t, const Mont*, uint64_t);

static uint64_t smooth_walk(uint64_t n, uint64_t x0, uint64_t B1, StepFn step,
                            uint64_t aux, uint64_t target) {
    Mont M;
    mont_setup(&M, n);
    size_t np;
    const uint32_t* P = fi_small_primes(&np);
    if (B1 >= FI_SMALL_PRIME_LIMIT) B1 = FI_SMALL_PRIME_LIMIT - 1;
    uint64_t x = to_mont(x0, &M), tg = to_mont(target, &M);
    aux = aux ? to_mont(aux, &M) : 0;
    for (size_t i = 0; i < np && P[i] <= B1; ) {
        size_t j0 = i;
        uint64_t saved = x;
        for (; i < np && P[i] <= B1 && i < j0 + 32; i++) {
            uint64_t q = P[i];
            while (q <= B1 / P[i]) q *= P[i];
            x = step(x, q, &M, aux);
        }
        uint64_t g = fi_u64_gcd(subm(x, tg, n), n);
        if (g == 1) continue;
        if (g < n) return g;
        /* Overshoot: replay prime by prime, then power by power. */
        x = saved;
        for (size_t j = j0; j < i; j++) {
            uint64_t p = P[j], q = p;
            while (q <= B1 / p) q *= p;
            uint64_t y = step(x, q, &M, aux);
            g = fi_u64_gcd(subm(y, tg, n), n);
            if (g == 1) { x = y; continue; }
            if (g < n) return g;
            for (uint64_t qq = p; qq <= q; qq *= p) {
                x = step(x, p, &M, aux);
                g = fi_u64_gcd(subm(x, tg, n), n);
                if (g != 1) return g < n ? g : 0;
            }
            return 0;
        }
        return 0;
    }
    return 0;
}

uint64_t fi_u64_pm1(uint64_t n, uint64_t B1) {
    if (!(n & 1)) return 2;
    return smooth_walk(n, 2, B1 ? B1 : 1000000, pm1_step, 0, 1);
}

uint64_t fi_u64_pp1(uint64_t n, uint64_t B1, uint64_t seed) {
    if (!(n & 1)) return 2;
    uint64_t v = 3 + seed % 1000;
    /* V_0 = 2 is the identity: p | n is found when V_E == 2 (mod p). */
    return smooth_walk(n, v % n, B1 ? B1 : 1000000, pp1_step, 2, 2);
}

/* ------------------------------------------------------------------ */
/* Complete factorisation of a word                                    */
/* ------------------------------------------------------------------ */

/* Split odd composite n (no factor below 59): micro-ECM above 2^40, rho
 * below, with SQUFOF and Lehman/rho as guaranteed backstops. */
static uint64_t u64_split(uint64_t n, uint64_t* rng) {
    uint64_t d;
    if (n < (UINT64_C(1) << 40)) {
        for (uint64_t c = 1; c < 64; c++)
            if ((d = fi_u64_rho(n, c, UINT64_C(1) << 22))) return d;
    }
    if ((d = fi_u64_ecm(n, rng, 128))) return d;
    if ((d = fi_u64_squfof(n, 0))) return d;
    for (uint64_t c = 1; c < 1024; c++)
        if ((d = fi_u64_rho(n, c, UINT64_C(1) << 26))) return d;
    return 0;
}

/* Output sink for a word factorisation: at most 64 prime powers. */
typedef struct { uint64_t p[64]; int64_t e[64]; int n; } WordFac;

static void wf_add(WordFac* w, uint64_t p, int64_t e) {
    for (int i = 0; i < w->n; i++) if (w->p[i] == p) { w->e[i] += e; return; }
    w->p[w->n] = p;
    w->e[w->n++] = e;
}

static void u64_factor_odd(uint64_t n, int64_t mult, WordFac* out, uint64_t* rng) {
    if (n == 1) return;
    if (fi_u64_is_prime(n)) { wf_add(out, n, mult); return; }
    uint64_t s = fi_u64_isqrt(n);
    if (s * s == n) { u64_factor_odd(s, 2 * mult, out, rng); return; }
    uint64_t d = u64_split(n, rng);
    if (!d) { wf_add(out, n, mult); return; }   /* unreachable in practice */
    u64_factor_odd(d, mult, out, rng);
    u64_factor_odd(n / d, mult, out, rng);
}

/* Odd primes below 2^10 for the inverse-multiply divisibility test:
 * n divisible by p  <=>  n * p^-1 (mod 2^64) <= (2^64 - 1) / p. */
static uint32_t tp[200];
static uint64_t tp_inv[200], tp_lim[200];
static size_t tp_n = 0;
static void tp_init(void) {
    for (uint32_t p = 3; p < 1024; p += 2) {
        int ok = 1;
        for (uint32_t q = 3; q * q <= p; q += 2) if (p % q == 0) { ok = 0; break; }
        if (!ok) continue;
        tp[tp_n] = p;
        tp_inv[tp_n] = inv64(p);
        tp_lim[tp_n] = UINT64_MAX / p;
        tp_n++;
    }
}

static void word_factor(uint64_t n, WordFac* out) {
    out->n = 0;
    if (n <= 1) return;
    if (!tp_n) tp_init();
    if (!sq_ready) sq_init();
    int z = ctz64(n);
    if (z) { wf_add(out, 2, z); n >>= z; }
    for (size_t i = 0; i < tp_n && n > 1; i++) {
        uint64_t p = tp[i];
        if (p * p > n) break;
        if (n * tp_inv[i] <= tp_lim[i]) {
            int64_t e = 0;
            do { n *= tp_inv[i]; e++; } while (n * tp_inv[i] <= tp_lim[i]);
            wf_add(out, p, e);
        }
    }
    if (n == 1) return;
    if (n < 1024u * 1024u) { wf_add(out, n, 1); return; }   /* no factor below 2^10 */
    uint64_t rng = n * UINT64_C(0x9E3779B97F4A7C15) | 1;
    u64_factor_odd(n, 1, out, &rng);
}

void fi_u64_factor(uint64_t n, int64_t mult, FiFactors* out) {
    WordFac w;
    word_factor(n, &w);
    for (int i = 0; i < w.n; i++) fi_factors_add_ui(out, w.p[i], w.e[i] * mult);
}

/* Allocation-free variant for FactorInteger's machine-integer fast path:
 * fills p[]/e[] (room for 64) sorted ascending, returns the count. */
int fi_u64_factor_array(uint64_t n, uint64_t* p, int64_t* e) {
    WordFac w;
    word_factor(n, &w);
    for (int i = 1; i < w.n; i++) {                  /* insertion sort */
        uint64_t pp = w.p[i];
        int64_t ee = w.e[i];
        int j = i;
        while (j > 0 && w.p[j - 1] > pp) { w.p[j] = w.p[j - 1]; w.e[j] = w.e[j - 1]; j--; }
        w.p[j] = pp; w.e[j] = ee;
    }
    for (int i = 0; i < w.n; i++) { p[i] = w.p[i]; e[i] = w.e[i]; }
    return w.n;
}
