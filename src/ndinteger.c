/* ndinteger.c — exact-integer NDArray kernels for the number-theoretic heads.
 * See ndinteger.h for what separates these from ndkernels.c's `double` set. */

#include "partitions.h"      /* partitions_p_i64 / partitions_q_i64 */
#include "ndinteger.h"
#include "ndarray.h"
#include "expr.h"
#include "pack.h"
#include "symtab.h"
#include "sym_names.h"
#include "checked_int.h"
#include "facint.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

/* ---------------------------------------------------------------------------
 *  int64 factorisation, by trial division.
 *
 *  The scalar builtins factor through GMP, which is the right choice for one
 *  number of unbounded size and the wrong one for 200000 of them: the per-call
 *  mpz_init/clear alone dominates. Stripping the factors of 2 and then dividing
 *  by the odd numbers needs no allocation at all and settles a 6-digit element
 *  in well under a microsecond.
 *
 *  It is also why there is a size ceiling. Trial division is O(sqrt(n)), so a
 *  single 10^18 element would cost 10^9 iterations and be far slower than the
 *  GMP path it replaced. NDI_FACTOR_MAX caps sqrt(n) at 10^6 iterations; past
 *  it the kernel declines, which abandons the whole array and lets the List
 *  path answer through GMP exactly as before. A ceiling that turns a fast path
 *  off is a performance decision; the ANSWER is identical either way.
 * -------------------------------------------------------------------------- */
#define NDI_FACTOR_MAX ((int64_t)1000000000000LL)   /* 10^12; sqrt = 10^6 */
#define NDI_MAX_FACTORS 64                          /* 2^64 > any int64 */

/* Factor |n| (n != 0) into distinct primes and exponents. Returns false when n
 * is 0 or beyond NDI_FACTOR_MAX. n == 1 gives the empty factorisation. */
static bool ndi_factor(int64_t n, int64_t* p, int* e, int* np) {
    if (n == 0) return false;
    int64_t m;
    if (ci_abs_i64(n, &m)) return false;            /* INT64_MIN has no |n| */
    if (m > NDI_FACTOR_MAX) return false;
    int k = 0;
    if ((m & 1) == 0) {
        p[k] = 2; e[k] = 0;
        while ((m & 1) == 0) { e[k]++; m >>= 1; }
        k++;
    }
    for (int64_t d = 3; d * d <= m; d += 2) {
        if (m % d) continue;
        if (k >= NDI_MAX_FACTORS) return false;
        p[k] = d; e[k] = 0; k++;
        while (m % d == 0) { e[k-1]++; m /= d; }
    }
    if (m > 1) {
        if (k >= NDI_MAX_FACTORS) return false;
        p[k] = m; e[k] = 1; k++;
    }
    *np = k;
    return true;
}

/* Exported for the scalar number-theory builtins (MultiplicativeOrder's
 * machine path): the same trial-division factorisation, same ceiling. */
bool ndint_factor_i64(int64_t n, int64_t* p, int* e, int* np) {
    return ndi_factor(n, p, e, np);
}

/* ---- unary integer kernels ---------------------------------------------- */

/* IntegerLength[n]: decimal digits of |n|. IntegerLength[0] is 0, matching
 * builtin_integerlength (and Mathematica). */
static bool ndk_IntegerLength_ii(int64_t n, int64_t* out) {
    int64_t m;
    if (ci_abs_i64(n, &m)) return false;
    int64_t d = 0;
    while (m > 0) { d++; m /= 10; }
    *out = d;
    return true;
}

/* BitLength[n]: binary bits of n. BitLength[0] is 0 and, for n < 0,
 * BitLength[n] = BitLength[BitNot[n]] with BitNot[n] = ~n = -n-1. Working on the
 * unsigned magnitude `~(uint64_t)n` for negatives is exact for every int64
 * (including INT64_MIN, whose complement is INT64_MAX), so this kernel — unlike
 * ndk_IntegerLength_ii, which negates — never overflows and always succeeds. */
static bool ndk_BitLength_ii(int64_t n, int64_t* out) {
    uint64_t m = (n >= 0) ? (uint64_t)n : ~(uint64_t)n;
    int64_t d = 0;
    while (m) { d++; m >>= 1; }
    *out = d;
    return true;
}

/* BitNot[n] = ~n = -1 - n. Total on int64 (INT64_MIN <-> INT64_MAX), so it
 * never declines; the scalar builtin and Compile's OP_BNOT_I compute the same
 * `~n`, so every surface agrees bit for bit. */
static bool ndk_BitNot_ii(int64_t n, int64_t* out) {
    *out = ~n;
    return true;
}

/* Whole-buffer twin: one branch-free loop the compiler vectorises, instead of
 * an indirect call per element. */
static bool ndk_BitNot_whole(const int64_t* in, size_t sz, int64_t* out) {
    for (size_t i = 0; i < sz; i++) out[i] = ~in[i];
    return true;
}

/* EulerPhi[n] = |n| * prod (1 - 1/p). EulerPhi[0] is 0 and the sign of n is
 * ignored, both matching the scalar builtin. Divide before multiplying so the
 * running value never exceeds |n| and cannot overflow. */
static int64_t ndi_phi_from_factors(int64_t m, const int64_t* p, int np) {
    int64_t phi = m;
    for (int i = 0; i < np; i++) phi = phi / p[i] * (p[i] - 1);
    return phi;
}

/* Shared by the scalar builtin and Compile's OP_PHI_I. */
bool ndint_euler_phi_i64(int64_t n, int64_t* out) {
    if (n == 0) { *out = 0; return true; }
    int64_t p[NDI_MAX_FACTORS]; int e[NDI_MAX_FACTORS]; int np;
    if (!ndi_factor(n, p, e, &np)) return false;
    int64_t m;
    if (ci_abs_i64(n, &m)) return false;
    *out = ndi_phi_from_factors(m, p, np);
    return true;
}

static bool ndk_EulerPhi_ii(int64_t n, int64_t* out) {
    return ndint_euler_phi_i64(n, out);
}

/* MoebiusMu[n]: 0 if any exponent >= 2, else (-1)^(distinct primes).
 * MoebiusMu[0] is 0 (0 is divisible by every square). Sign
 * ignored. Shared by the scalar builtin and Compile's OP_MU_I. */
static int64_t ndi_mu_from_exps(const int* e, int np) {
    for (int i = 0; i < np; i++) if (e[i] >= 2) return 0;
    return (np % 2 == 0) ? 1 : -1;
}

bool ndint_moebius_mu_i64(int64_t n, int64_t* out) {
    if (n == 0) { *out = 0; return true; }
    int64_t p[NDI_MAX_FACTORS]; int e[NDI_MAX_FACTORS]; int np;
    if (!ndi_factor(n, p, e, &np)) return false;
    *out = ndi_mu_from_exps(e, np);
    return true;
}

static bool ndk_MoebiusMu_ii(int64_t n, int64_t* out) {
    return ndint_moebius_mu_i64(n, out);
}

/* ---- binary integer kernels --------------------------------------------- */

static int64_t ndi_gcd_u(int64_t a, int64_t b) {
    while (b) { int64_t t = a % b; a = b; b = t; }
    return a;
}

/* GCD[a, b] over the non-negative representatives: GCD[-12, 18] is 6 and
 * GCD[0, 5] is 5, both matching the scalar builtin. */
static bool ndk_GCD_ii(int64_t a, int64_t b, int64_t* out) {
    int64_t x, y;
    if (ci_abs_i64(a, &x) || ci_abs_i64(b, &y)) return false;
    *out = ndi_gcd_u(x, y);
    return true;
}

/* BitXor[a, b]: two's-complement XOR. Total and exact on int64 (the XOR of two
 * int64 values is an int64), so it never declines. The scalar builtin's
 * all-machine-integer path computes the same `a ^ b`, and so does Compile's
 * OP_BXOR_I, so every surface agrees bit for bit. */
static bool ndk_BitXor_ii(int64_t a, int64_t b, int64_t* out) {
    *out = a ^ b;
    return true;
}

/* Whole-buffer twin: XOR against a broadcast scalar is symmetric, so the side
 * the buffer was on does not matter. No call per element -- the loop
 * vectorises. */
static bool ndk_BitXor_whole(const int64_t* in, size_t sz, int64_t scal,
                             bool arr_first, int64_t* out) {
    (void)arr_first;
    for (size_t i = 0; i < sz; i++) out[i] = in[i] ^ scal;
    return true;
}

/* Two int64 buffers of the same shape: element-wise XOR in one tight loop. */
static bool ndk_BitXor_whole2(const int64_t* a, const int64_t* b, size_t sz,
                              int64_t* out) {
    for (size_t i = 0; i < sz; i++) out[i] = a[i] ^ b[i];
    return true;
}

/* BitAnd[a, b]: two's-complement AND. Total and exact on int64, so it never
 * declines; the scalar builtin and Compile's OP_BAND_I compute the same
 * `a & b`, so every surface agrees bit for bit. */
static bool ndk_BitAnd_ii(int64_t a, int64_t b, int64_t* out) {
    *out = a & b;
    return true;
}

/* Whole-buffer twin: AND against a broadcast scalar is symmetric. */
static bool ndk_BitAnd_whole(const int64_t* in, size_t sz, int64_t scal,
                             bool arr_first, int64_t* out) {
    (void)arr_first;
    for (size_t i = 0; i < sz; i++) out[i] = in[i] & scal;
    return true;
}

/* Two int64 buffers of the same shape: element-wise AND in one tight loop. */
static bool ndk_BitAnd_whole2(const int64_t* a, const int64_t* b, size_t sz,
                              int64_t* out) {
    for (size_t i = 0; i < sz; i++) out[i] = a[i] & b[i];
    return true;
}

/* BitOr[a, b]: two's-complement OR. Total and exact on int64, so it never
 * declines; the scalar builtin and Compile's OP_BOR_I compute the same
 * `a | b`, so every surface agrees bit for bit. */
static bool ndk_BitOr_ii(int64_t a, int64_t b, int64_t* out) {
    *out = a | b;
    return true;
}

/* Whole-buffer twin: OR against a broadcast scalar is symmetric. */
static bool ndk_BitOr_whole(const int64_t* in, size_t sz, int64_t scal,
                            bool arr_first, int64_t* out) {
    (void)arr_first;
    for (size_t i = 0; i < sz; i++) out[i] = in[i] | scal;
    return true;
}

/* Two int64 buffers of the same shape: element-wise OR in one tight loop. */
static bool ndk_BitOr_whole2(const int64_t* a, const int64_t* b, size_t sz,
                             int64_t* out) {
    for (size_t i = 0; i < sz; i++) out[i] = a[i] | b[i];
    return true;
}

/* BitSet[n, k]: n with bit k set, through ndint_bitset_i64 -- the helper the
 * scalar builtin and Compile's OP_BSET_I share. Declines k < 0 (BitSet stays
 * unevaluated) and n >= 0 with k >= 63 (a bignum answer); the whole array
 * then takes the List path. */
static bool ndk_BitSet_ii(int64_t n, int64_t k, int64_t* out) {
    return ndint_bitset_i64(n, k, out);
}

/* BitSet with a per-element k: n | (1 << k) for every lane, branch-free, with
 * any k outside [0, 62] folded into a flag -- if one is seen the arm DECLINES
 * and the per-element loop answers it (a negative n with k >= 63 is still an
 * int64 there; anything else degrades to the List path). The masked shift
 * count keeps the discarded lanes defined behaviour.
 *
 * The loop is bound by the VARIABLE 64-bit shift. arm64's NEON shifts each
 * lane by its own count, so the plain loop vectorises there; baseline x86-64
 * (SSE2) has no per-lane shift at all and runs it scalar, ~1.5x behind NumPy,
 * which dispatches at run time. So on x86 the same body is also compiled for
 * AVX2 (vpsllvq) and chosen once, in ndinteger_init, from the CPU's own
 * feature bits. Both clones are the identical C body: bit-identical answers. */
#define NDI_BITSET_VAR_BODY(AEXPR)                                          \
    uint64_t bad = 0;                                                       \
    for (size_t i = 0; i < sz; i++) {                                       \
        uint64_t k = (uint64_t)ks[i];                                       \
        bad |= (uint64_t)(k > 62);                                          \
        out[i] = (int64_t)((uint64_t)(AEXPR) | (UINT64_C(1) << (k & 63)));  \
    }                                                                       \
    return bad == 0;

static bool ndi_bitset_sn(int64_t n, const int64_t* ks, size_t sz, int64_t* out) {
    NDI_BITSET_VAR_BODY(n)
}
static bool ndi_bitset_an(const int64_t* a, const int64_t* ks, size_t sz, int64_t* out) {
    NDI_BITSET_VAR_BODY(a[i])
}

#if defined(__GNUC__) && !defined(__clang__) && defined(__x86_64__) && !defined(__AVX2__)
#define NDI_BITSET_AVX2 1
__attribute__((target("avx2")))
static bool ndi_bitset_sn_avx2(int64_t n, const int64_t* ks, size_t sz, int64_t* out) {
    NDI_BITSET_VAR_BODY(n)
}
__attribute__((target("avx2")))
static bool ndi_bitset_an_avx2(const int64_t* a, const int64_t* ks, size_t sz, int64_t* out) {
    NDI_BITSET_VAR_BODY(a[i])
}
#endif

/* Selected once by ndinteger_init (before any evaluation, so no thread ever
 * sees them change); the portable clones until then. */
static bool (*ndi_bitset_sn_fn)(int64_t, const int64_t*, size_t, int64_t*) = ndi_bitset_sn;
static bool (*ndi_bitset_an_fn)(const int64_t*, const int64_t*, size_t, int64_t*) = ndi_bitset_an;

/* Whole-buffer twin. A buffer of n with a scalar k in [0, 62] is one OR with
 * a constant mask, which vectorises everywhere; a scalar n over a buffer of
 * k takes the per-element-k loop above. */
static bool ndk_BitSet_whole(const int64_t* in, size_t sz, int64_t scal,
                             bool arr_first, int64_t* out) {
    if (arr_first) {
        if (scal < 0 || scal > 62) return false;
        uint64_t mask = UINT64_C(1) << scal;
        for (size_t i = 0; i < sz; i++)
            out[i] = (int64_t)((uint64_t)in[i] | mask);
        return true;
    }
    return ndi_bitset_sn_fn(scal, in, sz, out);
}

/* Two int64 buffers of the same shape: n from the first, k from the second. */
static bool ndk_BitSet_whole2(const int64_t* a, const int64_t* b, size_t sz,
                              int64_t* out) {
    return ndi_bitset_an_fn(a, b, sz, out);
}

/* LCM[a, b] = |a|/gcd * |b|, in that order so the product is formed from the
 * already-reduced factor. LCM[0, n] is 0. Overflow declines. */
static bool ndk_LCM_ii(int64_t a, int64_t b, int64_t* out) {
    int64_t x, y;
    if (ci_abs_i64(a, &x) || ci_abs_i64(b, &y)) return false;
    if (x == 0 || y == 0) { *out = 0; return true; }
    int64_t g = ndi_gcd_u(x, y);
    return !ci_mul_i64(x / g, y, out);
}

/* DivisorSigma[k, n] = prod_i (1 + q_i + q_i^2 + ... + q_i^e_i),  q_i = p_i^k.
 *
 * The shared helper behind every machine surface of DivisorSigma: the element
 * kernel below, its sieve twin, the Compile VM's OP_DSIG_I and the scalar
 * builtin's machine path all reach this one function, so the four answers are
 * the same int64 by construction. q is raised once per prime and the geometric
 * sum is Horner's, so a prime of multiplicity e costs e multiply-adds rather
 * than e powerings. Declines k < 0 (sigma_-1[6] is 2 only by luck: sigma_-2[10]
 * is 13/10, a Rational, which no int64 holds) and any overflow. */
bool ndint_sigma_from_factors(int64_t k, const int64_t* p, const int* e, int np,
                              int64_t* out) {
    if (k < 0) return false;
    int64_t total = 1;
    for (int i = 0; i < np; i++) {
        int64_t q, term = 1;
        if (ci_powi_i64(p[i], k, &q)) return false;
        for (int j = 0; j < e[i]; j++) {
            if (ci_mul_i64(term, q, &term)) return false;
            if (ci_add_i64(term, 1, &term)) return false;
        }
        if (ci_mul_i64(total, term, &total)) return false;
    }
    *out = total;
    return true;
}

/* Argument order is the CALL order, which ndarray_map_binary composes for
 * either operand position -- so `k` is always the first argument and `n` the
 * second, whichever of the two was the array. Declines n == 0 (unevaluated in
 * the scalar builtin) as well as everything ndint_sigma_from_factors declines.
 * The sign of n is ignored, matching the builtin. */
bool ndint_divisor_sigma_i64(int64_t k, int64_t n, int64_t* out) {
    if (n == 0 || k < 0) return false;
    int64_t p[NDI_MAX_FACTORS]; int e[NDI_MAX_FACTORS]; int np;
    if (!ndi_factor(n, p, e, &np)) return false;
    return ndint_sigma_from_factors(k, p, e, np, out);
}

static bool ndk_DivisorSigma_ii(int64_t k, int64_t n, int64_t* out) {
    return ndint_divisor_sigma_i64(k, n, out);
}

/* ---- DivisorSigma over a whole buffer: one smallest-prime-factor sieve ----- *
 * Per element, ndi_factor is trial division, O(sqrt(n)): over Range[10^7] that
 * is ~10^3 divisions an element and the array cost grows as N^1.5 (18 s at
 * 10^7 against 0.66 s at 10^6). A smallest-prime-factor table to max|n| makes
 * every factorisation O(Omega(n)) table reads after one O(max log log max)
 * pass, so the array costs about what the sieve does.
 *
 * The table covers ODD m only, and stores the smallest prime factor of a
 * composite m (or 0 for a prime). That factor is at most sqrt(max) < 2^16 for
 * every max this ever sieves to, so a uint16 suffices: one byte per integer of
 * range, the same footprint as Prime's byte sieve, under the same ceilings. */
#define NDI_SPF_FREE  ((int64_t)4000000)    /* 2 MB, always worth building */
#define NDI_SPF_MAX   ((int64_t)50000000)   /* 25 MB, the hard ceiling     */

/* Plan a sieve for the buffer: the largest |element|, and whether a sieve to
 * it is worth building (see the measured break-even below). Declines an
 * INT64_MIN element (no |n|), and a zero one unless `allow_zero`. */
static bool ndi_spf_plan(const int64_t* in, size_t sz, bool allow_zero,
                         int64_t* nmax_out) {
    if (sz == 0) return false;
    int64_t nmax = 0;
    for (size_t i = 0; i < sz; i++) {
        int64_t v = in[i];
        if (v == INT64_MIN || (v == 0 && !allow_zero)) return false;
        if (v < 0) v = -v;
        if (v > nmax) nmax = v;
    }
    if (nmax > NDI_SPF_MAX) return false;
    /* Past _FREE the sieve must be paid for by the elements. Measured (Apple
     * M-series, 2026-10-19): the sieve costs ~6.5 ns per integer of range and
     * trial division ~1.5 ns per unit of sqrt(n) per element (10.6 us an
     * element near 5*10^7), so the sieve wins once sz >= ~4.3 sqrt(max) --
     * 30000 elements at the 5*10^7 ceiling, where it is a 0.3 s build. */
    if (nmax > NDI_SPF_FREE && (double)sz < 4.0 * sqrt((double)nmax))
        return false;
    *nmax_out = nmax;
    return true;
}

/* The smallest-prime-factor table over the odd numbers up to nmax (index
 * m >> 1), 0 for a prime. NULL on allocation failure. Caller frees. */
static uint16_t* ndi_spf_build(int64_t nmax) {
    size_t half = (size_t)(nmax / 2) + 1;
    uint16_t* spf = (uint16_t*)calloc(half, sizeof(uint16_t));
    if (!spf) return NULL;
    for (int64_t q = 3; q * q <= nmax; q += 2) {
        if (spf[q >> 1]) continue;                        /* q composite */
        for (int64_t j = q * q; j <= nmax; j += 2 * q)
            if (!spf[j >> 1]) spf[j >> 1] = (uint16_t)q;
    }
    return spf;
}

/* Factor m (1 <= m <= the table's nmax) off the table: distinct primes in
 * ascending order with their exponents. Walks the smallest-factor chain, one
 * exact division per prime factor (with multiplicity) and never a trial one. */
static int ndi_spf_factor(const uint16_t* spf, int64_t m, int64_t* p, int* e) {
    int np = 0;
    if ((m & 1) == 0) {
        p[np] = 2; e[np] = 0;
        while ((m & 1) == 0) { e[np]++; m >>= 1; }
        np++;
    }
    while (m > 1) {
        int64_t q = spf[m >> 1] ? (int64_t)spf[m >> 1] : m;
        m /= q;
        if (np > 0 && p[np - 1] == q) e[np - 1]++;
        else { p[np] = q; e[np] = 1; np++; }
    }
    return np;
}

static bool ndk_DivisorSigma_whole(const int64_t* in, size_t sz, int64_t scal,
                                   bool arr_first, int64_t* out) {
    /* Only DivisorSigma[k, arr]: an array of k is not a factoring workload. */
    if (arr_first || scal < 0) return false;
    int64_t nmax;
    if (!ndi_spf_plan(in, sz, false, &nmax)) return false;
    uint16_t* spf = ndi_spf_build(nmax);
    if (!spf) return false;
    bool ok = true;
    for (size_t i = 0; i < sz && ok; i++) {
        int64_t p[NDI_MAX_FACTORS]; int e[NDI_MAX_FACTORS];
        int np = ndi_spf_factor(spf, in[i] < 0 ? -in[i] : in[i], p, e);
        ok = ndint_sigma_from_factors(scal, p, e, np, &out[i]);
    }
    free(spf);
    return ok;
}

/* EulerPhi and MoebiusMu over a whole buffer, off the same sieve. Both are
 * total on int64 bar INT64_MIN (phi(n) <= |n|, mu in {-1, 0, 1}), so once the
 * plan accepts the buffer every element answers. */
static bool ndk_EulerPhi_whole(const int64_t* in, size_t sz, int64_t* out) {
    int64_t nmax;
    if (!ndi_spf_plan(in, sz, true, &nmax)) return false;
    uint16_t* spf = ndi_spf_build(nmax);
    if (!spf) return false;
    for (size_t i = 0; i < sz; i++) {
        int64_t m = in[i] < 0 ? -in[i] : in[i];
        if (m == 0) { out[i] = 0; continue; }
        int64_t p[NDI_MAX_FACTORS]; int e[NDI_MAX_FACTORS];
        int np = ndi_spf_factor(spf, m, p, e);
        out[i] = ndi_phi_from_factors(m, p, np);
    }
    free(spf);
    return true;
}

static bool ndk_MoebiusMu_whole(const int64_t* in, size_t sz, int64_t* out) {
    int64_t nmax;
    if (!ndi_spf_plan(in, sz, true, &nmax)) return false;
    uint16_t* spf = ndi_spf_build(nmax);
    if (!spf) return false;
    for (size_t i = 0; i < sz; i++) {
        int64_t m = in[i] < 0 ? -in[i] : in[i];
        if (m == 0) { out[i] = 0; continue; }
        int64_t p[NDI_MAX_FACTORS]; int e[NDI_MAX_FACTORS];
        int np = ndi_spf_factor(spf, m, p, e);
        out[i] = ndi_mu_from_exps(e, np);
    }
    free(spf);
    return true;
}

/* JacobiSymbol[a, n] in int64: the full Kronecker symbol, exactly what the
 * scalar builtin's mpz_kronecker computes (Cohen, Alg. 1.4.10), so it is
 * total -- every int64 pair has an answer in {-1, 0, 1}. Shared by the
 * element kernel and Compile's OP_JACOBI_I. */
int64_t ndint_kronecker_i64(int64_t a, int64_t n) {
    static const int tab2[8] = { 0, 1, 0, -1, 0, -1, 0, 1 };
    if (n == 0) return (a == 1 || a == -1) ? 1 : 0;
    if ((a & 1) == 0 && (n & 1) == 0) return 0;
    uint64_t un = n < 0 ? (uint64_t)0 - (uint64_t)n : (uint64_t)n;
    int k = 1;
    int v = 0;
    while ((un & 1) == 0) { un >>= 1; v++; }
    if (v & 1) k = tab2[(uint64_t)a & 7];
    if (n < 0 && a < 0) k = -k;
    /* un is now odd and positive: (a/un) depends only on a mod un. */
    uint64_t ua;
    if (a >= 0) ua = (uint64_t)a % un;
    else {
        uint64_t r = ((uint64_t)0 - (uint64_t)a) % un;
        ua = r ? un - r : 0;
    }
    while (ua != 0) {
        v = 0;
        while ((ua & 1) == 0) { ua >>= 1; v++; }
        if (v & 1) k *= tab2[un & 7];
        if (ua & un & 2) k = -k;                        /* reciprocity */
        uint64_t t = un % ua; un = ua; ua = t;
    }
    return un == 1 ? k : 0;
}

static bool ndk_JacobiSymbol_ii(int64_t a, int64_t n, int64_t* out) {
    *out = ndint_kronecker_i64(a, n);
    return true;
}

/* Integer-only descriptors: no `cplx`, no `real`, `real_closed` false. Every
 * non-integer branch of ndarray_map_unary / _binary tests `cplx` and declines,
 * so a Real buffer degrades to the List path — which is what these heads want
 * anyway (EulerPhi[3.] is not EulerPhi[3]). */
static const NDUnaryKernel NDKU_IntegerLength =
    { NULL, NULL, false, false, NULL, ndk_IntegerLength_ii, true, NULL };
static const NDUnaryKernel NDKU_BitLength =
    { NULL, NULL, false, false, NULL, ndk_BitLength_ii, true, NULL };
static const NDUnaryKernel NDKU_BitNot =
    { NULL, NULL, false, false, NULL, ndk_BitNot_ii, true, ndk_BitNot_whole };
static const NDUnaryKernel NDKU_EulerPhi =
    { NULL, NULL, false, false, NULL, ndk_EulerPhi_ii, true, ndk_EulerPhi_whole };
static const NDUnaryKernel NDKU_MoebiusMu =
    { NULL, NULL, false, false, NULL, ndk_MoebiusMu_ii, true, ndk_MoebiusMu_whole };
/* PartitionsP / PartitionsQ: an exact int64 table lookup, through the helper
 * the scalar builtins and Compile share. Past the int64 range (p beyond
 * n = 405) it declines and the whole array takes the List path to GMP. */
static bool ndk_PartitionsP_ii(int64_t n, int64_t* out) { return partitions_p_i64(n, out); }
static bool ndk_PartitionsQ_ii(int64_t n, int64_t* out) { return partitions_q_i64(n, out); }
static const NDUnaryKernel NDKU_PartitionsP =
    { NULL, NULL, false, false, NULL, ndk_PartitionsP_ii, true, NULL };
static const NDUnaryKernel NDKU_PartitionsQ =
    { NULL, NULL, false, false, NULL, ndk_PartitionsQ_ii, true, NULL };
static const NDBinaryKernel NDKB_GCD =
    { NULL, false, NULL, ndk_GCD_ii, true, NULL, NULL };
static const NDBinaryKernel NDKB_LCM =
    { NULL, false, NULL, ndk_LCM_ii, true, NULL, NULL };
static const NDBinaryKernel NDKB_BitXor =
    { NULL, false, NULL, ndk_BitXor_ii, true, ndk_BitXor_whole, ndk_BitXor_whole2 };
static const NDBinaryKernel NDKB_BitAnd =
    { NULL, false, NULL, ndk_BitAnd_ii, true, ndk_BitAnd_whole, ndk_BitAnd_whole2 };
static const NDBinaryKernel NDKB_BitOr =
    { NULL, false, NULL, ndk_BitOr_ii, true, ndk_BitOr_whole, ndk_BitOr_whole2 };
static const NDBinaryKernel NDKB_BitSet =
    { NULL, false, NULL, ndk_BitSet_ii, true, ndk_BitSet_whole, ndk_BitSet_whole2 };
/* IntegerExponent[n, b]: the largest e with b^e | n. Declines n == 0
 * (Infinity, not an int64) and b < 2 (::ibase in the scalar builtin). */
static bool ndk_IntegerExponent_ii(int64_t n, int64_t b, int64_t* out) {
    if (n == 0 || b < 2) return false;
    uint64_t m = n < 0 ? (uint64_t)0 - (uint64_t)n : (uint64_t)n;
    int64_t e = 0;
    while (m % (uint64_t)b == 0) { m /= (uint64_t)b; e++; }
    *out = e;
    return true;
}

static const NDBinaryKernel NDKB_IntegerExponent =
    { NULL, false, NULL, ndk_IntegerExponent_ii, true, NULL, NULL };
static const NDBinaryKernel NDKB_JacobiSymbol =
    { NULL, false, NULL, ndk_JacobiSymbol_ii, true, NULL, NULL };
static const NDBinaryKernel NDKB_DivisorSigma =
    { NULL, false, NULL, ndk_DivisorSigma_ii, true, ndk_DivisorSigma_whole, NULL };

/* ---------------------------------------------------------------------------
 *  Head-level buffer paths.
 * -------------------------------------------------------------------------- */

/* The sole argument, when it is a rank-1-or-higher int64 NDArray. */
static const Expr* ndi_int_arg(const Expr* res, size_t which, size_t argc) {
    if (!res || res->type != EXPR_FUNCTION) return NULL;
    if (res->data.function.arg_count != argc) return NULL;
    const Expr* a = res->data.function.args[which];
    if (!is_ndarray(a) || a->data.ndarray.dtype != NDT_INT64) return NULL;
    return a;
}

/* An exact non-negative int64 from an Integer argument. */
static bool ndi_exact_int(const Expr* e, int64_t* out) {
    if (!e || e->type != EXPR_INTEGER) return false;
    *out = e->data.integer;
    return true;
}

/* ---- Prime[arr] --------------------------------------------------------- *
 * One sieve for the whole array, not one prime count per element.
 *
 * Prime[n] alone is a Meissel/Lucy count, which is the right algorithm for a
 * single large index and the wrong one repeated 200000 times. Given the largest
 * index in the array, a sieve of Eratosthenes to Rosser's bound
 *     p_n < n (ln n + ln ln n)      (n >= 6)
 * produces every requested prime in one linear pass, and the gather is an
 * indexed read. The array as a whole then costs about what its largest single
 * element used to.
 *
 * Declines an index below 1 (Prime[0] is unevaluated, with a message, in the
 * scalar builtin) and an upper bound past what a byte sieve should allocate,
 * where the builtin's own path is the better one anyway. */
/* THE SIEVE IS NOT UNCONDITIONALLY BETTER, and the condition is the array's
 * LENGTH, not its values. One sieve to `limit` costs O(limit log log limit)
 * whatever the array holds, where the scalar path costs one Meissel count per
 * element -- so for a two-element array asking for a large index the builtin's
 * own path wins, and taking the sieve anyway would be a regression paid in
 * hundreds of megabytes.
 *
 *   under NDI_SIEVE_FREE   a few milliseconds and a few MB. Cheaper than a
 *                          single Meissel count at that size, so no amount of
 *                          elements is needed to justify it.
 *   up to NDI_SIEVE_MAX    only once there are enough elements to amortise it.
 *   past NDI_SIEVE_MAX     never; decline, and each element takes the Meissel
 *                          count exactly as before.
 */
#define NDI_SIEVE_FREE  ((int64_t)4000000)    /* ~4 MB, always worth it */
#define NDI_SIEVE_MAX   ((int64_t)50000000)   /* ~50 MB, the hard ceiling */
#define NDI_SIEVE_AMORT ((size_t)64)          /* elements needed past _FREE */

Expr* ndint_prime(Expr* res) {
    const Expr* a = ndi_int_arg(res, 0, 1);
    if (!a) return NULL;
    size_t sz = ndarray_size(a);
    if (sz == 0) return NULL;
    const void* in = a->data.ndarray.data;

    int64_t nmax = 0;
    for (size_t i = 0; i < sz; i++) {
        int64_t v = ndt_get_i(in, i, NDT_INT64);
        if (v < 1) return NULL;                    /* Prime[0] stays symbolic */
        if (v > nmax) nmax = v;
    }

    /* Rosser & Schoenfeld's upper bound on the n-th prime, with the small-n
     * cases (where ln ln n is <= 0) handled by a flat floor. */
    double dn = (double)nmax;
    int64_t limit;
    if (nmax < 6) limit = 13;
    else {
        double ln = log(dn), lnln = log(ln);
        limit = (int64_t)(dn * (ln + lnln)) + 16;
    }
    if (limit > NDI_SIEVE_MAX) return NULL;
    if (limit > NDI_SIEVE_FREE && sz < NDI_SIEVE_AMORT) return NULL;

    unsigned char* comp = calloc((size_t)limit + 1, 1);
    if (!comp) return NULL;
    for (int64_t i = 2; i * i <= limit; i++)
        if (!comp[i]) for (int64_t j = i * i; j <= limit; j += i) comp[j] = 1;

    /* Index the primes in order; primes[k] is the (k+1)-th prime. */
    int64_t* primes = malloc(sizeof(int64_t) * (size_t)nmax);
    if (!primes) { free(comp); return NULL; }
    int64_t found = 0;
    for (int64_t i = 2; i <= limit && found < nmax; i++)
        if (!comp[i]) primes[found++] = i;
    free(comp);
    if (found < nmax) { free(primes); return NULL; }   /* bound was too tight */

    int64_t* out = malloc(sizeof(int64_t) * sz);
    if (!out) { free(primes); return NULL; }
    for (size_t i = 0; i < sz; i++)
        out[i] = primes[ndt_get_i(in, i, NDT_INT64) - 1];
    free(primes);
    return expr_new_ndarray_like(a, a->data.ndarray.rank, a->data.ndarray.dims,
                                 out, NDT_INT64);
}

/* ---- PowerMod[arr, e, m] ------------------------------------------------ *
 * Ternary, so no NDBinaryKernel can carry it. Restricted to a NON-NEGATIVE
 * exponent: PowerMod[b, -1, m] is a modular inverse, which exists only when
 * gcd(b, m) is 1, and an array where some elements have an inverse and some do
 * not has no uniform answer — the List path handles that case element by
 * element, and answers unevaluated where it must.
 *
 * The modulus is bounded so that the squaring step cannot overflow int64: with
 * |m| <= 3037000499 (< 2^31.5), every intermediate product is under 2^63. */
#define NDI_POWERMOD_MAX ((int64_t)3037000499LL)

Expr* ndint_powermod(Expr* res) {
    const Expr* a = ndi_int_arg(res, 0, 3);
    if (!a) return NULL;
    int64_t e, m;
    if (!ndi_exact_int(res->data.function.args[1], &e) || e < 0) return NULL;
    if (!ndi_exact_int(res->data.function.args[2], &m)) return NULL;
    /* A negative modulus is left unevaluated by the builtin,
     * so the buffer declines and the List path reproduces it. */
    if (m <= 0) return NULL;
    int64_t am = m;
    if (am > NDI_POWERMOD_MAX) return NULL;

    size_t sz = ndarray_size(a);
    const void* in = a->data.ndarray.data;
    int64_t* out = malloc(sizeof(int64_t) * (sz ? sz : 1));
    if (!out) return NULL;
    for (size_t i = 0; i < sz; i++) {
        int64_t b = ndt_get_i(in, i, NDT_INT64) % am;
        if (b < 0) b += am;
        int64_t r = 1 % am, k = e;
        while (k > 0) {
            if (k & 1) r = (r * b) % am;
            k >>= 1;
            if (k > 0) b = (b * b) % am;
        }
        out[i] = r;
    }
    return expr_new_ndarray_like(a, a->data.ndarray.rank, a->data.ndarray.dims,
                                 out, NDT_INT64);
}

/* ---- IntegerDigits[arr] ------------------------------------------------- *
 * The result is ragged, so the win is entirely on the input side: read each
 * element straight out of the buffer rather than materialising the whole array
 * into Integers first and then reading each one back. Rank 1 only — a
 * higher-rank argument's result nests one level deeper per axis, which the
 * List path already builds correctly and is not the shape anyone measures. */
Expr* ndint_integerdigits(Expr* res) {
    const Expr* a = ndi_int_arg(res, 0, 1);
    if (!a || a->data.ndarray.rank != 1) return NULL;
    size_t sz = ndarray_size(a);
    const void* in = a->data.ndarray.data;

    for (size_t i = 0; i < sz; i++)                 /* INT64_MIN: let GMP have it */
        if (ndt_get_i(in, i, NDT_INT64) == INT64_MIN) return NULL;

    /* Expressions are refcounted, so every digit node is one of ten shared
     * Integers and every row shares one List head: each row costs a single
     * node, filled in place -- 10^6 rows of ~7 digits allocate 10^6 nodes
     * instead of ~8 * 10^6 (IntegerDigits[Range[10^6]]: 1.23 s -> see the
     * v0.473 changelog). */
    Expr* dig[10];
    for (int j = 0; j < 10; j++) dig[j] = expr_new_integer(j);
    Expr* head = expr_new_symbol(SYM_List);
    Expr* out = expr_new_function(expr_new_symbol(SYM_List), NULL, sz);
    if (!out) { for (int j = 0; j < 10; j++) expr_free(dig[j]); expr_free(head); return NULL; }
    for (size_t i = 0; i < sz; i++) {
        int64_t v = ndt_get_i(in, i, NDT_INT64);
        uint64_t m = v < 0 ? (uint64_t)0 - (uint64_t)v : (uint64_t)v;
        int d[20]; int nd = 0;
        if (m == 0) d[nd++] = 0;
        while (m > 0) { d[nd++] = (int)(m % 10); m /= 10; }
        Expr* row = expr_new_function(expr_copy(head), NULL, (size_t)nd);
        for (int j = 0; j < nd; j++) row->data.function.args[j] = expr_copy(dig[d[nd - 1 - j]]);
        out->data.function.args[i] = row;
    }
    for (int j = 0; j < 10; j++) expr_free(dig[j]);
    expr_free(head);
    return out;
}

/* ---- Positive / Negative / NonNegative / NonPositive -------------------- */

Expr* ndint_sign_predicate(Expr* res, NDSignPred which) {
    if (!res || res->type != EXPR_FUNCTION) return NULL;
    if (res->data.function.arg_count != 1) return NULL;
    const Expr* a = res->data.function.args[0];
    if (!is_ndarray(a)) return NULL;
    NDType dt = a->data.ndarray.dtype;
    /* A complex buffer has no order, and the scalar builtins leave Positive[I]
     * unevaluated rather than answering False -- so only the ordered dtypes. */
    if (ndt_is_complex(dt)) return NULL;
    /* Rank > 1 threads to a nested List of True/False, which the List path
     * already builds; the flat vector is the shape that is measured. */
    if (a->data.ndarray.rank != 1) return NULL;

    size_t sz = ndarray_size(a);
    if (sz == 0) return NULL;           /* empty: let the List path build {} */
    const void* in = a->data.ndarray.data;

    /* The answer is now a one-byte-per-element bool BUFFER, not 10^6 shared
     * True/False Expr nodes (performance.md §13 gap C.1). It always presents as a
     * List -- a PACKED bool List, observably identical to the plain List of
     * True/False the predicate returned before this change, whichever surface the
     * input came from. (Unlike Sin, which inherits its input's presentation, a
     * sign predicate has always answered with a List; keeping that is what lets
     * the visible surface agree with the packed one and preserves the contract
     * the tests pin.) */
    void* buf = NULL;
    int64_t dims[1] = { (int64_t)sz };
    Expr* out = ndbuild_open_like(a, 1, dims, NDT_BOOL, &buf);
    if (!out) return NULL;
    out->data.ndarray.present_as = NDA_HEAD_LIST;   /* always a List, never NDArray[...] */
    pack_g_any_created = true;                       /* a packed list now exists */
    uint8_t* ob = (uint8_t*)buf;
    for (size_t i = 0; i < sz; i++) {
        bool val;
        if (dt == NDT_INT64) {
            int64_t v = ndt_get_i(in, i, dt);
            val = (which == NDSP_POSITIVE)    ? (v > 0)
                : (which == NDSP_NEGATIVE)    ? (v < 0)
                : (which == NDSP_NONNEGATIVE) ? (v >= 0)
                                              : (v <= 0);
        } else {
            double re, im;
            ndt_get(in, i, dt, &re, &im);
            /* An indeterminate element is unordered: Positive[Indeterminate]
             * stays symbolic rather than answering False, so abandon the array
             * and let the List path reproduce that element for element. */
            if (re != re) { expr_free(out); return NULL; }
            val = (which == NDSP_POSITIVE)    ? (re > 0.0)
                : (which == NDSP_NEGATIVE)    ? (re < 0.0)
                : (which == NDSP_NONNEGATIVE) ? (re >= 0.0)
                                              : (re <= 0.0);
        }
        ob[i] = val ? 1 : 0;
    }
    return out;
}

/* ---- PrimePowerQ ---------------------------------------------------------- */

Expr* ndint_prime_power_predicate(Expr* res) {
    if (!res || res->type != EXPR_FUNCTION) return NULL;
    /* An option argument (GaussianIntegers -> ...) changes the question; let
     * the element-by-element path answer it. */
    if (res->data.function.arg_count != 1) return NULL;
    const Expr* a = res->data.function.args[0];
    if (!is_ndarray(a)) return NULL;
    /* Only an exact integer buffer has an answer here: a real or complex
     * element draws PrimePowerQ::exact, which the List path emits per element. */
    if (a->data.ndarray.dtype != NDT_INT64) return NULL;
    if (a->data.ndarray.rank != 1) return NULL;
    size_t sz = ndarray_size(a);
    if (sz == 0) return NULL;
    const void* in = a->data.ndarray.data;

    /* Same output contract as ndint_sign_predicate: a packed NDT_BOOL List. */
    void* buf = NULL;
    int64_t dims[1] = { (int64_t)sz };
    Expr* out = ndbuild_open_like(a, 1, dims, NDT_BOOL, &buf);
    if (!out) return NULL;
    out->data.ndarray.present_as = NDA_HEAD_LIST;
    pack_g_any_created = true;
    uint8_t* ob = (uint8_t*)buf;
    for (size_t i = 0; i < sz; i++) {
        int64_t v = ndt_get_i(in, i, NDT_INT64);
        /* Magnitude in uint64, so INT64_MIN = -2^63 is answered (True). */
        uint64_t m = v < 0 ? (uint64_t)0 - (uint64_t)v : (uint64_t)v;
        ob[i] = mth_prime_power_base_u64(m, NULL) ? 1 : 0;
    }
    return out;
}

/* ---- PrimeQ ---------------------------------------------------------------- *
 * PrimeQ[buffer]: a packed bool List, as PrimePowerQ. A dense buffer of small
 * values (max |n| <= 2^26 and at least max/64 elements -- Range-like) is
 * answered off one byte sieve; anything else per element by the deterministic
 * Miller-Rabin the scalar builtin uses (mth_is_prime_u64). Sign ignored, as in
 * the scalar PrimeQ. Only rank-1 int64 buffers; the caller degrades the rest. */
Expr* ndint_primeq(Expr* res) {
    if (!res || res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) return NULL;
    const Expr* a = res->data.function.args[0];
    if (!is_ndarray(a) || a->data.ndarray.dtype != NDT_INT64 || a->data.ndarray.rank != 1)
        return NULL;
    size_t sz = ndarray_size(a);
    if (sz == 0) return NULL;
    const int64_t* in = (const int64_t*)a->data.ndarray.data;

    uint64_t nmax = 0;
    for (size_t i = 0; i < sz; i++) {
        uint64_t m = in[i] < 0 ? (uint64_t)0 - (uint64_t)in[i] : (uint64_t)in[i];
        if (m > nmax) nmax = m;
    }
    uint8_t* comp = NULL;
    if (nmax <= ((uint64_t)1 << 26) && (uint64_t)sz >= nmax / 64) {
        comp = calloc((size_t)nmax + 1, 1);
        if (comp) {
            comp[0] = 1; if (nmax >= 1) comp[1] = 1;
            for (uint64_t p = 2; p * p <= nmax; p++)
                if (!comp[p]) for (uint64_t m = p * p; m <= nmax; m += p) comp[m] = 1;
        }
    }

    void* buf = NULL;
    int64_t dims[1] = { (int64_t)sz };
    Expr* out = ndbuild_open_like(a, 1, dims, NDT_BOOL, &buf);
    if (!out) { free(comp); return NULL; }
    out->data.ndarray.present_as = NDA_HEAD_LIST;
    pack_g_any_created = true;
    uint8_t* ob = (uint8_t*)buf;
    for (size_t i = 0; i < sz; i++) {
        uint64_t m = in[i] < 0 ? (uint64_t)0 - (uint64_t)in[i] : (uint64_t)in[i];
        ob[i] = comp ? !comp[m] : (mth_is_prime_u64(m) ? 1 : 0);
    }
    free(comp);
    return out;
}

/* ---- CoprimeQ / Divisible ---------------------------------------------- *
 * Binary integer predicates: read an int64 buffer against an exact Integer (on
 * either side), or two equal-shape int64 buffers, and answer with the same
 * packed NDT_BOOL List as the sign predicates and PrimePowerQ. Declines -- so
 * the List path answers element by element, with its messages -- a zero
 * divisor (Divisible::divz), INT64_MIN (no |n|), and every other shape. */
static int64_t ndi_gcd_abs(int64_t a, int64_t b) {
    uint64_t x = a < 0 ? (uint64_t)0 - (uint64_t)a : (uint64_t)a;
    uint64_t y = b < 0 ? (uint64_t)0 - (uint64_t)b : (uint64_t)b;
    while (y) { uint64_t t = x % y; x = y; y = t; }
    return (int64_t)x;                     /* x <= 2^63 only for INT64_MIN pairs, declined */
}

Expr* ndint_binary_predicate(Expr* res, NDBinPred which) {
    if (!res || res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) return NULL;
    const Expr* a0 = res->data.function.args[0];
    const Expr* a1 = res->data.function.args[1];
    bool nd0 = is_ndarray(a0), nd1 = is_ndarray(a1);
    if (!nd0 && !nd1) return NULL;
    const Expr* shape = nd0 ? a0 : a1;
    if (nd0 && a0->data.ndarray.dtype != NDT_INT64) return NULL;
    if (nd1 && a1->data.ndarray.dtype != NDT_INT64) return NULL;
    if (!nd0 && a0->type != EXPR_INTEGER) return NULL;
    if (!nd1 && a1->type != EXPR_INTEGER) return NULL;
    size_t sz = ndarray_size(shape);
    if (sz == 0) return NULL;
    if (nd0 && nd1) {
        if (a0->data.ndarray.rank != a1->data.ndarray.rank) return NULL;
        for (int d = 0; d < a0->data.ndarray.rank; d++)
            if (a0->data.ndarray.dims[d] != a1->data.ndarray.dims[d]) return NULL;
    }
    const int64_t* in0 = nd0 ? (const int64_t*)a0->data.ndarray.data : NULL;
    const int64_t* in1 = nd1 ? (const int64_t*)a1->data.ndarray.data : NULL;
    int64_t s0 = nd0 ? 0 : a0->data.integer, s1 = nd1 ? 0 : a1->data.integer;

    /* Validate first, so a decline costs no allocation. */
    for (size_t i = 0; i < sz; i++) {
        int64_t x = nd0 ? in0[i] : s0, y = nd1 ? in1[i] : s1;
        if (x == INT64_MIN || y == INT64_MIN) return NULL;
        if (which == NDBP_DIVISIBLE && y == 0) return NULL;
        if (!nd0 && !nd1) break;
    }

    void* buf = NULL;
    Expr* out = ndbuild_open_like(shape, shape->data.ndarray.rank,
                                  shape->data.ndarray.dims, NDT_BOOL, &buf);
    if (!out) return NULL;
    out->data.ndarray.present_as = NDA_HEAD_LIST;   /* a List of True/False, as the sign predicates */
    pack_g_any_created = true;
    uint8_t* ob = (uint8_t*)buf;
    for (size_t i = 0; i < sz; i++) {
        int64_t x = nd0 ? in0[i] : s0, y = nd1 ? in1[i] : s1;
        ob[i] = (which == NDBP_DIVISIBLE)
                    ? (y == -1 || x % y == 0)
                    : (ndi_gcd_abs(x, y) == 1);
    }
    return out;
}

/* ---- registration ------------------------------------------------------- */

void ndinteger_init(void) {
#ifdef NDI_BITSET_AVX2
    __builtin_cpu_init();
    if (__builtin_cpu_supports("avx2")) {
        ndi_bitset_sn_fn = ndi_bitset_sn_avx2;
        ndi_bitset_an_fn = ndi_bitset_an_avx2;
    }
#endif
    symtab_set_ndarray_unary_kernel("IntegerLength", &NDKU_IntegerLength);
    symtab_set_ndarray_unary_kernel("BitLength",     &NDKU_BitLength);
    symtab_set_ndarray_unary_kernel("BitNot",        &NDKU_BitNot);
    symtab_set_ndarray_unary_kernel("EulerPhi",      &NDKU_EulerPhi);
    symtab_set_ndarray_unary_kernel("MoebiusMu",     &NDKU_MoebiusMu);
    symtab_set_ndarray_unary_kernel("PartitionsP",   &NDKU_PartitionsP);
    symtab_set_ndarray_unary_kernel("PartitionsQ",   &NDKU_PartitionsQ);
    symtab_set_ndarray_binary_kernel("GCD",          &NDKB_GCD);
    symtab_set_ndarray_binary_kernel("BitXor",       &NDKB_BitXor);
    symtab_set_ndarray_binary_kernel("BitAnd",       &NDKB_BitAnd);
    symtab_set_ndarray_binary_kernel("BitOr",        &NDKB_BitOr);
    symtab_set_ndarray_binary_kernel("BitSet",       &NDKB_BitSet);
    symtab_set_ndarray_binary_kernel("LCM",          &NDKB_LCM);
    symtab_set_ndarray_binary_kernel("DivisorSigma", &NDKB_DivisorSigma);
    symtab_set_ndarray_binary_kernel("JacobiSymbol", &NDKB_JacobiSymbol);
    symtab_set_ndarray_binary_kernel("IntegerExponent", &NDKB_IntegerExponent);
}
