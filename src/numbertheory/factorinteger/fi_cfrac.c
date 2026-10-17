/*
 * fi_cfrac.c — the continued-fraction method (Morrison & Brillhart, 1975).
 *
 * The partial quotients of the continued fraction of sqrt(kN) give
 * A_{i-1}^2 == (-1)^i Q_i (mod N) with 0 < Q_i < 2 sqrt(kN), far smaller
 * than the residues of Dixon's method, so a useful fraction of them factor
 * over a base of small primes.  Ingredients:
 *   - Knuth–Schroeppel multiplier k, chosen to maximise the expected
 *     contribution of small primes to Q_i (primes with (kN/p) = 1 divide
 *     some Q_i; p | k always does);
 *   - factor base of the primes with (kN/p) != -1, sized by digit count;
 *   - batched trial division (one bignum remainder per word-sized product
 *     of base primes) with early abort when the cofactor after the first
 *     quarter of the base is still too large to plausibly finish;
 *   - the single large-prime variation (partials pair up in fi_gf2.c);
 *   - GF(2) elimination and square roots, shared (fi_gf2.c).
 *
 * fi_fb_factor() is shared with Dixon's method.
 */
#include "fi_internal.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

/* Factor v (> 0) over the factor base fb[1..nfb] (fb[0] is the sign column
 * and not a prime).  Writes the sorted (col, exp) list starting after any
 * entries already in col/exps (count in *cnt), and the large-prime cofactor
 * in *L (1 when v was smooth).  Returns 1 if v is smooth or a partial with
 * L <= lp_max, 0 otherwise.  `abort_bits`: after the first quarter of the
 * base, give up when the cofactor still has more bits than this. */
int fi_fb_factor(mpz_t v, const uint32_t* fb, size_t nfb, uint32_t* col, uint32_t* exps,
                 size_t* cnt, uint64_t* L, uint64_t lp_max, size_t abort_bits) {
    size_t k = *cnt;
    size_t quarter = nfb / 4 + 1;
    mp_bitcnt_t z = mpz_scan1(v, 0);
    size_t c = 1;
    if (fb[1] == 2) {
        if (z) { mpz_tdiv_q_2exp(v, v, z); col[k] = 1; exps[k++] = (uint32_t)z; }
        c = 2;
    }
    while (c <= nfb) {
        /* Batch base primes into a word-sized product. */
        uint64_t M = 1;
        size_t c0 = c;
        while (c <= nfb) {
            uint64_t h, l;
            fi_mul64(M, fb[c], &h, &l);
            if (h) break;
            M = l;
            c++;
        }
        uint64_t r = mpz_fdiv_ui(v, (unsigned long)M);
        for (size_t j = c0; j < c; j++) {
            if (r % fb[j]) continue;
            uint32_t e = 0;
            while (mpz_divisible_ui_p(v, fb[j])) { mpz_divexact_ui(v, v, fb[j]); e++; }
            if (e) { col[k] = (uint32_t)j; exps[k++] = e; }
        }
        if (mpz_cmp_ui(v, 1) == 0) break;
        if (c0 <= quarter && c > quarter && mpz_sizeinbase(v, 2) > abort_bits) { *cnt = k; return 0; }
    }
    *cnt = k;
    if (mpz_cmp_ui(v, 1) == 0) { *L = 1; return 1; }
    if (mpz_sizeinbase(v, 2) <= 64) {
        uint64_t rem = mpz_get_ui(v);
        if (rem <= lp_max) { *L = rem; return 1; }   /* < pmax^2: prime */
    }
    return 0;
}

/* Jacobi symbol (a/m) for odd m > 0, word arithmetic. */
static int jacobi_u32(uint32_t a, uint32_t m) {
    int t = 1;
    a %= m;
    while (a) {
        while (!(a & 1)) {
            a >>= 1;
            uint32_t r = m & 7;
            if (r == 3 || r == 5) t = -t;
        }
        uint32_t tmp = a; a = m; m = tmp;
        if ((a & 3) == 3 && (m & 3) == 3) t = -t;
        a %= m;
    }
    return m == 1 ? t : 0;
}

/* Knuth–Schroeppel multiplier (shared with the quadratic sieve): maximise
 * f(k) = -log(k)/2 + sum_p g(p, kN) log p over the odd primes p < 1000 and
 * p = 2, where g = 2/(p-1) when (kN/p) = 1 and 1/p when p | k.  (kN/p) is
 * computed as (k/p)(N/p), with (N/p) taken once per prime. */
unsigned fi_ks_multiplier(const mpz_t n) {
    static const unsigned K[] = { 1, 2, 3, 5, 6, 7, 10, 11, 13, 14, 15, 17, 19, 21, 22, 23,
                                  26, 29, 30, 31, 33, 34, 35, 37, 38, 39, 41, 42, 43, 46, 47,
                                  51, 53, 55, 57, 58, 59, 61, 62, 65, 66, 67, 69, 70, 71, 73 };
    size_t np;
    const uint32_t* P = fi_small_primes(&np);
    enum { MAXP = 200 };                 /* odd primes below ~1230 */
    int legn[MAXP];
    double w[MAXP], wk[MAXP];
    uint32_t pr[MAXP];
    size_t m = 0;
    for (size_t i = 1; i < np && m < MAXP; i++) {
        pr[m] = P[i];
        legn[m] = jacobi_u32((uint32_t)mpz_fdiv_ui(n, P[i]), P[i]);
        w[m] = 2.0 * log((double)P[i]) / (P[i] - 1);
        wk[m] = log((double)P[i]) / P[i];
        m++;
    }
    unsigned n8 = (unsigned)mpz_fdiv_ui(n, 8);
    unsigned best = 1;
    double bs = -1e9;
    for (size_t t = 0; t < sizeof K / sizeof K[0]; t++) {
        unsigned k = K[t];
        double s = -0.5 * log((double)k);
        unsigned m8 = (n8 * k) % 8;
        s += (m8 == 1) ? 2 * log(2.0) : (m8 == 5) ? log(2.0) : 0.5 * log(2.0);
        for (size_t i = 0; i < m; i++) {
            if (k % pr[i] == 0) { s += wk[i]; continue; }
            if (legn[i] * jacobi_u32(k, pr[i]) == 1) s += w[i];
        }
        if (s > bs) { bs = s; best = k; }
    }
    return best;
}

/* Factor-base size by digit count (interpolated). */
size_t fi_cfrac_fb_size(size_t digits) {
    static const struct { size_t d, f; } T[] = {
        { 10, 40 }, { 20, 100 }, { 25, 160 }, { 30, 260 }, { 35, 420 }, { 40, 650 },
        { 45, 1000 }, { 50, 1500 }, { 55, 2200 }, { 60, 3200 }, { 65, 4500 }, { 70, 5800 }
    };
    size_t nt = sizeof T / sizeof T[0];
    if (digits <= T[0].d) return T[0].f;
    for (size_t i = 1; i < nt; i++)
        if (digits <= T[i].d)
            return T[i - 1].f + (T[i].f - T[i - 1].f) * (digits - T[i - 1].d) / (T[i].d - T[i - 1].d);
    return T[nt - 1].f;
}

/* Build fb[1..] = primes p with (kn/p) != -1; returns count. */
size_t fi_build_fb(const mpz_t kn, size_t want, uint32_t* fb) {
    size_t np, c = 1;
    const uint32_t* P = fi_small_primes(&np);
    for (size_t i = 0; i < np && c <= want; i++) {
        uint32_t p = P[i];
        if (p == 2) { fb[c++] = 2; continue; }
        if (mpz_kronecker_ui(kn, p) != -1) fb[c++] = p;
    }
    return c - 1;
}

int fi_split_cfrac(mpz_t d, const mpz_t n, FiCtx* ctx) {
    size_t digits = mpz_sizeinbase(n, 10);
    unsigned k = fi_ks_multiplier(n);
    mpz_t kn, g, P, Q, Qp, q, t, A, Ap, v;
    mpz_inits(kn, g, P, Q, Qp, q, t, A, Ap, v, NULL);
    mpz_mul_ui(kn, n, k);
    mpz_sqrt(g, kn);

    size_t want = fi_cfrac_fb_size(digits);
    uint32_t* fb = malloc((want + 2) * sizeof(uint32_t));
    fb[0] = 0;
    size_t nfb = fi_build_fb(kn, want, fb);
    /* A base prime dividing n itself is a factor. */
    int found = 0;
    for (size_t i = 1; i <= nfb && !found; i++)
        if (mpz_divisible_ui_p(n, fb[i]) && mpz_cmp_ui(n, fb[i]) != 0) { mpz_set_ui(d, fb[i]); found = 1; }

    FiRels* R = fi_rels_new(fb, nfb + 1);
    uint32_t* col = malloc((nfb + 2) * sizeof(uint32_t));
    uint32_t* exps = malloc((nfb + 2) * sizeof(uint32_t));
    uint64_t pmax = fb[nfb];
    uint64_t lp_max = pmax * 128 < pmax * pmax ? pmax * 128 : pmax * pmax;
    size_t qbits = mpz_sizeinbase(g, 2) + 1;
    size_t abort_bits = (size_t)(0.70 * qbits) + 8;
    uint64_t budget = ctx->iterations ? ctx->iterations : (UINT64_C(1) << 32);

    /* i = 0: P0 = 0, Q0 = 1, q0 = g, A_{-1} = 1, A_0 = g. */
    mpz_set_ui(P, 0);
    mpz_set_ui(Qp, 1);              /* Q_{i-1} */
    mpz_set(q, g);
    mpz_set_ui(Ap, 1);              /* A_{i-2} */
    mpz_mod(A, g, n);               /* A_{i-1} */
    size_t need = nfb + 1 + 48;
    for (uint64_t i = 1; i < budget && !found; i++) {
        /* P_i = q_{i-1} Q_{i-1} - P_{i-1};  Q_i = (kN - P_i^2) / Q_{i-1}. */
        mpz_mul(t, q, Qp);
        mpz_sub(P, t, P);
        mpz_mul(t, P, P);
        mpz_sub(t, kn, t);
        mpz_divexact(Q, t, Qp);
        if (mpz_sgn(Q) == 0) break;
        /* Relation A_{i-1}^2 == (-1)^i Q_i (mod N). */
        mpz_set(v, Q);
        size_t cnt = 0;
        if (i & 1) { col[0] = 0; exps[0] = 1; cnt = 1; }
        uint64_t L;
        if (fi_fb_factor(v, fb, nfb, col, exps, &cnt, &L, lp_max, abort_bits)) {
            fi_rels_add(R, A, col, exps, cnt, L, n);
            if (fi_rels_count(R) >= need) {
                found = fi_rels_solve(R, n, d);
                if (!found) need += 32;      /* every dependency trivial: gather more */
            }
        }
        /* q_i = floor((g + P_i) / Q_i);  A_i = q_i A_{i-1} + A_{i-2}. */
        mpz_add(t, g, P);
        mpz_fdiv_q(q, t, Q);
        mpz_mul(t, q, A);
        mpz_add(t, t, Ap);
        mpz_set(Ap, A);
        mpz_mod(A, t, n);
        mpz_swap(Qp, Q);
        /* A square Q at an even step also gives a direct (SQUFOF-like) try. */
    }
    fi_rels_free(R);
    free(col); free(exps); free(fb);
    mpz_clears(kn, g, P, Q, Qp, q, t, A, Ap, v, NULL);
    return found;
}
