/*
 * fi_rbd_modacc.c — the division-free window accumulator of the exact
 * descent: a running product of candidates modulo an odd N, kept in
 * Montgomery-reduced form on fixed-size limb arrays.
 *
 * The batched descent only ever asks one question of its accumulator:
 * is gcd(prod, N) non-trivial?  Every Montgomery reduction (REDC) divides
 * by a power of R = 2^GMP_NUMB_BITS, and R is a unit modulo an odd N, so a
 * REDC changes the accumulator by a unit factor and leaves that gcd
 * untouched.  The accumulator therefore never converts in or out of
 * Montgomery form, and folding an s-limb factor x costs exactly
 *
 *     mpn_mul(acc, n, x, s)  +  s passes of mpn_addmul_1(N, n)
 *
 * — about 2 n s limb products and no division.  The mpz path this replaces
 * paid two schoolbook divisions (mpz_tdiv_r) per step on top of the
 * multiplication, which was ~45% of the whole search.  A factor wider than N
 * is first shrunk to n limbs by a partial REDC, again only a unit factor.
 *
 * Plain GMP mpn, no allocation, no global state: reentrant.
 */
#include "fi_rbd.h"

#include <stdlib.h>
#include <string.h>

#if GMP_NAIL_BITS != 0
#error "fi_rbd_modacc.c assumes a nail-free GMP"
#endif

int fi_rbd_modctx_init(fi_rbd_modctx_t* c, const mpz_t N) {
    c->N = NULL;
    c->n = 0;
    c->ninv = 0;
    if (mpz_sgn(N) <= 0 || mpz_even_p(N) || mpz_cmp_ui(N, 1) == 0) return 0;
    mp_size_t n = (mp_size_t)mpz_size(N);
    mp_limb_t* p = malloc((size_t)n * sizeof(mp_limb_t));
    if (!p) return 0;
    memcpy(p, mpz_limbs_read(N), (size_t)n * sizeof(mp_limb_t));
    /* Newton iteration for N0^{-1} mod R: x = N0 is right to 3 bits for an
     * odd N0, and each step doubles the number of correct bits. */
    mp_limb_t n0 = p[0], x = n0;
    for (int bits = 3; bits < GMP_NUMB_BITS; bits *= 2) x *= 2 - n0 * x;
    c->N = p;
    c->n = n;
    c->ninv = (mp_limb_t)0 - x;
    return 1;
}

void fi_rbd_modctx_clear(fi_rbd_modctx_t* c) {
    free(c->N);
    c->N = NULL;
    c->n = 0;
}

size_t fi_rbd_modacc_work(const fi_rbd_modctx_t* c, size_t s) {
    size_t n = (size_t)c->n;
    /* shrink buffer (s + 1) + product buffer (2n + 1) */
    return s + 2 * n + 4;
}

void fi_rbd_modacc_one(const fi_rbd_modctx_t* c, mp_limb_t* acc) {
    memset(acc, 0, (size_t)c->n * sizeof(mp_limb_t));
    acc[0] = 1;
}

/* REDC the low k limbs out of t[0..tn] (t[tn] must be a zero guard limb,
 * k + n <= tn): t <- (t + m N) / R^k with m chosen limb by limb so that the
 * low k limbs vanish.  The result is left at t + k, tn + 1 - k limbs. */
static void redc_partial(mp_limb_t* t, mp_size_t tn, mp_size_t k, const fi_rbd_modctx_t* c) {
    const mp_limb_t* N = c->N;
    mp_size_t n = c->n;
    for (mp_size_t i = 0; i < k; i++) {
        mp_limb_t q = t[i] * c->ninv;
        mp_limb_t cy = mpn_addmul_1(t + i, N, n, q);
        mpn_add_1(t + i + n, t + i + n, tn + 1 - (i + n), cy);
    }
}

void fi_rbd_modacc_fold(const fi_rbd_modctx_t* c, mp_limb_t* acc, const mp_limb_t* x,
                        mp_size_t s, mp_limb_t* work) {
    const mp_limb_t* N = c->N;
    mp_size_t n = c->n;
    while (s > 0 && x[s - 1] == 0) s--;
    if (s == 0) {                           /* a zero factor: acc = 0, gcd = N */
        memset(acc, 0, (size_t)n * sizeof(mp_limb_t));
        return;
    }
    mp_limb_t* t = work;
    if (s > n) {
        /* x < R^s  =>  (x + m N) / R^(s-n) < R^n + N: n limbs plus a top
         * bit, which one subtraction of N clears. */
        mp_limb_t* w = work;                 /* w[0..s]; the product goes after it */
        memcpy(w, x, (size_t)s * sizeof(mp_limb_t));
        w[s] = 0;
        redc_partial(w, s, s - n, c);
        t = w + s + 1;
        mp_limb_t* r = w + (s - n);
        if (r[n]) mpn_sub_n(r, r, N, n);
        x = r;
        s = n;
        while (s > 0 && x[s - 1] == 0) s--;
        if (s == 0) { memset(acc, 0, (size_t)n * sizeof(mp_limb_t)); return; }
    }
    /* acc < N, x < R^s  =>  acc x < N R^s, and the REDC by s limbs leaves a
     * value below 2N in n + 1 limbs. */
    mpn_mul(t, acc, n, x, s);
    t[n + s] = 0;
    redc_partial(t, n + s, s, c);
    mp_limb_t* r = t + s;
    if (r[n] || mpn_cmp(r, N, n) >= 0) mpn_sub_n(acc, r, N, n);
    else                               memcpy(acc, r, (size_t)n * sizeof(mp_limb_t));
}
