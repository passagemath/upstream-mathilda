/*
 * fi_rbd_bases.c — enumeration of rational bases for Blake rational base
 * descent.
 *
 *   Stern-Brocot  breadth-first over the right half of the Stern-Brocot tree
 *                 (the rationals > 1), pruned at numerator H.  Every
 *                 rational appears exactly once and already in lowest terms,
 *                 so no coprimality test is needed, and the order is coarse to
 *                 fine: each level refines the mediant intervals of the last
 *                 ([RBD] Section 4.1, [LRBD] Section 6).  Mediants only grow
 *                 numerators, so a node with a > H prunes its whole subtree.
 *   Sum           the reference order of [RBD] Appendix A: a + b ascending,
 *                 b ascending, coprime pairs only.
 *
 * Both apply the perfect-power filter: a/b = (x/y)^k is redundant with x/y.
 * fi_rbd_simplest_between gives the least-height rational in an interval,
 * which realises a grid point of the lattice pass P2.
 */
#include "fi_rbd.h"

#include <math.h>
#include <stdlib.h>

/* x^k <= lim?  Writes x^k to *out when it does not exceed lim. */
static int pow_le(uint64_t x, unsigned k, uint64_t lim, uint64_t* out) {
    uint64_t r = 1;
    for (unsigned i = 0; i < k; i++) {
        if (x != 0 && r > lim / x) return 0;
        r *= x;
    }
    *out = r;
    return 1;
}

/* The integer k-th root of n when n is a perfect k-th power, else 0. */
static uint64_t exact_root(uint64_t n, unsigned k) {
    if (n <= 1) return n;
    uint64_t r = (uint64_t)(pow((double)n, 1.0 / (double)k) + 0.5);
    for (uint64_t x = r > 1 ? r - 1 : 1; x <= r + 1; x++) {
        uint64_t p;
        if (pow_le(x, k, n, &p) && p == n) return x;
    }
    return 0;
}

int fi_rbd_is_primitive(uint64_t a, uint64_t b) {
    static const unsigned PK[] = { 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47,
                                   53, 59, 61 };
    if (a < 4) return 1;
    unsigned bits = 0;
    for (uint64_t t = a; t; t >>= 1) bits++;
    /* a = x^k, b = y^k for some k >= 2 iff it holds for a prime k. */
    for (size_t i = 0; i < sizeof PK / sizeof PK[0] && PK[i] < bits; i++) {
        unsigned k = PK[i];
        if (b > 1 && !exact_root(b, k)) continue;
        if (exact_root(a, k)) return 0;
    }
    return 1;
}

/* ------------------------------------------------------------------ */
/* Stern-Brocot                                                        */
/* ------------------------------------------------------------------ */

static int sb_push(fi_rbd_sb_t* it, uint64_t pl, uint64_t ql, uint64_t pr, uint64_t qr) {
    if (it->len == it->cap) {
        size_t nc = it->cap ? 2 * it->cap : 1024;
        uint64_t* nq = malloc(nc * 4 * sizeof(uint64_t));
        if (!nq) return 0;
        for (size_t i = 0; i < it->len; i++) {
            size_t src = ((it->head + i) % it->cap) * 4;
            for (int c = 0; c < 4; c++) nq[i * 4 + c] = it->q[src + c];
        }
        free(it->q);
        it->q = nq;
        it->cap = nc;
        it->head = 0;
    }
    size_t at = ((it->head + it->len) % it->cap) * 4;
    it->q[at] = pl; it->q[at + 1] = ql; it->q[at + 2] = pr; it->q[at + 3] = qr;
    it->len++;
    return 1;
}

void fi_rbd_sb_init(fi_rbd_sb_t* it, uint64_t H) {
    it->q = NULL;
    it->head = it->len = it->cap = 0;
    it->H = H;
    sb_push(it, 1, 1, 1, 0);               /* (1/1, 1/0): the rationals > 1 */
}

void fi_rbd_sb_clear(fi_rbd_sb_t* it) {
    free(it->q);
    it->q = NULL;
    it->head = it->len = it->cap = 0;
}

size_t fi_rbd_sb_next(fi_rbd_sb_t* it, fi_rbd_base_t* out, size_t max) {
    size_t n = 0;
    while (n < max && it->len) {
        size_t at = it->head * 4;
        uint64_t pl = it->q[at], ql = it->q[at + 1], pr = it->q[at + 2], qr = it->q[at + 3];
        it->head = (it->head + 1) % it->cap;
        it->len--;
        uint64_t a = pl + pr, b = ql + qr;  /* the mediant */
        if (a > it->H || a < pl) continue;  /* prune (and guard wrap-around) */
        if (fi_rbd_is_primitive(a, b)) { out[n].a = a; out[n].b = b; n++; }
        sb_push(it, pl, ql, a, b);
        sb_push(it, a, b, pr, qr);
    }
    return n;
}

/* ------------------------------------------------------------------ */
/* Sum order                                                           */
/* ------------------------------------------------------------------ */

static uint64_t gcd64(uint64_t u, uint64_t v) {
    while (v) { uint64_t t = u % v; u = v; v = t; }
    return u;
}

void fi_rbd_sum_init(fi_rbd_sum_t* it, uint64_t H) {
    it->s = 3;
    it->b = 1;
    it->H = H;
}

size_t fi_rbd_sum_next(fi_rbd_sum_t* it, fi_rbd_base_t* out, size_t max) {
    size_t n = 0;
    while (n < max) {
        if (2 * it->b >= it->s) {          /* next sum */
            it->s++;
            it->b = 1;
            /* The least numerator at sum s is s/2 + 1 (a > b): past H, done. */
            if (it->H && it->s / 2 + 1 > it->H) break;
            continue;
        }
        uint64_t b = it->b++, a = it->s - b;
        if (it->H && a > it->H) continue;
        if (gcd64(a, b) != 1 || !fi_rbd_is_primitive(a, b)) continue;
        out[n].a = a;
        out[n].b = b;
        n++;
    }
    return n;
}

/* ------------------------------------------------------------------ */
/* Simplest rational in an interval                                    */
/* ------------------------------------------------------------------ */

/* Continued-fraction descent: with x = lo, y = hi and f = floor(x), the
 * answer is x itself if x is an integer, f + 1 if f + 1 <= y, and otherwise
 * f + 1 / simplest(1/(y - f), 1/(x - f)).  The iteration keeps the
 * accumulated convergent as a Moebius map (p0 + p1 z) / (q0 + q1 z). */
void fi_rbd_simplest_between(mpz_t p, mpz_t q, mpq_srcptr lo, mpq_srcptr hi) {
    mpq_t x, y, t;
    mpz_t f, p0, p1, q0, q1, u;
    mpq_inits(x, y, t, NULL);
    mpz_inits(f, p0, p1, q0, q1, u, NULL);
    mpq_set(x, lo);
    mpq_set(y, hi);
    mpz_set_ui(p0, 0); mpz_set_ui(p1, 1);  /* z -> z */
    mpz_set_ui(q0, 1); mpz_set_ui(q1, 0);
    for (;;) {
        mpz_fdiv_q(f, mpq_numref(x), mpq_denref(x));
        /* z = f + 1/z': fold the integer part into the map. */
        int done = 0;
        mpz_t zn, zd;                      /* the final z = zn/zd when done */
        mpz_inits(zn, zd, NULL);
        if (mpz_cmp_ui(mpq_denref(x), 1) == 0) {
            mpz_set(zn, mpq_numref(x)); mpz_set_ui(zd, 1); done = 1;
        } else {
            mpq_set_z(t, f);
            mpz_add_ui(mpq_numref(t), mpq_numref(t), 1);   /* f + 1 */
            if (mpq_cmp(t, y) <= 0) { mpz_add_ui(zn, f, 1); mpz_set_ui(zd, 1); done = 1; }
        }
        if (done) {
            /* p/q = (p0 zd + p1 zn) / (q0 zd + q1 zn) */
            mpz_mul(p, p0, zd); mpz_addmul(p, p1, zn);
            mpz_mul(q, q0, zd); mpz_addmul(q, q1, zn);
            mpz_gcd(u, p, q);
            if (mpz_cmp_ui(u, 1) > 0) { mpz_divexact(p, p, u); mpz_divexact(q, q, u); }
            if (mpz_sgn(q) < 0) { mpz_neg(p, p); mpz_neg(q, q); }
            mpz_clears(zn, zd, NULL);
            break;
        }
        mpz_clears(zn, zd, NULL);
        /* Substitute z = f + 1/w:
         *   (p0 + p1 z) / (q0 + q1 z) = (p1 + (p0 + p1 f) w) / (q1 + (q0 + q1 f) w),
         * so (p0, p1) <- (p1, p0 + p1 f), and likewise for q. */
        mpz_addmul(p0, p1, f);
        mpz_swap(p0, p1);
        mpz_addmul(q0, q1, f);
        mpz_swap(q0, q1);
        /* (x, y) <- (1/(y - f), 1/(x - f)) */
        mpq_set_z(t, f);
        mpq_sub(x, x, t);
        mpq_sub(y, y, t);
        mpq_inv(x, x);
        mpq_inv(y, y);
        mpq_swap(x, y);
    }
    mpq_clears(x, y, t, NULL);
    mpz_clears(f, p0, p1, q0, q1, u, NULL);
}
