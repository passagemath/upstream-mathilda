/*
 * fi_rbd.c — Blake rational base descent.
 *
 * Moved verbatim from facint.c (only the option plumbing changed: the
 * "Base" -> a/b suboption now arrives in FiCtx).  Its optimisation is
 * deliberately deferred to a later session.
 *
 * Idea: assume a factor p ~ c (a/b)^k for a rational base a/b; descend
 * Q <- floor(b Q / a) from n and probe gcd(Q +- j, n) for small j.
 */
#include "fi_internal.h"

/* Try a single (a, b) pair.  Returns 1 if a non-trivial factor is set in
 * `factor`, 0 otherwise.  Requires a > b > 0. */
static int rbd_try_pair(mpz_t factor, const mpz_t n, const mpz_t a, const mpz_t b) {
    if (mpz_sgn(b) <= 0 || mpz_cmp(a, b) <= 0) return 0;

    mpz_t diff, search_limit_z, tmp;
    mpz_inits(diff, search_limit_z, tmp, NULL);
    mpz_sub(diff, a, b);

    /* search_limit ~ (2a - b - 1)/(a - b) + 2: small for a/b far from 1,
     * larger when the descent step is gentle. */
    mpz_add(tmp, a, diff);
    mpz_sub_ui(tmp, tmp, 1);
    mpz_fdiv_q(search_limit_z, tmp, diff);
    mpz_add_ui(search_limit_z, search_limit_z, 2);

    if (!mpz_fits_ulong_p(search_limit_z)) {
        mpz_clears(diff, search_limit_z, tmp, NULL);
        return 0;
    }
    unsigned long search_limit = mpz_get_ui(search_limit_z);

    mpz_t Q, g;
    mpz_init_set(Q, n);
    mpz_init(g);

    int found = 0;
    while (1) {
        mpz_mul(tmp, b, Q);
        mpz_fdiv_q(Q, tmp, a);

        if (mpz_cmp_ui(Q, 2) < 0) break;

        for (unsigned long j = 0; j <= search_limit; j++) {
            mpz_sub_ui(tmp, Q, j);
            mpz_gcd(g, tmp, n);
            if (mpz_cmp_ui(g, 1) > 0 && mpz_cmp(g, n) < 0) {
                mpz_set(factor, g);
                found = 1;
                goto done;
            }

            if (j > 0) {
                mpz_add_ui(tmp, Q, j);
                mpz_gcd(g, tmp, n);
                if (mpz_cmp_ui(g, 1) > 0 && mpz_cmp(g, n) < 0) {
                    mpz_set(factor, g);
                    found = 1;
                    goto done;
                }
            }
        }
    }

done:
    mpz_clears(diff, search_limit_z, tmp, Q, g, NULL);
    return found;
}

int fi_split_rbd(mpz_t factor, const mpz_t n, FiCtx* c) {
    mpz_set_ui(factor, 0);
    if (mpz_cmp_ui(n, 4) < 0) return 0;
    if (mpz_even_p(n)) { mpz_set_ui(factor, 2); return 1; }
    if (mpz_perfect_square_p(n)) { mpz_sqrt(factor, n); return 1; }

    if (c->rbd_have_base) return rbd_try_pair(factor, n, c->rbd_a, c->rbd_b);

    mpz_t a, b;
    mpz_inits(a, b, NULL);

    /* No explicit base: walk j = 3, 4, ... and try every coprime partition
     * a + b = j with a > b >= 1.  Mirrors the reference Python
     * factor_rational_base() search, but bounded by a pair budget so a hard
     * composite cannot hang the REPL (the reference runs until interrupted). */
    const unsigned long RBD_MAX_PAIRS = 200000UL;
    unsigned long pairs_tried = 0;
    int found = 0;
    for (unsigned long j = 3; pairs_tried < RBD_MAX_PAIRS && !found; j++) {
        for (unsigned long bb = 1; 2 * bb < j; bb++) {
            unsigned long aa = j - bb;
            unsigned long u = aa, v = bb;
            while (v) { unsigned long t = u % v; u = v; v = t; }
            if (u != 1) continue;

            mpz_set_ui(a, aa);
            mpz_set_ui(b, bb);
            if (rbd_try_pair(factor, n, a, b)) { found = 1; break; }
            if (++pairs_tried >= RBD_MAX_PAIRS) break;
        }
    }
    mpz_clears(a, b, NULL);
    return found;
}
