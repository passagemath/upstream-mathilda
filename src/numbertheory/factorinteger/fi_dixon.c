/*
 * fi_dixon.c — Dixon's random-squares method.
 *
 * Dixon draws random x and keeps those whose square residue z == x^2
 * (mod N) factors over a base of small primes; a GF(2) dependency among the
 * exponent vectors gives X^2 == Y^2 (mod N).  Textbook Dixon takes x
 * uniformly, so z is as large as N and almost never smooth.  Here each x is
 * drawn as floor(sqrt(kN)) + j with random k < 2^16 and random j < 2^20, so
 * z = x^2 - kN is only O(j sqrt(kN)) — the residues are square-root sized,
 * which is what makes the method usable at all, while the choices stay
 * random and every residue is trial divided (no sieve: that refinement is
 * the quadratic sieve).  Since k varies, every prime can divide some z, so
 * the base is all primes up to the bound.
 *
 * Shares batched factor-base trial division with early abort, the single
 * large-prime variation, linear algebra and square roots with CFRAC.
 */
#include "fi_internal.h"

#include <stdlib.h>

int    fi_fb_factor(mpz_t v, const uint32_t* fb, size_t nfb, uint32_t* col, uint32_t* exps,
                    size_t* cnt, uint64_t* L, uint64_t lp_max, size_t abort_bits);
size_t fi_cfrac_fb_size(size_t digits);

int fi_split_dixon(mpz_t d, const mpz_t n, FiCtx* ctx) {
    size_t digits = mpz_sizeinbase(n, 10);
    /* All primes, so a base about 1.5x CFRAC's covers the same QR primes. */
    size_t want = fi_cfrac_fb_size(digits) * 3 / 2;
    size_t np;
    const uint32_t* P = fi_small_primes(&np);
    if (want > np - 1) want = np - 1;
    uint32_t* fb = malloc((want + 2) * sizeof(uint32_t));
    fb[0] = 0;
    for (size_t i = 0; i < want; i++) fb[i + 1] = P[i];
    size_t nfb = want;

    int found = 0;
    for (size_t i = 1; i <= nfb && !found; i++)
        if (mpz_divisible_ui_p(n, fb[i]) && mpz_cmp_ui(n, fb[i]) != 0) { mpz_set_ui(d, fb[i]); found = 1; }

    FiRels* R = fi_rels_new(fb, nfb + 1);
    uint32_t* col = malloc((nfb + 2) * sizeof(uint32_t));
    uint32_t* exps = malloc((nfb + 2) * sizeof(uint32_t));
    uint64_t pmax = fb[nfb];
    uint64_t lp_max = pmax * 128;
    mpz_t kn, x, z, s;
    mpz_inits(kn, x, z, s, NULL);
    size_t zbits = mpz_sizeinbase(n, 2) / 2 + 40;
    size_t abort_bits = (size_t)(0.70 * zbits) + 8;
    uint64_t budget = ctx->iterations ? ctx->iterations : (UINT64_C(1) << 30);
    size_t need = nfb + 48;
    for (uint64_t it = 0; it < budget && !found; it++) {
        uint64_t r = fi_rand(ctx);
        unsigned long k = 1 + (unsigned long)(r & 0xffff);
        unsigned long j = 1 + (unsigned long)((r >> 16) & 0xfffff);
        mpz_mul_ui(kn, n, k);
        mpz_sqrt(s, kn);
        mpz_add_ui(x, s, j);
        mpz_mul(z, x, x);
        mpz_sub(z, z, kn);                  /* z = x^2 - kN > 0, z == x^2 (mod N) */
        size_t cnt = 0;
        uint64_t L;
        if (!fi_fb_factor(z, fb, nfb, col, exps, &cnt, &L, lp_max, abort_bits)) continue;
        if (L > 1 && mpz_divisible_ui_p(n, (unsigned long)L) && mpz_cmp_ui(n, (unsigned long)L) > 0) {
            mpz_set_ui(d, (unsigned long)L);    /* a prime of n as the large prime */
            found = 1;
            break;
        }
        mpz_mod(x, x, n);
        fi_rels_add(R, x, col, exps, cnt, L, n);
        if (fi_rels_count(R) >= need) {
            found = fi_rels_solve(R, n, d);
            if (!found) need += 32;
        }
    }
    mpz_clears(kn, x, z, s, NULL);
    fi_rels_free(R);
    free(col); free(exps); free(fb);
    return found;
}
