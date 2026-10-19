/* nt_gaussian.c -- shared Gaussian-integer & divisor enumeration (is_gaussian_integer, df_*, divisors_ordinary, divisors_gaussian).
 * Split from numbertheory.c; see numbertheory.h and
 * numbertheory_internal.h for the subsystem layout. */

#include "numbertheory.h"
#include "numbertheory_internal.h"
#include "arithmetic.h"
#include "eval.h"
#include "sym_names.h"
#include "internal.h"
#include "print.h"
#include "symtab.h"
#include "attr.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <gmp.h>

/* True when e is an exact Gaussian integer Complex[a, b] with both real and
 * imaginary parts integer-like (machine int or BigInt). */
bool is_gaussian_integer(Expr* e) {
    Expr *re, *im;
    if (!is_complex(e, &re, &im)) return false;
    return expr_is_integer_like(re) && expr_is_integer_like(im);
}
/* ---- Divisors ----------------------------------------------------------- */

/* qsort comparator for an array of mpz_t (ascending). */
static int df_mpz_cmp_qsort(const void* a, const void* b) {
    return mpz_cmp(*(const mpz_t*)a, *(const mpz_t*)b);
}

/* A Gaussian integer re + im i, used so divisor pairs can be qsort-ed as a
 * unit.  Relocating an mpz_t by memcpy (as qsort does) is safe: only the
 * struct moves, the heap-allocated limbs it points at do not. */
typedef struct { mpz_t re; mpz_t im; } df_gpair;

/* Sort Gaussian divisors by (Re, Im) ascending — matches Mathematica's
 * canonical ordering of the resulting list. */
static int df_gpair_cmp(const void* x, const void* y) {
    const df_gpair* p = (const df_gpair*)x;
    const df_gpair* q = (const df_gpair*)y;
    int c = mpz_cmp(p->re, q->re);
    if (c != 0) return c;
    return mpz_cmp(p->im, q->im);
}

/* Read a FactorInteger result {{p1, e1}, ...} (borrowed, not freed) into
 * distinct primes and exponents.  On success returns true and sets
 * *out_primes (an array of *out_n init'd mpz_t) and *out_exps (matching
 * exponents); the caller frees both arrays and clears each prime.  A unit
 * entry {1, 1} contributes nothing, so FactorInteger[1] reads as the empty
 * factorisation (n == 0, arrays NULL).  Returns false on any structural
 * problem AND on any base that is not prime -- see the comment in the loop.
 *
 * Split out of df_factor_mpz so that the decline guard can be tested on a
 * crafted incomplete factorisation: FactorInteger's Automatic method now
 * finishes every composite a test can afford to hand it, so no natural input
 * reaches this path in reasonable time. */
bool df_factor_list_read(const Expr* fact, mpz_t** out_primes,
                         unsigned long** out_exps, size_t* out_n) {
    if (!fact || fact->type != EXPR_FUNCTION ||
        fact->data.function.head->type != EXPR_SYMBOL ||
        fact->data.function.head->data.symbol.name != SYM_List)
        return false;

    size_t n = fact->data.function.arg_count;
    mpz_t* primes = (n > 0) ? (mpz_t*)malloc(n * sizeof(mpz_t)) : NULL;
    unsigned long* exps = (n > 0) ? (unsigned long*)malloc(n * sizeof(unsigned long)) : NULL;
    bool ok = true;
    size_t got = 0;
    for (size_t i = 0; i < n; i++) {
        Expr* pair = fact->data.function.args[i];
        if (pair->type != EXPR_FUNCTION || pair->data.function.arg_count != 2 ||
            !expr_is_integer_like(pair->data.function.args[0]) ||
            pair->data.function.args[1]->type != EXPR_INTEGER ||
            pair->data.function.args[1]->data.integer <= 0) {
            ok = false; break;
        }
        /* FactorInteger[1] is {{1, 1}} (Mathematica): a unit entry is not a
         * prime and contributes nothing. */
        expr_to_mpz(pair->data.function.args[0], primes[got]);  /* inits primes[got] */
        if (mpz_cmpabs_ui(primes[got], 1) == 0) { mpz_clear(primes[got]); continue; }
        exps[got] = (unsigned long)pair->data.function.args[1]->data.integer;
        got++;
        /* Every base must actually BE prime.  FactorInteger's Automatic method
         * is bounded, and on a hard composite that survives trial division,
         * Pollard rho and ECM it returns the cofactor unfactored with exponent
         * 1 (it warns, FactorInteger::nofac, but it still answers).  Reading
         * that as a prime factor is how MoebiusMu, PrimeNu, PrimeOmega and
         * LiouvilleLambda came to return confident WRONG numbers rather than
         * declining:
         *
         *     m = NextPrime[10^40] NextPrime[10^41]   (* a 82-digit semiprime *)
         *     PrimeQ[m]       False     (* Mathilda knows it is composite *)
         *     PrimeNu[m]      1         (* true: 2 *)
         *     PrimeOmega[m]   1         (* true: 2 *)
         *     MoebiusMu[m]    -1        (* true: 1 *)
         *
         * facint.c's own facint_factor_complete already states the contract --
         * "the caller must then DECLINE rather than trust a possibly-composite
         * prime" -- but this helper goes through internal_factorinteger and so
         * never saw the flag.  40 Miller-Rabin rounds, the same bound that
         * function uses.  A unit base (+-1) is not a prime and not a failure.
         * Every caller already returns NULL on false, which leaves the head
         * unevaluated: an honest "I cannot" in place of a wrong integer. */
        if (mpz_probab_prime_p(primes[got - 1], 40) == 0) {
            ok = false; break;
        }
    }

    if (!ok) {
        for (size_t i = 0; i < got; i++) mpz_clear(primes[i]);
        free(primes);
        free(exps);
        return false;
    }
    if (got == 0) { free(primes); free(exps); primes = NULL; exps = NULL; }
    *out_primes = primes;
    *out_exps = exps;
    *out_n = got;
    return true;
}

/* Factor a positive integer m into distinct primes and exponents by delegating
 * to internal_factorinteger (FactorInteger) and reading the result with
 * df_factor_list_read (same outputs and ownership).  Returns false on any
 * structural problem or an incomplete factorisation. */
bool df_factor_mpz(const mpz_t m, mpz_t** out_primes,
                          unsigned long** out_exps, size_t* out_n) {
    /* internal_factorinteger takes ownership of the args array contents, so the
     * constructed integer must NOT be freed here. */
    Expr* args[1] = { expr_bigint_normalize(expr_new_bigint_from_mpz(m)) };
    Expr* fact = internal_factorinteger(args, 1);
    bool ok = df_factor_list_read(fact, out_primes, out_exps, out_n);
    if (fact) expr_free(fact);
    return ok;
}

/* Build the ascending list of positive divisors of nabs (>= 2 in practice).
 * Enumerates the divisor lattice with a mixed-radix exponent odometer over the
 * prime factorisation.  Returns a List Expr, or NULL if factoring fails. */
Expr* divisors_ordinary(const mpz_t nabs) {
    mpz_t* primes;
    unsigned long* exps;
    size_t np;
    if (!df_factor_mpz(nabs, &primes, &exps, &np)) return NULL;

    /* Number of divisors = prod(exps[i] + 1).  Bail (leave unevaluated) if this
     * product overflows size_t — e.g. Divisors[100!] has ~10^28 divisors, which
     * is intractable to materialise as a list. */
    size_t ndiv = 1;
    bool overflow = false;
    for (size_t i = 0; i < np; i++) {
        size_t f = exps[i] + 1;
        if (ndiv > SIZE_MAX / f) { overflow = true; break; }
        ndiv *= f;
    }
    if (overflow) {
        for (size_t i = 0; i < np; i++) mpz_clear(primes[i]);
        free(primes);
        free(exps);
        return NULL;
    }

    mpz_t* divs = (mpz_t*)malloc(ndiv * sizeof(mpz_t));
    if (!divs) {
        for (size_t i = 0; i < np; i++) mpz_clear(primes[i]);
        free(primes);
        free(exps);
        return NULL;
    }
    unsigned long* cur = (unsigned long*)calloc(np, sizeof(unsigned long));
    mpz_t acc, pw;
    mpz_inits(acc, pw, NULL);
    for (size_t d = 0; d < ndiv; d++) {
        mpz_set_ui(acc, 1);
        for (size_t i = 0; i < np; i++) {
            mpz_pow_ui(pw, primes[i], cur[i]);
            mpz_mul(acc, acc, pw);
        }
        mpz_init_set(divs[d], acc);
        /* advance mixed-radix odometer: digit i ranges 0..exps[i] */
        for (size_t i = 0; i < np; i++) {
            if (++cur[i] <= exps[i]) break;
            cur[i] = 0;
        }
    }
    mpz_clears(acc, pw, NULL);
    free(cur);

    qsort(divs, ndiv, sizeof(mpz_t), df_mpz_cmp_qsort);

    Expr** items = (Expr**)malloc(ndiv * sizeof(Expr*));
    for (size_t d = 0; d < ndiv; d++)
        items[d] = expr_bigint_normalize(expr_new_bigint_from_mpz(divs[d]));
    Expr* out = expr_new_function(expr_new_symbol(SYM_List), items, ndiv);
    free(items);

    for (size_t d = 0; d < ndiv; d++) mpz_clear(divs[d]);
    free(divs);
    for (size_t i = 0; i < np; i++) mpz_clear(primes[i]);
    free(primes);
    free(exps);
    return out;
}

/* Solve p = c^2 + d^2 (c, d > 0) for a prime p ≡ 1 (mod 4) via Cornacchia's
 * algorithm.  Finds a square root of -1 modulo p (t^((p-1)/4) for a quadratic
 * non-residue t), then Euclidean-reduces (p, root) until the remainder drops
 * below sqrt(p).  Sets c, d (caller pre-inits) and returns true on success. */
static bool df_two_squares(const mpz_t p, mpz_t c, mpz_t d) {
    mpz_t t, x, e, a, b, r, s, sq;
    mpz_inits(t, x, e, a, b, r, s, sq, NULL);

    bool found = false;
    for (unsigned long tv = 2; tv < 1000000UL; tv++) {
        mpz_set_ui(t, tv);
        if (mpz_legendre(t, p) == -1) { found = true; break; }
    }
    bool ok = false;
    if (found) {
        mpz_sub_ui(e, p, 1);
        mpz_fdiv_q_ui(e, e, 4);
        mpz_powm(x, t, e, p);           /* x^2 ≡ -1 (mod p) */

        mpz_set(a, p);
        mpz_set(b, x);
        mpz_mul(sq, b, b);
        while (mpz_cmp(sq, p) >= 0) {    /* reduce until b^2 < p */
            mpz_mod(r, a, b);
            mpz_set(a, b);
            mpz_set(b, r);
            mpz_mul(sq, b, b);
        }
        mpz_set(c, b);
        mpz_mul(sq, b, b);
        mpz_sub(s, p, sq);              /* s = p - c^2 */
        mpz_sqrt(d, s);
        mpz_mul(sq, d, d);
        ok = (mpz_cmp(sq, s) == 0);     /* confirm s is a perfect square */
    }

    mpz_clears(t, x, e, a, b, r, s, sq, NULL);
    return ok;
}

/* If the Gaussian integer g = u + v i divides z = a + b i exactly, overwrite
 * (a, b) with the quotient and return true; otherwise leave (a, b) unchanged
 * and return false.  z/g = z conj(g) / N(g). */
static bool df_gauss_divide(mpz_t a, mpz_t b, const mpz_t u, const mpz_t v) {
    mpz_t nrm, re, im, t;
    mpz_inits(nrm, re, im, t, NULL);
    mpz_mul(nrm, u, u); mpz_mul(t, v, v); mpz_add(nrm, nrm, t);   /* N(g) */
    mpz_mul(re, a, u); mpz_mul(t, b, v); mpz_add(re, re, t);      /* a u + b v */
    mpz_mul(im, b, u); mpz_mul(t, a, v); mpz_sub(im, im, t);      /* b u - a v */
    bool ok = mpz_divisible_p(re, nrm) && mpz_divisible_p(im, nrm);
    if (ok) {
        mpz_divexact(a, re, nrm);
        mpz_divexact(b, im, nrm);
    }
    mpz_clears(nrm, re, im, t, NULL);
    return ok;
}

/* Rotate (re, im) by multiples of i into the canonical first-quadrant
 * associate: the unique unit-multiple with Re > 0 and Im >= 0. */
void df_normalize_quadrant(mpz_t re, mpz_t im) {
    mpz_t t;
    mpz_init(t);
    for (int k = 0; k < 4; k++) {
        if (mpz_sgn(re) > 0 && mpz_sgn(im) >= 0) break;
        mpz_neg(t, im);    /* multiply by i: (re, im) -> (-im, re) */
        mpz_set(im, re);
        mpz_set(re, t);
    }
    mpz_clear(t);
}

/* Extract a Gaussian-integer expression into two mpz_t parts a, b (both are
 * INITIALIZED here — the caller must NOT pre-init them).  Accepts a rational
 * integer (machine int / BigInt; im 0) or Complex[a, b] with integer-like
 * parts.  Returns false (leaving a, b uninitialised) for anything else. */
bool df_to_gaussian(Expr* e, mpz_t a, mpz_t b) {
    if (expr_is_integer_like(e)) {
        expr_to_mpz(e, a);          /* inits a */
        mpz_init_set_ui(b, 0);
        return true;
    }
    Expr *r, *i;
    if (is_complex(e, &r, &i) && expr_is_integer_like(r) && expr_is_integer_like(i)) {
        expr_to_mpz(r, a);          /* inits a */
        expr_to_mpz(i, b);          /* inits b */
        return true;
    }
    return false;
}

/* Factor the Gaussian integer z = a_in + b_in i into Gaussian prime powers
 * (gu[i] + gv[i] i)^ge[i].  The primes are NOT associate-normalised (a split
 * prime's conjugate appears as c - d i).  On success returns true and sets
 * *out_gu, *out_gv (arrays of init'd mpz_t), *out_ge and *out_gn; the caller
 * clears each prime and frees all three arrays.  A unit (N(z) == 1) factors as
 * the empty product (*out_gn == 0).  Returns false — allocating nothing — for
 * z == 0 or any factoring failure. */
bool df_gaussian_prime_factor(const mpz_t a_in, const mpz_t b_in,
                                     mpz_t** out_gu, mpz_t** out_gv,
                                     unsigned long** out_ge, size_t* out_gn) {
    mpz_t N, t;
    mpz_inits(N, t, NULL);
    mpz_mul(N, a_in, a_in);
    mpz_mul(t, b_in, b_in);
    mpz_add(N, N, t);
    mpz_clear(t);

    if (mpz_sgn(N) == 0) { mpz_clear(N); return false; }   /* z = 0 */

    mpz_t* primes;
    unsigned long* exps;
    size_t np;
    bool fok = df_factor_mpz(N, &primes, &exps, &np);
    mpz_clear(N);
    if (!fok) return false;

    /* Gaussian prime factors of z: (gu[i] + gv[i] i)^ge[i]. */
    size_t cap = 2 * np + 1;
    mpz_t* gu = (mpz_t*)malloc(cap * sizeof(mpz_t));
    mpz_t* gv = (mpz_t*)malloc(cap * sizeof(mpz_t));
    unsigned long* ge = (unsigned long*)malloc(cap * sizeof(unsigned long));
    size_t gn = 0;
    bool ok = true;

    for (size_t i = 0; i < np && ok; i++) {
        if (mpz_cmp_ui(primes[i], 2) == 0) {
            /* 2 ramifies: (1 + i) with exponent v_2(N). */
            mpz_init_set_ui(gu[gn], 1);
            mpz_init_set_ui(gv[gn], 1);
            ge[gn] = exps[i];
            gn++;
        } else if (mpz_fdiv_ui(primes[i], 4) == 3) {
            /* p ≡ 3 (mod 4) is inert; v_p(N) is even, exponent in z is half. */
            mpz_init_set(gu[gn], primes[i]);
            mpz_init_set_ui(gv[gn], 0);
            ge[gn] = exps[i] / 2;
            gn++;
        } else {
            /* p ≡ 1 (mod 4) splits as (c + d i)(c - d i). */
            mpz_t c, d;
            mpz_inits(c, d, NULL);
            if (!df_two_squares(primes[i], c, d)) {
                mpz_clears(c, d, NULL);
                ok = false;
                break;
            }
            /* Exponent of c + d i in z by repeated exact division. */
            mpz_t wa, wb;
            mpz_init_set(wa, a_in);
            mpz_init_set(wb, b_in);
            unsigned long epi = 0;
            while (df_gauss_divide(wa, wb, c, d)) epi++;
            mpz_clears(wa, wb, NULL);
            unsigned long econj = exps[i] - epi;
            if (epi > 0) {
                mpz_init_set(gu[gn], c);
                mpz_init_set(gv[gn], d);
                ge[gn] = epi;
                gn++;
            }
            if (econj > 0) {
                mpz_init_set(gu[gn], c);
                mpz_init(gv[gn]);
                mpz_neg(gv[gn], d);
                ge[gn] = econj;
                gn++;
            }
            mpz_clears(c, d, NULL);
        }
    }

    for (size_t i = 0; i < np; i++) mpz_clear(primes[i]);
    free(primes);
    free(exps);

    if (!ok) {
        for (size_t i = 0; i < gn; i++) mpz_clears(gu[i], gv[i], NULL);
        free(gu); free(gv); free(ge);
        return false;
    }
    *out_gu = gu;
    *out_gv = gv;
    *out_ge = ge;
    *out_gn = gn;
    return true;
}

/* Build the list of Gaussian divisors of z = a_in + b_in i, one representative
 * per associate class (first quadrant), sorted by (Re, Im).  Factors z into
 * Gaussian primes, enumerates the divisor lattice, and normalises each product
 * to its first-quadrant associate.  Returns a List Expr (just {1} for a unit),
 * or NULL on z == 0 / factoring failure. */
Expr* divisors_gaussian(const mpz_t a_in, const mpz_t b_in) {
    mpz_t* gu;
    mpz_t* gv;
    unsigned long* ge;
    size_t gn;
    if (!df_gaussian_prime_factor(a_in, b_in, &gu, &gv, &ge, &gn)) return NULL;

    size_t ndiv = 1;
    bool ovf = false;
    for (size_t i = 0; i < gn; i++) {
        size_t f = ge[i] + 1;
        if (ndiv > SIZE_MAX / f) { ovf = true; break; }
        ndiv *= f;
    }
    df_gpair* divs = ovf ? NULL : (df_gpair*)malloc(ndiv * sizeof(df_gpair));
    if (!divs) {
        for (size_t i = 0; i < gn; i++) mpz_clears(gu[i], gv[i], NULL);
        free(gu); free(gv); free(ge);
        return NULL;
    }
    unsigned long* cur = (unsigned long*)calloc(gn, sizeof(unsigned long));
    mpz_t re, im, pu, pv, nu, nv, tt;
    mpz_inits(re, im, pu, pv, nu, nv, tt, NULL);
    for (size_t d = 0; d < ndiv; d++) {
        mpz_set_ui(re, 1);
        mpz_set_ui(im, 0);
        for (size_t i = 0; i < gn; i++) {
            /* (pu + pv i) = (gu[i] + gv[i] i)^cur[i] */
            mpz_set_ui(pu, 1);
            mpz_set_ui(pv, 0);
            for (unsigned long k = 0; k < cur[i]; k++) {
                mpz_mul(nu, pu, gu[i]); mpz_mul(tt, pv, gv[i]); mpz_sub(nu, nu, tt);
                mpz_mul(nv, pu, gv[i]); mpz_mul(tt, pv, gu[i]); mpz_add(nv, nv, tt);
                mpz_set(pu, nu); mpz_set(pv, nv);
            }
            /* (re + im i) *= (pu + pv i) */
            mpz_mul(nu, re, pu); mpz_mul(tt, im, pv); mpz_sub(nu, nu, tt);
            mpz_mul(nv, re, pv); mpz_mul(tt, im, pu); mpz_add(nv, nv, tt);
            mpz_set(re, nu); mpz_set(im, nv);
        }
        df_normalize_quadrant(re, im);
        mpz_init_set(divs[d].re, re);
        mpz_init_set(divs[d].im, im);
        for (size_t i = 0; i < gn; i++) {
            if (++cur[i] <= ge[i]) break;
            cur[i] = 0;
        }
    }
    mpz_clears(re, im, pu, pv, nu, nv, tt, NULL);
    free(cur);
    for (size_t i = 0; i < gn; i++) mpz_clears(gu[i], gv[i], NULL);
    free(gu); free(gv); free(ge);

    qsort(divs, ndiv, sizeof(df_gpair), df_gpair_cmp);

    Expr** items = (Expr**)malloc(ndiv * sizeof(Expr*));
    for (size_t d = 0; d < ndiv; d++) {
        if (mpz_sgn(divs[d].im) == 0) {
            items[d] = expr_bigint_normalize(expr_new_bigint_from_mpz(divs[d].re));
        } else {
            Expr* cargs[2] = {
                expr_bigint_normalize(expr_new_bigint_from_mpz(divs[d].re)),
                expr_bigint_normalize(expr_new_bigint_from_mpz(divs[d].im))
            };
            items[d] = expr_new_function(expr_new_symbol(SYM_Complex), cargs, 2);
        }
    }
    Expr* out = expr_new_function(expr_new_symbol(SYM_List), items, ndiv);
    free(items);
    for (size_t d = 0; d < ndiv; d++) mpz_clears(divs[d].re, divs[d].im, NULL);
    free(divs);
    return out;
}

/* ---- Gaussian GCD / extended GCD (GCD, LCM, ExtendedGCD, CoprimeQ) -------- *
 * The Euclidean algorithm in Z[i]: the quotient of x by y is x conj(y) / N(y)
 * rounded component-wise to the NEAREST integer, ties to even (Round's
 * convention) -- and it starts from the operand of LARGER norm (swapping the
 * cofactors back). Together those give the reference ExtendedGCD
 * coefficients on 399 of 400 random pairs (the remaining one is a norm tie
 * between remainders, where both answers are valid Bezout certificates);
 * round-half-up matched 319, and no norm ordering 390. The gcd is normalised
 * to its first-quadrant associate (Re > 0, Im >= 0), which matches
 * the reference on all 400. */

/* q = num / den rounded to nearest, ties to even (den > 0). */
static void gi_round_even(mpz_t q, const mpz_t num, const mpz_t den) {
    mpz_t r, r2;
    mpz_inits(r, r2, NULL);
    mpz_fdiv_qr(q, r, num, den);                 /* num = q den + r, 0 <= r < den */
    mpz_mul_2exp(r2, r, 1);
    int c = mpz_cmp(r2, den);
    if (c > 0 || (c == 0 && mpz_odd_p(q))) mpz_add_ui(q, q, 1);
    mpz_clears(r, r2, NULL);
}

/* (dr, di) = (xr + xi i)(yr + yi i). Outputs may alias nothing. */
static void gi_mul(mpz_t dr, mpz_t di, const mpz_t xr, const mpz_t xi,
                   const mpz_t yr, const mpz_t yi) {
    mpz_t t;
    mpz_init(t);
    mpz_mul(dr, xr, yr); mpz_mul(t, xi, yi); mpz_sub(dr, dr, t);
    mpz_mul(di, xr, yi); mpz_mul(t, xi, yr); mpz_add(di, di, t);
    mpz_clear(t);
}

/* Multiply (re, im) in place by i^k. */
static void gi_rotate(mpz_t re, mpz_t im, int k) {
    mpz_t t;
    mpz_init(t);
    for (int j = 0; j < (k & 3); j++) {          /* (re, im) * i = (-im, re) */
        mpz_neg(t, im);
        mpz_set(im, re);
        mpz_set(re, t);
    }
    mpz_clear(t);
}

/* The power k of i that takes (re, im) != 0 to its first-quadrant associate. */
static int gi_quadrant_unit(const mpz_t re, const mpz_t im) {
    int sr = mpz_sgn(re), si = mpz_sgn(im);
    if (sr > 0 && si >= 0) return 0;
    if (sr <= 0 && si > 0) return 3;             /* * (-i) */
    if (sr < 0 && si <= 0) return 2;             /* * (-1) */
    return 1;                                    /* * i    */
}

static void gi_xgcd_ordered(mpz_t gr, mpz_t gi, mpz_t sr, mpz_t si, mpz_t tr, mpz_t ti,
                            const mpz_t ar, const mpz_t ai, const mpz_t br, const mpz_t bi);

void gi_xgcd(mpz_t gr, mpz_t gi, mpz_t sr, mpz_t si, mpz_t tr, mpz_t ti,
             const mpz_t ar, const mpz_t ai, const mpz_t br, const mpz_t bi) {
    mpz_t na, nb;
    mpz_inits(na, nb, NULL);
    mpz_mul(na, ar, ar); mpz_addmul(na, ai, ai);
    mpz_mul(nb, br, br); mpz_addmul(nb, bi, bi);
    if (mpz_cmp(na, nb) < 0)
        gi_xgcd_ordered(gr, gi, tr, ti, sr, si, br, bi, ar, ai);
    else
        gi_xgcd_ordered(gr, gi, sr, si, tr, ti, ar, ai, br, bi);
    mpz_clears(na, nb, NULL);
}

static void gi_xgcd_ordered(mpz_t gr, mpz_t gi, mpz_t sr, mpz_t si, mpz_t tr, mpz_t ti,
                            const mpz_t ar, const mpz_t ai, const mpz_t br, const mpz_t bi) {
    mpz_t r0r, r0i, r1r, r1i, s0r, s0i, s1r, s1i, t0r, t0i, t1r, t1i;
    mpz_t nb, num, t, qr, qi, pr, pi;
    mpz_inits(r0r, r0i, r1r, r1i, s0r, s0i, s1r, s1i, t0r, t0i, t1r, t1i,
              nb, num, t, qr, qi, pr, pi, NULL);
    mpz_set(r0r, ar); mpz_set(r0i, ai); mpz_set(r1r, br); mpz_set(r1i, bi);
    mpz_set_ui(s0r, 1);                          /* s = (1, 0), t = (0, 1) */
    mpz_set_ui(t1r, 1);
    while (mpz_sgn(r1r) != 0 || mpz_sgn(r1i) != 0) {
        mpz_mul(nb, r1r, r1r); mpz_mul(t, r1i, r1i); mpz_add(nb, nb, t);
        mpz_mul(num, r0r, r1r); mpz_mul(t, r0i, r1i); mpz_add(num, num, t);
        gi_round_even(qr, num, nb);              /* Re(r0 conj(r1)) / N(r1) */
        mpz_mul(num, r0i, r1r); mpz_mul(t, r0r, r1i); mpz_sub(num, num, t);
        gi_round_even(qi, num, nb);              /* Im(r0 conj(r1)) / N(r1) */
        /* (r0, r1) <- (r1, r0 - q r1), and likewise for s and t. */
        gi_mul(pr, pi, qr, qi, r1r, r1i);
        mpz_sub(pr, r0r, pr); mpz_sub(pi, r0i, pi);
        mpz_swap(r0r, r1r); mpz_swap(r0i, r1i); mpz_swap(r1r, pr); mpz_swap(r1i, pi);
        gi_mul(pr, pi, qr, qi, s1r, s1i);
        mpz_sub(pr, s0r, pr); mpz_sub(pi, s0i, pi);
        mpz_swap(s0r, s1r); mpz_swap(s0i, s1i); mpz_swap(s1r, pr); mpz_swap(s1i, pi);
        gi_mul(pr, pi, qr, qi, t1r, t1i);
        mpz_sub(pr, t0r, pr); mpz_sub(pi, t0i, pi);
        mpz_swap(t0r, t1r); mpz_swap(t0i, t1i); mpz_swap(t1r, pr); mpz_swap(t1i, pi);
    }
    if (mpz_sgn(r0r) != 0 || mpz_sgn(r0i) != 0) {
        int k = gi_quadrant_unit(r0r, r0i);
        gi_rotate(r0r, r0i, k);
        gi_rotate(s0r, s0i, k);
        gi_rotate(t0r, t0i, k);
    }
    mpz_set(gr, r0r); mpz_set(gi, r0i);
    mpz_set(sr, s0r); mpz_set(si, s0i);
    mpz_set(tr, t0r); mpz_set(ti, t0i);
    mpz_clears(r0r, r0i, r1r, r1i, s0r, s0i, s1r, s1i, t0r, t0i, t1r, t1i,
               nb, num, t, qr, qi, pr, pi, NULL);
}

void gi_gcd(mpz_t gr, mpz_t gi, const mpz_t ar, const mpz_t ai,
            const mpz_t br, const mpz_t bi) {
    mpz_t sr, si, tr, ti;
    mpz_inits(sr, si, tr, ti, NULL);
    gi_xgcd(gr, gi, sr, si, tr, ti, ar, ai, br, bi);
    mpz_clears(sr, si, tr, ti, NULL);
}

void gi_normalize(mpz_t re, mpz_t im, int* unit_power) {
    int k = 0;
    if (mpz_sgn(re) != 0 || mpz_sgn(im) != 0) {
        k = gi_quadrant_unit(re, im);
        gi_rotate(re, im, k);
    }
    if (unit_power) *unit_power = k;
}

void gi_rotate_by(mpz_t re, mpz_t im, int unit_power) {
    gi_rotate(re, im, unit_power);
}

void gi_lcm(mpz_t lr, mpz_t li, const mpz_t ar, const mpz_t ai,
            const mpz_t br, const mpz_t bi) {
    if ((mpz_sgn(ar) == 0 && mpz_sgn(ai) == 0) || (mpz_sgn(br) == 0 && mpz_sgn(bi) == 0)) {
        mpz_set_ui(lr, 0); mpz_set_ui(li, 0);
        return;
    }
    mpz_t gr, gi_, pr, pi;
    mpz_inits(gr, gi_, pr, pi, NULL);
    gi_gcd(gr, gi_, ar, ai, br, bi);
    gi_mul(pr, pi, ar, ai, br, bi);
    /* exact division of the product by the gcd: p conj(g) / N(g) */
    mpz_t nrm, t, qr, qi;
    mpz_inits(nrm, t, qr, qi, NULL);
    mpz_mul(nrm, gr, gr); mpz_mul(t, gi_, gi_); mpz_add(nrm, nrm, t);
    mpz_mul(qr, pr, gr); mpz_mul(t, pi, gi_); mpz_add(qr, qr, t);
    mpz_mul(qi, pi, gr); mpz_mul(t, pr, gi_); mpz_sub(qi, qi, t);
    mpz_divexact(lr, qr, nrm);
    mpz_divexact(li, qi, nrm);
    gi_normalize(lr, li, NULL);
    mpz_clears(gr, gi_, pr, pi, nrm, t, qr, qi, NULL);
}

/* A Gaussian rational (re + im i) / den, den > 0, from an Integer, Rational,
 * or Complex with rational parts. re, im, den must already be init'd. */
bool gq_from_expr(Expr* e, mpz_t re, mpz_t im, mpz_t den) {
    if (is_rational_like(e)) {
        rational_like_to_mpz_pair(e, re, den);
        mpz_set_ui(im, 0);
        return true;
    }
    Expr *r, *i;
    if (!is_complex(e, &r, &i) || !is_rational_like(r) || !is_rational_like(i))
        return false;
    mpz_t rn, rd, in_, id;
    mpz_inits(rn, rd, in_, id, NULL);
    rational_like_to_mpz_pair(r, rn, rd);
    rational_like_to_mpz_pair(i, in_, id);
    mpz_lcm(den, rd, id);
    mpz_divexact(rd, den, rd);                   /* scale factors */
    mpz_divexact(id, den, id);
    mpz_mul(re, rn, rd);
    mpz_mul(im, in_, id);
    mpz_clears(rn, rd, in_, id, NULL);
    return true;
}

/* (re + im i) / den as a canonical number: Integer, Rational or Complex. */
Expr* gq_to_expr(const mpz_t re, const mpz_t im, const mpz_t den) {
    Expr* r = make_rational_mpz(re, den);
    if (mpz_sgn(im) == 0) return r;
    Expr* i = make_rational_mpz(im, den);
    Expr* args[2] = { r, i };
    return expr_new_function(expr_new_symbol(SYM_Complex), args, 2);
}

/* True when e is an exact non-real Complex number with rational parts. */
bool is_gaussian_rational(Expr* e) {
    Expr *r, *i;
    return is_complex(e, &r, &i) && is_rational_like(r) && is_rational_like(i);
}

/* x mod m in Z[i], the nearest-quotient convention: x - m Round[x / m], the quotient
 * rounded component-wise to nearest, ties to even (Mod[3 + 4 I, 2] is -1).
 * m must be nonzero. Writes (rr, ri); outputs init'd by the caller and may not
 * alias the inputs. */
void gi_mod(mpz_t rr, mpz_t ri, const mpz_t xr, const mpz_t xi,
            const mpz_t mr, const mpz_t mi) {
    mpz_t nm, num, t, qr, qi;
    mpz_inits(nm, num, t, qr, qi, NULL);
    mpz_mul(nm, mr, mr); mpz_addmul(nm, mi, mi);
    mpz_mul(num, xr, mr); mpz_addmul(num, xi, mi);         /* Re(x conj(m)) */
    gi_round_even(qr, num, nm);
    mpz_mul(num, xi, mr); mpz_submul(num, xr, mi);         /* Im(x conj(m)) */
    gi_round_even(qi, num, nm);
    gi_mul(rr, ri, qr, qi, mr, mi);
    mpz_sub(rr, xr, rr);
    mpz_sub(ri, xi, ri);
    mpz_clears(nm, num, t, qr, qi, NULL);
}

/* a^e mod m in Z[i] for an integer e (m != 0), by square-and-multiply with a
 * gi_mod reduction after every product -- which is what reproduces
 * the reference PowerMod on 390 of 400 random Gaussian cases (all 250 with
 * e >= 0). A negative e inverts a first through gi_xgcd. Returns false when a
 * is not invertible modulo m. Outputs init'd by the caller. */
bool gi_powmod(mpz_t rr, mpz_t ri, const mpz_t ar, const mpz_t ai,
               const mpz_t e_in, const mpz_t mr, const mpz_t mi) {
    mpz_t br, bi, e, xr, xi, tr, ti;
    mpz_inits(br, bi, e, xr, xi, tr, ti, NULL);
    bool ok = true;
    mpz_set(e, e_in);
    if (mpz_sgn(e) < 0) {
        mpz_t gr, gim, sr, si, ur, ui;
        mpz_inits(gr, gim, sr, si, ur, ui, NULL);
        gi_xgcd(gr, gim, sr, si, ur, ui, ar, ai, mr, mi);
        ok = mpz_cmp_ui(gr, 1) == 0 && mpz_sgn(gim) == 0;   /* gcd a unit */
        if (ok) gi_mod(br, bi, sr, si, mr, mi);             /* s a == 1 (mod m) */
        mpz_clears(gr, gim, sr, si, ur, ui, NULL);
        mpz_neg(e, e);
    } else {
        gi_mod(br, bi, ar, ai, mr, mi);
    }
    if (ok) {
        mpz_set_ui(xr, 1); mpz_set_ui(xi, 0);
        size_t bits = mpz_sizeinbase(e, 2);
        for (size_t k = 0; k < bits && mpz_sgn(e) != 0; k++) {
            if (mpz_tstbit(e, k)) {
                gi_mul(tr, ti, xr, xi, br, bi);
                gi_mod(xr, xi, tr, ti, mr, mi);
            }
            if (k + 1 < bits) {
                gi_mul(tr, ti, br, bi, br, bi);
                gi_mod(br, bi, tr, ti, mr, mi);
            }
        }
        gi_mod(rr, ri, xr, xi, mr, mi);
    }
    mpz_clears(br, bi, e, xr, xi, tr, ti, NULL);
    return ok;
}

/* ---- machine Gaussian arithmetic (Mod, PowerMod fast paths) -------------- *
 * The same nearest-quotient reduction as gi_mod, in int64, for a modulus whose
 * components are below 2^16 and operands below 2^40: every intermediate --
 * x conj(m) for x a reduced residue product (components < 2^36) or an operand
 * (< 2^41) -- stays under 2^58. Bit-identical to the mpz path it shortcuts. */
#define GI64_MOD_MAX ((int64_t)1 << 16)
#define GI64_ARG_MAX ((int64_t)1 << 40)

static int64_t gi64_round_even(int64_t num, int64_t den) {   /* den > 0 */
    int64_t q = num / den, r = num % den;
    if (r < 0) { q -= 1; r += den; }                          /* floored */
    if (2 * r > den || (2 * r == den && (q & 1))) q += 1;
    return q;
}

static void gi64_mod(int64_t xr, int64_t xi, int64_t mr, int64_t mi,
                     int64_t* rr, int64_t* ri) {
    int64_t nm = mr * mr + mi * mi;
    int64_t qr = gi64_round_even(xr * mr + xi * mi, nm);
    int64_t qi = gi64_round_even(xi * mr - xr * mi, nm);
    *rr = xr - (qr * mr - qi * mi);
    *ri = xi - (qr * mi + qi * mr);
}

static bool gi64_fits(const Expr* e, int64_t lim, int64_t* re, int64_t* im) {
    Expr *r, *i;
    if (e->type == EXPR_INTEGER) {
        *re = e->data.integer; *im = 0;
    } else if (is_complex((Expr*)e, &r, &i) &&
               r->type == EXPR_INTEGER && i->type == EXPR_INTEGER) {
        *re = r->data.integer; *im = i->data.integer;
    } else {
        return false;
    }
    return *re > -lim && *re < lim && *im > -lim && *im < lim;
}

static Expr* gi64_expr(int64_t re, int64_t im) {
    if (im == 0) return expr_new_integer(re);
    Expr* args[2] = { expr_new_integer(re), expr_new_integer(im) };
    return expr_new_function(expr_new_symbol(SYM_Complex), args, 2);
}

Expr* nt_gaussian_mod(Expr* m, Expr* n, Expr* d) {
    int64_t mr, mi, nr, ni, dr = 0, di = 0;
    if (gi64_fits(n, GI64_MOD_MAX, &nr, &ni) && gi64_fits(m, GI64_ARG_MAX, &mr, &mi) &&
        (!d || gi64_fits(d, GI64_ARG_MAX, &dr, &di)) && (nr != 0 || ni != 0)) {
        int64_t rr, ri;
        gi64_mod(mr - dr, mi - di, nr, ni, &rr, &ri);
        return gi64_expr(rr + dr, ri + di);
    }
    mpz_t a, b, c, e, f, g, x, y;
    if (!df_to_gaussian(m, a, b)) return NULL;
    if (!df_to_gaussian(n, c, e)) { mpz_clears(a, b, NULL); return NULL; }
    if (d) {
        if (!df_to_gaussian(d, f, g)) { mpz_clears(a, b, c, e, NULL); return NULL; }
    } else {
        mpz_inits(f, g, NULL);
    }
    Expr* out = NULL;
    if (mpz_sgn(c) != 0 || mpz_sgn(e) != 0) {
        mpz_inits(x, y, NULL);
        mpz_sub(a, a, f); mpz_sub(b, b, g);
        gi_mod(x, y, a, b, c, e);
        mpz_add(x, x, f); mpz_add(y, y, g);
        mpz_t one;
        mpz_init_set_ui(one, 1);
        out = gq_to_expr(x, y, one);
        mpz_clears(x, y, one, NULL);
    }
    mpz_clears(a, b, c, e, f, g, NULL);
    return out;
}

/* a^e mod m for small Gaussian a, m and a machine e >= 0; false to decline. */
bool gi64_powmod(Expr* a, int64_t e, Expr* m, Expr** out) {
    int64_t ar, ai, mr, mi;
    if (e < 0 || !gi64_fits(m, GI64_MOD_MAX, &mr, &mi) || (mr == 0 && mi == 0) ||
        !gi64_fits(a, GI64_ARG_MAX, &ar, &ai)) return false;
    int64_t br, bi, xr = 1, xi = 0, tr, ti;
    gi64_mod(ar, ai, mr, mi, &br, &bi);
    while (e > 0) {
        if (e & 1) {
            tr = xr * br - xi * bi; ti = xr * bi + xi * br;
            gi64_mod(tr, ti, mr, mi, &xr, &xi);
        }
        e >>= 1;
        if (e > 0) {
            tr = br * br - bi * bi; ti = 2 * br * bi;
            gi64_mod(tr, ti, mr, mi, &br, &bi);
        }
    }
    gi64_mod(xr, xi, mr, mi, &xr, &xi);
    *out = gi64_expr(xr, xi);
    return true;
}
