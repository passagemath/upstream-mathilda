/* coprimeq.c -- CoprimeQ[] (incl. Gaussian-integer coprimality).
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
#include "ndarray.h"
#include "ndinteger.h" /* ndint_binary_predicate */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <gmp.h>

/* ---- CoprimeQ ----------------------------------------------------------- */

/* Extract a Gaussian-integer expression into two initialized mpz_t parts.
 * Accepts a rational integer (machine int or BigInt; imaginary part 0) or an
 * exact Gaussian integer Complex[a, b] with integer-like parts.  Returns true
 * on success (re/im are set); false for anything else (rationals, reals,
 * symbols, ...).  re and im must already be mpz_init'd by the caller. */
static bool coprimeq_to_gaussian(Expr* e, mpz_t re, mpz_t im) {
    /* re/im are already init'd by the caller (see the contract above); the
     * integer branch's mpz_set_ui(im, 0) relies on that.  expr_to_mpz re-inits
     * its output, so calling it directly on re/im would leak their existing
     * limbs -- extract through a temporary and mpz_set into the targets. */
    mpz_t tmp;
    if (expr_is_integer_like(e)) {
        expr_to_mpz(e, tmp);          /* inits tmp */
        mpz_set(re, tmp);
        mpz_clear(tmp);
        mpz_set_ui(im, 0);
        return true;
    }
    Expr *r, *i;
    if (is_complex(e, &r, &i) && expr_is_integer_like(r) && expr_is_integer_like(i)) {
        expr_to_mpz(r, tmp);
        mpz_set(re, tmp);
        mpz_clear(tmp);
        expr_to_mpz(i, tmp);
        mpz_set(im, tmp);
        mpz_clear(tmp);
        return true;
    }
    return false;
}

/* True iff two Gaussian integers are coprime in Z[i]: their gcd is a unit,
 * which after first-quadrant normalisation means exactly 1. */
static bool gaussian_pair_coprime(const mpz_t a1, const mpz_t a2,
                                  const mpz_t b1, const mpz_t b2) {
    mpz_t gr, gi;
    mpz_inits(gr, gi, NULL);
    gi_gcd(gr, gi, a1, a2, b1, b2);
    bool coprime = mpz_cmp_ui(gr, 1) == 0 && mpz_sgn(gi) == 0;
    mpz_clears(gr, gi, NULL);
    return coprime;
}

/* CoprimeQ[n1, n2, ...] yields True when the arguments are pairwise relatively
 * prime, and False otherwise. Over the integers a pair is coprime when their
 * GCD is 1; when any argument is an exact Gaussian integer the test is made in
 * Z[i] (the gcd is a unit). As a *Q predicate it always answers True or False:
 * anything that is not an integer or Gaussian integer -- a rational, a real, a
 * symbol, a rule -- makes it False. With one argument it asks whether n is a
 * unit (CoprimeQ[1], CoprimeQ[-1], CoprimeQ[I] are True; CoprimeQ[3] and
 * CoprimeQ[0] are False); CoprimeQ[] is False. All as in Mathematica 14, which
 * has no GaussianIntegers option here: CoprimeQ[2, 5, GaussianIntegers -> True]
 * is False because the rule is not an integer. Listable threading is done by
 * the evaluator first. */
Expr* builtin_coprimeq(Expr* res) {
    if (res->type != EXPR_FUNCTION) return expr_new_symbol(SYM_False);
    size_t nnums = res->data.function.arg_count;
    Expr** nums = res->data.function.args;

    if (nnums == 0) return expr_new_symbol(SYM_False);
    bool gaussian = false;

    /* A packed / visible int64 buffer: one pass, answering with a packed bool
     * List; any shape the kernel declines degrades element by element. */
    for (size_t i = 0; i < nnums; i++) {
        if (is_ndarray(nums[i])) {
            Expr* nd = ndint_binary_predicate(res, NDBP_COPRIME);
            return nd ? nd : ndarray_delist_and_reeval(res);
        }
    }

    /* One argument: is it a unit of Z (or Z[i])? */
    if (nnums == 1) {
        mpz_t a, b;
        if (!df_to_gaussian(nums[0], a, b)) return expr_new_symbol(SYM_False);
        mpz_mul(a, a, a); mpz_addmul(a, b, b);           /* norm */
        bool unit = mpz_cmp_ui(a, 1) == 0;
        mpz_clears(a, b, NULL);
        return expr_new_symbol(unit ? SYM_True : SYM_False);
    }

    /* If any operand is an exact Gaussian integer (non-real), switch to Z[i]. */
    if (!gaussian) {
        for (size_t i = 0; i < nnums; i++) {
            if (is_gaussian_integer(nums[i])) { gaussian = true; break; }
        }
    }

    bool result = true;

    if (gaussian) {
        /* Extract every operand as a Gaussian integer up front; any operand
         * that is not integer-like or a Gaussian integer is not manifestly
         * coprime -> False. */
        mpz_t* re = (mpz_t*)malloc(nnums * sizeof(mpz_t));
        mpz_t* im = (mpz_t*)malloc(nnums * sizeof(mpz_t));
        bool ok = true;
        size_t got = 0;
        for (size_t i = 0; i < nnums; i++) {
            mpz_init(re[i]);
            mpz_init(im[i]);
            got = i + 1;
            if (!coprimeq_to_gaussian(nums[i], re[i], im[i])) { ok = false; break; }
        }
        if (!ok) {
            result = false;
        } else {
            for (size_t i = 0; i < nnums && result; i++) {
                for (size_t j = i + 1; j < nnums && result; j++) {
                    if (!gaussian_pair_coprime(re[i], im[i], re[j], im[j])) result = false;
                }
            }
        }
        for (size_t i = 0; i < got; i++) mpz_clears(re[i], im[i], NULL);
        free(re);
        free(im);
    } else {
        /* Ordinary integers.  Every operand must be integer-like. */
        for (size_t i = 0; i < nnums; i++) {
            if (!expr_is_integer_like(nums[i])) { result = false; break; }
        }
        if (result) {
            mpz_t* v = (mpz_t*)malloc(nnums * sizeof(mpz_t));
            for (size_t i = 0; i < nnums; i++) {
                expr_to_mpz(nums[i], v[i]);   /* inits v[i] */
            }
            mpz_t g;
            mpz_init(g);
            for (size_t i = 0; i < nnums && result; i++) {
                for (size_t j = i + 1; j < nnums && result; j++) {
                    mpz_gcd(g, v[i], v[j]);
                    if (mpz_cmp_ui(g, 1) != 0) result = false;
                }
            }
            mpz_clear(g);
            for (size_t i = 0; i < nnums; i++) mpz_clear(v[i]);
            free(v);
        }
    }

    return expr_new_symbol(result ? SYM_True : SYM_False);
}
