/* divisors.c -- Divisors[].
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
#include "message.h"   /* mth_message: Quiet/Check funnel */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <gmp.h>

/* Emit the wrong-argument-count diagnostic for Divisors and return NULL so the
 * call is left unevaluated. */
static Expr* divisors_emit_argx(size_t npos) {
    mth_message_gated(g_arith_warnings_muted, "Divisors", "argx",
                      "Divisors called with %zu arguments; 1 argument is expected.", npos);
    return NULL;
}

/* Divisors[n] gives the ascending list of positive integers dividing n (sign
 * ignored).  Divisors[n, GaussianIntegers -> True] (or a non-real Gaussian
 * input) returns the divisors in Z[i], one first-quadrant representative per
 * associate class, sorted by (Re, Im).  Machine ints and BigInts are handled
 * uniformly through GMP.  Divisors is Listable, so list arguments are threaded
 * by the evaluator before this builtin runs.  Non-integer and zero n stay
 * unevaluated; Divisors[] draws Divisors::argx, a second non-rule argument
 * Divisors::nonopt, and a bad option ::optx / ::opttf. */
Expr* builtin_divisors(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    Expr** args = res->data.function.args;

    /* One positional argument, then options (GaussianIntegers only), with
     * the standard argx / nonopt / optx / opttf messages. */
    if (argc == 0 || (argc >= 1 && args[0]->type == EXPR_FUNCTION &&
                      args[0]->data.function.head->type == EXPR_SYMBOL &&
                      (args[0]->data.function.head->data.symbol.name == SYM_Rule ||
                       args[0]->data.function.head->data.symbol.name == SYM_RuleDelayed)))
        return divisors_emit_argx(0);
    bool gaussian = false, gaussian_set = false;
    if (!nt_parse_gaussian_options(res, 1, "Divisors", &gaussian, &gaussian_set))
        return NULL;
    Expr* n = args[0];

    /* Explicit option wins; otherwise auto-enable Z[i] for non-real input. */
    bool use_gaussian = gaussian_set ? gaussian : is_gaussian_integer(n);

    if (use_gaussian) {
        mpz_t a, b;
        if (!df_to_gaussian(n, a, b)) return NULL;
        Expr* out = divisors_gaussian(a, b);
        mpz_clears(a, b, NULL);
        return out;
    }

    /* Ordinary integer path. */
    if (!expr_is_integer_like(n)) return NULL;
    mpz_t v;
    expr_to_mpz(n, v);                 /* inits v */
    if (mpz_sgn(v) == 0) { mpz_clear(v); return NULL; }   /* Divisors[0] */
    mpz_abs(v, v);

    Expr* out;
    if (mpz_cmp_ui(v, 1) == 0) {
        Expr* items[1] = { expr_new_integer(1) };
        out = expr_new_function(expr_new_symbol(SYM_List), items, 1);
    } else {
        out = divisors_ordinary(v);
    }
    mpz_clear(v);
    return out;
}
