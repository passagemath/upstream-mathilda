/* multinomial.c -- Multinomial[n1, n2, ...].
 * Split one builtin per file; see numbertheory.h and numbertheory_internal.h
 * for the subsystem layout. */

#include "numbertheory.h"
#include "numbertheory_internal.h"
#include "eval.h"
#include "sym_names.h"
#include "symtab.h"
#include <stdlib.h>
#include <stdbool.h>

/* An explicit number: Integer, BigInt, Real, MPFR, Rational or Complex.
 * Pi and E are numeric quantities but not numbers, and count as symbolic. */
static bool multinomial_is_number(const Expr* e) {
    switch (e->type) {
        case EXPR_INTEGER:
        case EXPR_BIGINT:
        case EXPR_REAL:
#ifdef USE_MPFR
        case EXPR_MPFR:
#endif
            return true;
        case EXPR_FUNCTION: {
            const Expr* h = e->data.function.head;
            return h->type == EXPR_SYMBOL &&
                   (h->data.symbol.name == SYM_Rational || h->data.symbol.name == SYM_Complex);
        }
        default:
            return false;
    }
}

static bool multinomial_is_exact_zero(const Expr* e) {
    return e->type == EXPR_INTEGER && e->data.integer == 0;
}

/* Multinomial[live...] -- the call with its exact zeros dropped. */
static Expr* multinomial_rebuild(Expr** live, size_t kept) {
    Expr** copies = malloc(sizeof(Expr*) * kept);
    for (size_t i = 0; i < kept; i++) copies[i] = expr_copy(live[i]);
    Expr* out = expr_new_function(expr_new_symbol(SYM_Multinomial), copies, kept);
    free(copies);
    return out;
}

static bool multinomial_has_head(const Expr* e, const char* head) {
    return e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL &&
           e->data.function.head->data.symbol.name == head;
}

/* Multinomial[n1, ..., nk] = (n1 + ... + nk)! / (n1! ... nk!), written as the
 * product of binomials over the partial sums,
 *
 *     Binomial[n1 + n2, n2] Binomial[n1 + n2 + n3, n3] ... ,
 *
 * so every case inherits Binomial's closed forms and its pole conventions
 * (Multinomial[-2, 1] = Binomial[-1, 1] = -1).  The arguments arrive sorted
 * (Orderless), numbers first, and the reference's rules are:
 *   - exact zeros drop out: Multinomial[a, b, 0] is Multinomial[a, b];
 *   - no argument, or one, gives 1;
 *   - two arguments always give Binomial[n1 + n2, n2];
 *   - three or more give the product when at most one argument is not an
 *     explicit number and at most one of the binomials stays unevaluated:
 *     Multinomial[3/2, 2, a] is 35/8 Binomial[7/2 + a, a], while
 *     Multinomial[2, a, b] and Multinomial[1/5, 1/3, 1/2] stay as they are. */
Expr* builtin_multinomial(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    Expr** args = res->data.function.args;

    size_t kept = 0, symbolic = 0;
    for (size_t i = 0; i < argc; i++) {
        if (multinomial_is_exact_zero(args[i])) continue;
        kept++;
        if (!multinomial_is_number(args[i])) symbolic++;
    }
    if (kept <= 1) return expr_new_integer(1);

    Expr** live = malloc(sizeof(Expr*) * kept);
    size_t k = 0;
    for (size_t i = 0; i < argc; i++)
        if (!multinomial_is_exact_zero(args[i])) live[k++] = args[i];

    if (kept == 2) {
        Expr* sum = expr_new_function(expr_new_symbol(SYM_Plus),
            (Expr*[]){ expr_copy(live[0]), expr_copy(live[1]) }, 2);
        Expr* out = eval_and_free(expr_new_function(expr_new_symbol(SYM_Binomial),
            (Expr*[]){ sum, expr_copy(live[1]) }, 2));
        free(live);
        return out;
    }

    if (symbolic > 1) {
        Expr* out = kept < argc ? multinomial_rebuild(live, kept) : NULL;   /* only the zeros drop */
        free(live);
        return out;
    }

    /* Binomial factors over the partial sums; keep the product only if at
     * most one of them stays unevaluated. */
    Expr** factors = malloc(sizeof(Expr*) * (kept - 1));
    Expr* partial = expr_copy(live[0]);
    size_t unevaluated = 0;
    for (size_t i = 1; i < kept; i++) {
        partial = eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus),
            (Expr*[]){ partial, expr_copy(live[i]) }, 2));
        factors[i - 1] = eval_and_free(expr_new_function(expr_new_symbol(SYM_Binomial),
            (Expr*[]){ expr_copy(partial), expr_copy(live[i]) }, 2));
        if (multinomial_has_head(factors[i - 1], SYM_Binomial)) unevaluated++;
    }
    expr_free(partial);

    Expr* out = NULL;
    if (unevaluated <= 1) {
        out = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times), factors, kept - 1));
    } else {
        for (size_t i = 0; i + 1 < kept; i++) expr_free(factors[i]);
        if (kept < argc) out = multinomial_rebuild(live, kept);
    }
    free(factors);
    free(live);
    return out;
}
