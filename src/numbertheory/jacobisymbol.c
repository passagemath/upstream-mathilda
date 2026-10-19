/* jacobisymbol.c -- JacobiSymbol[].
 * Split from numbertheory.c; see numbertheory.h and
 * numbertheory_internal.h for the subsystem layout. */

#include "numbertheory.h"
#include "numbertheory_internal.h"
#include "arithmetic.h"
#include "expr.h"
#include "message.h"   /* mth_message: Quiet/Check funnel */
#include "ndinteger.h" /* ndint_kronecker_i64 */
#include <stdio.h>
#include <gmp.h>

/* ---- JacobiSymbol ------------------------------------------------------- */

/* Emit the wrong-argument-count diagnostic for JacobiSymbol and return NULL so
 * the call is left unevaluated. */
static Expr* jacobisymbol_emit_argrx(size_t npos) {
    /* The expected tag is argr for one argument, argrx otherwise. */
    if (npos == 1)
        mth_message_gated(g_arith_warnings_muted, "JacobiSymbol", "argr",
            "JacobiSymbol called with 1 argument; 2 arguments are expected.");
    else
        mth_message_gated(g_arith_warnings_muted, "JacobiSymbol", "argrx",
            "JacobiSymbol called with %zu arguments; 2 arguments are expected.", npos);
    return NULL;
}

/* JacobiSymbol[n, m] gives the Jacobi symbol (n/m).
 *
 * Following the Wolfram Language, this is the full Kronecker-symbol
 * generalisation of the Jacobi symbol: m need not be odd or positive, and n
 * may be negative.  For prime m it reduces to the Legendre symbol, which is
 * +-1 according to whether n is a quadratic residue modulo m (and 0 when m
 * divides n).  GMP's mpz_kronecker computes exactly this and is O((log m)^2),
 * so machine integers and arbitrary-precision bigints are handled uniformly.
 *
 * Symbolic arguments leave the call unevaluated; Listable threading over
 * lists is applied by the evaluator before this builtin is reached. */
Expr* builtin_jacobisymbol(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc != 2) return jacobisymbol_emit_argrx(argc);

    Expr* n_expr = res->data.function.args[0];
    Expr* m_expr = res->data.function.args[1];
    if (!expr_is_integer_like(n_expr) || !expr_is_integer_like(m_expr))
        return NULL;

    /* Machine pair: the int64 Kronecker symbol the NDArray kernel and
     * Compile's OP_JACOBI_I use. */
    if (n_expr->type == EXPR_INTEGER && m_expr->type == EXPR_INTEGER)
        return expr_new_integer(ndint_kronecker_i64(n_expr->data.integer,
                                                    m_expr->data.integer));

    mpz_t n, m;
    expr_to_mpz(n_expr, n);
    expr_to_mpz(m_expr, m);
    int k = mpz_kronecker(n, m);   /* full Kronecker: handles even/neg/zero m */
    mpz_clears(n, m, NULL);
    return expr_new_integer(k);    /* always -1, 0, or 1 */
}
