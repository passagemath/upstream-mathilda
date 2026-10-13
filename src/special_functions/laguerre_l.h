/* Mathilda -- LaguerreL, the (generalized) Laguerre polynomial.
 *
 *   LaguerreL[n, x]      gives the Laguerre polynomial L_n(x).
 *   LaguerreL[n, a, x]   gives the generalized Laguerre polynomial L_n^(a)(x),
 *                        with L_n[n, x] == L_n[n, 0, x].
 *
 * See laguerre_l.c for the evaluation strategy (exact integer-order polynomial,
 * stable numeric recurrence, and the confluent-hypergeometric form for a
 * non-integer / complex order).
 */
#ifndef MATHILDA_LAGUERRE_L_H
#define MATHILDA_LAGUERRE_L_H

#include "expr.h"

/* Builtin entry point. Takes ownership of `res` per the builtin contract. */
Expr* builtin_laguerre_l(Expr* res);

/* Register LaguerreL and its attributes. Called from core_init(). */
void laguerre_l_init(void);

#endif /* MATHILDA_LAGUERRE_L_H */
