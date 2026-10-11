/* Mathilda -- ChebyshevT, the Chebyshev polynomial of the first kind.
 *
 *   ChebyshevT[n, x]   gives T_n(x).
 *
 * See chebyshev_t.c for the evaluation strategy (exact integer-order polynomial,
 * stable numeric recurrence, and the closed forms for non-integer order). The
 * argument-shape predicates and polynomial/result assembly are shared with
 * ChebyshevU via chebyshev_common.{c,h}.
 */
#ifndef MATHILDA_CHEBYSHEV_T_H
#define MATHILDA_CHEBYSHEV_T_H

#include "expr.h"

/* Builtin entry point. Takes ownership of `res` per the builtin contract. */
Expr* builtin_chebyshev_t(Expr* res);

/* Register ChebyshevT and its attributes. Called from core_init(). */
void chebyshev_t_init(void);

#endif /* MATHILDA_CHEBYSHEV_T_H */
