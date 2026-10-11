/* Mathilda -- ChebyshevU, the Chebyshev polynomial of the second kind.
 *
 *   ChebyshevU[n, x]   gives U_n(x).
 *
 * See chebyshev_u.c for the evaluation strategy (exact integer-order polynomial,
 * stable numeric recurrence, and the hypergeometric form for non-integer order).
 */
#ifndef MATHILDA_CHEBYSHEV_U_H
#define MATHILDA_CHEBYSHEV_U_H

#include "expr.h"

/* Builtin entry point. Takes ownership of `res` per the builtin contract. */
Expr* builtin_chebyshev_u(Expr* res);

/* Register ChebyshevU and its attributes. Called from core_init(). */
void chebyshev_u_init(void);

#endif /* MATHILDA_CHEBYSHEV_U_H */
