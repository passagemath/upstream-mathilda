/* Mathilda -- HypergeometricU, the Tricomi confluent hypergeometric function.
 *
 *   HypergeometricU[a, b, z]   gives U(a, b, z), the second solution of the
 *                              confluent equation z w'' + (b - z) w' - a w = 0.
 *
 * Numeric evaluation layers (mirroring the Hypergeometric1F1 / WhittakerM design):
 *   a real machine argument uses the shared double helper sf_machine_hyperu
 *   (which also backs the NDArray/packed buffer loop and the Compile[] VM); an
 *   MPFR/complex/z<=0 argument uses FLINT's rigorous acb_hypgeom_u through
 *   flint_num_hypergeometric_u. Exact reductions: U(0,b,z) = 1 and
 *   U(a,a+1,z) = z^(-a). Everything else stays symbolic.
 */
#ifndef MATHILDA_HYPERGEOMETRICU_H
#define MATHILDA_HYPERGEOMETRICU_H

#include "expr.h"

/* Builtin entry point. Takes ownership of `res` per the builtin contract. */
Expr* builtin_hypergeometric_u(Expr* res);

/* Register HypergeometricU and its attributes. Called from core_init(). */
void hypergeometricu_init(void);

#endif /* MATHILDA_HYPERGEOMETRICU_H */
