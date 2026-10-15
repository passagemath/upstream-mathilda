/* Mathilda -- WhittakerM, the Whittaker M function.
 *
 *   WhittakerM[k, m, z]   gives the Whittaker function M_{k,m}(z).
 *
 * Related to the Kummer confluent hypergeometric function by
 *   M_{k,m}(z) = e^(-z/2) z^(m+1/2) Hypergeometric1F1[m-k+1/2, 2m+1, z].
 *
 * See whittaker.c for the evaluation strategy (special values, the shared
 * machine helper sf_machine_whittaker_m for the real-machine fast path /
 * NDArray buffer / Compile VM, and the confluent-hypergeometric template for
 * MPFR precision and complex arguments).
 */
#ifndef MATHILDA_WHITTAKER_H
#define MATHILDA_WHITTAKER_H

#include "expr.h"

/* Builtin entry point. Takes ownership of `res` per the builtin contract. */
Expr* builtin_whittaker_m(Expr* res);

/* Register WhittakerM and its attributes. Called from core_init(). */
void whittaker_init(void);

#endif /* MATHILDA_WHITTAKER_H */
