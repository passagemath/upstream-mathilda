/* Mathilda -- HermiteH, the physicists' Hermite polynomial H_n(x).
 *
 *   HermiteH[n, x]
 *
 * See hermite_h.c for the evaluation strategy. The builtin takes ownership of
 * `res` per the standard builtin contract; hermite_h_init() registers it and
 * its attributes (Listable, NumericFunction, Protected). The docstring lives in
 * info.c (info_init); the derivative rule lives in calculus/deriv.c.
 */
#ifndef MATHILDA_HERMITE_H_H
#define MATHILDA_HERMITE_H_H

#include "expr.h"

/* HermiteH[n, x]: the two-argument builtin entry point. */
Expr* builtin_hermite_h(Expr* res);

/* Register the builtin and its attributes. */
void hermite_h_init(void);

#endif /* MATHILDA_HERMITE_H_H */
