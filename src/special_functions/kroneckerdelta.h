#ifndef KRONECKERDELTA_H
#define KRONECKERDELTA_H

#include "expr.h"

/* KroneckerDelta[n1, n2, ...] -- 1 when all ni are equal, 0 when any two are
 * provably unequal, unevaluated otherwise. */
Expr* builtin_kroneckerdelta(Expr* res);
void  kroneckerdelta_init(void);

#endif /* KRONECKERDELTA_H */
