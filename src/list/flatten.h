#ifndef FLATTEN_H
#define FLATTEN_H

#include "expr.h"
#include <stdbool.h>

Expr* builtin_flatten(Expr* res);

/* Is the level argument of Flatten[list, n] Infinity (every level)? */
bool flatten_level_is_infinity(const Expr* lv);

/* Flatten::flev for a call whose level argument (position 2) is a negative
 * integer, with the call as written. Shared with the SparseArray kernel
 * (src/sparse/sparse_struct.c) so both surfaces report it the same way. */
void flatten_flev_message(const Expr* call);

#endif /* FLATTEN_H */
