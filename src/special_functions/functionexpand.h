/*
 * functionexpand.h -- FunctionExpand: MeijerG to named functions, PolyLog and
 * Log at roots of unity to closed forms. See functionexpand.c.
 */
#ifndef MATHILDA_FUNCTIONEXPAND_H
#define MATHILDA_FUNCTIONEXPAND_H

#include "expr.h"

Expr* builtin_function_expand(Expr* res);
void  functionexpand_init(void);

/* The final pass of a definite Integrate: rewrite the MeijerG nodes of `r` to
 * named functions (Refining a resulting Piecewise under `assumptions`, which
 * may be NULL), and PolyLog / Log values at roots of unity to Zeta / Catalan /
 * Pi closed forms -- each family accepted only when it leaves a genuinely
 * closed result (see functionexpand.c). Takes ownership of `r`; returns the
 * owned result (r itself when nothing applied). */
Expr* fe_reduce_definite_result(Expr* r, const Expr* assumptions);

#endif /* MATHILDA_FUNCTIONEXPAND_H */
