/*
 * functionexpand.h -- FunctionExpand (minimal: MeijerG -> Slater expansion).
 * See functionexpand.c.
 */
#ifndef MATHILDA_FUNCTIONEXPAND_H
#define MATHILDA_FUNCTIONEXPAND_H

#include "expr.h"

Expr* builtin_function_expand(Expr* res);
void  functionexpand_init(void);

#endif /* MATHILDA_FUNCTIONEXPAND_H */
