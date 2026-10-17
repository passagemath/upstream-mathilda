#ifndef MATHILDA_DIRICHLET_L_H
#define MATHILDA_DIRICHLET_L_H

#include "expr.h"

/* DirichletL[k, j, s]: the Dirichlet L-function of DirichletCharacter[k, j, .]. */
Expr* builtin_dirichletl(Expr* res);
void dirichletl_init(void);

#endif
