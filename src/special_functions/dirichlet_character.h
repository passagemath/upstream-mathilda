#ifndef MATHILDA_DIRICHLET_CHARACTER_H
#define MATHILDA_DIRICHLET_CHARACTER_H

#include "expr.h"

/* DirichletCharacter[k, j, n]: the Dirichlet character chi_{k,j}(n). */
Expr* builtin_dirichletcharacter(Expr* res);
void dirichletcharacter_init(void);

#endif
