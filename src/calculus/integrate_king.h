/*
 * integrate_king.h -- Definite integration by the general King reflection
 * (reflection-averaging) technique.
 *
 * For a finite interval [a,b] the substitution x -> a+b-x is an exact
 * measure-preserving bijection, so
 *     Integrate[f(x), {x,a,b}] = (1/2) Integrate[f(x) + f(a+b-x), {x,a,b}].
 * When the reflected sum h(x) = f(x) + f(a+b-x) collapses to a CONSTANT c the
 * value is immediate:  c (b-a)/2.  This is the classic trick behind
 *     Integrate[1/(1 + Tan[x]^s), {x,0,Pi/2}] = Pi/4
 * (f(x)+f(Pi/2-x) = 1 for every real s), which has no elementary antiderivative.
 *
 * The constant c is read off as h at the interval midpoint (Simplify[h /. x->(a+b)/2]);
 * because the symbolic zero-test is defeated by irrational exponents / cofunction
 * identities (Simplify cannot prove 1/(1+Tan^s)+1/(1+Cot^s) == 1), constancy is
 * screened NUMERICALLY across several interior points and the final value is
 * cross-checked against NIntegrate.  If h is not in fact constant the candidate
 * c (b-a)/2 fails the NIntegrate check and the stage declines -- so a mis-detected
 * constancy can never produce a wrong answer.
 *
 * This complements integrate_kingreflect_try (integrate_beta.c), which handles
 * the distinct f = x g(x), g invariant, shape where h is NOT constant.
 *
 * Reachable via Method -> "King" / "Reflection" and the explicit entry point
 *   Integrate`King[f, {x, a, b}].
 */

#ifndef MATHILDA_INTEGRATE_KING_H
#define MATHILDA_INTEGRATE_KING_H

#include "expr.h"

/* General King reflection on a finite interval [a,b].  Borrowed args;
 * assumptions may be NULL.  Returns a fresh value or NULL to fall through. */
Expr* integrate_king_try(Expr* f, Expr* x, Expr* a, Expr* b, Expr* assumptions);

/* Explicit builtin.  Strict: NULL on any non-applicable input. */
Expr* builtin_integrate_king(Expr* res);

/* Register the builtin + attributes + docstring. */
void integrate_king_init(void);

#endif /* MATHILDA_INTEGRATE_KING_H */
