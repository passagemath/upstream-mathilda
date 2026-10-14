/*
 * integrate_feynparam.h -- Definite integration of the parameter / Frullani-
 * logarithm family by differentiation under the integral sign.
 *
 * Closes the [0,1] integrals
 *     Integrate[(Sum_i c_i x^(e_i)) / Log[x], {x, 0, 1}] = Sum_i c_i Log[e_i + 1]
 * valid when Sum_i c_i = 0 (so the 1/Log[x] singularity at x = 1 is integrable)
 * and every e_i + 1 > 0.  The archetype is the "Feynman on a present parameter"
 * integral
 *     Integrate[(x^a - 1)/Log[x], {x, 0, 1}] = Log[a + 1]   (a > -1),
 * obtained by differentiating under the integral sign w.r.t. a -- d/da of the
 * integrand is x^a (the Log cancels), Integrate[x^a,{x,0,1}] = 1/(a+1), and
 * integrating back from a = 0 (where the integrand vanishes) gives Log[a+1].
 * By linearity this extends to any x-power combination whose coefficients sum to
 * zero, e.g. Integrate[(x^a - x^b)/Log[x], {x,0,1}] = Log[(a+1)/(b+1)].
 *
 * The inner x-power integral is taken DIRECTLY as 1/(e+1), so the method does no
 * heavy re-entrant Integrate (it never forms the Beta/Gamma shape the simplifier
 * cannot reduce, and stays clear of the long-running parametric cascade).  Each
 * e_i + 1 > 0 is proved from Assumptions, and the closed form is numerically
 * cross-checked against NIntegrate at an in-domain sample before it is returned;
 * the stage declines (returns NULL) otherwise.
 *
 * Reachable via Method -> "FeynmanParameter" and the explicit entry point
 *   Integrate`FeynmanParameter[f, {x, 0, 1}].
 */

#ifndef MATHILDA_INTEGRATE_FEYNPARAM_H
#define MATHILDA_INTEGRATE_FEYNPARAM_H

#include "expr.h"

/* Parameter / Frullani-logarithm family on [0,1].  Borrowed args; assumptions
 * may be NULL.  Returns a fresh value or NULL to fall through. */
Expr* integrate_feynparam_try(Expr* f, Expr* x, Expr* a, Expr* b, Expr* assumptions);

/* Explicit builtin.  Strict: NULL on any non-applicable input. */
Expr* builtin_integrate_feynparam(Expr* res);

/* Register the builtin + attributes + docstring. */
void integrate_feynparam_init(void);

#endif /* MATHILDA_INTEGRATE_FEYNPARAM_H */
