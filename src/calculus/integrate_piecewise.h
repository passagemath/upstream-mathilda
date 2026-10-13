/*
 * integrate_piecewise.h -- Definite integration of piecewise integrands.
 *
 * Handles integrands built from the branch-defining heads
 *   Abs, Sign, HeavisideTheta, UnitStep, Min, Max, Piecewise
 * (with arguments depending on the integration variable) over a real interval
 * [a, b].  The method:
 *   1. Collect the breakpoints -- the x-values where any branch switches
 *      (root of the Abs/Sign/step argument, u==v for Min/Max, the condition
 *      boundaries of a Piecewise).
 *   2. Keep the breakpoints strictly interior to (a, b), decided
 *      assumption-aware with Refine, and sort them.
 *   3. On each sub-interval the integrand reduces to a single smooth branch
 *      (selected by the sign of each argument at an interior test point);
 *      integrate that branch over the sub-interval with the ordinary cascade.
 *   4. Sum the pieces.
 *
 * Correctness: a branch is selected only when its governing inequality is
 * *proved* under the assumptions (Refine === True); any undecidable breakpoint
 * position, branch sign, or non-closing piece yields NULL so the integral is
 * left unevaluated rather than guessed.
 *
 * Reachable as:
 *   Integrate[f, {x, a, b}]                          (automatic dispatch)
 *   Integrate[f, {x, a, b}, Method -> "Piecewise"]
 *   Integrate`Piecewise[f, {x, a, b}]                (explicit entry point)
 */

#ifndef MATHILDA_INTEGRATE_PIECEWISE_H
#define MATHILDA_INTEGRATE_PIECEWISE_H

#include "expr.h"

/*
 * Core entry point for the Integrate dispatcher's definite path.
 *   f, x, a, b, assumptions are borrowed (not consumed); assumptions may be NULL.
 * Returns a freshly-allocated Expr* (the definite value), or NULL to leave the
 * integral unevaluated (no piecewise head, undecidable breakpoint position or
 * branch sign, or a piece that does not close).
 */
Expr* integrate_piecewise_try(Expr* f, Expr* x, Expr* a, Expr* b,
                              Expr* assumptions);

/* `Integrate`Piecewise[f, {x, a, b}]` builtin.  Strict: NULL on any
 * non-applicable input. */
Expr* builtin_integrate_piecewise(Expr* res);

/* Register the builtin + attributes + docstring. */
void integrate_piecewise_init(void);

#endif /* MATHILDA_INTEGRATE_PIECEWISE_H */
