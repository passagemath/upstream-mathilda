/*
 * integrate_series.h -- Definite integration by series expansion, term-by-term
 * integration, and resummation (the "beat it with a series" technique).
 *
 * Closes log-weighted rational integrals on [0,1] whose antiderivative is a
 * non-elementary PolyLog and whose denominator is a root-of-unity cyclotomic
 * factor 1 +/- x^d.  It is the ALTERNATING sibling of integrate_ratlog01_try
 * (integrate_beta.c): that stage closes the 1 - x^m family via HurwitzZeta but
 * cannot touch 1 + x^d, whose resummation needs Catalan / DirichletBeta /
 * Dirichlet-eta -- all of which the existing Sum[] engine already delivers.
 *
 * Method (interval exactly [0,1], integrand R(x) Log[x]^k with R rational, k>=1):
 *   1. Strip the Log[x]^k weight (k>=1 required: it is the whole point, and it
 *      forces the resummed series to be a p-series with p = k+1 >= 2, i.e.
 *      ABSOLUTELY convergent -- the correctness guarantee for interchanging the
 *      sum and the integral).
 *   2. Put R over a common denominator N(x)/Q(x), Q(0) != 0, and find the
 *      smallest d (and sign eps) with Q | (1 - eps x^d) by an EXACT polynomial
 *      remainder.  A near-miss leaves a nonzero remainder and is declined, so a
 *      false denominator match is impossible.
 *   3. 1/(1 - eps x^d) = Sum[eps^n x^(d n)], so with N = Quot (1 - eps x^d) + S
 *      (deg S < d) the value is
 *        Sum_m [x^m in Quot]  c_m (-1)^k k! / (m+1)^(k+1)          (finite)
 *      + Sum_j [x^j in S]     s_j (-1)^k k! Sum[eps^n/(d n+j+1)^(k+1),{n,0,Inf}]
 *      using the exact identity  Integrate[x^m Log[x]^k,{x,0,1}]
 *      = (-1)^k k!/(m+1)^(k+1).  Each inner Sum is handed verbatim to Sum[].
 *   4. The closed form is numerically cross-checked against NIntegrate before it
 *      is returned; any disagreement (or a Sum[] that does not close) declines.
 *
 * Correct-by-construction on the rational cyclotomic class; declines everything
 * else (non-rational factor, non-[0,1] interval, pole at 0, denominator not a
 * root-of-unity factor within Dmax, or an unclosed resummation).
 *
 * Reachable via Method -> "SeriesExpansion" / "SeriesTermByTerm" and the
 * explicit entry point  Integrate`SeriesExpansion[f, {x, 0, 1}].
 */

#ifndef MATHILDA_INTEGRATE_SERIES_H
#define MATHILDA_INTEGRATE_SERIES_H

#include "expr.h"

/* Series term-by-term definite integration on [0,1] of R(x) Log[x]^k.  Borrowed
 * args; assumptions may be NULL.  Returns a fresh value or NULL to fall through. */
Expr* integrate_series_try(Expr* f, Expr* x, Expr* a, Expr* b, Expr* assumptions);

/* Explicit builtin.  Strict: NULL on any non-applicable input. */
Expr* builtin_integrate_series(Expr* res);

/* Register the builtin + attributes + docstring. */
void integrate_series_init(void);

#endif /* MATHILDA_INTEGRATE_SERIES_H */
