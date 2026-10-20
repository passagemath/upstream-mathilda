#ifndef STIELTJESGAMMA_H
#define STIELTJESGAMMA_H

#include "expr.h"

#include "numeric.h"

/* StieltjesGamma[n] -- the n-th Stieltjes constant gamma_n, the coefficients
 * in the Laurent expansion of the Riemann zeta function about s = 1:
 *
 *   zeta(s) = 1/(s-1) + Sum_{n>=0} ((-1)^n / n!) gamma_n (s-1)^n,
 *
 * and StieltjesGamma[n, a], the same for the Hurwitz zeta function.
 * gamma_0 = EulerGamma. An exact call stays symbolic (except
 * StieltjesGamma[0] -> EulerGamma); an inexact a evaluates numerically. It is
 * produced by the series expansion of Zeta about s = 1.
 *
 * Attributes: Listable, NHoldFirst, Protected. */
Expr* builtin_stieltjesgamma(Expr* res);

/* N's route to the value: `call` is StieltjesGamma[n] or StieltjesGamma[n, a]
 * with an exact order n (NHoldFirst keeps it exact); returns the numeric value
 * at `spec`, or NULL. Called from numericalize. */
Expr* numeric_stieltjes_n(const Expr* call, NumericSpec spec);

void stieltjesgamma_init(void);

#endif /* STIELTJESGAMMA_H */
