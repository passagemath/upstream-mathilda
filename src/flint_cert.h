/* flint_cert.h
 *
 * Rigorous arb ball / interval evaluation of an exact univariate expression,
 * exported from flint_num_bridge.c (the evaluator behind flint_cert_parts_sign)
 * for callers that already speak arb.  Only available when the build has FLINT
 * and MPFR; callers guard their use with the same condition.
 */
#ifndef FLINT_CERT_H
#define FLINT_CERT_H

#if defined(USE_FLINT) && defined(USE_MPFR)
#include <flint/arb.h>
#include "expr.h"

/* Enclose the value of the exact expression `e`, with the symbol `var`
 * (interned name, or NULL for a constant) bound to the real ball `xball`, in
 * the real ball `out` at working precision `prec` bits.  The supported
 * language is that of flint_cert_parts_sign (rationals, Pi/E/..., Plus, Times,
 * Power, Exp, Log, Sqrt, trig/hyperbolic and their inverses, Gamma, Erf, Abs,
 * ...).  Returns 1 when the enclosure is finite and certified REAL (exact-zero
 * imaginary part, e.g. no Sqrt of a ball reaching below 0), 0 otherwise. */
int flint_cert_eval_real(const Expr* e, const char* var, const arb_t xball,
                         arb_t out, slong prec);
#endif

#endif /* FLINT_CERT_H */
