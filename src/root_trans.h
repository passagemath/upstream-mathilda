/* root_trans.h
 *
 * Transcendental Root objects and certified real-root isolation.
 *
 *   Root[{f, x0}]
 *       The real root of the univariate (generally transcendental) function
 *       f -- a pure Function in #1, e.g. -Sin[#1] + #1/2 & -- that lies near
 *       the approximate value x0.  It is a NUMERIC constant: N refines it to
 *       any precision by interval Newton in arb, which also CERTIFIES that f
 *       has exactly one simple root in a small window around x0 (otherwise the
 *       object is invalid: Root::invrt).  This is the form Reduce emits for the
 *       breakpoints of an exp-trig-polynomial inequality, as Mathematica does.
 *
 * Isolation (root_trans_isolate) finds ALL real roots of an exact univariate
 * expression on a closed bounded interval rigorously: arb interval evaluation
 * of f excludes root-free subintervals, a derivative that excludes zero
 * certifies monotonicity (hence at most one root, decided by the endpoint
 * signs), and a root that is an exact rational (x == 0 of Sin[x] - x/2) is
 * recognised by exact evaluation, including a multiple one (its multiplicity k
 * is certified by f^(k) excluding zero, which by Rolle allows no other root).
 * Anything it cannot certify -- an irrational multiple root, two roots closer
 * than the resolution, a point outside the supported language -- DECLINES.
 *
 * Everything here needs FLINT + MPFR; without them every entry declines.
 */
#ifndef ROOT_TRANS_H
#define ROOT_TRANS_H

#include "expr.h"
#include "numeric.h"
#include <stdbool.h>

/* Structural test: e is Root[{Function[body], x0}] with x0 a real number and
 * body built only from numeric constants, #1 and NumericFunction heads. */
bool root_trans_is_object(const Expr* e);

/* N[Root[{f, x0}]] at the requested precision (certified interval Newton).
 * Owned Real / MPFR, or NULL when the root cannot be certified near x0. */
Expr* root_trans_numericalize(const Expr* e, NumericSpec spec);

/* Builtin-side canonicalisation of a transcendental Root: returns a new Root
 * whose machine x0 is the double nearest the certified root, or NULL when e
 * is already canonical (or not of this form).  *invalid is set when e has the
 * form but no simple root can be certified near x0 (Root::invrt). */
Expr* root_trans_canonicalize(const Expr* e, bool* invalid);

/* Build Root[{Function[f /. x -> #1], approx}] (f is borrowed). */
Expr* root_trans_make(const Expr* f, const Expr* x, double approx);

/* Prove that f (exact, univariate in the symbol x) has NO real root with
 * |x| >= R and return that R as an owned non-negative Integer.  Works when f
 * expands to a polynomial whose coefficients may carry BOUNDED factors (Sin,
 * Cos, Tanh, Erf, Sech, ArcTan of real polynomials) below a top-degree
 * coefficient that is a pure nonzero constant (Cauchy-type bound).  Returns
 * NULL when no such bound can be proved (Sin[x], x Sin[x] - 1, Exp[x] - x). */
Expr* root_trans_dominance_bound(const Expr* f, const Expr* x);

/* Append every real root of f in the closed interval [lo, hi] (exact numeric
 * Exprs, borrowed) to (*out, *n, *cap) as owned values: an exact Integer /
 * Rational when the root is rational, else a Root[{f, approx}] object.  Roots
 * found marginally outside [lo, hi] (the endpoints are enclosed by balls) may
 * also be appended; each is still a genuine root.  Returns false -- with the
 * appended entries freed and *n restored -- when the roots cannot all be
 * certified. */
bool root_trans_isolate(const Expr* f, const Expr* x, const Expr* lo,
                        const Expr* hi, Expr*** out, int* n, int* cap);

#endif /* ROOT_TRANS_H */
