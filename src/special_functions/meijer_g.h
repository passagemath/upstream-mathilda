/*
 * meijer_g.h -- MeijerG, the Meijer G-function, for Mathilda.
 *
 *   MeijerG[{{a1..an},{an+1..ap}}, {{b1..bm},{bm+1..bq}}, z]
 *   MeijerG[alist, blist, z, r]          (generalized form; r = 1 by default)
 *
 * The builtin (meijer_g.c) owns the exact and symbolic surface: argument
 * validation and the hdiv / argt messages, canonical parameter order, the
 * values at z = 0, and the automatic reduction to simpler functions. The
 * numerics live in sf_machine_meijerg.c (machine) and flint_num_bridge.c (arb).
 */
#ifndef MATHILDA_MEIJER_G_H
#define MATHILDA_MEIJER_G_H

#include <stdbool.h>
#include "expr.h"

Expr* builtin_meijer_g(Expr* res);
void  meijer_g_init(void);

/* Slater's theorem as an Expr: the expansion of MeijerG `g` (a well-formed
 * 3- or 4-argument MeijerG call) as a sum of HypergeometricPFQ terms --
 * the expansion in z for p < q (or p == q, the shorter of the two), in 1/z
 * for p > q. With `named`, each term's pFq is rewritten to Bessel / Exp /
 * power / elliptic forms where one exists. Returns NULL (nothing owned) when
 * there is no such expansion: the log case (parameters an integer apart),
 * p == q with m + n == p (the two expansions are not continuations of each
 * other), or a malformed call. The result is evaluated and owned by the
 * caller. Used by FunctionExpand and Series. */
Expr* meijerg_slater_expand(const Expr* g, bool named);

/* The leading-power data the series code needs, in the variable z of `g`:
 * the truncated expansion  sum_h  z^{e_h} * sum_{k<=order} c_{h,k} z^k
 * about 0 (`at_infinity` false, p <= q) or in powers of 1/z about Infinity
 * (`at_infinity` true, p >= q), with the log case resolved by perturbation
 * (the result then carries Log[z] and PolyGamma terms). Returns an evaluated,
 * owned Expr that is a finite sum (Normal form, no O[] term), or NULL. */
Expr* meijerg_series_terms(const Expr* g, bool at_infinity, long order);

#endif /* MATHILDA_MEIJER_G_H */
