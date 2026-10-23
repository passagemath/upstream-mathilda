#ifndef PIECEWISE_H
#define PIECEWISE_H

#include "expr.h"

Expr* builtin_floor(Expr* res);
Expr* builtin_ceiling(Expr* res);
Expr* builtin_round(Expr* res);
Expr* builtin_integerpart(Expr* res);
Expr* builtin_fractionalpart(Expr* res);
Expr* builtin_unitstep(Expr* res);
Expr* builtin_ramp(Expr* res);
Expr* builtin_unitbox(Expr* res);
Expr* builtin_unitize(Expr* res);

/* Unitize's machine helpers -- the ONE implementation of the machine-number
 * arithmetic behind Unitize, shared by the scalar builtin (machine Integer /
 * Real arguments) and the packed / NDArray kernels in src/ndkernels.c
 * (NDKU_Unitize, NDKB_Unitize -- and through them SparseArray's value-buffer
 * threading), so every surface answers bit for bit alike. Each writes the
 * exact 0/1 to *out and returns true, or returns false to DECLINE (a NaN
 * argument, a threshold that is not a finite positive number); the caller then
 * falls back to the exact path, which answers or leaves the call unevaluated.
 *
 *   unitize_machine_r(x)          x != 0  (-0. counts as zero)
 *   unitize_machine_i(x)          x != 0
 *   unitize_machine_tol_r(x, dx)  |x| >= dx, dx > 0 and finite
 *   unitize_machine_tol_i(x, dx)  |x| >= dx, dx > 0, written as
 *                                 x >= dx || x <= -dx so INT64_MIN cannot
 *                                 overflow an Abs. */
bool unitize_machine_r(double x, int64_t* out);
bool unitize_machine_i(int64_t x, int64_t* out);
bool unitize_machine_tol_r(double x, double dx, int64_t* out);
bool unitize_machine_tol_i(int64_t x, int64_t dx, int64_t* out);

void piecewise_init(void);

#endif // PIECEWISE_H
