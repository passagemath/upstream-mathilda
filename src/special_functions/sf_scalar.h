/* The machine-scalar gate shared by the special-function builtins.
 *
 * A special-function builtin sees `F[a1, ..., an]`. When the call is a
 * MACHINE-precision numeric evaluation, the answer must come from the same
 * raw-double kernel (sf_machine.c and its complex siblings) that the NDArray
 * element-wise path and the Compile[] VM use:
 *   - It is 10-100x faster than computing in MPFR / an `ncpx` pair and then
 *     rounding.
 *   - It keeps the scalar answer bit-identical to the packed / NDArray /
 *     Compile answers (see CLAUDE.md, "the machine-scalar path reuses the
 *     kernel's helper").
 * The test for "this is a machine-precision call" was re-derived differently in
 * every module. This header is the one definition.
 *
 * A call is machine-scalar when every argument is one of:
 *   - a machine Real
 *   - an exact Integer / BigInt / Rational (these impose no precision of their own)
 *   - a Complex[re, im] whose parts are each one of the above
 * and AT LEAST ONE argument is a machine Real or a Complex with a machine-Real
 * part. An MPFR leaf anywhere makes it an arbitrary-precision call (false), and
 * so does an all-exact call (an exact call stays exact). Symbols, including
 * Infinity and Pi, also decline: the special-value ladders own them.
 *
 * Usage, at the TOP of a builtin's numeric branch, before any strcmp-on-name
 * special-value ladder:
 *
 *     SFMachineArgs ma;
 *     if (sf_machine_scalar_args(res, &ma)) {
 *         if (!ma.any_complex) {
 *             double v;
 *             if (sf_machine_foo(ma.re[0], &v)) return expr_new_real(v);
 *         } else {
 *             double vr, vi;
 *             if (sf_machine_foo_c(ma.re[0], ma.im[0], &vr, &vi))
 *                 return sf_machine_complex_result(vr, vi);
 *         }
 *         // kernel declined (domain / overflow): fall through unchanged
 *     }
 */
#ifndef MATHILDA_SF_SCALAR_H
#define MATHILDA_SF_SCALAR_H

#include <stdbool.h>
#include "expr.h"

#define SF_MACHINE_MAXARGS 6

typedef struct {
    int    n;                         /* argument count */
    double re[SF_MACHINE_MAXARGS];    /* real part of each argument */
    double im[SF_MACHINE_MAXARGS];    /* imaginary part (0 for a real argument) */
    bool   is_complex[SF_MACHINE_MAXARGS]; /* argument i was a Complex[...] */
    bool   is_exact[SF_MACHINE_MAXARGS];   /* argument i was exact (Integer / Rational) */
    bool   any_complex;               /* some argument was a Complex[...] */
} SFMachineArgs;

/* Classify the arguments of the call `res` (an EXPR_FUNCTION). Returns true
 * iff it is a machine-scalar call (see above), and fills `out`. Allocation-free
 * and cheap: one switch per leaf. */
bool sf_machine_scalar_args(const Expr* res, SFMachineArgs* out);

/* Read one leaf as a machine real, accepting exactly the leaf kinds listed
 * above. `*inexact` is set when the leaf was a machine Real. */
bool sf_machine_leaf(const Expr* e, double* out, bool* inexact);

/* Build a machine-complex result as Complex[re, im]. This is the shape the
 * evaluator's own complex arithmetic produces, so a Complex with a zero
 * imaginary part stays Complex: a machine-complex input yields a machine-complex
 * output, as in Mathematica. */
Expr* sf_machine_complex_result(double re, double im);

#endif /* MATHILDA_SF_SCALAR_H */
