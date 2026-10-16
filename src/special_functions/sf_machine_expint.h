/* Machine-precision (double) kernels for the exponential-integral family:
 * SinIntegral, CosIntegral, SinhIntegral, CoshIntegral, ExpIntegralEi, E_1,
 * LogIntegral and ExpIntegralE[n, x].
 *
 * ONE implementation serves every surface. The scalar builtins
 * (sinintegral.c, ...) call these for a machine-precision argument, the NDArray
 * element-wise kernels (ndkernels.c) map them over a buffer, and Compile[]
 * reaches them through the same kernel registry -- so a machine answer is
 * bit-identical whichever surface produced it (CLAUDE.md: "the machine-scalar
 * path reuses the kernel's helper").
 *
 * Real entry points (declared in sf_machine.h, defined in sf_machine_expint.c):
 *   sf_machine_si / _ci / _shi / _chi / _ei / _e1 / _li
 * Complex entry points (raw re/im doubles, no <complex.h>):
 *   sf_machine_*_c below.
 *
 * Contract: return false (decline) wherever the double answer is unusable --
 * a pole (Ei(0), li(1), Ci(0)), a non-finite argument, an overflow or a result
 * that would underflow to a subnormal. The caller then falls back to the
 * MPFR/arb path, which owns those cases. A real kernel also declines where the
 * true value is complex (Ci/Chi/li of a negative argument): the complex kernels
 * cover those.
 *
 * Accuracy (tools/sf_accuracy/expint.py against mpmath at 30 digits): see the
 * per-head notes in docs/spec/builtins/special_functions.md.
 */
#ifndef MATHILDA_SF_MACHINE_EXPINT_H
#define MATHILDA_SF_MACHINE_EXPINT_H

#include <stdbool.h>

/* Complex kernels, principal branches matching the interpreter (Mathematica):
 *   Si, Shi    entire
 *   Ci, Chi    cut along (-inf, 0]; the upper side (Im z >= +0) on the cut
 *   Ei         cut along (-inf, 0]: Ei(z) = -E1(-z) + i pi sgn(Im z)
 *   li         li(z) = Ei(log z) with the principal log
 *   E1         cut along (-inf, 0]                                        */
bool sf_machine_si_c(double re, double im, double* ore, double* oim);
bool sf_machine_ci_c(double re, double im, double* ore, double* oim);
bool sf_machine_shi_c(double re, double im, double* ore, double* oim);
bool sf_machine_chi_c(double re, double im, double* ore, double* oim);
bool sf_machine_ei_c(double re, double im, double* ore, double* oim);
bool sf_machine_e1_c(double re, double im, double* ore, double* oim);
bool sf_machine_li_c(double re, double im, double* ore, double* oim);

/* E_n(x), the generalised exponential integral ExpIntegralE[n, x], for an
 * integer order n >= 0 and real x > 0. */
bool sf_machine_expn(double n, double x, double* out);

#endif /* MATHILDA_SF_MACHINE_EXPINT_H */
