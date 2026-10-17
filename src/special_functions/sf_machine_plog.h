/* Machine-precision (double) kernels for ProductLog (Lambert W, principal
 * branch) and the integer-degree Legendre functions LegendreP / LegendreQ,
 * real and complex.
 *
 * Every entry point is shared bit-for-bit by the scalar builtin (the
 * sf_machine_scalar_args gate in productlog.c / legendre.c), the NDArray
 * element-wise kernels (src/ndkernels.c) and the Compile[] VM, per CLAUDE.md
 * ("the machine-scalar path reuses the kernel's helper").  Each returns FALSE
 * ("no usable machine value") outside its regime, and the caller falls back to
 * the MPFR / symbolic path unchanged.  Complex arithmetic is spelled out on
 * re/im doubles: no <complex.h>.
 */
#ifndef MATHILDA_SF_MACHINE_PLOG_H
#define MATHILDA_SF_MACHINE_PLOG_H

#include <stdbool.h>

/* W0(z) for a complex z = zr + i zi (principal branch).  A real z >= -1/e
 * (zi == 0) answers through the real kernel sf_machine_productlog with *wi = 0,
 * so the complex arm agrees bit-for-bit with the real one on the real axis.  A
 * real z < -1/e, or any z off the axis, takes the complex Halley iteration.
 * Declines for a non-finite input/result or an iterate that fails the
 * w e^w = z / W0-region verification. */
bool sf_machine_productlog_c(double zr, double zi, double* wr, double* wi);

/* LegendreP[n, x] for an integer-valued degree n (any sign: P_{-n-1} = P_n) and
 * x = xr + i xi.  xi == 0 answers through the real recurrence with *ri = 0. */
bool sf_machine_legendre_p_c(double n, double xr, double xi, double* rr, double* ri);

/* LegendreQ[n, x] (type 1, Mathematica's default) for an integer degree n >= 0
 * and x = xr + i xi, with
 *     Q_0(x) = (Log[1 + x] - Log[1 - x]) / 2      (principal logs)
 * and the P-recurrence upward from Q_0, Q_1 = x Q_0 - 1.  On the real axis the
 * convention is the closed form's (legendre.c legendre_q_closed): a real x with
 * |x| > 1 has Log[1 - x] = log(x - 1) + i Pi for x > 1 and Log[1 + x] =
 * log(-1 - x) + i Pi for x < -1, i.e. a negative-zero imaginary part is read as
 * +0.  xi == 0 with |xr| < 1 answers through the real kernel (*ri = 0).
 * Declines at x = +-1 (the logarithmic singularity). */
bool sf_machine_legendre_q_c(double n, double xr, double xi, double* rr, double* ri);

#endif /* MATHILDA_SF_MACHINE_PLOG_H */
