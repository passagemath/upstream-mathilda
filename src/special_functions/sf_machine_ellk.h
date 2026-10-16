/* Machine-precision (double) kernels for the COMPLETE elliptic integrals
 * K(m) and E(m), parameter convention m = k^2.
 *
 * One kernel per (function, real|complex), shared bit-for-bit by every surface:
 * the scalar builtins (EllipticK[m], EllipticE[m] in elliptic.c), the NDArray
 * element-wise path (src/ndkernels.c) and the Compile[] VM (which calls the
 * same registered kernels).
 *
 * Contract (the engine's standing one): return false for "no usable machine
 * value" -- a pole, a non-finite argument or result, or (for the REAL entry
 * points) an argument whose value is genuinely complex (m > 1). The caller then
 * falls back to its slower exact path; nothing is invented.
 */
#ifndef MATHILDA_SF_MACHINE_ELLK_H
#define MATHILDA_SF_MACHINE_ELLK_H

#include <stdbool.h>

/* K(m) for real m < 1. Declines m >= 1 (pole at 1, complex above). */
bool sf_machine_ellk_k(double m, double* out);
/* E(m) for real m <= 1 (E(1) = 1). Declines m > 1 (complex). */
bool sf_machine_ellk_e(double m, double* out);

/* K(m), E(m) for complex m = mre + i mim, principal branch (cut on the real
 * axis m >= 1, where the value is that of the limit from below, Im m -> 0-,
 * which is the Mathematica / mpmath convention: K(2) = 1.31103 - 1.31103 i).
 * A zero imaginary part of either sign is treated as +0. Declines m = 1. */
bool sf_machine_ellk_k_c(double mre, double mim, double* rre, double* rim);
bool sf_machine_ellk_e_c(double mre, double mim, double* rre, double* rim);

#endif /* MATHILDA_SF_MACHINE_ELLK_H */
