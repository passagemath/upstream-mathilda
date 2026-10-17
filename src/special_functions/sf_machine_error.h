/* Machine-precision kernels for the error-function family:
 *   Erf, Erfc, Erfi (complex), FresnelC / FresnelS (real and complex),
 *   Erfi (real), Sinc (complex), InverseErf / InverseErfc (real).
 *
 * Every kernel works on raw `double` real/imaginary parts (no <complex.h>)
 * and follows the kernel contract of sf_machine.h: it returns true and writes
 * the value, or returns false ("no usable machine value": overflow, a domain
 * the kernel does not cover) so the caller falls back to its arbitrary-precision
 * path unchanged.
 *
 * The SAME helper serves every surface -- the scalar builtin, the NDArray /
 * packed element-wise map (src/ndkernels.c) and Compile[] (through the
 * registered kernel) -- so the five surfaces agree bit for bit.
 *
 * Algorithms and their sources are documented at each kernel in
 * sf_machine_error.c. */
#ifndef MATHILDA_SF_MACHINE_ERROR_H
#define MATHILDA_SF_MACHINE_ERROR_H

#include <stdbool.h>

/* Real erf and erfc (all real x; erfc underflows to 0 past x = 27.39). */
bool sf_erf_r (double x, double* out);
bool sf_erfc_r(double x, double* out);

/* Complex error functions of z = x + i y. */
bool sf_erf_c (double x, double y, double* re, double* im);
bool sf_erfc_c(double x, double y, double* re, double* im);
bool sf_erfi_c(double x, double y, double* re, double* im);

/* Real erfi(x) = -i erf(i x). */
bool sf_erfi_r(double x, double* out);

/* Fresnel integrals, Mathematica's pi/2 convention:
 *   C(z) = int_0^z cos(pi t^2/2) dt,  S(z) = int_0^z sin(pi t^2/2) dt. */
bool sf_fresnel_c_r(double x, double* out);
bool sf_fresnel_s_r(double x, double* out);
bool sf_fresnel_c_c(double x, double y, double* re, double* im);
bool sf_fresnel_s_c(double x, double y, double* re, double* im);

/* Sinc(z) = sin(z)/z, Sinc(0) = 1. */
bool sf_sinc_r(double x, double* out);
bool sf_sinc_c(double x, double y, double* re, double* im);

/* Inverse error functions on the real line. inverse_erf needs -1 < p < 1,
 * inverse_erfc needs 0 < q < 2; the endpoints (+-Infinity) and anything
 * outside decline. */
bool sf_inverse_erf_r (double p, double* out);
bool sf_inverse_erfc_r(double q, double* out);

#endif /* MATHILDA_SF_MACHINE_ERROR_H */
