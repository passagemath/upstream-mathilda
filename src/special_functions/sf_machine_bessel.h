/* Machine-precision (double) kernels for the Bessel and Airy functions.
 *
 * One home for the raw-double Bessel J/Y/I/K and Airy Ai/Bi/Ai'/Bi' kernels,
 * real and complex. Each entry point is shared, bit for bit, by every surface
 * that evaluates the head at machine precision:
 *   - the scalar builtin (bessel.c / airyai.c / airybi.c), via the
 *     sf_machine_scalar_args gate;
 *   - the NDArray / packed element-wise path (ndkernels.c);
 *   - the Compile[] VM, which calls the same NDBinaryKernel / NDUnaryKernel.
 *
 * Contract, as for sf_machine.h: a kernel returns FALSE when it cannot deliver
 * a machine answer -- the true value is complex (a REAL entry point on the
 * branch cut), infinite (a pole), overflows or underflows the double range, or
 * an algorithm's convergence test fails. The caller then falls back to the MPFR
 * / arb path, which handles those cases with an extended exponent. A kernel
 * never returns a value it has not computed to machine accuracy.
 *
 * Complex values are raw (re, im) double pairs; no <complex.h> (portability,
 * and libgcc's __divdc3 is slow).
 *
 * Orders are REAL (machine or exact); a complex order is left to MPFR.
 */
#ifndef MATHILDA_SF_MACHINE_BESSEL_H
#define MATHILDA_SF_MACHINE_BESSEL_H

#include <stdbool.h>

/* Real order nu, real argument x. Declines where the value is complex
 * (x < 0 with non-integer nu; x < 0 for Y and K) or infinite. */
bool sf_bessel_j_r(double nu, double x, double* out);
bool sf_bessel_y_r(double nu, double x, double* out);
bool sf_bessel_i_r(double nu, double x, double* out);
bool sf_bessel_k_r(double nu, double x, double* out);

/* Real order nu, complex argument z = zr + i zi, principal branch (cut along
 * the negative real axis; zi == 0 with zr < 0 is taken on the upper side,
 * arg z = pi, as Mathematica does). */
bool sf_bessel_j_c(double nu, double zr, double zi, double* rr, double* ri);
bool sf_bessel_y_c(double nu, double zr, double zi, double* rr, double* ri);
bool sf_bessel_i_c(double nu, double zr, double zi, double* rr, double* ri);
bool sf_bessel_k_c(double nu, double zr, double zi, double* rr, double* ri);

/* Airy functions of a real argument. */
bool sf_airy_ai_r(double x, double* out);
bool sf_airy_bi_r(double x, double* out);
bool sf_airy_aip_r(double x, double* out);
bool sf_airy_bip_r(double x, double* out);

/* Airy functions of a complex argument (entire: no branch cut). */
bool sf_airy_ai_c(double zr, double zi, double* rr, double* ri);
bool sf_airy_bi_c(double zr, double zi, double* rr, double* ri);
bool sf_airy_aip_c(double zr, double zi, double* rr, double* ri);
bool sf_airy_bip_c(double zr, double zi, double* rr, double* ri);

/* ---- the shared dispatch every surface calls ----------------------------
 *
 * Kind selectors for the Bessel dispatch. */
enum { SF_BESSEL_J = 0, SF_BESSEL_Y = 1, SF_BESSEL_I = 2, SF_BESSEL_K = 3 };
enum { SF_AIRY_AI = 0, SF_AIRY_AIP = 1, SF_AIRY_BI = 2, SF_AIRY_BIP = 3 };

/* The ND / Compile kernel, NDBinaryKernel.cplx shape: order (nre, nim),
 * argument (zre, zim). A complex order declines. A real argument (zim == 0)
 * takes the REAL kernel and declines where the value is not real -- the
 * array path's real-closed contract; the builtin then answers that element
 * through the scalar path. A complex argument takes the complex kernel. */
bool sf_bessel_nd(int kind, double nre, double nim, double zre, double zim,
                  double* rr, double* ri);

/* The same for Airy (NDUnaryKernel.cplx shape). */
bool sf_airy_nd(int which, double zre, double zim, double* rr, double* ri);

/* The scalar builtin's dispatch: real order nu, argument (zr, zi). Tries the
 * real kernel when zi == 0 and, if that declines or zi != 0, the complex one.
 * *is_real reports which answered (the caller builds a Real only when the
 * input was real and the real kernel answered). Bit-identical to
 * sf_bessel_nd wherever sf_bessel_nd answers. */
bool sf_bessel_scalar(int kind, double nu, double zr, double zi,
                      double* rr, double* ri, bool* is_real);
bool sf_airy_scalar(int which, double zr, double zi,
                    double* rr, double* ri, bool* is_real);

#endif /* MATHILDA_SF_MACHINE_BESSEL_H */
