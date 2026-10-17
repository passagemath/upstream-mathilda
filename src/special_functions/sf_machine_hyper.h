/* Machine-precision kernels for the hypergeometric family: 0F1, 1F1, 2F1, pFq,
 * HypergeometricU and the Whittaker functions -- real `double` and raw (re, im)
 * double-complex.
 *
 * These are THE machine kernels for the family: the scalar builtins (for
 * machine-precision arguments), the NDArray / packed element-wise paths and the
 * Compile[] VM all call them, so the five surfaces agree bit-for-bit.
 *
 * The REAL entry points keep the names declared in sf_machine.h
 * (sf_machine_hyper0f1 / _hyper1f1 / _hyper2f1 / _pfq); they are defined in
 * sf_machine_hyper.c. This header adds the complex twins and the
 * constant-parameter buffer driver.
 *
 * Contract, common to every kernel: return true and write a finite result, or
 * return false ("decline") -- for a pole, a branch cut the real kernel cannot
 * represent, an overflow, or a regime where the double algorithms cannot certify
 * ~machine accuracy (measured cancellation too large, series too long). A decline
 * is never an error: the caller falls through to the arbitrary-precision path.
 *
 * Algorithms (see sf_machine_hyper.c for the citations):
 *   - every series measures its own cancellation (peak term / |sum|) and re-sums
 *     in double-double arithmetic when the double sum would lose bits;
 *   - 0F1: series, Miller backward recurrence for J_nu (z < 0), Hankel
 *     asymptotic for I_nu (large z > 0);
 *   - 1F1: Kummer's transformation for Re z < 0, DLMF 13.7.2 asymptotic for
 *     large |z|, elementary closed forms (a = 0, a = b, (a,b) = (1,2));
 *   - 2F1: Pfaff (z/(z-1)), 1-z and 1/(1-z) transformations with the
 *     degenerate (c-a-b or b-a integer) logarithmic limit formulas
 *     (Abramowitz & Stegun 15.3.6, 15.3.10-15.3.12), Gauss's sum at z = 1.
 */
#ifndef MATHILDA_SF_MACHINE_HYPER_H
#define MATHILDA_SF_MACHINE_HYPER_H

#include <stdbool.h>
#include <stddef.h>

/* pFq with complex parameters and argument. Arrays ar/ai (length p) and br/bi
 * (length q) carry the parameters' real/imaginary parts. Dispatches the
 * (p, q) = (0,1), (1,1), (2,1) shapes to the specialised kernels below, so
 * HypergeometricPFQ[{a}, {b}, z] and Hypergeometric1F1[a, b, z] agree bit-for-bit. */
bool sf_machine_pfq_c(const double* ar, const double* ai, size_t p,
                      const double* br, const double* bi, size_t q,
                      double zr, double zi, double* rr, double* ri);

/* Buffer driver for REAL parameters and a real buffer: out[i] = pFq(a; b; x[i]).
 * Returns how many leading elements were computed: n on success, otherwise the
 * index of the first element the kernel declined (the caller then takes the
 * element-wise route for the whole buffer). Bit-identical to calling
 * sf_machine_pfq per element -- it only hoists the parameter-dependent work (the
 * 2F1 coefficient tables) out of the loop. */
size_t sf_machine_pfq_vec(const double* a, size_t p, const double* b, size_t q,
                          const double* x, double* out, size_t n);

/* 0F1(; b; z), complex b and z. */
bool sf_machine_hyper0f1_c(double br, double bi, double zr, double zi,
                           double* rr, double* ri);

/* 1F1(a; b; z), complex a, b, z. */
bool sf_machine_hyper1f1_c(double ar, double ai, double br, double bi,
                           double zr, double zi, double* rr, double* ri);

/* 2F1(a, b; c; z), complex parameters and argument. The transformations need
 * Gamma of the parameters, so they are applied for REAL a, b, c only; complex
 * parameters are summed directly inside |z| <= 0.9 and decline elsewhere.
 * On the cut z > 1 (zi == 0) the value is the limit from below (zi = -0), the
 * convention of Mathematica and mpmath; pass zi = -0.0 explicitly to select it. */
bool sf_machine_hyper2f1_c(double ar, double ai, double br, double bi,
                           double cr, double ci, double zr, double zi,
                           double* rr, double* ri);

/* Real U(a, b, z) for z > 0, a > 0 by double-exponential quadrature of its
 * integral representation -- the moderate-z regime of sf_machine_hyperu, used
 * when the asymptotic series and the connection formulas cannot certify the
 * result. Declines when the integrand is too sharply peaked for the rule. */
bool sf_machine_hyperu_quad(double a, double b, double z, double* out);

/* HypergeometricU(a, b, z) for complex z (real a, b). */
bool sf_machine_hyperu_c(double a, double b, double zr, double zi,
                         double* rr, double* ri);

/* WhittakerM / WhittakerW with complex z (real k, m), principal branch of
 * z^(m+1/2). */
bool sf_machine_whittaker_m_c(double k, double m, double zr, double zi,
                              double* rr, double* ri);
bool sf_machine_whittaker_w_c(double k, double m, double zr, double zi,
                              double* rr, double* ri);

#endif /* MATHILDA_SF_MACHINE_HYPER_H */
