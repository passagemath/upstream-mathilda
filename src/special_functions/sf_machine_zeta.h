/* Machine-precision (raw double / double re-im) kernels for the zeta family:
 * Zeta, HurwitzZeta (and two-argument Zeta), LerchPhi and PolyLog.
 *
 * These are the ONE implementation behind every machine surface of those heads:
 * the scalar builtin (machine-precision arguments), the packed / NDArray
 * element-wise kernels (src/ndkernels.c) and the Compile[] VM, which reaches the
 * same registry. Keeping a single helper is what makes the surfaces agree
 * bit-for-bit (CLAUDE.md, "the machine-scalar path reuses the kernel's helper").
 *
 * Contract, as for every sf_machine_* kernel: return true and write a finite
 * result, or return false ("decline") for an argument outside the kernel's
 * domain, at a pole, or where the double result would overflow / lose accuracy.
 * The caller then falls through to the arbitrary-precision path unchanged.
 *
 * No <complex.h>: complex values travel as (re, im) double pairs.
 */
#ifndef MATHILDA_SF_MACHINE_ZETA_H
#define MATHILDA_SF_MACHINE_ZETA_H

#include <stdbool.h>

/* Riemann zeta(s), real s != 1. */
bool sf_zeta_r(double s, double* out);
/* Riemann zeta(s), complex s. */
bool sf_zeta_c(double sr, double si, double* vr, double* vi);

/* Hurwitz zeta(s, a) = Sum_{k>=0} (k+a)^-s (principal powers), real s != 1,
 * real a > 0. */
bool sf_hurwitz_r(double s, double a, double* out);
/* Hurwitz zeta, complex s and complex a with Re a > 0. */
bool sf_hurwitz_c(double sr, double si, double ar, double ai,
                  double* vr, double* vi);

/* Zeta[s, a]: Hurwitz zeta in the symmetric convention -- for a real a <= 0
 * the terms are |k+a|^-s and the k + a = 0 term is dropped (0^0 = 1 at s = 0).
 * Equal to sf_hurwitz_* for Re a > 0 / complex a. */
bool sf_zeta2_r(double s, double a, double* out);
bool sf_zeta2_c(double sr, double si, double ar, double ai, double* vr, double* vi);

/* LerchPhi(z, s, a) = Sum_{k>=0} z^k / (k+a)^s, real z in (-1, 1) or z = 1,
 * real s; a <= 0 in the symmetric convention. */
bool sf_lerch_r(double z, double s, double a, double* out);
/* LerchPhi(z, s, a) = Sum_{k>=0} z^k / (k+a)^s, complex z (|z| < 1, or z = 1),
 * complex s, complex a with Re a > 0. A wholly real call delegates to the real
 * kernel sf_machine_lerchphi. */
bool sf_lerch_c(double zr, double zi, double sr, double si, double ar, double ai,
                double* vr, double* vi);

/* Dilogarithm Li_2(x), real x <= 1 (x > 1 is complex: use sf_dilog_c). */
bool sf_dilog_r(double x, double* out);
/* Dilogarithm Li_2(z), complex z; a real z > 1 takes the value continuous from
 * below the cut [1, Infinity), as Mathematica does (Im Li_2(x) = -Pi Log x). */
bool sf_dilog_c(double zr, double zi, double* vr, double* vi);

/* PolyLog Li_s(x), real order s, real x in [-1, 1] (and every real x <= 1 for
 * s = 2). */
bool sf_polylog_r(double s, double x, double* out);
/* PolyLog Li_s(z), real order s, complex z. Same cut convention as sf_dilog_c. */
bool sf_polylog_c(double sr, double si, double zr, double zi,
                  double* vr, double* vi);

#endif /* MATHILDA_SF_MACHINE_ZETA_H */
