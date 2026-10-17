/* Machine-precision (raw double re/im) kernels for the Gamma family:
 * Gamma, LogGamma, PolyGamma, Beta and Pochhammer.
 *
 * These are the ONE implementation behind every machine surface of those heads:
 * the scalar builtin (via the sf_machine_scalar_args gate), the NDArray /
 * packed element-wise kernels registered in ndkernels.c, and the Compile[] VM,
 * which reaches the same registered kernels. Sharing the helper is what makes
 * scalar == List == packed == NDArray == Compile bit-for-bit (CLAUDE.md).
 *
 * Contract (the engine's usual one): every entry point returns FALSE when it
 * cannot produce a finite machine answer it trusts -- a pole, an argument outside
 * the regime the algorithm covers, overflow -- and the caller then falls through
 * to its MPFR / symbolic path unchanged. A kernel never invents a value.
 *
 * Complex arithmetic is spelled out on (re, im) pairs: no <complex.h>, whose
 * libgcc __divdc3 / cpow are slow and whose availability varies (CLAUDE.md §10).
 */
#ifndef MATHILDA_SF_MACHINE_GAMMA_H
#define MATHILDA_SF_MACHINE_GAMMA_H

#include <stdbool.h>

/* Digamma psi(x) for real x (poles at x = 0, -1, -2, ... decline). */
bool sfg_digamma(double x, double* out);

/* PolyGamma psi^(n)(x) for integer 0 <= n <= SFG_POLYGAMMA_MAX_N and real x.
 * n >= 1 requires x > 0 (the negative axis is left to the MPFR path). */
#define SFG_POLYGAMMA_MAX_N 40
bool sfg_polygamma(int n, double x, double* out);

/* PolyGamma psi^(n)(z) for a complex z = zr + i zi. n = 0 covers the whole plane
 * bar the poles (reflection for Re z < 0); n >= 1 requires Re z >= 0. */
bool sfg_polygamma_c(int n, double zr, double zi, double* rr, double* ri);

/* Gamma(z) and the principal (continuous) LogGamma(z) for complex z. LogGamma is
 * analytic off the negative real axis and agrees with Log[Gamma[x]] for x > 0,
 * Mathematica's convention (Hare 1997); it is NOT the principal log of Gamma. */
bool sfg_gamma_c(double zr, double zi, double* rr, double* ri);
bool sfg_loggamma_c(double zr, double zi, double* rr, double* ri);

/* Beta(a, b) = Gamma(a) Gamma(b) / Gamma(a + b), real and complex. */
bool sfg_beta(double a, double b, double* out);
bool sfg_beta_c(double ar, double ai, double br, double bi, double* rr, double* ri);

/* Pochhammer (a)_n = Gamma(a + n) / Gamma(a), real and complex. An integer n is a
 * direct product (n > 0) or the reciprocal product (n < 0). */
bool sfg_pochhammer(double a, double n, double* out);
bool sfg_pochhammer_c(double ar, double ai, double nr, double ni, double* rr, double* ri);

#endif /* MATHILDA_SF_MACHINE_GAMMA_H */
