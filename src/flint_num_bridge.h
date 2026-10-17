/*
 * flint_num_bridge.h
 * ------------------
 * Boundary between Mathilda's numeric Expr scalars and FLINT's arbitrary-
 * precision complex ball arithmetic (arb/acb), for the transcendental special
 * functions. FLINT's acb_hypgeom / acb_dirichlet kernels are rigorous (correct
 * rounding, full complex plane, standard branch cuts), so they back the FLINT`
 * context wrappers exposed at the REPL and, over time, transparent numeric
 * fast paths for the hand-rolled MPFR kernels.
 *
 * Each entry point takes numeric Expr arguments (Integer / BigInt / Rational /
 * Real / MPFR / Complex of those) and returns a freshly-owned numeric Expr
 * (EXPR_MPFR, or Complex[...] when the imaginary part is non-zero), using the
 * working precision implied by the arguments (machine 53-bit floor). It returns
 * NULL — leaving the caller's FLINT`f[...] unevaluated — when an argument is
 * non-numeric / symbolic, or the result is a pole / non-finite, or FLINT is
 * absent (USE_FLINT undefined). No argument is mutated.
 */
#ifndef FLINT_NUM_BRIDGE_H
#define FLINT_NUM_BRIDGE_H

#include "expr.h"

#ifdef __cplusplus
extern "C" {
#endif

Expr* flint_num_zeta(const Expr* s);                       /* Zeta[s]            */
Expr* flint_num_hurwitz_zeta(const Expr* s, const Expr* a);/* HurwitzZeta[s, a]  */
Expr* flint_num_polygamma(const Expr* n, const Expr* z);   /* PolyGamma[n, z]    */
Expr* flint_num_polylog(const Expr* s, const Expr* z);     /* PolyLog[s, z]      */
Expr* flint_num_airy_ai(const Expr* z);                    /* AiryAi[z]          */
Expr* flint_num_airy_ai_prime(const Expr* z);              /* AiryAiPrime[z]     */
Expr* flint_num_airy_bi(const Expr* z);                    /* AiryBi[z]          */
Expr* flint_num_airy_bi_prime(const Expr* z);              /* AiryBiPrime[z]     */
Expr* flint_num_lerchphi(const Expr* z, const Expr* s, const Expr* a); /* LerchPhi[z,s,a] */
Expr* flint_num_stieltjes(const Expr* n, const Expr* a);   /* StieltjesGamma[n,a]*/

/* The Legendre elliptic integrals, in the PARAMETER convention m = k^2 that
 * Mathematica uses (not the modulus k) -- which is also Arb's, so no argument
 * translation happens here:
 *
 *   EllipticK[m]           = F(Pi/2 | m)
 *   EllipticE[m]           = E(Pi/2 | m)
 *   EllipticF[phi, m]      = Int_0^phi dt / Sqrt(1 - m Sin[t]^2)
 *   EllipticE[phi, m]      = Int_0^phi Sqrt(1 - m Sin[t]^2) dt
 *   EllipticPi[n, m]       = Pi(n; Pi/2 | m)
 *   EllipticPi[n, phi, m]  = Int_0^phi dt / ((1 - n Sin[t]^2) Sqrt(1 - m Sin[t]^2))
 *
 * Arb defines the incomplete forms on -Pi/2 <= Re phi <= Pi/2 and extends them
 * by the quasi-period (F(phi + k Pi | m) = F(phi | m) + 2 k K(m), and likewise
 * for E and Pi), so a phi of any size or complexity is handled here rather than
 * by the caller. EllipticPi with n > 1 crosses the pole at Sin[t]^2 = 1/n and is
 * the Cauchy principal value, which acb_elliptic_pi_inc computes -- a complex
 * value on the principal branch, matching Mathematica. */
Expr* flint_num_elliptic_k(const Expr* m);                              /* EllipticK[m]        */
Expr* flint_num_elliptic_e(const Expr* m);                              /* EllipticE[m]        */
Expr* flint_num_elliptic_f(const Expr* phi, const Expr* m);             /* EllipticF[phi, m]   */
Expr* flint_num_elliptic_e_inc(const Expr* phi, const Expr* m);         /* EllipticE[phi, m]   */
Expr* flint_num_elliptic_pi(const Expr* n, const Expr* m);              /* EllipticPi[n, m]    */
Expr* flint_num_elliptic_pi_inc(const Expr* n, const Expr* phi,
                                const Expr* m);                          /* EllipticPi[n,phi,m] */

/* HypergeometricU[a, b, z] -- the Tricomi confluent hypergeometric U, via FLINT's
 * rigorous acb_hypgeom_u (correct rounding, full complex plane, standard branch
 * cut on z from -infinity to 0, integer b handled natively). Used for the MPFR /
 * complex / z<=0 numerics of HypergeometricU and (through the confluent template)
 * WhittakerW. */
Expr* flint_num_hypergeometric_u(const Expr* a, const Expr* b, const Expr* z);

/* StruveH[nu, z] -- the Struve function H_nu(z). FLINT/arb ships no Struve
 * kernel, so this is composed from the rigorous regularized generalized
 * hypergeometric: H_nu(z) = (z/2)^(nu+1) * 1F2~(1; 3/2, nu+3/2; -z^2/4). The
 * regularized form folds the 1/Gamma factors into the sum, so it is entire in
 * nu (no poles) and valid over the whole complex plane. A cancellation
 * allowance (~1.44*|z| bits) is pre-added to the working precision so a large
 * |z| gets its bits up front rather than by blind doubling. */
Expr* flint_num_struveh(const Expr* nu, const Expr* z);     /* StruveH[nu, z]    */

/* Hypergeometric2F1[a, b, c, z] at arbitrary precision (acb_hypgeom_2f1, with
 * its transformations near z = 1). Used by HypergeometricPFQ's MPFR path. */
Expr* flint_num_hyp2f1(const Expr* a, const Expr* b, const Expr* c, const Expr* z);

/* MeijerG[{{a_1..a_n},{a_n+1..a_p}}, {{b_1..b_m},{b_m+1..b_q}}, z, r] via
 * Slater's theorem over arb's regularized pFq (the dual expansion in 1/w for
 * p > q or p == q, |w| > 1; w = z^(1/r)), with the log case (integer-differing
 * b's) handled by perturbation. `a` holds a_1..a_p, `b` b_1..b_q; `r` may be
 * NULL (r = 1). The output precision is the minimum inexact precision of the
 * arguments (machine floor). NULL for non-numeric input, z = 0, a value that
 * does not exist, or p == q >= 3 with |w| = 1 (no continuation implemented). */
Expr* flint_num_meijerg(const Expr* const* a, size_t n, size_t p,
                        const Expr* const* b, size_t m, size_t q,
                        const Expr* z, const Expr* r);

/* The same engine on machine inputs, rounded to double: the fallback behind
 * src/special_functions/sf_machine_meijerg.c when the machine Slater sum
 * declines (log case, cancellation), so every machine surface still answers.
 * ai / bi may be NULL. *is_c is false when the imaginary part is exactly zero.
 * Returns false without FLINT or when the value does not exist. */
bool flint_mg_eval_d(const double* ar, const double* ai, size_t n, size_t p,
                     const double* br, const double* bi, size_t m, size_t q,
                     double r, double zr, double zi,
                     double* rr, double* ri, bool* is_c);

/* Laurent expansion about eps = 0 of
 *
 *   F(eps) = prod_{i < nf} Gamma(x_i + s_i eps)^{e_i},   s_i, e_i in {+1, -1},
 *
 * as F = eps^v (f_0 + f_1 eps + ... + f_{T-1} eps^{T-1} + O(eps^T)). x_i is
 * read EXACTLY from x[4i .. 4i+3] = (re.hi, re.lo, im.hi, im.lo); an x_i with
 * lo parts and imaginary part zero and an integer hi <= 0 is a pole of Gamma
 * (a zero of 1/Gamma) and moves v by -e_i. The coefficients come back as
 * double-doubles, out[4k .. 4k+3] = (re.hi, re.lo, im.hi, im.lo) of f_k,
 * computed in arb to >= 110 certified bits. This is the z-independent start of
 * sf_machine_meijerg.c's log-case residue table. False without FLINT, or when
 * the accuracy is not reached (the caller then has no table and declines). */
bool flint_mg_laurent(size_t nf, const double* x, const int* s, const int* e,
                      size_t T, long* v, double* out);

/* Registers the FLINT` context numeric builtins (Zeta / HurwitzZeta /
 * PolyGamma / StieltjesGamma). Called from core_init(). No-op without FLINT. */
void flint_num_bridge_init(void);

#ifdef __cplusplus
}
#endif

#endif /* FLINT_NUM_BRIDGE_H */
