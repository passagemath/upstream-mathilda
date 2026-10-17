/* sf_machine_dirichlet.h -- machine-precision Dirichlet L-functions.
 *
 * The raw-double kernel behind every machine surface of DirichletL[k, j, s]:
 * the scalar builtin, the packed / NDArray buffer loop, and Compile[] (through
 * the n-ary ND kernel ndk_DirichletL_n). One helper, so the surfaces agree
 * bit-for-bit.
 *
 * Method, for the primitive character chi* mod f inducing chi mod k:
 *   L(s, chi) = L(s, chi*) prod_{p | k, p !| f} (1 - chi*(p) p^-s).
 *   f = 1     L(s, chi*) = Zeta[s] (the shared sf_zeta_* kernels).
 *   Re s >= 1/2
 *             Euler-Maclaurin on L itself: the direct sum over n <= M f plus,
 *             per residue class a, the E-M tail of f^-s Zeta[s, M + a/f]. The
 *             integral terms are summed as
 *               -M^(1-s) sum_a chi(a) log1p(a/(fM)) phi((1-s) log1p(a/(fM))),
 *             phi(u) = expm1(u)/u, which is exact because sum_a chi(a) = 0: no
 *             cancellation as s -> 1, and s = 1 itself is just phi(0) = 1.
 *   Re s < 1/2
 *             the functional equation of the primitive character,
 *               L(s, chi*) = eps (f/Pi)^(1/2-s) Gamma((1-s+d)/2)/Gamma((s+d)/2)
 *                            L(1-s, conj chi*),
 *             with the root number eps from the character table; exact 0 at
 *             the trivial zeros s = -d - 2m.
 * Declines (returns false; the caller falls through to arbitrary precision):
 * a pole (principal character at s = 1), overflow, a conductor beyond the
 * table limit, or a direct sum that would exceed ~2*10^7 terms.
 */
#ifndef MATHILDA_SF_MACHINE_DIRICHLET_H
#define MATHILDA_SF_MACHINE_DIRICHLET_H

#include <stdbool.h>
#include "dirichlet.h"

/* L(s, chi) for the table T at complex s. */
bool sf_dirichlet_l_tab_c(const DchTable* T, double sr, double si, double* vr, double* vi);
/* L(s, chi) at real s for a REAL character (false otherwise). */
bool sf_dirichlet_l_tab_r(const DchTable* T, double s, double* v);

/* The same with (k, j) given as doubles, as the ND / Compile surfaces pass
 * them: false unless k, j are integers in range. */
bool sf_dirichlet_l_c(double k, double j, double sr, double si, double* vr, double* vi);
bool sf_dirichlet_l_r(double k, double j, double s, double* v);

#endif /* MATHILDA_SF_MACHINE_DIRICHLET_H */
