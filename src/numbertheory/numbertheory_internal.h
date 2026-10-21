#ifndef NUMBERTHEORY_INTERNAL_H
#define NUMBERTHEORY_INTERNAL_H

/* Cross-cutting helpers shared between the per-builtin translation units of
 * src/numbertheory/.  Everything declared here is internal to the
 * number-theory subsystem and must NOT leak into any other compilation unit
 * (mirrors the eigen_internal.h / ludecomp_internal.h pattern in
 * src/linalg/).  Each helper is defined in exactly one .c file as noted. */

#include "numbertheory.h"   /* Expr */
#include <stdbool.h>
#include <stdint.h>
#include <stddef.h>
#include <gmp.h>

/* Maximum distinct prime factors collected by pr_collect_distinct_primes;
 * the caller-supplied mpz_t arrays must have this capacity.  Shared by
 * primitiveroot.c and multiplicativeorder.c. */
#define PR_MAX_DISTINCT_PRIMES 256

/* --- nt_util.c : generic GCD/LCM/rational helpers ---------------------- */
Expr*   single_arg_abs_or_copy(Expr* arg);
void    rational_like_to_mpz_pair(const Expr* e, mpz_t num, mpz_t den);
bool    rational_like_needs_bigint(const Expr* e);
Expr*   mpz_pair_to_rational_expr(const mpz_t num_in, const mpz_t den_in);
int64_t lcm_checked(int64_t a, int64_t b, bool *overflow);
/* Shared diagnostics (nt_util.c), routed through the Quiet/Check funnel and
 * gated on the arithmetic probe mute. Both return NULL. */
Expr*   nt_first_inexact_arg(Expr* res);
Expr*   nt_emit_exact(const char* head, Expr* arg, Expr* res);
Expr*   nt_emit_int(const char* head, size_t pos, Expr* res);

/* --- primitiveroot.c (shared with multiplicativeorder.c) -------------- */
/* Gather the distinct prime factors of m_in into out[0..*count); out must
 * have PR_MAX_DISTINCT_PRIMES capacity and *count must be 0 on entry. */
bool    pr_collect_distinct_primes(const mpz_t m_in, mpz_t* out, size_t* count);

/* --- nt_dlog.c : discrete logarithm (MultiplicativeOrder's 3-arg form) - */
/* Smallest m >= 1 with k^m == targets[i] (mod n) over all i, written to
 * m_out.  Pre: n >= 1, k a unit mod n of order d = prod q[j]^e[j] (q prime);
 * cyclic says whether (Z/nZ)^* is cyclic.  Pohlig-Hellman with BSGS / Pollard
 * rho per prime; every answer is verified.  Returns false when no target is a
 * power of k, or on allocation failure / an expired TimeConstrained. */
bool    nt_dlog_min(const mpz_t k, const mpz_t n, const mpz_t d,
                    const mpz_t* q, const unsigned long* e, size_t nq,
                    int cyclic, const mpz_t* targets, size_t nt, mpz_t m_out);

/* --- nt_gaussian.c : Gaussian-integer & divisor enumeration ----------- */
bool    is_gaussian_integer(Expr* e);
bool    df_factor_mpz(const mpz_t m, mpz_t** out_primes,
                      unsigned long** out_exps, size_t* out_n);
bool    df_factor_list_read(const Expr* fact, mpz_t** out_primes,
                            unsigned long** out_exps, size_t* out_n);
bool    df_to_gaussian(Expr* e, mpz_t a, mpz_t b);
bool    df_gaussian_prime_factor(const mpz_t a_in, const mpz_t b_in,
                                 mpz_t** out_gu, mpz_t** out_gv,
                                 unsigned long** out_ge, size_t* out_gn);
void    df_normalize_quadrant(mpz_t re, mpz_t im);
Expr*   divisors_ordinary(const mpz_t nabs);
Expr*   divisors_gaussian(const mpz_t a_in, const mpz_t b_in);

/* Gaussian GCD arithmetic (nt_gaussian.c). Every output is init'd by the
 * caller. gcd / lcm / the xgcd gcd are first-quadrant normalised (Re > 0,
 * Im >= 0); gi_xgcd gives g = s a + t b with the canonical coefficients
 * (nearest-rounding Euclid, ties to even). gi_normalize reports the power of
 * i it multiplied by; gi_rotate_by multiplies by i^k. */
void    gi_xgcd(mpz_t gr, mpz_t gi, mpz_t sr, mpz_t si, mpz_t tr, mpz_t ti,
                const mpz_t ar, const mpz_t ai, const mpz_t br, const mpz_t bi);
void    gi_gcd(mpz_t gr, mpz_t gi, const mpz_t ar, const mpz_t ai,
               const mpz_t br, const mpz_t bi);
void    gi_lcm(mpz_t lr, mpz_t li, const mpz_t ar, const mpz_t ai,
               const mpz_t br, const mpz_t bi);
void    gi_normalize(mpz_t re, mpz_t im, int* unit_power);
void    gi_rotate_by(mpz_t re, mpz_t im, int unit_power);
/* x mod m = x - m Round[x/m] (ties to even), and a^e mod m for integer e
 * (false when e < 0 and a is not invertible). Outputs init'd by the caller. */
void    gi_mod(mpz_t rr, mpz_t ri, const mpz_t xr, const mpz_t xi,
               const mpz_t mr, const mpz_t mi);
bool    gi64_powmod(Expr* a, int64_t e, Expr* m, Expr** out);
bool    gi_powmod(mpz_t rr, mpz_t ri, const mpz_t ar, const mpz_t ai,
                  const mpz_t e, const mpz_t mr, const mpz_t mi);
/* Gaussian rationals (re + im i) / den with den > 0 (outputs init'd). */
bool    gq_from_expr(Expr* e, mpz_t re, mpz_t im, mpz_t den);
Expr*   gq_to_expr(const mpz_t re, const mpz_t im, const mpz_t den);
bool    is_gaussian_rational(Expr* e);

#endif
