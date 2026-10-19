#ifndef FACINT_H
#define FACINT_H

#include "expr.h"
#include <gmp.h>
#include <stdbool.h>
#include <stdint.h>

Expr* builtin_primeq(Expr* res);
Expr* builtin_nextprime(Expr* res);
Expr* builtin_factorinteger(Expr* res);
Expr* builtin_eulerphi(Expr* res);
Expr* builtin_primepowerq(Expr* res);

/* THE PrimeQ primality test, shared by every caller that must agree with
 * PrimeQ (RandomPrime draws its candidates through it).
 *
 *   mth_is_prime_u64  deterministic Miller-Rabin in Montgomery form
 *                     (fi_u64_is_prime): exact for every n < 2^64.
 *   mth_is_prime_mpz  dispatches to mth_is_prime_u64 when 0 <= n < 2^64,
 *                     else mpz_probab_prime_p(n, 25) (BPSW + 25 MR rounds).
 *
 * Both return 1 for prime, 0 otherwise; n < 2 (including negatives) is 0 --
 * callers wanting PrimeQ's |n| semantics take the absolute value first. */
int mth_is_prime_u64(uint64_t n);
int mth_is_prime_mpz(const mpz_t n);

/* THE PrimePowerQ test: is n = p^k for a prime p and k >= 1?  On true,
 * writes the prime p to *p_out / p_out (either may be NULL).  n < 2 is
 * false -- callers wanting |n| semantics take the absolute value first.
 *
 *   mth_prime_power_base_u64  trial division by the odd primes < 64, then
 *                             one mth_is_prime_u64, then exact k-th roots
 *                             for k in {2, 3, 5, 7}; no GMP.
 *   mth_prime_power_base_mpz  < 2^64 dispatches to the u64 arm; else trial
 *                             division by the primes < 2^10 (a hit is one
 *                             mpz_remove), reduction to the primitive root
 *                             (FLINT fmpz_is_perfect_power, or GMP
 *                             mpz_perfect_power_p + mpz_root without
 *                             FLINT), then one mth_is_prime_mpz on it. */
bool mth_prime_power_base_u64(uint64_t n, uint64_t* p_out);
bool mth_prime_power_base_mpz(const mpz_t n, mpz_t p_out);

/* Complete prime factorisation of |n| for callers that need an exact,
 * proven-complete factor list (the Thue solver's maximal-order gate).
 * Fills ascending primes[]/exps[] (primes mpz_init'd on success; caller
 * mpz_clears), sets *count.  Returns true iff COMPLETE (every base a
 * strong probable prime and count <= cap); false => caller must DECLINE. */
bool facint_factor_complete(const mpz_t n, mpz_t* primes, int64_t* exps,
                            int cap, int* count);

void facint_init(void);

#endif
