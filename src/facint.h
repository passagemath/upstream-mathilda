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

/* Complete prime factorisation of |n| for callers that need an exact,
 * proven-complete factor list (the Thue solver's maximal-order gate).
 * Fills ascending primes[]/exps[] (primes mpz_init'd on success; caller
 * mpz_clears), sets *count.  Returns true iff COMPLETE (every base a
 * strong probable prime and count <= cap); false => caller must DECLINE. */
bool facint_factor_complete(const mpz_t n, mpz_t* primes, int64_t* exps,
                            int cap, int* count);

void facint_init(void);

#endif
