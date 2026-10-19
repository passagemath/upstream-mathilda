/*
 * facint.c — PrimeQ, NextPrime, EulerPhi, and the facint_factor_complete
 * API.  Integer factorisation itself (FactorInteger and every
 * FactorInteger`Method) lives in src/numbertheory/factorinteger/.
 */
#include "facint.h"
#include "fi_internal.h"
#include "arithmetic.h"
#include "expr.h"
#include "eval.h"
#include "symtab.h"
#include "attr.h"
#include "sym_names.h"
#include "gmp_compat.h"   /* mpz_prevprime fallback for GMP < 6.3.0 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
/* See facint.h. The 64-bit arm skips GMP entirely: no limb allocation, one
 * Montgomery setup and at most seven witnesses, versus mpz_probab_prime_p's
 * trial division + BPSW on an mpz. Both are exact below 2^64 (BPSW has no
 * counterexample there; the witness sets are proven), so PrimeQ's answers are
 * unchanged -- only its cost on machine integers drops. */
int mth_is_prime_u64(uint64_t n) {
    return fi_u64_is_prime(n);
}

int mth_is_prime_mpz(const mpz_t n) {
    if (mpz_sgn(n) <= 0) return 0;
    if (mpz_sizeinbase(n, 2) <= 64) {
        /* mpz_export, not mpz_get_ui: unsigned long is 32 bits on LLP64. */
        uint64_t v = 0;
        mpz_export(&v, NULL, -1, sizeof v, 0, 0, n);
        return fi_u64_is_prime(v);
    }
    return mpz_probab_prime_p(n, 25) != 0;
}

/* Gaussian-integer primality test for a complex number a + b*I with
 * integer parts. Sets *out_is_gprime to 1/0. Inputs `a` and `b` may carry
 * any sign; they are taken in absolute value internally. */
static void gaussian_prime_test(const mpz_t a_in, const mpz_t b_in, int* out_is_gprime) {
    mpz_t a, b;
    mpz_init(a);
    mpz_init(b);
    mpz_abs(a, a_in);
    mpz_abs(b, b_in);
    int is_gprime = 0;
    if (mpz_sgn(a) == 0 && mpz_sgn(b) == 0) {
        is_gprime = 0;  /* 0 + 0i is not prime */
    } else if (mpz_sgn(b) == 0) {
        /* Pure real: Gaussian prime iff |a| is prime and |a| ≡ 3 (mod 4). */
        if (mpz_probab_prime_p(a, 25) && mpz_fdiv_ui(a, 4) == 3) {
            is_gprime = 1;
        }
    } else if (mpz_sgn(a) == 0) {
        /* Pure imaginary: Gaussian prime iff |b| is prime and |b| ≡ 3 (mod 4). */
        if (mpz_probab_prime_p(b, 25) && mpz_fdiv_ui(b, 4) == 3) {
            is_gprime = 1;
        }
    } else {
        /* Both nonzero: Gaussian prime iff a^2 + b^2 is prime. */
        mpz_t norm;
        mpz_init(norm);
        mpz_mul(norm, a, a);
        mpz_addmul(norm, b, b);
        is_gprime = mpz_probab_prime_p(norm, 25);
        mpz_clear(norm);
    }
    mpz_clear(a);
    mpz_clear(b);
    *out_is_gprime = is_gprime;
}

/* Parse PrimeQ option arguments after the first. Returns true on success,
 * filling *out_gaussian. Returns false if any option is unrecognized or
 * malformed (caller should return False from the predicate in that case
 * to honour the *Q always-boolean contract). */
static bool primeq_parse_options(const Expr* res, size_t first_opt,
                                 bool* out_gaussian) {
    *out_gaussian = false;
    size_t argc = res->data.function.arg_count;
    for (size_t i = first_opt; i < argc; i++) {
        Expr* a = res->data.function.args[i];
        if (a->type != EXPR_FUNCTION ||
            a->data.function.head->type != EXPR_SYMBOL ||
            a->data.function.head->data.symbol.name != SYM_Rule ||
            a->data.function.arg_count != 2) {
            return false;
        }
        Expr* key = a->data.function.args[0];
        Expr* val = a->data.function.args[1];
        if (key->type != EXPR_SYMBOL) return false;
        if (key->data.symbol.name == SYM_GaussianIntegers) {
            if (val->type != EXPR_SYMBOL) return false;
            if (strcmp(val->data.symbol.name, "True") == 0) *out_gaussian = true;
            else if (strcmp(val->data.symbol.name, "False") == 0) *out_gaussian = false;
            else return false;
        } else {
            return false;  /* unknown option */
        }
    }
    return true;
}

Expr* builtin_primeq(Expr* res) {
    /* *Q predicates must always evaluate to True or False — never remain
     * unevaluated. Zero-arg or non-function res (the latter unreachable
     * via the evaluator, kept defensive) is therefore not prime. */
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count < 1) {
        return expr_new_symbol(SYM_False);
    }
    Expr* arg = res->data.function.args[0];

    /* Parse optional GaussianIntegers->True/False. A malformed option list
     * makes the call evaluate to False (per *Q always-boolean contract). */
    bool gaussian = false;
    if (res->data.function.arg_count > 1) {
        if (!primeq_parse_options(res, 1, &gaussian)) {
            return expr_new_symbol(SYM_False);
        }
    }

    if (arg->type == EXPR_INTEGER || arg->type == EXPR_BIGINT) {
        mpz_t n;
        expr_to_mpz(arg, n);
        int is_prime = 0;
        if (gaussian) {
            /* A rational integer n is a Gaussian prime iff |n| is a prime
             * in Z AND |n| ≡ 3 (mod 4). (n=2 = -i(1+i)^2, and primes ≡ 1
             * mod 4 split as (a+bi)(a-bi).) */
            mpz_t zero;
            mpz_init_set_ui(zero, 0);
            gaussian_prime_test(n, zero, &is_prime);
            mpz_clear(zero);
        } else {
            mpz_abs(n, n);
            if (mpz_sgn(n) > 0) {
                is_prime = mth_is_prime_mpz(n);
            }
        }
        mpz_clear(n);
        return expr_new_symbol(is_prime ? "True" : "False");
    }

    if (arg->type == EXPR_FUNCTION && arg->data.function.head->type == EXPR_SYMBOL) {
        const char* head = arg->data.function.head->data.symbol.name;
        /* Gaussian prime test for Complex[a, b] with integer parts. Works
         * in both option-on and option-off modes since membership in Z[i]
         * is the relevant universe in either case. */
        if (head == SYM_Complex && arg->data.function.arg_count == 2) {
            Expr* re_expr = arg->data.function.args[0];
            Expr* im_expr = arg->data.function.args[1];
            bool re_int = (re_expr->type == EXPR_INTEGER || re_expr->type == EXPR_BIGINT);
            bool im_int = (im_expr->type == EXPR_INTEGER || im_expr->type == EXPR_BIGINT);
            if (!re_int || !im_int) {
                return expr_new_symbol(SYM_False);
            }
            mpz_t a, b;
            expr_to_mpz(re_expr, a);
            expr_to_mpz(im_expr, b);
            int is_gprime = 0;
            gaussian_prime_test(a, b, &is_gprime);
            mpz_clear(a);
            mpz_clear(b);
            return expr_new_symbol(is_gprime ? "True" : "False");
        }
    }

    /* Everything else (symbols, reals, rationals, strings, symbolic
     * functions, etc.) is not a (Gaussian) prime. *Q predicates must
     * always return True or False — never remain unevaluated. */
    return expr_new_symbol(SYM_False);
}

Expr* builtin_nextprime(Expr* res) {
    if (res->type != EXPR_FUNCTION || (res->data.function.arg_count != 1 && res->data.function.arg_count != 2)) return NULL;

    Expr* x_expr = res->data.function.args[0];
    int64_t k = 1;
    if (res->data.function.arg_count == 2) {
        if (res->data.function.args[1]->type != EXPR_INTEGER) return NULL;
        k = res->data.function.args[1]->data.integer;
    }

    if (k == 0) return expr_copy(x_expr);

    mpz_t current;
    if (x_expr->type == EXPR_INTEGER || x_expr->type == EXPR_BIGINT) {
        expr_to_mpz(x_expr, current);
    } else if (x_expr->type == EXPR_REAL) {
        mpz_init_set_d(current, floor(x_expr->data.real));
    } else {
        int64_t n, d;
        if (is_rational(x_expr, &n, &d)) {
            int64_t start_val = n / d;
            if (n < 0 && n % d != 0) start_val--; 
            mpz_init_set_si(current, start_val);
        } else {
            return NULL;
        }
    }

    if (k > 0) {
        for (int64_t i = 0; i < k; i++) {
            mpz_nextprime(current, current);
        }
    } else {
        int64_t count = -k;
        for (int64_t i = 0; i < count; i++) {
            if (mpz_cmp_ui(current, 2) <= 0) {
                mpz_clear(current);
                return NULL;
            }
            int found = mathilda_mpz_prevprime(current, current);
            if (!found) {
                mpz_clear(current);
                return NULL;
            }
        }
    }

    Expr* result = expr_bigint_normalize(expr_new_bigint_from_mpz(current));
    mpz_clear(current);
    return result;
}

/* --------------------------------------------------------------------- *
 *  Public complete-factorisation helper (number-field maximal-order gate *
 *  src/numbertheory/numberfield.c, Dirichlet characters).                *
 *                                                                        *
 *  Fully factor |n| with the Automatic method.  The distinct primes are  *
 *  written in ascending order to primes[] (each mpz_init'd on success -- *
 *  caller must mpz_clear) and their exponents to exps[], with *count set *
 *  to the number of distinct primes.                                     *
 *                                                                        *
 *  Returns true iff the factorisation is COMPLETE (no composite survived *
 *  the search) and the distinct-prime count fits in `cap`.  Returns      *
 *  false -- writing nothing, *count = 0 -- otherwise; the caller must    *
 *  then DECLINE rather than trust a possibly-composite "prime".          *
 *  n == 0 returns false; |n| == 1 returns true with *count = 0.          */
bool facint_factor_complete(const mpz_t n, mpz_t* primes, int64_t* exps,
                            int cap, int* count) {
    *count = 0;
    if (mpz_sgn(n) == 0) return false;
    if (mpz_cmpabs_ui(n, 1) == 0) return true;
    FiFactors f;
    fi_factors_init(&f);
    bool complete = fi_factor_complete_mpz(n, &f) && f.n <= (size_t)cap;
    if (complete) {
        for (size_t i = 0; i < f.n; i++) {
            mpz_init_set(primes[i], f.v[i].p);
            exps[i] = f.v[i].e;
        }
        *count = (int)f.n;
    }
    fi_factors_clear(&f);
    return complete;
}

Expr* builtin_eulerphi(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) return NULL;
    Expr* arg = res->data.function.args[0];
    if (!expr_is_integer_like(arg)) return NULL;

    mpz_t n;
    expr_to_mpz(arg, n);
    if (mpz_sgn(n) == 0) { mpz_clear(n); return expr_new_integer(0); }
    mpz_abs(n, n);                                          /* phi(-n) = phi(n) */
    if (mpz_cmp_ui(n, 1) == 0) { mpz_clear(n); return expr_new_integer(1); }

    FiFactors f;
    fi_factors_init(&f);
    if (!fi_factor_complete_mpz(n, &f)) {
        /* A composite survived the search: phi would be wrong, decline. */
        fi_factors_clear(&f);
        mpz_clear(n);
        return NULL;
    }

    /* phi(n) = prod p^(e-1) (p - 1), applied as (phi / p) * (p - 1). */
    mpz_t phi, pm1;
    mpz_init_set(phi, n);
    mpz_init(pm1);
    for (size_t i = 0; i < f.n; i++) {
        mpz_divexact(phi, phi, f.v[i].p);
        mpz_sub_ui(pm1, f.v[i].p, 1);
        mpz_mul(phi, phi, pm1);
    }
    fi_factors_clear(&f);
    mpz_clears(pm1, n, NULL);
    Expr* out = expr_bigint_normalize(expr_new_bigint_from_mpz(phi));
    mpz_clear(phi);
    return out;
}

void fi_init(void);

void facint_init(void) {
    symtab_add_builtin("PrimeQ", builtin_primeq);
    symtab_add_builtin("NextPrime", builtin_nextprime);
    symtab_add_builtin("EulerPhi", builtin_eulerphi);

    symtab_get_def("PrimeQ")->attributes |= (ATTR_PROTECTED | ATTR_LISTABLE);
    /* PrimePi is registered in numbertheory_init() (src/numbertheory/prime.c). */
    symtab_get_def("EulerPhi")->attributes |= (ATTR_PROTECTED | ATTR_LISTABLE);
    symtab_get_def("NextPrime")->attributes |= (ATTR_PROTECTED | ATTR_LISTABLE);

    /* FactorInteger and FactorInteger`<Method> (src/numbertheory/factorinteger/). */
    fi_init();
}
