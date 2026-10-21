/* multiplicativeorder.c -- MultiplicativeOrder[].
 * Split from numbertheory.c; see numbertheory.h and
 * numbertheory_internal.h for the subsystem layout. */

#include "numbertheory.h"
#include "numbertheory_internal.h"
#include "arithmetic.h"
#include "eval.h"
#include "sym_names.h"
#include "internal.h"
#include "print.h"
#include "symtab.h"
#include "attr.h"
#include "ndinteger.h" /* int64 phi / factor helpers */
#include "message.h"   /* mth_message: Quiet/Check funnel */
#include "facint.h"    /* facint_factor_complete, mth_is_prime_mpz */
#include "core.h"      /* tc_check_deadline */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <gmp.h>

/* =====================================================================
 * MultiplicativeOrder
 *
 *   MultiplicativeOrder[k, n]
 *       smallest positive m with k^m == 1 (mod n).
 *   MultiplicativeOrder[k, n, {r1, r2, ...}]
 *       smallest positive m with k^m == r_i (mod n) for some i.
 *
 * Returns unevaluated when:
 *   - n is 0 (no group);
 *   - gcd(k, n) != 1 (k is not a unit mod n, no finite order exists);
 *   - in the 3-arg form, no power of k lands in the residue set.
 * The 3-arg form is a discrete logarithm, solved by Pohlig-Hellman with
 * baby-step giant-step / Pollard rho per prime of ord(k) (nt_dlog.c) -- no
 * search cap, so any order whose primes are tractable answers.
 * All arithmetic is on GMP mpz_t so bignum k and n share the same code
 * path with machine integers.  Negative n is treated as |n|; negative or
 * out-of-range k is reduced modulo n.
 * ===================================================================*/

/* The prime factorisation of a group order: primes q[0..n) with exponents
 * e[].  Entries are mpz_init'd as they are added; mo_factors_clear frees. */
typedef struct {
    mpz_t         q[PR_MAX_DISTINCT_PRIMES];
    unsigned long e[PR_MAX_DISTINCT_PRIMES];
    size_t        n;
} MoFactors;

static void mo_factors_clear(MoFactors* f) {
    for (size_t i = 0; i < f->n; i++) mpz_clear(f->q[i]);
    f->n = 0;
}

/* Multiply the recorded factorisation by p^e, merging a repeated prime. */
static bool mo_factors_add(MoFactors* f, const mpz_t p, unsigned long e) {
    if (e == 0) return true;
    for (size_t i = 0; i < f->n; i++)
        if (mpz_cmp(f->q[i], p) == 0) { f->e[i] += e; return true; }
    if (f->n >= PR_MAX_DISTINCT_PRIMES) return false;
    mpz_init_set(f->q[f->n], p);
    f->e[f->n++] = e;
    return true;
}

/* Complete factorisation of m >= 1 (FactorInteger's engine), merged into f. */
static bool mo_factor_into(MoFactors* f, const mpz_t m) {
    mpz_t pr[PR_MAX_DISTINCT_PRIMES];
    int64_t ex[PR_MAX_DISTINCT_PRIMES];
    int cnt = 0;
    if (!facint_factor_complete(m, pr, ex, PR_MAX_DISTINCT_PRIMES, &cnt)) return false;
    bool ok = true;
    for (int i = 0; i < cnt; i++) {
        if (ok) ok = mo_factors_add(f, pr[i], (unsigned long)ex[i]);
        mpz_clear(pr[i]);
    }
    return ok;
}

/* Multiplicative order of k_in modulo n_in.  Pre: order_out is mpz_init'd
 * by the caller; n_in / k_in are read-only.  Returns true on success.
 * Returns false (without setting order_out) when n is zero, k is not a
 * unit modulo n, or n cannot be factored.
 *
 * Strategy: any element order divides phi(n).  Factor n, build phi(n) =
 * prod p^(a-1) (p - 1) together with its factorisation (each p - 1 is
 * factored on its own -- far easier than factoring phi whole), then for each
 * prime q | phi divide it out of the candidate while k^(order/q) stays 1.
 * The residual is exactly the order.  When fac is non-NULL it receives the
 * factorisation of the order, and *cyclic whether (Z/nZ)^* is cyclic
 * (n = 1, 2, 4, p^a, 2 p^a) -- what the discrete-log engine needs. */
static bool mo_order_mpz(const mpz_t k_in, const mpz_t n_in, mpz_t order_out,
                         MoFactors* fac, int* cyclic) {
    if (mpz_sgn(n_in) == 0) return false;
    mpz_t n;
    mpz_init_set(n, n_in);
    mpz_abs(n, n);                               /* allow negative n */
    if (mpz_cmp_ui(n, 1) == 0) {
        mpz_set_ui(order_out, 1);
        if (fac) fac->n = 0;
        if (cyclic) *cyclic = 1;
        mpz_clear(n);
        return true;
    }

    mpz_t k, g;
    mpz_inits(k, g, NULL);
    mpz_mod(k, k_in, n);                         /* k in [0, n-1] */
    mpz_gcd(g, k, n);
    if (mpz_cmp_ui(g, 1) != 0) {
        mpz_clears(k, g, n, NULL);
        return false;                            /* not coprime */
    }

    mpz_t np[PR_MAX_DISTINCT_PRIMES];
    int64_t ne[PR_MAX_DISTINCT_PRIMES];
    int nn = 0;
    /* A prime modulus (the cryptographic case) skips FactorInteger's trial
     * division: one primality test settles it. */
    if (mth_is_prime_mpz(n)) {
        mpz_init_set(np[0], n);
        ne[0] = 1;
        nn = 1;
    } else if (!facint_factor_complete(n, np, ne, PR_MAX_DISTINCT_PRIMES, &nn)) {
        mpz_clears(k, g, n, NULL);
        return false;
    }

    /* phi(n) and its factorisation. */
    MoFactors phi_f;
    phi_f.n = 0;
    mpz_t phi, pm1;
    mpz_init_set_ui(phi, 1);
    mpz_init(pm1);
    bool ok = true;
    for (int i = 0; i < nn && ok; i++) {
        mpz_sub_ui(pm1, np[i], 1);
        mpz_mul(phi, phi, pm1);
        mpz_pow_ui(g, np[i], (unsigned long)(ne[i] - 1));
        mpz_mul(phi, phi, g);
        ok = mo_factors_add(&phi_f, np[i], (unsigned long)(ne[i] - 1)) &&
             (mpz_cmp_ui(pm1, 1) <= 0 || mo_factor_into(&phi_f, pm1));
    }
    if (cyclic) {
        bool two = nn > 0 && mpz_cmp_ui(np[0], 2) == 0;   /* ascending */
        *cyclic = (nn == 1 && (!two || ne[0] <= 2)) ||
                  (nn == 2 && two && ne[0] == 1);
    }
    for (int i = 0; i < nn; i++) mpz_clear(np[i]);
    if (!ok) {
        mo_factors_clear(&phi_f);
        mpz_clears(k, g, n, phi, pm1, NULL);
        return false;
    }

    /* Strip each prime of phi while k^(order/q) == 1. */
    mpz_t order, q, pow;
    mpz_inits(order, q, pow, NULL);
    mpz_set(order, phi);
    if (fac) fac->n = 0;
    for (size_t i = 0; i < phi_f.n; i++) {
        unsigned long left = phi_f.e[i];
        while (left > 0) {
            mpz_divexact(q, order, phi_f.q[i]);
            mpz_powm(pow, k, q, n);
            if (mpz_cmp_ui(pow, 1) != 0) break;
            mpz_set(order, q);
            left--;
        }
        if (fac) mo_factors_add(fac, phi_f.q[i], left);   /* fits: <= phi_f.n */
    }

    mpz_set(order_out, order);
    mo_factors_clear(&phi_f);
    mpz_clears(order, q, pow, k, g, n, phi, pm1, NULL);
    return true;
}

static Expr* mo_emit_argt(size_t argc) {
    /* The expected tag is argtu for one argument, argt otherwise. */
    if (argc == 1)
        mth_message("MultiplicativeOrder", "argtu", "MultiplicativeOrder called with 1 argument; 2 or 3 arguments are expected.");
    else
        mth_message("MultiplicativeOrder", "argt", "MultiplicativeOrder called with %zu arguments; 2 or 3 arguments are expected.", argc);
    return NULL;
}

/* The order of k modulo n for a machine 1 < n < 2^32 with gcd(k, n) == 1,
 * entirely in uint64 (every product of two residues fits). Same algorithm as
 * mo_order_mpz -- strip each prime of phi(n) while k^(order/p) stays 1 --
 * with phi and its factorisation from the int64 helpers the NDArray kernels
 * use. Returns 0 when it declines (the GMP path then answers). */
#define MO_MACHINE_MAX ((int64_t)1 << 32)

static uint64_t mo_powmod_u64(uint64_t b, uint64_t e, uint64_t m) {
    uint64_t r = 1 % m;
    b %= m;
    while (e) {
        if (e & 1) r = r * b % m;
        e >>= 1;
        if (e) b = b * b % m;
    }
    return r;
}

static int64_t mo_order_machine(int64_t k, int64_t n) {
    if (n <= 1 || n >= MO_MACHINE_MAX) return 0;
    int64_t kk = k % n;
    if (kk < 0) kk += n;
    uint64_t a = (uint64_t)kk, b = (uint64_t)n;
    while (b) { uint64_t t = a % b; a = b; b = t; }
    if (a != 1) return 0;                        /* not a unit: GMP path decides */
    int64_t phi;
    if (!ndint_euler_phi_i64(n, &phi)) return 0;
    int64_t p[NDINT_MAX_FACTORS]; int e[NDINT_MAX_FACTORS]; int np;
    if (!ndint_factor_i64(phi, p, e, &np)) return 0;
    uint64_t order = (uint64_t)phi;
    for (int i = 0; i < np; i++)
        while (order % (uint64_t)p[i] == 0 &&
               mo_powmod_u64((uint64_t)kk, order / (uint64_t)p[i], (uint64_t)n) == 1)
            order /= (uint64_t)p[i];
    return (int64_t)order;
}

Expr* builtin_multiplicativeorder(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 2 || argc > 3) return mo_emit_argt(argc);

    Expr* k_expr = res->data.function.args[0];
    Expr* n_expr = res->data.function.args[1];

    /* Integer-only contract: non-integer numerics (Real, Complex, Rational)
     * and symbolic args flow through unevaluated with no diagnostic, matching
     * MultiplicativeOrder[10., 21] -> MultiplicativeOrder[10., 21]. */
    if (!expr_is_integer_like(k_expr) || !expr_is_integer_like(n_expr)) return NULL;

    /* The common call: two machine integers, a modulus below 2^32. */
    if (argc == 2 && k_expr->type == EXPR_INTEGER && n_expr->type == EXPR_INTEGER) {
        int64_t o = mo_order_machine(k_expr->data.integer, n_expr->data.integer);
        if (o > 0) return expr_new_integer(o);
    }

    mpz_t k, n;
    expr_to_mpz(k_expr, k);
    expr_to_mpz(n_expr, n);
    /* A zero or negative modulus stays unevaluated
     * (MultiplicativeOrder[2, -7] is not 3). */
    if (mpz_sgn(n) <= 0) { mpz_clears(k, n, NULL); return NULL; }

    mpz_t order;
    mpz_init(order);
    MoFactors fac;
    fac.n = 0;
    int cyclic = 0;
    if (!mo_order_mpz(k, n, order, argc == 3 ? &fac : NULL, &cyclic)) {
        mpz_clears(k, n, order, NULL);
        return NULL;                             /* not coprime / failure */
    }

    if (argc == 2) {
        Expr* out = expr_bigint_normalize(expr_new_bigint_from_mpz(order));
        mpz_clears(k, n, order, NULL);
        return out;
    }

    /* 3-arg form: third argument must be a List of integer residues.  A
     * non-List or any non-integer element leaves the call unevaluated. */
    Expr* list = res->data.function.args[2];
    size_t lcount = 0;
    if (list->type == EXPR_FUNCTION &&
        list->data.function.head->type == EXPR_SYMBOL &&
        list->data.function.head->data.symbol.name == SYM_List)
        lcount = list->data.function.arg_count;
    mpz_t* residues = lcount ? (mpz_t*)malloc(sizeof(mpz_t) * lcount) : NULL;
    size_t nres = 0;
    bool list_ok = residues != NULL;
    for (size_t i = 0; list_ok && i < lcount; i++) {
        Expr* r = list->data.function.args[i];
        if (!expr_is_integer_like(r)) { list_ok = false; break; }
        expr_to_mpz(r, residues[nres]);   /* inits residues[nres] */
        nres++;
    }

    Expr* out = NULL;
    mpz_t m;
    mpz_init(m);
    /* The engine mallocs its tables outside GMP's guarded allocator, so defer
     * TimeConstrained's async jump across it (it polls the deadline itself);
     * tc_run_guarded restores the count if a GMP safepoint unwinds us. */
    tc_async_defer_push();
    bool hit = list_ok && nt_dlog_min(k, n, order, (const mpz_t*)fac.q, fac.e,
                                      fac.n, cyclic, (const mpz_t*)residues,
                                      nres, m);
    tc_async_defer_pop();
    if (hit) out = expr_bigint_normalize(expr_new_bigint_from_mpz(m));

    for (size_t i = 0; i < nres; i++) mpz_clear(residues[i]);
    free(residues);
    mo_factors_clear(&fac);
    mpz_clears(k, n, order, m, NULL);
    /* The engine returns cleanly when a TimeConstrained deadline passes;
     * with everything freed, take the timeout here rather than leave the
     * call unevaluated as if no target were a power of k. */
    if (!out) tc_check_deadline();
    return out;
}
