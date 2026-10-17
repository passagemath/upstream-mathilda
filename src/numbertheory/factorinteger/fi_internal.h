/*
 * fi_internal.h — shared internals of the FactorInteger subsystem.
 *
 * Layout (one concern per file, all under src/numbertheory/factorinteger/):
 *
 *   fi_core.c      FiFactors list, FiCtx, deterministic RNG, the recursive
 *                  driver fi_factor() and the Automatic tier table.
 *   fi_builtin.c   FactorInteger[...] and every FactorInteger`Method[...]
 *                  entry point: argument/option parsing, messages, Gaussian
 *                  integers, output construction, registration + docstrings.
 *   fi_u64.c       word-size engine: Montgomery arithmetic, deterministic
 *                  primality, rho / SQUFOF / Lehman / micro-ECM on uint64_t.
 *   fi_u128.c      two-word Montgomery arithmetic and Pollard rho for
 *                  2^64 <= n < 2^127.
 *   fi_trial.c     batched trial division (one bignum pass per word-sized
 *                  product of primes) and the shared small-prime table.
 *   fi_rho.c       Pollard rho (Brent) on bignums.
 *   fi_fermat.c    Fermat (residue-sieved) and Lehman.
 *   fi_squfof.c    Shanks' square forms, racing multipliers.
 *   fi_ecm.c       ECM, Pollard p-1, Williams p+1 (GMP-ECM, with native
 *                  fallbacks when the build has no GMP-ECM).
 *   fi_qs.c        self-initialising quadratic sieve (FLINT qsieve or native).
 *   fi_cfrac.c     continued-fraction factoring (Morrison–Brillhart).
 *   fi_dixon.c     Dixon's random-squares method.
 *   fi_gf2.c       GF(2) nullspace solver shared by CFRAC / Dixon / QS.
 *   fi_rbd.c       Blake rational base descent (moved verbatim).
 *
 * A *splitter* takes an odd composite n that is not a perfect power and
 * either writes a non-trivial divisor 1 < d < n and returns 1, or returns 0
 * when its budget is exhausted.  Splitters never loop unboundedly: every one
 * carries an iteration/curve/relation budget, overridable per call through
 * the method suboptions in FiCtx.
 */
#ifndef MATHILDA_FI_INTERNAL_H
#define MATHILDA_FI_INTERNAL_H

#include <gmp.h>
#include <stddef.h>
#include <stdint.h>

/* ------------------------------------------------------------------ */
/* Methods                                                             */
/* ------------------------------------------------------------------ */

typedef enum {
    FI_AUTOMATIC = 0,
    FI_TRIAL,
    FI_RHO,
    FI_FERMAT,
    FI_LEHMAN,
    FI_SQUFOF,
    FI_CFRAC,
    FI_DIXON,
    FI_QS,
    FI_PM1,
    FI_PP1,
    FI_ECM,
    FI_RBD,
    FI_NMETHODS
} FiMethod;

/* ------------------------------------------------------------------ */
/* Factor list                                                         */
/* ------------------------------------------------------------------ */

typedef struct {
    mpz_t   p;      /* prime (or, for an incomplete result, a composite) */
    int64_t e;      /* exponent; negative for denominator primes */
} FiEntry;

typedef struct {
    FiEntry* v;
    size_t   n, cap;
} FiFactors;

void fi_factors_init(FiFactors* f);
void fi_factors_clear(FiFactors* f);
/* Add p^e, merging with an existing entry for the same p. */
void fi_factors_add(FiFactors* f, const mpz_t p, int64_t e);
void fi_factors_add_ui(FiFactors* f, uint64_t p, int64_t e);
/* Sort ascending by p. */
void fi_factors_sort(FiFactors* f);

/* ------------------------------------------------------------------ */
/* Context                                                             */
/* ------------------------------------------------------------------ */

typedef struct {
    FiMethod method;

    /* FactorInteger[n, k]: at most k entries in the answer (the last may be
     * an unfactored remainder).  0 = complete factorisation. */
    int64_t  k_limit;

    /* Method suboptions; 0 means "the method's default". */
    uint64_t bound;       /* TrialDivision "Bound"                       */
    uint64_t iterations;  /* Rho/Fermat/Lehman/SQUFOF/CFRAC/Dixon budget */
    uint64_t starts;      /* PollardRho "Starts"                         */
    uint64_t curves;      /* ECM "Curves"                                */
    double   B1, B2;      /* ECM / p-1 / p+1 stage bounds                */
    uint64_t sigma;       /* ECM "Sigma" (first curve; then sigma+1, ...)*/
    int      rbd_have_base;
    mpz_t    rbd_a, rbd_b;/* BlakeRationalBaseDescent "Base" -> a/b      */

    /* State. */
    int      incomplete;  /* some composite was returned unfactored      */
    uint64_t rng;         /* xorshift64* state, seeded from n            */
    int      warned;      /* set: suppress nofac (internal callers)      */
} FiCtx;

void     fi_ctx_init(FiCtx* c, FiMethod m);
void     fi_ctx_clear(FiCtx* c);
void     fi_ctx_seed(FiCtx* c, const mpz_t n);
uint64_t fi_rand(FiCtx* c);

/* Probable-prime test used throughout (BPSW; deterministic below 2^64). */
int fi_is_prime(const mpz_t n);

/* Factor n >= 1 completely (subject to c->k_limit) with c->method, adding
 * prime powers to out.  A composite the method could not split is added
 * with exponent equal to its multiplicity and c->incomplete is set. */
void fi_factor(const mpz_t n, FiFactors* out, FiCtx* c);

/* Complete factorisation of |n| (Automatic, silent), sorted; returns 0 if a
 * composite survived the search (the caller must then decline). */
int fi_factor_complete_mpz(const mpz_t n, FiFactors* out);

/* ------------------------------------------------------------------ */
/* Word-size engine (fi_u64.c)                                         */
/* ------------------------------------------------------------------ */

int      fi_u64_is_prime(uint64_t n);
/* Complete factorisation of n >= 1 (multiplicity mult). */
void     fi_u64_factor(uint64_t n, int64_t mult, FiFactors* out);
/* Allocation-free: sorted prime powers of n into p[]/e[] (room for 64). */
int      fi_u64_factor_array(uint64_t n, uint64_t* p, int64_t* e);
/* Splitters on odd composite n; return a non-trivial divisor or 0. */
uint64_t fi_u64_rho(uint64_t n, uint64_t c, uint64_t max_iters);
uint64_t fi_u64_squfof(uint64_t n, uint64_t max_iters);
uint64_t fi_u64_lehman(uint64_t n);
uint64_t fi_u64_fermat(uint64_t n, uint64_t max_iters);
uint64_t fi_u64_ecm(uint64_t n, uint64_t* rng, int max_curves);
uint64_t fi_u64_pm1(uint64_t n, uint64_t B1);
uint64_t fi_u64_pp1(uint64_t n, uint64_t B1, uint64_t seed);
uint64_t fi_u64_gcd(uint64_t a, uint64_t b);
uint64_t fi_u64_isqrt(uint64_t n);
/* Prime bitmap up to *limit (65536), built once. */
const unsigned char* fi_prime_bitmap(uint32_t* limit);
/* Lowest-level helpers shared with fi_u128.c. */
void     fi_mul64(uint64_t a, uint64_t b, uint64_t* hi, uint64_t* lo);

/* ------------------------------------------------------------------ */
/* Two-word engine (fi_u128.c)                                         */
/* ------------------------------------------------------------------ */

/* Pollard–Brent rho on 2^64 <= n < 2^126 using two-limb Montgomery
 * arithmetic.  Returns 1 and sets d on success. */
int fi_u128_rho(mpz_t d, const mpz_t n, uint64_t c, uint64_t max_iters);
/* Two-limb micro-ECM ("tinyecm") for 2^64 <= n < 2^126: up to max_curves
 * Suyama curves at stage-1 bound B1 (B2 = 25 B1).  Returns 1 with d set. */
int fi_u128_ecm(mpz_t d, const mpz_t n, uint64_t* rng, uint32_t B1, int max_curves);

/* ------------------------------------------------------------------ */
/* Trial division (fi_trial.c)                                         */
/* ------------------------------------------------------------------ */

/* Primes below FI_SMALL_PRIME_LIMIT, generated once. */
#define FI_SMALL_PRIME_LIMIT (1u << 20)
const uint32_t* fi_small_primes(size_t* count);

/* Divide out every prime p <= bound from n (in place), adding p^e * mult to
 * out.  Returns the number of distinct primes found. */
size_t fi_trial_divide(mpz_t n, uint64_t bound, int64_t mult, FiFactors* out);

/* ------------------------------------------------------------------ */
/* Bignum splitters                                                    */
/* ------------------------------------------------------------------ */

int fi_split_rho(mpz_t d, const mpz_t n, FiCtx* c, uint64_t max_iters, uint64_t starts);
int fi_split_fermat(mpz_t d, const mpz_t n, FiCtx* c, uint64_t max_iters);
int fi_split_lehman(mpz_t d, const mpz_t n, FiCtx* c, uint64_t max_iters);
int fi_split_squfof(mpz_t d, const mpz_t n, FiCtx* c, uint64_t max_iters);
int fi_split_cfrac(mpz_t d, const mpz_t n, FiCtx* c);
int fi_split_dixon(mpz_t d, const mpz_t n, FiCtx* c);
int fi_split_rbd(mpz_t d, const mpz_t n, FiCtx* c);

int fi_split_qs(mpz_t d, const mpz_t n, FiCtx* c);
/* Knuth–Schroeppel multiplier for CFRAC / QS. */
unsigned fi_ks_multiplier(const mpz_t n);

/* ECM family.  B1/curves of 0 take the method's defaults. */
int fi_split_ecm(mpz_t d, const mpz_t n, FiCtx* c, double B1, uint64_t curves);
int fi_split_pm1(mpz_t d, const mpz_t n, FiCtx* c, double B1, double B2);
int fi_split_pp1(mpz_t d, const mpz_t n, FiCtx* c, double B1, double B2);
int fi_split_pp1_seeds(mpz_t d, const mpz_t n, FiCtx* c, double B1, double B2, int nseeds);
/* Automatic's ECM tier: run curves at the B1 level that targets factors of
 * `digits` digits; returns 1 with d set on success. */
int fi_ecm_level(mpz_t d, const mpz_t n, FiCtx* c, int digits);
/* Expected curve count and B1 for a digit level (GMP-ECM table). */
void fi_ecm_level_params(int digits, double* B1, uint64_t* curves);
/* 1 when this build has GMP-ECM. */
int fi_have_gmp_ecm(void);

/* ------------------------------------------------------------------ */
/* GF(2) linear algebra (fi_gf2.c)                                     */
/* ------------------------------------------------------------------ */

/* Sparse relation matrix over GF(2): `nrows` relations, each a list of
 * column indices (factor-base primes with odd exponent).  Finds up to
 * `max_deps` linear dependencies among the rows; each dependency is a bit
 * vector over rows, returned in deps[k * words(nrows) ...].  Returns the
 * number of dependencies found. */
typedef struct {
    uint32_t  nrows, ncols;
    uint32_t* start;     /* row r's columns are col[start[r] .. start[r+1]) */
    uint32_t* col;
} FiSparseMat;

size_t fi_gf2_dependencies(const FiSparseMat* m, uint64_t* deps, size_t max_deps);
static inline size_t fi_gf2_words(size_t nrows) { return (nrows + 63) / 64; }

/* Relation store for the congruence-of-squares methods.  Column 0 is the
 * sign (-1); column c >= 1 is the prime fb[c].  A relation records
 * X^2 == prod fb[col]^exp * L (mod n), col sorted ascending; L = 1 for a
 * full relation, else one large prime (partials pair up automatically). */
typedef struct FiRels FiRels;
FiRels* fi_rels_new(const uint32_t* fb, size_t ncols);
void    fi_rels_free(FiRels* R);
int     fi_rels_add(FiRels* R, const mpz_t X, const uint32_t* col, const uint32_t* exps,
                    size_t cnt, uint64_t L, const mpz_t n);
size_t  fi_rels_count(const FiRels* R);
size_t  fi_rels_partials(const FiRels* R);
/* Columns that occur with an odd exponent in some full relation: once the
 * relation count exceeds this, dependencies must exist. */
size_t  fi_rels_active(const FiRels* R);
/* Linear algebra + square roots; returns 1 with a proper divisor in d. */
int     fi_rels_solve(FiRels* R, const mpz_t n, mpz_t d);

/* ------------------------------------------------------------------ */
/* Messages                                                            */
/* ------------------------------------------------------------------ */

/* FactorInteger::nofac for composite m (unless c->warned). */
void fi_warn_incomplete(FiCtx* c, const mpz_t m);

/* ------------------------------------------------------------------ */
/* Bit helpers (GCC builtins where available, portable otherwise)      */
/* ------------------------------------------------------------------ */

static inline int fi_ctz64(uint64_t x) {
#if defined(__GNUC__)
    return __builtin_ctzll(x);
#else
    int c = 0;
    while (!(x & 1)) { x >>= 1; c++; }
    return c;
#endif
}

static inline int fi_parity64(uint64_t x) {
#if defined(__GNUC__)
    return __builtin_parityll(x);
#else
    x ^= x >> 32; x ^= x >> 16; x ^= x >> 8; x ^= x >> 4; x ^= x >> 2; x ^= x >> 1;
    return (int)(x & 1);
#endif
}

/* ------------------------------------------------------------------ */
/* Tuning hooks                                                        */
/* ------------------------------------------------------------------ */

/* Build with -DFI_TUNING to let environment variables override the tuned
 * constants (FI_QS_FB, FI_QS_BLK, FI_QS_LP, FI_QS_SLACK, FI_QS_SMALL,
 * FI_QS_MED, FI_TINY_B1, FI_TINY_CURVES, FI_DENSE_MAX) for parameter
 * sweeps.  Compiled out otherwise. */
#ifdef FI_TUNING
#include <stdlib.h>
#define FI_TUNE_U(name, var) do { const char* s_ = getenv(name); if (s_) (var) = (unsigned)atol(s_); } while (0)
#define FI_TUNE_D(name, var) do { const char* s_ = getenv(name); if (s_) (var) = atof(s_); } while (0)
#else
#define FI_TUNE_U(name, var) ((void)0)
#define FI_TUNE_D(name, var) ((void)0)
#endif

#endif /* MATHILDA_FI_INTERNAL_H */
