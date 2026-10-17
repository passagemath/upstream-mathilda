/*
 * fi_builtin.c — FactorInteger[...] and the FactorInteger`Method[...] entry
 * points: argument and option parsing, rational and Gaussian inputs, partial
 * factorisation, messages, output, registration and docstrings.
 *
 *   FactorInteger[n]                       complete factorisation
 *   FactorInteger[n, k]                    at most k entries; the last may be
 *                                          an unfactored remainder
 *   FactorInteger[n, Automatic]            same as FactorInteger[n]
 *   FactorInteger[n, Method -> m]          m a method name, Automatic, or
 *                                          {"Name", "Sub" -> v, ...}
 *   FactorInteger[n, GaussianIntegers -> True]   factor over Z[i]
 *   FactorInteger`Name[n, ...]             = FactorInteger[n, Method -> {"Name", ...}]
 *
 * n may be an integer, a rational (the denominator's primes get negative
 * exponents), or a Gaussian integer / Gaussian rational (factored over Z[i]
 * automatically, unit first, primes in the first quadrant sorted by real
 * then imaginary part — Mathematica's normalisation).
 */
#include "fi_internal.h"

#include "arithmetic.h"
#include "attr.h"
#include "eval.h"
#include "expr.h"
#include "message.h"
#include "print.h"
#include "sym_names.h"
#include "symtab.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ------------------------------------------------------------------ */
/* Method table                                                        */
/* ------------------------------------------------------------------ */

enum { SUB_BOUND = 1, SUB_ITER = 2, SUB_STARTS = 4, SUB_B1 = 8, SUB_B2 = 16,
       SUB_CURVES = 32, SUB_SIGMA = 64, SUB_BASE = 128 };

typedef struct {
    const char* name;       /* canonical Method string = FactorInteger`name */
    FiMethod    m;
    unsigned    subs;       /* accepted suboptions */
    const char* doc;        /* FactorInteger`name docstring */
} MethodInfo;

static const MethodInfo METHODS[] = {
    { "TrialDivision", FI_TRIAL, SUB_BOUND,
      "FactorInteger`TrialDivision[n] factors n by trial division only.\n"
      "\tPrimes are tested in word-sized batches (one bignum remainder per\n"
      "\tbatch).  Default bound: sqrt(n) when that is below 10^7, else 10^7;\n"
      "\ta cofactor left above the bound is returned with FactorInteger::nofac.\n"
      "\tSuboption: \"Bound\" -> B (at most 2^40)." },
    { "PollardRho", FI_RHO, SUB_ITER | SUB_STARTS,
      "FactorInteger`PollardRho[n] factors n with Brent's variant of Pollard rho\n"
      "\t(batched gcds, backtracking; Montgomery arithmetic on one and two words).\n"
      "\tFinds a factor p in about sqrt(p) steps.\n"
      "\tSuboptions: \"Iterations\" -> steps per start (default 2^24 on bignums,\n"
      "\t2^26 on words), \"Starts\" -> number of polynomials y^2 + c (default 4/16)." },
    { "Fermat", FI_FERMAT, SUB_ITER,
      "FactorInteger`Fermat[n] factors n by Fermat's method: it finds n = a^2 - b^2\n"
      "\tquickly when n has two factors close to sqrt(n).  Candidates are\n"
      "\tfiltered by quadratic-residue tables modulo nine small moduli.\n"
      "\tSuboption: \"Iterations\" -> number of a values tried (default 2^26)." },
    { "Lehman", FI_LEHMAN, SUB_ITER,
      "FactorInteger`Lehman[n] factors n by Lehman's O(n^(1/3)) method: trial\n"
      "\tdivision to n^(1/3), then Fermat on 4 k n for k <= n^(1/3).  Complete\n"
      "\tbelow 2^42; on larger n it runs within \"Iterations\" (default 2^24)." },
    { "ShanksSquareForms", FI_SQUFOF, SUB_ITER,
      "FactorInteger`ShanksSquareForms[n] factors n by Shanks' square forms\n"
      "\t(SQUFOF), trying the 16 Gower-Wagstaff multipliers in turn.  Expected\n"
      "\twork O(n^(1/4)); runs in word arithmetic for k n < 2^124 and declines\n"
      "\tbeyond.  Suboption: \"Iterations\" -> forward steps per multiplier." },
    { "CFRAC", FI_CFRAC, SUB_ITER,
      "FactorInteger`CFRAC[n] factors n by the continued-fraction method of\n"
      "\tMorrison and Brillhart: Knuth-Schroeppel multiplier, factor base sized\n"
      "\tby digit count, early-abort trial division, single large primes, and\n"
      "\tGF(2) elimination.  Suboption: \"Iterations\" -> continued-fraction\n"
      "\tsteps allowed (default 2^32)." },
    { "Dixon", FI_DIXON, SUB_ITER,
      "FactorInteger`Dixon[n] factors n by Dixon's random-squares method, with\n"
      "\tthe squares drawn near sqrt(k n) so the residues are square-root\n"
      "\tsized, early-abort trial division, single large primes and GF(2)\n"
      "\telimination.  Suboption: \"Iterations\" -> squares tried (default 2^30)." },
    { "QuadraticSieve", FI_QS, 0,
      "FactorInteger`QuadraticSieve[n] factors n by the self-initialising\n"
      "\tquadratic sieve, the fastest general method for 30-100 digit\n"
      "\tcomposites with no small factors.  Below 2^64 the word engine is used." },
    { "PollardPMinus1", FI_PM1, SUB_B1 | SUB_B2,
      "FactorInteger`PollardPMinus1[n] factors n by Pollard's p-1 method: it\n"
      "\tfinds a prime p of n when p - 1 is B1-smooth apart from one prime up to\n"
      "\tB2.  Uses GMP-ECM's fast stage 2 when available.\n"
      "\tSuboptions: \"B1\" (default 10^6), \"B2\" (default: GMP-ECM's choice)." },
    { "WilliamsPPlus1", FI_PP1, SUB_B1 | SUB_B2,
      "FactorInteger`WilliamsPPlus1[n] factors n by Williams' p+1 method: it\n"
      "\tfinds p when p + 1 is smooth (three seeds are run, since one seed\n"
      "\tworks only when its discriminant is a non-residue mod p).\n"
      "\tSuboptions: \"B1\" (default 10^6), \"B2\"." },
    { "ECM", FI_ECM, SUB_B1 | SUB_B2 | SUB_CURVES | SUB_SIGMA,
      "FactorInteger`ECM[n] factors n by Lenstra's elliptic-curve method (GMP-ECM),\n"
      "\twhose running time depends on the size of the factor found, not of n.\n"
      "\tBy default it climbs GMP-ECM's optimal levels (B1 = 2000, 11000, 5e4,\n"
      "\t2.5e5, 1e6, ...) with the expected curve count at each.\n"
      "\tSuboptions: \"B1\", \"B2\", \"Curves\", \"Sigma\" (first curve parameter)." },
    { "BlakeRationalBaseDescent", FI_RBD, SUB_BASE,
      "FactorInteger`BlakeRationalBaseDescent[n] factors n by rational base\n"
      "\tdescent: assuming a factor near c (a/b)^k it descends Q <- floor(b Q/a)\n"
      "\tand probes gcd(Q +- j, n).  Suboption: \"Base\" -> a/b (default: a\n"
      "\tbounded search over coprime a > b)." },
};
#define NMETHODS (sizeof METHODS / sizeof METHODS[0])

static const struct { const char* alias; FiMethod m; } ALIASES[] = {
    { "Automatic", FI_AUTOMATIC },
    { "PollardP-1", FI_PM1 }, { "PollardPM1", FI_PM1 }, { "P-1", FI_PM1 },
    { "WilliamsP+1", FI_PP1 }, { "WilliamsPP1", FI_PP1 }, { "P+1", FI_PP1 },
    { "SQUFOF", FI_SQUFOF },
    { "QS", FI_QS }, { "SIQS", FI_QS }, { "MPQS", FI_QS },
    { "ContinuedFraction", FI_CFRAC },
    { "EllipticCurve", FI_ECM },
};

/* Symbol aliases registered in addition to the canonical names. */
static const struct { const char* sym; FiMethod m; } SYM_ALIASES[] = {
    { "SQUFOF", FI_SQUFOF }, { "QS", FI_QS }, { "SIQS", FI_QS },
};

static const MethodInfo* method_info(FiMethod m) {
    for (size_t i = 0; i < NMETHODS; i++) if (METHODS[i].m == m) return &METHODS[i];
    return NULL;
}

static int lookup_method(const char* s, FiMethod* out) {
    for (size_t i = 0; i < NMETHODS; i++)
        if (strcmp(METHODS[i].name, s) == 0) { *out = METHODS[i].m; return 1; }
    for (size_t i = 0; i < sizeof ALIASES / sizeof ALIASES[0]; i++)
        if (strcmp(ALIASES[i].alias, s) == 0) { *out = ALIASES[i].m; return 1; }
    return 0;
}

/* ------------------------------------------------------------------ */
/* Messages                                                            */
/* ------------------------------------------------------------------ */

/* One message per unfactored composite (a call can leave several). */
void fi_warn_incomplete(FiCtx* c, const mpz_t m) {
    if (c->warned) return;              /* internal callers: silent */
    char* s = mpz_get_str(NULL, 10, m);
    mth_message_gated(g_arith_warnings_muted, "FactorInteger", "nofac",
        "%s is composite but no factor was found within the search bounds; "
        "it is returned unfactored with exponent 1.", s ? s : "The argument");
    free(s);
}

/* ------------------------------------------------------------------ */
/* Expression helpers                                                  */
/* ------------------------------------------------------------------ */

static int is_head(const Expr* e, const char* sym, size_t argc) {
    return e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL &&
           e->data.function.head->data.symbol.name == sym && e->data.function.arg_count == argc;
}

static int is_rule(const Expr* e) {
    return is_head(e, SYM_Rule, 2) || is_head(e, SYM_RuleDelayed, 2);
}

/* n as a rational p/q (q > 0).  Returns 0 if n is not an exact rational. */
static int get_rational(const Expr* e, mpz_t p, mpz_t q) {
    if (expr_is_integer_like(e)) { expr_to_mpz((Expr*)e, p); mpz_init_set_ui(q, 1); return 1; }
    if (is_head(e, SYM_Rational, 2) && expr_is_integer_like(e->data.function.args[0]) &&
        expr_is_integer_like(e->data.function.args[1])) {
        expr_to_mpz(e->data.function.args[0], p);
        expr_to_mpz(e->data.function.args[1], q);
        if (mpz_sgn(q) < 0) { mpz_neg(p, p); mpz_neg(q, q); }
        return 1;
    }
    return 0;
}

static Expr* mpz_expr(const mpz_t z) { return expr_bigint_normalize(expr_new_bigint_from_mpz(z)); }

static Expr* pair_expr(Expr* a, int64_t e) {
    Expr* args[2] = { a, expr_new_integer(e) };
    return expr_new_function(expr_new_symbol(SYM_List), args, 2);
}

/* ------------------------------------------------------------------ */
/* Option parsing                                                      */
/* ------------------------------------------------------------------ */

typedef struct {
    int gaussian;
    int have_k;
    int64_t k;
} CallOpts;

/* Read a numeric suboption value: positive integer or real. */
static int sub_number(const Expr* v, double* out) {
    if (v->type == EXPR_INTEGER) { *out = (double)v->data.integer; return v->data.integer > 0; }
    if (v->type == EXPR_REAL) { *out = v->data.real; return v->data.real > 0; }
    if (v->type == EXPR_BIGINT) { mpz_t z; expr_to_mpz((Expr*)v, z); *out = mpz_get_d(z); mpz_clear(z); return *out > 0; }
    return 0;
}

/* Apply one method suboption ("Name" -> value).  Returns 0 (after a message)
 * on an unknown or malformed suboption. */
static int apply_sub(FiCtx* c, const Expr* key, const Expr* val) {
    const char* k = key->type == EXPR_STRING ? key->data.string :
                    key->type == EXPR_SYMBOL ? key->data.symbol.name : NULL;
    const MethodInfo* mi = method_info(c->method);
    unsigned allowed = mi ? mi->subs : 0;
    const char* mname = mi ? mi->name : "Automatic";
    double x = 0;
    if (!k) goto bad;
    if (strcmp(k, "Base") == 0 && (allowed & SUB_BASE)) {
        mpz_t p, q;
        if (!get_rational(val, p, q)) goto bad;
        mpz_set(c->rbd_a, p);
        mpz_set(c->rbd_b, q);
        c->rbd_have_base = 1;
        mpz_clears(p, q, NULL);
        return 1;
    }
    if (!sub_number(val, &x)) goto bad;
    if      (strcmp(k, "Bound") == 0 && (allowed & SUB_BOUND))       c->bound = (uint64_t)x;
    else if (strcmp(k, "Iterations") == 0 && (allowed & SUB_ITER))   c->iterations = (uint64_t)x;
    else if (strcmp(k, "Starts") == 0 && (allowed & SUB_STARTS))     c->starts = (uint64_t)x;
    else if (strcmp(k, "Curves") == 0 && (allowed & SUB_CURVES))     c->curves = (uint64_t)x;
    else if (strcmp(k, "B1") == 0 && (allowed & SUB_B1))             c->B1 = x;
    else if (strcmp(k, "B2") == 0 && (allowed & SUB_B2))             c->B2 = x;
    else if (strcmp(k, "Sigma") == 0 && (allowed & SUB_SIGMA))       c->sigma = (uint64_t)x;
    else goto bad;
    return 1;
bad:
    mth_message("FactorInteger", "mopt",
        "%s is not a valid suboption setting for Method -> \"%s\".",
        k ? k : "The suboption", mname);
    return 0;
}

/* Method -> spec.  Returns 0 (after a message) when spec is not recognised. */
static int apply_method(FiCtx* c, const Expr* spec) {
    if (spec->type == EXPR_SYMBOL && spec->data.symbol.name == SYM_Automatic) {
        c->method = FI_AUTOMATIC;
        return 1;
    }
    if (spec->type == EXPR_STRING) {
        if (lookup_method(spec->data.string, &c->method)) return 1;
        mth_message("FactorInteger", "method",
            "Method -> \"%s\" is not a known factoring method; see ?FactorInteger.",
            spec->data.string);
        return 0;
    }
    if (spec->type == EXPR_FUNCTION && spec->data.function.head->type == EXPR_SYMBOL &&
        spec->data.function.head->data.symbol.name == SYM_List && spec->data.function.arg_count >= 1) {
        if (!apply_method(c, spec->data.function.args[0])) return 0;
        for (size_t i = 1; i < spec->data.function.arg_count; i++) {
            const Expr* r = spec->data.function.args[i];
            if (!is_rule(r) || !apply_sub(c, r->data.function.args[0], r->data.function.args[1]))
                return 0;
        }
        return 1;
    }
    mth_message("FactorInteger", "method",
        "The Method setting is not a method name, Automatic, or {name, suboptions...}.");
    return 0;
}

/* Parse args[first..] of FactorInteger (or of a FactorInteger`X call when
 * `named` is set: bare suboption rules allowed, Method not).  Returns 0 to
 * leave the call unevaluated. */
static int parse_args(const Expr* res, size_t first, FiCtx* c, CallOpts* o, int named) {
    size_t argc = res->data.function.arg_count;
    for (size_t i = first; i < argc; i++) {
        const Expr* a = res->data.function.args[i];
        if (a->type == EXPR_INTEGER && i == first) {
            if (a->data.integer < 1) return 0;          /* FactorInteger[n, 0]: unevaluated */
            o->have_k = 1;
            o->k = a->data.integer;
            continue;
        }
        if (a->type == EXPR_SYMBOL && a->data.symbol.name == SYM_Automatic && i == first) continue;
        if (!is_rule(a)) return 0;
        const Expr* lhs = a->data.function.args[0];
        const Expr* rhs = a->data.function.args[1];
        if (lhs->type == EXPR_SYMBOL && lhs->data.symbol.name == SYM_GaussianIntegers) {
            if (rhs->type != EXPR_SYMBOL) return 0;
            if (rhs->data.symbol.name == SYM_True) o->gaussian = 1;
            else if (rhs->data.symbol.name == SYM_False) o->gaussian = 0;
            else return 0;
            continue;
        }
        if (!named && lhs->type == EXPR_SYMBOL && lhs->data.symbol.name == SYM_Method) {
            if (!apply_method(c, rhs)) return 0;
            continue;
        }
        if (named && lhs->type == EXPR_STRING) {
            if (!apply_sub(c, lhs, rhs)) return 0;
            continue;
        }
        if (lhs->type == EXPR_SYMBOL) {
            mth_message("FactorInteger", "optx", "Unknown option %s in FactorInteger.",
                        lhs->data.symbol.name);
        }
        return 0;
    }
    return 1;
}

/* ------------------------------------------------------------------ */
/* Rational factorisation                                              */
/* ------------------------------------------------------------------ */

/* Factor p/q into out (q's primes negative). */
static void factor_rational(const mpz_t p, const mpz_t q, FiFactors* out, FiCtx* c) {
    mpz_t a;
    mpz_init(a);
    mpz_abs(a, p);
    fi_factor(a, out, c);
    if (mpz_cmp_ui(q, 1) > 0) {
        FiFactors den;
        fi_factors_init(&den);
        int64_t k = c->k_limit;
        c->k_limit = 0;                 /* partial applies to the whole answer */
        fi_factor(q, &den, c);
        c->k_limit = k;
        for (size_t i = 0; i < den.n; i++) fi_factors_add(out, den.v[i].p, -den.v[i].e);
        fi_factors_clear(&den);
    }
    mpz_clear(a);
}

/* Build the output list for a rational input. */
static Expr* rational_result(const mpz_t p, const mpz_t q, FiFactors* f, int64_t k) {
    fi_factors_sort(f);
    int neg = mpz_sgn(p) < 0;
    size_t total = f->n + (size_t)neg;
    Expr** items;
    size_t ni = 0;
    if (k > 0 && total > (size_t)k) {
        /* Partial: keep the k-1 smallest primes, fold the rest (and the sign)
         * into one remainder entry {r, 1}. */
        size_t keep = (size_t)k - 1;
        mpq_t r;
        mpq_init(r);
        mpq_set_si(r, neg ? -1 : 1, 1);
        mpz_t pw;
        mpz_init(pw);
        for (size_t i = keep; i < f->n; i++) {
            int64_t e = f->v[i].e;
            mpz_pow_ui(pw, f->v[i].p, (unsigned long)(e < 0 ? -e : e));
            if (e > 0) mpz_mul(mpq_numref(r), mpq_numref(r), pw);
            else mpz_mul(mpq_denref(r), mpq_denref(r), pw);
        }
        mpq_canonicalize(r);
        items = malloc((keep + 1) * sizeof(Expr*));
        Expr* rem = make_rational_mpz(mpq_numref(r), mpq_denref(r));
        /* Sort the remainder into place by value. */
        int placed = 0;
        for (size_t i = 0; i < keep; i++) {
            if (!placed && mpq_cmp_z(r, f->v[i].p) < 0) { items[ni++] = pair_expr(rem, 1); placed = 1; }
            items[ni++] = pair_expr(mpz_expr(f->v[i].p), f->v[i].e);
        }
        if (!placed) items[ni++] = pair_expr(rem, 1);
        mpz_clear(pw);
        mpq_clear(r);
    } else {
        items = malloc((total + 1) * sizeof(Expr*));
        if (neg) items[ni++] = pair_expr(expr_new_integer(-1), 1);
        for (size_t i = 0; i < f->n; i++) items[ni++] = pair_expr(mpz_expr(f->v[i].p), f->v[i].e);
    }
    (void)q;
    Expr* out = expr_new_function(expr_new_symbol(SYM_List), items, ni);
    free(items);
    return out;
}

/* ------------------------------------------------------------------ */
/* Gaussian integers                                                   */
/* ------------------------------------------------------------------ */

typedef struct { mpz_t re, im; int64_t e; } GEntry;
typedef struct { GEntry* v; size_t n, cap; } GList;

static void gl_add(GList* g, const mpz_t re, const mpz_t im, int64_t e) {
    for (size_t i = 0; i < g->n; i++)
        if (mpz_cmp(g->v[i].re, re) == 0 && mpz_cmp(g->v[i].im, im) == 0) { g->v[i].e += e; return; }
    if (g->n == g->cap) { g->cap = g->cap ? 2 * g->cap : 8; g->v = realloc(g->v, g->cap * sizeof(GEntry)); }
    mpz_init_set(g->v[g->n].re, re);
    mpz_init_set(g->v[g->n].im, im);
    g->v[g->n].e = e;
    g->n++;
}
static void gl_clear(GList* g) {
    for (size_t i = 0; i < g->n; i++) mpz_clears(g->v[i].re, g->v[i].im, NULL);
    free(g->v);
}

/* Rotate (re, im) by powers of i into the first quadrant: re > 0, im >= 0. */
static void normalize_quadrant(mpz_t re, mpz_t im) {
    mpz_t t;
    mpz_init(t);
    for (int k = 0; k < 4 && !(mpz_sgn(re) > 0 && mpz_sgn(im) >= 0); k++) {
        mpz_neg(t, im); mpz_set(im, re); mpz_set(re, t);   /* times i */
    }
    mpz_clear(t);
}

/* p == 1 (mod 4) prime: write p = c^2 + d^2 (Cornacchia / Hermite–Serret). */
static void two_squares(const mpz_t p, mpz_t c, mpz_t d) {
    mpz_t t, x, e, a, b, r;
    mpz_inits(t, x, e, a, b, r, NULL);
    for (unsigned long tv = 2;; tv++) {
        mpz_set_ui(t, tv);
        if (mpz_legendre(t, p) == -1) break;
    }
    mpz_sub_ui(e, p, 1);
    mpz_fdiv_q_2exp(e, e, 2);
    mpz_powm(x, t, e, p);                   /* x^2 == -1 (mod p) */
    mpz_set(a, p);
    mpz_set(b, x);
    mpz_mul(r, b, b);
    while (mpz_cmp(r, p) > 0) {
        mpz_mod(r, a, b);
        mpz_set(a, b);
        mpz_set(b, r);
        mpz_mul(r, b, b);
    }
    mpz_set(c, b);
    mpz_sub(r, p, r);
    mpz_sqrt(d, r);
    mpz_clears(t, x, e, a, b, r, NULL);
}

/* If (u + v i) divides (a + b i), replace (a, b) by the quotient. */
static int gdiv(mpz_t a, mpz_t b, const mpz_t u, const mpz_t v) {
    mpz_t N, re, im, t;
    mpz_inits(N, re, im, t, NULL);
    mpz_mul(N, u, u); mpz_addmul(N, v, v);
    mpz_mul(re, a, u); mpz_addmul(re, b, v);
    mpz_mul(im, b, u); mpz_submul(im, a, v);
    int ok = mpz_divisible_p(re, N) && mpz_divisible_p(im, N);
    if (ok) { mpz_divexact(a, re, N); mpz_divexact(b, im, N); }
    mpz_clears(N, re, im, t, NULL);
    return ok;
}

/* Factor the nonzero Gaussian integer a + b i into first-quadrant primes;
 * (ua, ub) receives the unit.  Returns 0 if the norm could not be fully
 * factored (the call is then left unevaluated). */
static int gaussian_factor(const mpz_t a0, const mpz_t b0, GList* out, mpz_t ua, mpz_t ub, FiCtx* c) {
    mpz_t N, a, b, cc, dd, cr, ci;
    mpz_inits(N, a, b, cc, dd, cr, ci, NULL);
    mpz_set(a, a0);
    mpz_set(b, b0);
    mpz_mul(N, a, a);
    mpz_addmul(N, b, b);
    FiFactors nf;
    fi_factors_init(&nf);
    int64_t kl = c->k_limit;
    c->k_limit = 0;
    fi_factor(N, &nf, c);
    c->k_limit = kl;
    int ok = !c->incomplete;
    for (size_t i = 0; i < nf.n && ok; i++) {
        const mpz_t* p = (const mpz_t*)&nf.v[i].p;
        int64_t e = nf.v[i].e;
        if (mpz_cmp_ui(*p, 2) == 0) {
            mpz_set_ui(cr, 1); mpz_set_ui(ci, 1);
            for (int64_t k = 0; k < e; k++) gdiv(a, b, cr, ci);
            gl_add(out, cr, ci, e);
        } else if (mpz_fdiv_ui(*p, 4) == 3) {
            mpz_set(cr, *p); mpz_set_ui(ci, 0);
            for (int64_t k = 0; k < e / 2; k++) gdiv(a, b, cr, ci);
            gl_add(out, cr, ci, e / 2);
        } else {
            two_squares(*p, cc, dd);
            /* pi = c + d i and its conjugate, each taken in its first-quadrant
             * form BEFORE dividing, so the leftover unit is relative to the
             * primes actually reported. */
            int64_t ep = 0, eq = 0;
            mpz_set(cr, cc); mpz_set(ci, dd);
            normalize_quadrant(cr, ci);
            while (ep < e && gdiv(a, b, cr, ci)) ep++;
            if (ep) gl_add(out, cr, ci, ep);
            mpz_set(cr, cc); mpz_neg(ci, dd);
            normalize_quadrant(cr, ci);
            while (eq < e - ep && gdiv(a, b, cr, ci)) eq++;
            if (eq) gl_add(out, cr, ci, eq);
        }
    }
    mpz_set(ua, a);                         /* what remains is a unit */
    mpz_set(ub, b);
    fi_factors_clear(&nf);
    mpz_clears(N, a, b, cc, dd, cr, ci, NULL);
    return ok;
}

static int gcmp(const void* x, const void* y) {
    const GEntry* p = x;
    const GEntry* q = y;
    int c = mpz_cmp(p->re, q->re);
    return c ? c : mpz_cmp(p->im, q->im);
}

static Expr* gaussian_expr(const mpz_t re, const mpz_t im) {
    if (mpz_sgn(im) == 0) return mpz_expr(re);
    return make_complex(mpz_expr(re), mpz_expr(im));
}

/* Factor z = (pr/qr) + (pi/qi) i over Z[i]. */
static Expr* gaussian_result(const mpz_t pr, const mpz_t qr, const mpz_t pi, const mpz_t qi, FiCtx* c) {
    mpz_t m, a, b, ua, ub, va, vb, t;
    mpz_inits(m, a, b, ua, ub, va, vb, t, NULL);
    mpz_lcm(m, qr, qi);
    mpz_divexact(t, m, qr); mpz_mul(a, pr, t);
    mpz_divexact(t, m, qi); mpz_mul(b, pi, t);
    Expr* out = NULL;
    if (mpz_sgn(a) == 0 && mpz_sgn(b) == 0) {
        Expr* items[1] = { pair_expr(expr_new_integer(0), 1) };
        out = expr_new_function(expr_new_symbol(SYM_List), items, 1);
        goto done;
    }
    GList g = { NULL, 0, 0 };
    int ok = gaussian_factor(a, b, &g, ua, ub, c);
    mpz_set_ui(va, 1); mpz_set_ui(vb, 0);
    if (ok && mpz_cmp_ui(m, 1) > 0) {
        GList h = { NULL, 0, 0 };
        mpz_set_ui(t, 0);
        ok = gaussian_factor(m, t, &h, va, vb, c);
        for (size_t i = 0; i < h.n; i++) gl_add(&g, h.v[i].re, h.v[i].im, -h.v[i].e);
        gl_clear(&h);
    }
    if (ok) {
        /* unit = u / v = u * conj(v) for units. */
        mpz_t wr, wi;
        mpz_inits(wr, wi, NULL);
        mpz_mul(wr, ua, va); mpz_addmul(wr, ub, vb);
        mpz_mul(wi, ub, va); mpz_submul(wi, ua, vb);
        size_t kept = 0;
        for (size_t i = 0; i < g.n; i++) if (g.v[i].e) g.v[kept++] = g.v[i]; else mpz_clears(g.v[i].re, g.v[i].im, NULL);
        g.n = kept;
        qsort(g.v, g.n, sizeof(GEntry), gcmp);
        int unit_one = mpz_cmp_ui(wr, 1) == 0 && mpz_sgn(wi) == 0;
        size_t ni = 0;
        Expr** items = malloc((g.n + 1) * sizeof(Expr*));
        if (!unit_one || g.n == 0) items[ni++] = pair_expr(gaussian_expr(wr, wi), 1);
        for (size_t i = 0; i < g.n; i++) items[ni++] = pair_expr(gaussian_expr(g.v[i].re, g.v[i].im), g.v[i].e);
        out = expr_new_function(expr_new_symbol(SYM_List), items, ni);
        free(items);
        mpz_clears(wr, wi, NULL);
    }
    gl_clear(&g);
done:
    mpz_clears(m, a, b, ua, ub, va, vb, t, NULL);
    return out;
}

/* ------------------------------------------------------------------ */
/* Entry points                                                        */
/* ------------------------------------------------------------------ */

static Expr* factorinteger_core(Expr* res, FiCtx* c, int named) {
    size_t argc = res->data.function.arg_count;
    if (argc < 1) return NULL;
    Expr* n = res->data.function.args[0];
    CallOpts o = { 0, 0, 0 };
    if (!parse_args(res, 1, c, &o, named)) return NULL;
    c->k_limit = o.have_k ? o.k : 0;

    mpz_t p, q;
    if (get_rational(n, p, q)) {
        Expr* out;
        if (o.gaussian) {
            mpz_t zi, oi;
            mpz_init_set_ui(zi, 0);
            mpz_init_set_ui(oi, 1);
            out = gaussian_result(p, q, zi, oi, c);
            mpz_clears(zi, oi, NULL);
        } else if (mpz_sgn(p) == 0) {
            Expr* items[1] = { pair_expr(expr_new_integer(0), 1) };
            out = expr_new_function(expr_new_symbol(SYM_List), items, 1);
        } else if (mpz_cmpabs_ui(p, 1) == 0 && mpz_cmp_ui(q, 1) == 0) {
            Expr* items[1] = { pair_expr(mpz_expr(p), 1) };
            out = expr_new_function(expr_new_symbol(SYM_List), items, 1);
        } else {
            FiFactors f;
            fi_factors_init(&f);
            factor_rational(p, q, &f, c);
            out = rational_result(p, q, &f, c->k_limit);
            fi_factors_clear(&f);
        }
        mpz_clears(p, q, NULL);
        return out;
    }
    Expr *re, *im;
    if (is_complex(n, &re, &im)) {
        mpz_t pr, qr, pi, qi;
        if (!get_rational(re, pr, qr)) goto inexact;
        if (!get_rational(im, pi, qi)) { mpz_clears(pr, qr, NULL); goto inexact; }
        Expr* out = gaussian_result(pr, qr, pi, qi, c);
        mpz_clears(pr, qr, pi, qi, NULL);
        return out;
    }
inexact:
    /* Mathematica: FactorInteger::exact for an inexact or symbolic argument
     * (lists thread through Listable before reaching here). */
    if (n->type == EXPR_REAL || n->type == EXPR_SYMBOL || is_complex(n, &re, &im)) {
        if (!(n->type == EXPR_SYMBOL && n->data.symbol.name == SYM_List)) {
            char* s = expr_to_string(n);
            mth_message("FactorInteger", "exact", "Argument %s in FactorInteger is not an exact number.",
                        s ? s : "");
            free(s);
        }
    }
    return NULL;
}

Expr* builtin_factorinteger(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    /* Fast path: FactorInteger[n] for a machine integer n >= 2 — no context,
     * no bignums, no heap until the result list (most calls are this). */
    if (res->data.function.arg_count == 1 && res->data.function.args[0]->type == EXPR_INTEGER &&
        res->data.function.args[0]->data.integer >= 2) {
        uint64_t p[64];
        int64_t e[64];
        int k = fi_u64_factor_array((uint64_t)res->data.function.args[0]->data.integer, p, e);
        Expr* items[64];
        for (int i = 0; i < k; i++) items[i] = pair_expr(expr_new_integer((int64_t)p[i]), e[i]);
        return expr_new_function(expr_new_symbol(SYM_List), items, (size_t)k);
    }
    FiCtx c;
    fi_ctx_init(&c, FI_AUTOMATIC);
    Expr* out = factorinteger_core(res, &c, 0);
    fi_ctx_clear(&c);
    return out;
}

/* FactorInteger`Name[n, ...]: the method comes from the head's name. */
static Expr* builtin_factorinteger_method(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.head->type != EXPR_SYMBOL) return NULL;
    const char* full = res->data.function.head->data.symbol.name;
    const char* tick = strrchr(full, '`');
    FiMethod m;
    if (!tick || !lookup_method(tick + 1, &m)) return NULL;
    FiCtx c;
    fi_ctx_init(&c, m);
    Expr* out = factorinteger_core(res, &c, 1);
    fi_ctx_clear(&c);
    return out;
}

/* ------------------------------------------------------------------ */
/* Public helpers for other modules                                    */
/* ------------------------------------------------------------------ */

/* Complete factorisation of |n| with the Automatic method (see facint.h). */
int fi_factor_complete_mpz(const mpz_t n, FiFactors* out) {
    FiCtx c;
    fi_ctx_init(&c, FI_AUTOMATIC);
    c.warned = 1;                         /* callers decline on incompleteness */
    mpz_t a;
    mpz_init(a);
    mpz_abs(a, n);
    fi_factor(a, out, &c);
    mpz_clear(a);
    int ok = !c.incomplete;
    fi_ctx_clear(&c);
    fi_factors_sort(out);
    return ok;
}

/* ------------------------------------------------------------------ */
/* Registration                                                        */
/* ------------------------------------------------------------------ */

static const char FACTORINTEGER_DOC[] =
    "FactorInteger[n] gives a list of the prime factors of the integer n, together\n"
    "\twith their exponents: {{p1, e1}, {p2, e2}, ...}.  For a negative n the first\n"
    "\tentry is {-1, 1}; for a rational n the denominator's primes have negative\n"
    "\texponents; a Gaussian integer or Gaussian rational is factored over Z[i]\n"
    "\t(unit first, primes in the first quadrant).\n"
    "FactorInteger[n, k] does a partial factorisation: at most k entries, the last\n"
    "\tof which may be an unfactored remainder.\n"
    "FactorInteger[n, Automatic] is FactorInteger[n].\n"
    "Options:\n"
    "\tGaussianIntegers -> True factors over the Gaussian integers.\n"
    "\tMethod -> m chooses the algorithm; m is a name below, Automatic, or\n"
    "\t{name, \"Sub\" -> v, ...} with suboptions.  Each method is also callable\n"
    "\tdirectly as FactorInteger`name[n, \"Sub\" -> v, ...] (see ?FactorInteger`name).\n"
    "Methods:\n"
    "\tAutomatic                  word engine below 2^64 (micro-ECM, rho, SQUFOF);\n"
    "\t                           above: trial division to 2^16, rho, micro-ECM,\n"
    "\t                           Fermat probe, p-1, p+1, then ECM levels, switching\n"
    "\t                           to the quadratic sieve for composites up to 100 digits\n"
    "\t\"TrialDivision\"            batched trial division        (\"Bound\")\n"
    "\t\"PollardRho\"               Pollard-Brent rho             (\"Iterations\", \"Starts\")\n"
    "\t\"Fermat\"                   Fermat, residue-sieved        (\"Iterations\")\n"
    "\t\"Lehman\"                   Lehman, O(n^(1/3))            (\"Iterations\")\n"
    "\t\"ShanksSquareForms\"        SQUFOF; alias \"SQUFOF\"        (\"Iterations\")\n"
    "\t\"CFRAC\"                    continued fractions           (\"Iterations\")\n"
    "\t\"Dixon\"                    random squares                (\"Iterations\")\n"
    "\t\"QuadraticSieve\"           SIQS; aliases \"QS\", \"SIQS\"\n"
    "\t\"PollardPMinus1\"           p-1; alias \"PollardP-1\"       (\"B1\", \"B2\")\n"
    "\t\"WilliamsPPlus1\"           p+1; alias \"WilliamsP+1\"      (\"B1\", \"B2\")\n"
    "\t\"ECM\"                      elliptic curves (GMP-ECM)     (\"B1\", \"B2\", \"Curves\", \"Sigma\")\n"
    "\t\"BlakeRationalBaseDescent\" rational base descent         (\"Base\")\n"
    "A composite no factor of which is found within a method's bounds is returned\n"
    "\twith exponent 1 and FactorInteger::nofac.";

static void reg(const char* name, const char* doc) {
    symtab_add_builtin(name, builtin_factorinteger_method);
    symtab_get_def(name)->attributes |= ATTR_PROTECTED | ATTR_LISTABLE;
    symtab_set_docstring(name, doc);
}

void fi_init(void) {
    symtab_add_builtin("FactorInteger", builtin_factorinteger);
    symtab_get_def("FactorInteger")->attributes |= ATTR_PROTECTED | ATTR_LISTABLE;
    symtab_set_docstring("FactorInteger", FACTORINTEGER_DOC);
    char buf[96];
    for (size_t i = 0; i < NMETHODS; i++) {
        snprintf(buf, sizeof buf, "FactorInteger`%s", METHODS[i].name);
        reg(buf, METHODS[i].doc);
    }
    for (size_t i = 0; i < sizeof SYM_ALIASES / sizeof SYM_ALIASES[0]; i++) {
        snprintf(buf, sizeof buf, "FactorInteger`%s", SYM_ALIASES[i].sym);
        reg(buf, method_info(SYM_ALIASES[i].m)->doc);
    }
}
