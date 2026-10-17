/*
 * fi_core.c — factor list, context, and the recursive driver.
 *
 * fi_factor() runs a work stack of (composite, multiplicity, stage) items.
 * Every method shares the same pipeline:
 *
 *   1. strip powers of two (and, for Automatic / TrialDivision, every prime
 *      up to the trial bound, once, at the top);
 *   2. for each work item: word-sized → the u64 engine (Automatic only; a
 *      named method keeps its own algorithm at every size), prime → output,
 *      perfect power → its primitive root with the multiplicity scaled;
 *   3. otherwise hand it to the method's splitter, and push both parts.
 *
 * A composite the splitter cannot crack is emitted with its multiplicity,
 * FactorInteger::nofac is issued once, and c->incomplete is set.
 *
 * Automatic's splitter is a tier table, resumed per item through `stage`:
 * cheap Pollard rho, one Pollard p-1, then ECM in rising digit levels, with a
 * switch to the quadratic sieve once the ECM effort spent reaches 2/7 of the
 * composite's digit count (yafu's crossover rule) and the composite is small
 * enough for QS to be the faster finisher.  A part produced by a split keeps
 * its parent's stage: the parent's search already ran over it.
 */
#include "fi_internal.h"

#include <stdlib.h>
#include <string.h>

/* ------------------------------------------------------------------ */
/* Factor list                                                         */
/* ------------------------------------------------------------------ */

void fi_factors_init(FiFactors* f) { f->v = NULL; f->n = f->cap = 0; }

void fi_factors_clear(FiFactors* f) {
    for (size_t i = 0; i < f->n; i++) mpz_clear(f->v[i].p);
    free(f->v);
    f->v = NULL;
    f->n = f->cap = 0;
}

static FiEntry* grow(FiFactors* f) {
    if (f->n == f->cap) {
        f->cap = f->cap ? 2 * f->cap : 16;
        f->v = realloc(f->v, f->cap * sizeof(FiEntry));
    }
    return &f->v[f->n++];
}

void fi_factors_add(FiFactors* f, const mpz_t p, int64_t e) {
    for (size_t i = 0; i < f->n; i++)
        if (mpz_cmp(f->v[i].p, p) == 0) { f->v[i].e += e; return; }
    FiEntry* x = grow(f);
    mpz_init_set(x->p, p);
    x->e = e;
}

void fi_factors_add_ui(FiFactors* f, uint64_t p, int64_t e) {
    for (size_t i = 0; i < f->n; i++)
        if (mpz_cmp_ui(f->v[i].p, (unsigned long)p) == 0) { f->v[i].e += e; return; }
    FiEntry* x = grow(f);
    mpz_init_set_ui(x->p, (unsigned long)p);
    x->e = e;
}

static int cmp_entry(const void* a, const void* b) {
    return mpz_cmp(((const FiEntry*)a)->p, ((const FiEntry*)b)->p);
}
void fi_factors_sort(FiFactors* f) {
    if (f->n > 1) qsort(f->v, f->n, sizeof(FiEntry), cmp_entry);
}

/* ------------------------------------------------------------------ */
/* Context                                                             */
/* ------------------------------------------------------------------ */

void fi_ctx_init(FiCtx* c, FiMethod m) {
    memset(c, 0, sizeof *c);
    c->method = m;
    mpz_init(c->rbd_a);
    mpz_init(c->rbd_b);
    c->rng = UINT64_C(0x9E3779B97F4A7C15);
}

void fi_ctx_clear(FiCtx* c) {
    mpz_clear(c->rbd_a);
    mpz_clear(c->rbd_b);
}

/* The randomised methods (rho starts, ECM sigmas, QS) draw from a generator
 * seeded by n itself, so a given input always takes the same path: the
 * answer, the time and any nofac depend on n alone, never on what the
 * process factored before (MATHILDA_DIVERGENCES F9). */
void fi_ctx_seed(FiCtx* c, const mpz_t n) {
    uint64_t h = UINT64_C(0x9E3779B97F4A7C15);
    size_t sz = mpz_size(n);
    for (size_t i = 0; i < sz; i++) {
        h ^= (uint64_t)mpz_getlimbn(n, (mp_size_t)i);
        h *= UINT64_C(0xff51afd7ed558ccd);
        h ^= h >> 33;
    }
    c->rng = h ? h : 1;
}

uint64_t fi_rand(FiCtx* c) {
    uint64_t x = c->rng;
    x ^= x >> 12; x ^= x << 25; x ^= x >> 27;
    c->rng = x;
    return x * UINT64_C(0x2545F4914F6CDD1D);
}

int fi_is_prime(const mpz_t n) {
    if (mpz_sizeinbase(n, 2) <= 64) return fi_u64_is_prime(mpz_get_ui(n));
    return mpz_probab_prime_p(n, 25) != 0;
}

/* ------------------------------------------------------------------ */
/* Work stack                                                          */
/* ------------------------------------------------------------------ */

typedef struct {
    mpz_t   m;
    int64_t mult;
    int     stage;      /* Automatic: next tier to run */
    int     qs_tried;   /* Automatic: QS already attempted on this item */
} Work;

typedef struct { Work* v; size_t n, cap; } WorkStack;

static void push(WorkStack* s, const mpz_t m, int64_t mult, int stage) {
    if (s->n == s->cap) {
        s->cap = s->cap ? 2 * s->cap : 16;
        s->v = realloc(s->v, s->cap * sizeof(Work));
    }
    Work* w = &s->v[s->n++];
    mpz_init_set(w->m, m);
    w->mult = mult;
    w->stage = stage;
    w->qs_tried = 0;
}

/* Reduce m to its primitive root: m = b^e with b not a perfect power.
 * Returns e (1 when m is not a perfect power).  Restarting the exponent scan
 * after every hit accumulates the full exponent (p^6 -> p via 2 then 3). */
static int64_t primitive_root(mpz_t b, const mpz_t m) {
    mpz_set(b, m);
    if (!mpz_perfect_power_p(b)) return 1;
    int64_t e = 1;
    mpz_t r;
    mpz_init(r);
    for (int reduced = 1; reduced; ) {
        reduced = 0;
        unsigned long emax = (unsigned long)mpz_sizeinbase(b, 2);
        for (unsigned long k = 2; k <= emax; k++) {
            if (mpz_root(r, b, k)) {
                mpz_set(b, r);
                e *= (int64_t)k;
                reduced = 1;
                break;
            }
        }
    }
    mpz_clear(r);
    return e;
}

/* ------------------------------------------------------------------ */
/* Named-method splitters                                              */
/* ------------------------------------------------------------------ */

#define DEF(v, d) ((v) ? (v) : (d))

static int split_named(mpz_t d, const mpz_t m, FiCtx* c) {
    int word = mpz_sizeinbase(m, 2) <= 64;
    uint64_t w = word ? mpz_get_ui(m) : 0, r = 0;
    switch (c->method) {
    case FI_RHO:
        if (word) {
            uint64_t starts = DEF(c->starts, 16), it = DEF(c->iterations, UINT64_C(1) << 26);
            for (uint64_t k = 0; k < starts && !r; k++) r = fi_u64_rho(w, 1 + k, it);
            break;
        }
        return fi_split_rho(d, m, c, c->iterations, c->starts);
    case FI_FERMAT:
        if (word && !(w >> 62)) { r = fi_u64_fermat(w, DEF(c->iterations, UINT64_C(1) << 24)); break; }
        return fi_split_fermat(d, m, c, c->iterations);
    case FI_LEHMAN:
        if (word && !(w >> 42)) { r = fi_u64_lehman(w); break; }
        return fi_split_lehman(d, m, c, c->iterations);
    case FI_SQUFOF:
        if (word && !(w >> 62)) { r = fi_u64_squfof(w, c->iterations); if (r) break; }
        return fi_split_squfof(d, m, c, c->iterations);
    case FI_CFRAC:  return fi_split_cfrac(d, m, c);
    case FI_DIXON:  return fi_split_dixon(d, m, c);
    case FI_QS:     return fi_split_qs(d, m, c);
    case FI_PM1:    return fi_split_pm1(d, m, c, c->B1, c->B2);
    case FI_PP1:    return fi_split_pp1(d, m, c, c->B1, c->B2);
    case FI_ECM:    return fi_split_ecm(d, m, c, c->B1, c->curves);
    case FI_RBD:    return fi_split_rbd(d, m, c);
    default:        return 0;
    }
    if (r > 1 && r < w) { mpz_set_ui(d, (unsigned long)r); return 1; }
    return 0;
}

/* ------------------------------------------------------------------ */
/* Automatic tier table                                                */
/* ------------------------------------------------------------------ */

/* Trial-division bound applied once at the top of an Automatic run. */
#define AUTO_TRIAL_BOUND   (1u << 16)
/* Cheap rho budget: finds factors up to ~2^24 before p-1 / ECM. */
#define AUTO_RHO_ITERS     (1u << 13)
/* Cheap Fermat probe: catches two factors within ~n^(1/4) * 2^7 of each
 * other (a common construction, and the one rho/ECM/QS are worst at). */
#define AUTO_FERMAT_ITERS  (1u << 14)
/* One Pollard p-1 and one Williams p+1 pass (GMP-ECM picks stage-2 B2). */
#define AUTO_PM1_B1        150000.0
#define AUTO_PP1_B1        100000.0
/* QS is the finisher for composites up to this many digits. */
#define AUTO_QS_MAX_DIGITS 100
/* ECM digit levels: 15, 20, ..., stopping after AUTO_ECM_MAX. */
#define AUTO_ECM_FIRST     15
#define AUTO_ECM_STEP      5
#define AUTO_ECM_MAX       65

enum { ST_RHO = 0, ST_TINY, ST_FERMAT, ST_PM1, ST_PP1, ST_ECM0 };

/* Run tier `s` once on m.  Stage s >= ST_ECM0 is ECM level
 * AUTO_ECM_FIRST + (s - ST_ECM0) * STEP.  The p-+1 bounds scale with the
 * size of m so a small composite is not charged a large-number search. */
static int run_tier(int s, mpz_t d, const mpz_t m, size_t bits, int digits, FiCtx* c) {
    switch (s) {
    case ST_RHO:
        return bits < 126 ? fi_u128_rho(d, m, 1, AUTO_RHO_ITERS)
                          : fi_split_rho(d, m, c, AUTO_RHO_ITERS, 1);
    case ST_TINY: {
        /* Two-limb micro-ECM, a light pass: 8 curves at B1 = 150 catch the
         * 9..12 digit factors rho misses for ~0.25 ms.  Measured: a heavier
         * pass (B1 = 300..1000) loses to QS on balanced 96-bit semiprimes. */
        if (bits > 126) return 0;
        uint32_t B1 = 150;
        int curves = 8;
        FI_TUNE_U("FI_TINY_B1", B1);
        FI_TUNE_U("FI_TINY_CURVES", curves);
        return fi_u128_ecm(d, m, &c->rng, B1, curves);
    }
    case ST_FERMAT: return fi_split_fermat(d, m, c, AUTO_FERMAT_ITERS);
    case ST_PM1: {
        double b1 = 2500.0 * digits;
        return fi_split_pm1(d, m, c, b1 < AUTO_PM1_B1 ? b1 : AUTO_PM1_B1, 0);
    }
    case ST_PP1: {
        /* Each seed succeeds with probability 1/2 (it tests p+1 or p-1
         * depending on a Legendre symbol), so composites big enough that ECM
         * or QS will cost seconds get three seeds; smaller ones one probe. */
        double b1 = 1500.0 * digits;
        return fi_split_pp1_seeds(d, m, c, b1 < AUTO_PP1_B1 ? b1 : AUTO_PP1_B1, 0, digits >= 60 ? 3 : 1);
    }
    default:        return fi_ecm_level(d, m, c, AUTO_ECM_FIRST + (s - ST_ECM0) * AUTO_ECM_STEP);
    }
}

/* Search effort spent once tiers [0, s) have run, as a digit count: the
 * size of factor they would have found with good probability.  Calibrated
 * so that QS (whose cost at 40 digits is about one p-1 pass at B1 = 10^5)
 * takes over directly after the Fermat probe below ~45 digits, after p-1
 * below ~52, after p+1 below ~56, and after ECM level L below 3.5 L. */
static int effort_digits(int s) {
    switch (s) {
    case ST_RHO:    return 0;
    case ST_TINY:   return 7;
    case ST_FERMAT: return 7;
    case ST_PM1:    return 13;
    case ST_PP1:    return 15;
    case ST_ECM0:   return 16;
    default: {
        int done = AUTO_ECM_FIRST + (s - ST_ECM0 - 1) * AUTO_ECM_STEP;
        return done > 16 ? done : 16;
    }
    }
}

/* Returns 1 with d set, 0 when the whole table is exhausted.  The stage only
 * advances when a tier FAILS: after a split both parts re-enter at the tier
 * that succeeded, so a composite rich in small factors keeps peeling them
 * off with cheap rho instead of falling through to ECM / QS.  Before every
 * tier past Fermat, QS takes over once the effort spent reaches 2/7 of the
 * composite's digit count (yafu's crossover): beyond that point a factor is
 * more cheaply found by sieving than by searching further. */
static int split_auto(mpz_t d, Work* w, FiCtx* c) {
    const mpz_t* mp = (const mpz_t*)&w->m;
    size_t bits = mpz_sizeinbase(*mp, 2);
    int digits = (int)(bits * 0.30103) + 1;
    for (;;) {
        int s = w->stage;
        if (s >= ST_PM1 && digits <= AUTO_QS_MAX_DIGITS && !w->qs_tried &&
            7 * effort_digits(s) >= 2 * digits) {
            w->qs_tried = 1;
            if (fi_split_qs(d, *mp, c)) return 1;
            /* QS failed (should not happen): keep searching with ECM. */
        }
        if (s >= ST_ECM0 && AUTO_ECM_FIRST + (s - ST_ECM0) * AUTO_ECM_STEP > AUTO_ECM_MAX) return 0;
        if (run_tier(s, d, *mp, bits, digits, c)) return 1;
        w->stage++;
    }
}

/* ------------------------------------------------------------------ */
/* Driver                                                              */
/* ------------------------------------------------------------------ */

void fi_warn_incomplete(FiCtx* c, const mpz_t m);

static size_t distinct_primes(const FiFactors* out) { return out->n; }

void fi_factor(const mpz_t n_in, FiFactors* out, FiCtx* c) {
    if (mpz_cmp_ui(n_in, 1) <= 0) return;
    fi_ctx_seed(c, n_in);

    /* FactorInteger[n, 1]: the answer is n itself. */
    if (c->k_limit == 1) { fi_factors_add(out, n_in, 1); return; }

    mpz_t n;
    mpz_init_set(n, n_in);

    /* Shared pre-step: powers of two; Automatic and TrialDivision also take
     * every prime up to their bound here, once. */
    if (c->method == FI_AUTOMATIC && mpz_sizeinbase(n, 2) > 64) {
        fi_trial_divide(n, AUTO_TRIAL_BOUND, 1, out);
    } else if (c->method == FI_TRIAL) {
        uint64_t B = c->bound;
        if (!B) {
            /* Default: complete trial division where it is cheap (sqrt(n)
             * below 10^7), else primes up to 10^7. */
            mpz_t s; mpz_init(s); mpz_sqrt(s, n);
            B = (mpz_cmp_ui(s, 10000000u) < 0) ? mpz_get_ui(s) : 10000000u;
            mpz_clear(s);
        }
        if (B > (UINT64_C(1) << 40)) B = UINT64_C(1) << 40;
        fi_trial_divide(n, B, 1, out);
        if (mpz_cmp_ui(n, 1) > 0) {
            if (fi_is_prime(n)) fi_factors_add(out, n, 1);
            else {
                mpz_t b; mpz_init(b);
                int64_t e = primitive_root(b, n);
                if (e > 1 && fi_is_prime(b)) fi_factors_add(out, b, e);
                else { fi_factors_add(out, n, 1); c->incomplete = 1; fi_warn_incomplete(c, n); }
                mpz_clear(b);
            }
        }
        mpz_clear(n);
        return;
    } else {
        mp_bitcnt_t z = mpz_scan1(n, 0);
        if (z) { mpz_tdiv_q_2exp(n, n, z); fi_factors_add_ui(out, 2, (int64_t)z); }
    }

    WorkStack st = { NULL, 0, 0 };
    if (mpz_cmp_ui(n, 1) > 0) push(&st, n, 1, 0);
    mpz_clear(n);

    mpz_t d, q, b;
    mpz_inits(d, q, b, NULL);
    mpz_t rest;
    mpz_init_set_ui(rest, 1);
    int truncated = 0;

    while (st.n) {
        Work w = st.v[--st.n];            /* take ownership of w.m */

        if (c->k_limit > 1 && distinct_primes(out) + 1 >= (size_t)c->k_limit &&
            !(st.n == 0 && !truncated && fi_is_prime(w.m))) {
            /* Partial factorisation: everything still pending becomes the
             * final entry, the unfactored remainder (a lone prime is still
             * emitted as itself).  The caller merges any overshoot. */
            mpz_t pw; mpz_init(pw);
            mpz_pow_ui(pw, w.m, (unsigned long)w.mult);
            mpz_mul(rest, rest, pw);
            mpz_clear(pw);
            mpz_clear(w.m);
            truncated = 1;
            continue;
        }

        if (c->method == FI_AUTOMATIC && mpz_sizeinbase(w.m, 2) <= 64) {
            fi_u64_factor(mpz_get_ui(w.m), w.mult, out);
            mpz_clear(w.m);
            continue;
        }
        if (fi_is_prime(w.m)) { fi_factors_add(out, w.m, w.mult); mpz_clear(w.m); continue; }
        int64_t e = primitive_root(b, w.m);
        if (e > 1) {
            push(&st, b, w.mult * e, w.stage);
            mpz_clear(w.m);
            continue;
        }

        int r;
        if (c->method == FI_AUTOMATIC) r = split_auto(d, &w, c);
        else r = split_named(d, w.m, c);

        if (r == 1 && mpz_cmp_ui(d, 1) > 0 && mpz_cmp(d, w.m) < 0 && mpz_divisible_p(w.m, d)) {
            mpz_divexact(q, w.m, d);
            /* Push the larger part first so the smaller is processed first:
             * small factors come out early, which is what FactorInteger[n,k]
             * wants. */
            if (mpz_cmp(d, q) < 0) { push(&st, q, w.mult, w.stage); push(&st, d, w.mult, w.stage); }
            else                   { push(&st, d, w.mult, w.stage); push(&st, q, w.mult, w.stage); }
        } else {
            fi_factors_add(out, w.m, w.mult);
            c->incomplete = 1;
            fi_warn_incomplete(c, w.m);
        }
        mpz_clear(w.m);
    }
    if (truncated && mpz_cmp_ui(rest, 1) > 0) fi_factors_add(out, rest, 1);
    free(st.v);
    mpz_clears(d, q, b, rest, NULL);
}
