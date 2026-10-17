/*
 * fi_rbd.c — FactorInteger`BlakeRationalBaseDescent: the Mathilda glue
 * around the reentrant core of fi_rbd.h.
 *
 * Variants (the "Variant" suboption; each also has an alias head
 * FactorInteger`BlakeRationalBaseDescent<Variant>):
 *
 *   Exact        scaled exact descent ([RBD] Section 5.1) with deferred
 *                testing and batched gcds.  The default.
 *   Classic      the paper's reference algorithm R1-R4, symmetric window.
 *   AllHits      the reference descent run to the end, every hit recorded
 *                and refined into pairwise-coprime divisors ([RBD] Section
 *                6, truncation divisibility); the pieces are cached so the
 *                driver's later calls on the cofactors split instantly.
 *   Lattice      lattice-augmented descent ([LRBD] L1-L4): the
 *                Howgrave-Graham detector at every depth of the band
 *                X <~ Q <~ sqrt(cmax N) + X.
 *   TwoPass      [LRBD] Section 6: P1 = Exact over the enumerated bases,
 *                then P2 = Lattice over a coarse grid of "BaseRange" (or
 *                over the same bases when no range is given).
 *   Coppersmith  known a/b, c and n: one detector call at floor(c a^n/b^n)
 *                ([LRBD] Remark "the fixed-base case is Coppersmith").
 *
 * Bases come from "Base" (one rational or a list), else from the
 * Stern-Brocot enumeration to height "Height" (or the a + b order), within
 * a budget of "Bases".  The search over many bases (or, for a single
 * lattice base, over the depths of the band) runs on fi_nthreads() workers
 * claiming chunks of indices; a hit lowers the claim's end, and the hit with
 * the LOWEST index wins, so the answer — and which base found it — is the
 * same for every thread count.
 */
#include "fi_internal.h"
#include "fi_rbd.h"

#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

/* ------------------------------------------------------------------ */
/* Options                                                             */
/* ------------------------------------------------------------------ */

#define RBD_DEFAULT_HEIGHT    1000
#define RBD_DEFAULT_BASES     200000
#define RBD_DEFAULT_BETA      0.5
#define RBD_DEFAULT_DIM       12
#define RBD_P2_HEIGHT         16
/* Bases generated per parallel region: short enough that a TimeConstrained
 * abort is honoured promptly, long enough to amortise thread start-up. */
#define RBD_BLOCK             2048

void fi_rbd_opts_init(FiRbdOpts* o) {
    memset(o, 0, sizeof *o);
    mpz_inits(o->cmax, o->mult, o->tol, o->coef, NULL);
    mpq_inits(o->lo, o->hi, NULL);
}

void fi_rbd_opts_reset_bases(FiRbdOpts* o) {
    for (size_t i = 0; i < o->nbases; i++) { mpz_clear(o->ba[i]); mpz_clear(o->bb[i]); }
    free(o->ba);
    free(o->bb);
    o->ba = o->bb = NULL;
    o->nbases = 0;
}

static void drop_cache(FiRbdOpts* o) {
    if (o->cache) fi_rbd_free_mpz_array(o->cache, o->ncache);
    o->cache = NULL;
    o->ncache = 0;
}

void fi_rbd_opts_clear(FiRbdOpts* o) {
    fi_rbd_opts_reset_bases(o);
    drop_cache(o);
    mpz_clears(o->cmax, o->mult, o->tol, o->coef, NULL);
    mpq_clears(o->lo, o->hi, NULL);
}

void fi_rbd_opts_add_base(FiRbdOpts* o, const mpz_t p, const mpz_t q) {
    mpz_t* na = realloc(o->ba, (o->nbases + 1) * sizeof(mpz_t));
    if (!na) return;
    o->ba = na;
    mpz_t* nb = realloc(o->bb, (o->nbases + 1) * sizeof(mpz_t));
    if (!nb) return;
    o->bb = nb;
    /* a/b > 1 in lowest terms: a base below 1 describes the same powers
     * inverted, so it is flipped. */
    mpz_init(o->ba[o->nbases]);
    mpz_init(o->bb[o->nbases]);
    mpz_t g;
    mpz_init(g);
    mpz_gcd(g, p, q);
    if (mpz_cmp(p, q) > 0) { mpz_divexact(o->ba[o->nbases], p, g); mpz_divexact(o->bb[o->nbases], q, g); }
    else                   { mpz_divexact(o->ba[o->nbases], q, g); mpz_divexact(o->bb[o->nbases], p, g); }
    mpz_clear(g);
    o->nbases++;
}

/* ------------------------------------------------------------------ */
/* Per-call setup                                                      */
/* ------------------------------------------------------------------ */

typedef struct {
    const FiRbdOpts* o;
    mpz_srcptr n;
    int        variant;     /* the kernel run per base (P1/P2 resolved) */
    mpz_t      T;           /* deferred-testing ceiling (exact window) */
    int        have_T;
    mpz_t      X;           /* lattice tolerance */
    int        t, d;
    int        lattice_ok;
} Setup;

static void setup_init(Setup* S, const mpz_t n, const FiRbdOpts* o, int variant) {
    S->o = o;
    S->n = n;
    S->variant = variant;
    mpz_inits(S->T, S->X, NULL);
    mpz_srcptr mult = o->mult_set ? o->mult : NULL;
    S->have_T = !(o->cmax_set && o->cmax_inf);
    if (S->have_T) {
        mpz_t c;
        mpz_init(c);
        if (o->cmax_set) mpz_set(c, o->cmax);
        else             mpz_ui_pow_ui(c, 2, 64);
        fi_rbd_ceiling(S->T, n, mult, c, 4);
        mpz_clear(c);
    }
    /* Lattice parameters: a given tolerance picks the least dimension that
     * covers it (unless t, d are both given); otherwise the dimension (12 by
     * default) picks the multiplicity and the largest tolerance its LLL
     * guarantee covers, less two bits of margin. */
    S->lattice_ok = 0;
    S->t = o->t;
    S->d = o->d;
    double L = (double)mpz_sizeinbase(n, 2);
    double beta = o->beta > 0 ? o->beta : RBD_DEFAULT_BETA;
    if (o->tol_set) {
        mpz_set(S->X, o->tol);
        if (S->t > 0 && S->d > S->t) S->lattice_ok = 1;
        else {
            double lx = (double)mpz_sizeinbase(o->tol, 2);
            int t = 0, d = 0;
            if (fi_rbd_hg_params(&t, &d, NULL, L, beta, lx)) {
                S->t = t; S->d = d; S->lattice_ok = 1;
            } else {
                /* Beyond any guarantee: still try the requested size at the
                 * paper's (6, 12) — the LLL bound is pessimistic. */
                S->t = S->t > 0 ? S->t : 6;
                S->d = S->d > S->t ? S->d : 12;
                S->lattice_ok = 1;
            }
        }
    } else {
        int t = 0, d = S->d > 1 ? S->d : RBD_DEFAULT_DIM;
        double lx = 0;
        if (fi_rbd_hg_params(&t, &d, &lx, L, beta, -1.0)) {
            if (S->t <= 0 || S->t >= d) S->t = t;
            S->d = d;
            long bits = (long)floor(lx) - 2;
            if (bits < 1) bits = 1;
            mpz_ui_pow_ui(S->X, 2, (unsigned long)bits);
            S->lattice_ok = 1;
        }
    }
}

static void setup_clear(Setup* S) { mpz_clears(S->T, S->X, NULL); }

/* ------------------------------------------------------------------ */
/* One base                                                            */
/* ------------------------------------------------------------------ */

typedef struct {
    fi_rbd_scratch_t s;
    fi_rbd_hits_t    hits;
    mpz_t            x0, lo, hi;
} Work;

static void work_init(Work* w) {
    fi_rbd_scratch_init(&w->s);
    fi_rbd_hits_init(&w->hits);
    mpz_inits(w->x0, w->lo, w->hi, NULL);
}

static void work_clear(Work* w) {
    fi_rbd_scratch_clear(&w->s);
    fi_rbd_hits_clear(&w->hits);
    mpz_clears(w->x0, w->lo, w->hi, NULL);
}

/* Run the variant's kernel on base a/b.  For AllHits, w->hits holds the
 * record afterwards (and g its first hit). */
static int try_base(mpz_t g, const Setup* S, const mpz_t a, const mpz_t b, Work* w) {
    fi_rbd_opts_t o;
    memset(&o, 0, sizeof o);
    o.mult = S->o->mult_set ? S->o->mult : NULL;
    switch (S->variant) {
    case FI_RBD_EXACT:
        o.exact = 1;
        o.batch = 1;
        o.ceiling = S->have_T ? S->T : NULL;
        return fi_rbd_descend(g, NULL, NULL, S->n, a, b, &o, NULL, &w->s);
    case FI_RBD_CLASSIC:
        return fi_rbd_descend(g, NULL, NULL, S->n, a, b, &o, NULL, &w->s);
    case FI_RBD_ALLHITS:
        o.batch = 1;
        o.record_all = 1;
        fi_rbd_hits_clear(&w->hits);
        return fi_rbd_descend(g, NULL, NULL, S->n, a, b, &o, &w->hits, &w->s) > 0;
    default:
        return 0;
    }
}

/* ------------------------------------------------------------------ */
/* Threaded search over indexed items                                  */
/* ------------------------------------------------------------------ */

/* A search over items 0..count (bases, or depths of one band).  Workers
 * claim chunks; the first hit in a worker caps the claim at its index, so
 * nothing past it is handed out, and every index below the final minimum
 * is examined by someone — the minimum is thread-count independent. */
typedef struct Search Search;
typedef int (*ItemFn)(Search* sh, size_t idx, mpz_t g, Work* w);

struct Search {
    const Setup* S;
    FiClaim      claim;
    size_t       chunk;
    ItemFn       item;
    /* Items: bases (a[i], b[i]) or band depths Q[i]. */
    mpz_t*       a;
    mpz_t*       b;
    mpz_t*       Q;
};

typedef struct {
    Search* sh;
    Work    w;
    size_t  best;           /* lowest hit index, or SIZE_MAX */
    mpz_t   g;
    mpz_t   tmp;
    fi_rbd_hits_t hits;     /* AllHits: the record of the best base */
} Job;

static void search_worker(void* arg) {
    Job* J = (Job*)arg;
    Search* sh = J->sh;
    size_t lo, hi;
    while (fi_claim_next(&sh->claim, sh->chunk, &lo, &hi)) {
        for (size_t i = lo; i < hi; i++) {
            if (i >= J->best || i >= fi_claim_end(&sh->claim)) break;
            if (sh->item(sh, i, J->tmp, &J->w)) {
                J->best = i;
                mpz_set(J->g, J->tmp);
                if (sh->S->variant == FI_RBD_ALLHITS) {
                    fi_rbd_hits_t t = J->hits;
                    J->hits = J->w.hits;
                    J->w.hits = t;
                }
                fi_claim_cap(&sh->claim, i);
                break;
            }
        }
    }
}

static int item_base(Search* sh, size_t i, mpz_t g, Work* w) {
    return try_base(g, sh->S, sh->a[i], sh->b[i], w);
}

static int item_depth(Search* sh, size_t i, mpz_t g, Work* w) {
    return fi_rbd_hg_detect(g, w->x0, sh->S->n, sh->Q[i], sh->S->X, sh->S->t, sh->S->d);
}

/* Run items [0, count) on the workers; returns the winning index (or
 * SIZE_MAX) with g, and for AllHits the winning record in *hits_out. */
static size_t search_run(Search* sh, size_t count, mpz_t g, fi_rbd_hits_t* hits_out) {
    int nt = fi_nthreads();
    if ((size_t)nt > count) nt = (int)(count ? count : 1);
    Job jobs[FI_MAXTHREADS];
    for (int i = 0; i < nt; i++) {
        jobs[i].sh = sh;
        work_init(&jobs[i].w);
        jobs[i].best = (size_t)-1;
        mpz_init(jobs[i].g);
        mpz_init(jobs[i].tmp);
        fi_rbd_hits_init(&jobs[i].hits);
    }
    fi_claim_reset(&sh->claim, count);
    fi_run_parallel(search_worker, jobs, sizeof(Job), nt);
    size_t best = (size_t)-1;
    int bj = -1;
    for (int i = 0; i < nt; i++)
        if (jobs[i].best < best) { best = jobs[i].best; bj = i; }
    if (bj >= 0) {
        mpz_set(g, jobs[bj].g);
        if (hits_out) { fi_rbd_hits_t t = *hits_out; *hits_out = jobs[bj].hits; jobs[bj].hits = t; }
    }
    for (int i = 0; i < nt; i++) {
        work_clear(&jobs[i].w);
        mpz_clears(jobs[i].g, jobs[i].tmp, NULL);
        fi_rbd_hits_clear(&jobs[i].hits);
    }
    return best;
}

/* ------------------------------------------------------------------ */
/* Base sources                                                        */
/* ------------------------------------------------------------------ */

typedef struct {
    int           kind;     /* 0 explicit list, 1 Stern-Brocot, 2 sum, 3 grid */
    size_t        next;     /* explicit / grid: next index */
    size_t        emitted, budget;
    fi_rbd_sb_t   sb;
    fi_rbd_sum_t  sum;
    fi_rbd_base_t* wb;      /* word scratch */
    /* Grid: ascending fixed-point positions r = R / 2^prec. */
    mpz_t         R, Rhi, step_num;
    mpq_t         qlo, qhi;
    mpz_t         pa, pb, last_a, last_b;
    unsigned long prec;
    double        lstep_base;  /* log2 of the relative spacing, less log2 n */
    double        L;           /* log2 N^(1-beta) / c */
} Source;

static void source_init(Source* src, const Setup* S, int grid) {
    const FiRbdOpts* o = S->o;
    memset(src, 0, sizeof *src);
    src->budget = o->max_bases ? (size_t)o->max_bases : RBD_DEFAULT_BASES;
    src->wb = malloc(RBD_BLOCK * sizeof(fi_rbd_base_t));
    mpz_inits(src->R, src->Rhi, src->step_num, src->pa, src->pb, src->last_a, src->last_b, NULL);
    mpq_inits(src->qlo, src->qhi, NULL);
    if (o->nbases && !(grid && o->range_set)) { src->kind = 0; return; }
    if (grid && o->range_set) {
        /* [LRBD] P2: relative spacing X / (2 c n N^beta) with n the depth at
         * which a factor of size N^(1-beta) is aligned for the base r:
         * n = log(N^(1-beta) / c) / log r.  c is "MaxCoefficient" when given
         * (the teeth narrow as c grows), else 1. */
        src->kind = 3;
        double beta = o->beta > 0 ? o->beta : RBD_DEFAULT_BETA;
        double Ln = (double)mpz_sizeinbase(S->n, 2);
        double lc = 0;
        if (o->cmax_set && !o->cmax_inf) lc = (double)mpz_sizeinbase(o->cmax, 2);
        double lX = (double)mpz_sizeinbase(S->X, 2) - 1;
        src->lstep_base = lX - 1.0 - lc - beta * Ln;
        src->L = (1.0 - beta) * Ln - lc;
        double minstep = src->lstep_base - log2(src->L > 1 ? src->L : 1) - 8;
        src->prec = (unsigned long)(-minstep > 64 ? -minstep + 64 : 128);
        /* R = ceil(lo 2^prec), Rhi = floor(hi 2^prec) */
        mpz_mul_2exp(src->R, mpq_numref(o->lo), src->prec);
        mpz_cdiv_q(src->R, src->R, mpq_denref(o->lo));
        mpz_mul_2exp(src->Rhi, mpq_numref(o->hi), src->prec);
        mpz_fdiv_q(src->Rhi, src->Rhi, mpq_denref(o->hi));
        return;
    }
    /* The lattice pass of TwoPass defaults to a much lower height: a
     * lattice band costs one LLL reduction per depth, and a base near 1 has
     * O(a log N) depths. */
    uint64_t H = o->height ? o->height : grid ? RBD_P2_HEIGHT : RBD_DEFAULT_HEIGHT;
    if (o->order == 1) { src->kind = 2; fi_rbd_sum_init(&src->sum, o->height); }
    else               { src->kind = 1; fi_rbd_sb_init(&src->sb, H); }
}

static void source_clear(Source* src) {
    if (src->kind == 1) fi_rbd_sb_clear(&src->sb);
    free(src->wb);
    mpz_clears(src->R, src->Rhi, src->step_num, src->pa, src->pb, src->last_a, src->last_b, NULL);
    mpq_clears(src->qlo, src->qhi, NULL);
}

/* z = v for a 64-bit v (unsigned long is 32 bits on LLP64). */
static void set_u64(mpz_t z, uint64_t v) {
    if (v <= 0xffffffffu) { mpz_set_ui(z, (unsigned long)v); return; }
    mpz_set_ui(z, (unsigned long)(v >> 32));
    mpz_mul_2exp(z, z, 32);
    mpz_add_ui(z, z, (unsigned long)(v & 0xffffffffu));
}

/* Fill a[0..max), b[0..max) (initialised) with the next bases; returns the
 * count, 0 when the source or the budget is exhausted. */
static size_t source_next(Source* src, const Setup* S, mpz_t* a, mpz_t* b, size_t max) {
    const FiRbdOpts* o = S->o;
    size_t left = src->budget - src->emitted;
    if (max > left) max = left;
    size_t n = 0;
    switch (src->kind) {
    case 0:
        while (n < max && src->next < o->nbases) {
            mpz_set(a[n], o->ba[src->next]);
            mpz_set(b[n], o->bb[src->next]);
            src->next++;
            n++;
        }
        break;
    case 1:
    case 2: {
        size_t got = src->kind == 1 ? fi_rbd_sb_next(&src->sb, src->wb, max)
                                    : fi_rbd_sum_next(&src->sum, src->wb, max);
        for (size_t i = 0; i < got; i++) {
            set_u64(a[i], src->wb[i].a);
            set_u64(b[i], src->wb[i].b);
        }
        n = got;
        break;
    }
    case 3:
        while (n < max && mpz_cmp(src->R, src->Rhi) <= 0) {
            /* Local spacing delta = 2^(lstep_base - log2 n(r)). */
            signed long e;
            double m = mpz_get_d_2exp(&e, src->R);
            double lr = log2(m) + (double)(e - (signed long)src->prec);   /* log2 r */
            double nn = lr > 0 ? src->L / lr : src->L;
            if (nn < 1) nn = 1;
            double ld = src->lstep_base - log2(nn);
            /* half-width w = R * delta / 2 in fixed point */
            long sh = (long)floor(ld) - 1;
            mpz_t w;
            mpz_init(w);
            if (sh >= 0) mpz_mul_2exp(w, src->R, (unsigned long)sh);
            else         mpz_fdiv_q_2exp(w, src->R, (unsigned long)(-sh));
            if (mpz_sgn(w) == 0) mpz_set_ui(w, 1);
            mpz_sub(mpq_numref(src->qlo), src->R, w);
            mpz_set_ui(mpq_denref(src->qlo), 1);
            mpz_mul_2exp(mpq_denref(src->qlo), mpq_denref(src->qlo), src->prec);
            mpq_canonicalize(src->qlo);
            mpz_add(mpq_numref(src->qhi), src->R, w);
            mpz_set_ui(mpq_denref(src->qhi), 1);
            mpz_mul_2exp(mpq_denref(src->qhi), mpq_denref(src->qhi), src->prec);
            mpq_canonicalize(src->qhi);
            fi_rbd_simplest_between(src->pa, src->pb, src->qlo, src->qhi);
            mpz_mul_2exp(w, w, 1);
            mpz_add(src->R, src->R, w);           /* next grid point */
            mpz_clear(w);
            if (mpz_cmp(src->pa, src->pb) <= 0) continue;
            if (mpz_cmp(src->pa, src->last_a) == 0 && mpz_cmp(src->pb, src->last_b) == 0) continue;
            mpz_set(src->last_a, src->pa);
            mpz_set(src->last_b, src->pb);
            mpz_set(a[n], src->pa);
            mpz_set(b[n], src->pb);
            n++;
        }
        break;
    }
    src->emitted += n;
    return n;
}

/* The lattice descent of one base, parallel over the depths of its band.
 * Each region covers two depths per thread, so a TimeConstrained abort is
 * taken within a few LLL reductions; regions run shallowest first and the
 * shallowest success wins, as in the serial descent. */
static int lattice_band_search(mpz_t d, const Setup* S, const mpz_t a, const mpz_t b) {
    if (!S->lattice_ok) return 0;
    Search bs;
    memset(&bs, 0, sizeof bs);
    bs.S = S;
    bs.item = item_depth;
    bs.chunk = 1;
    fi_claim_init(&bs.claim);
    mpz_t lo, hi;
    mpz_inits(lo, hi, NULL);
    mpz_set(lo, S->X);
    if (mpz_cmp_ui(lo, 2) < 0) mpz_set_ui(lo, 2);
    if (S->have_T) mpz_add(hi, S->T, S->X);
    else           mpz_set(hi, S->n);
    mpz_t* Qs;
    uint64_t* ks;
    size_t nb = fi_rbd_band(&Qs, &ks, S->n, S->o->mult_set ? S->o->mult : NULL, a, b, lo, hi);
    size_t step = 2 * (size_t)fi_nthreads();
    int found = 0;
    for (size_t off = 0; off < nb && !found; off += step) {
        bs.Q = Qs + off;
        size_t cnt = nb - off < step ? nb - off : step;
        found = search_run(&bs, cnt, d, NULL) != (size_t)-1;
    }
    fi_rbd_free_mpz_array(Qs, nb);
    free(ks);
    mpz_clears(lo, hi, NULL);
    fi_claim_destroy(&bs.claim);
    return found;
}

/* Search the bases of `src` with the setup's kernel.  Returns 1 with d set.
 * The descents over many bases run in parallel regions of RBD_BLOCK bases;
 * lattice bases run in order, each parallel over its band. */
static int search_bases(mpz_t d, Setup* S, Source* src, fi_rbd_hits_t* hits) {
    Search sh;
    memset(&sh, 0, sizeof sh);
    sh.S = S;
    sh.item = item_base;
    sh.chunk = 8;
    fi_claim_init(&sh.claim);
    sh.a = malloc(RBD_BLOCK * sizeof(mpz_t));
    sh.b = malloc(RBD_BLOCK * sizeof(mpz_t));
    for (size_t i = 0; i < RBD_BLOCK; i++) { mpz_init(sh.a[i]); mpz_init(sh.b[i]); }
    int found = 0;
    size_t cnt;
    while (!found && (cnt = source_next(src, S, sh.a, sh.b, RBD_BLOCK)) > 0) {
        if (S->variant == FI_RBD_LATTICE) {
            for (size_t i = 0; i < cnt && !found; i++) found = lattice_band_search(d, S, sh.a[i], sh.b[i]);
            continue;
        }
        found = search_run(&sh, cnt, d, hits) != (size_t)-1;
    }
    for (size_t i = 0; i < RBD_BLOCK; i++) { mpz_clear(sh.a[i]); mpz_clear(sh.b[i]); }
    free(sh.a);
    free(sh.b);
    fi_claim_destroy(&sh.claim);
    return found;
}

/* ------------------------------------------------------------------ */
/* Splitter                                                            */
/* ------------------------------------------------------------------ */

/* A cached piece that splits n. */
static int from_cache(mpz_t d, const mpz_t n, const FiRbdOpts* o) {
    for (size_t i = 0; i < o->ncache; i++) {
        mpz_gcd(d, n, o->cache[i]);
        if (mpz_cmp_ui(d, 1) > 0 && mpz_cmp(d, n) < 0) return 1;
    }
    return 0;
}

/* Store the coprime refinement of {n, hits} and return its least proper
 * piece in d. */
static int cache_hits(mpz_t d, const mpz_t n, FiRbdOpts* o, const fi_rbd_hits_t* hits) {
    size_t np;
    mpz_t* pcs = fi_rbd_coprime_base(hits, n, &np);
    if (!pcs) return 0;
    size_t old = o->ncache;
    mpz_t* nc = realloc(o->cache, (old + np) * sizeof(mpz_t));
    if (!nc) { fi_rbd_free_mpz_array(pcs, np); return 0; }
    o->cache = nc;
    for (size_t i = 0; i < np; i++) { mpz_init_set(o->cache[old + i], pcs[i]); }
    o->ncache = old + np;
    fi_rbd_free_mpz_array(pcs, np);
    return from_cache(d, n, o);
}

/* Coppersmith: P ~ floor(c a^n / b^n), one detector call. */
static int coppersmith(mpz_t d, const mpz_t n, const FiRbdOpts* o) {
    if (o->nbases < 1 || !o->coef_set || o->expo == 0) return 0;
    mpz_t Q, t, X, x0;
    mpz_inits(Q, t, X, x0, NULL);
    mpz_pow_ui(Q, o->ba[0], (unsigned long)o->expo);
    mpz_mul(Q, Q, o->coef);
    mpz_pow_ui(t, o->bb[0], (unsigned long)o->expo);
    mpz_fdiv_q(Q, Q, t);
    int found = 0;
    double L = (double)mpz_sizeinbase(n, 2);
    double lq = (double)mpz_sizeinbase(Q, 2);
    if (mpz_sgn(Q) > 0 && mpz_cmp(Q, n) < 0 && lq < L) {
        /* The factor near Q has size about N^(lq / L); a little margin for
         * the bit-length estimate. */
        double beta = o->beta > 0 ? o->beta : (lq - 1.0) / L;
        int tt = o->t, dd = o->d;
        if (o->tol_set) {
            mpz_set(X, o->tol);
            if (!(tt > 0 && dd > tt) &&
                !fi_rbd_hg_params(&tt, &dd, NULL, L, beta, (double)mpz_sizeinbase(X, 2))) {
                tt = 6; dd = 12;
            }
        } else {
            double lx = 0;
            dd = dd > 1 ? dd : RBD_DEFAULT_DIM;
            if (fi_rbd_hg_params(&tt, &dd, &lx, L, beta, -1.0)) {
                long bits = (long)floor(lx) - 2;
                mpz_ui_pow_ui(X, 2, (unsigned long)(bits > 1 ? bits : 1));
            } else {
                mpz_set_ui(X, 4);
                tt = 1; dd = 2;
            }
        }
        found = fi_rbd_hg_detect(d, x0, n, Q, X, tt, dd);
    }
    mpz_clears(Q, t, X, x0, NULL);
    return found;
}

/* The cofactor fallback: Exact over the Stern-Brocot bases of height at
 * most RBD_COFACTOR_HEIGHT.  Deliberately light: a cofactor that the
 * configured search could not split is usually unstructured (the leftover
 * cyclotomic parts of an AllHits run, say), and the full default budget would
 * only spend its time confirming that. */
#define RBD_COFACTOR_HEIGHT 64
#define RBD_COFACTOR_BASES  2000
static int default_search(mpz_t d, const mpz_t n) {
    FiRbdOpts def;
    fi_rbd_opts_init(&def);
    def.height = RBD_COFACTOR_HEIGHT;
    def.max_bases = RBD_COFACTOR_BASES;
    Setup S;
    setup_init(&S, n, &def, FI_RBD_EXACT);
    Source src;
    source_init(&src, &S, 0);
    int found = search_bases(d, &S, &src, NULL);
    source_clear(&src);
    setup_clear(&S);
    fi_rbd_opts_clear(&def);
    return found;
}

int fi_split_rbd(mpz_t d, const mpz_t n, FiCtx* c) {
    mpz_set_ui(d, 0);
    if (mpz_cmp_ui(n, 4) < 0) return 0;
    FiRbdOpts* o = &c->rbd;
    if (from_cache(d, n, o)) return 1;
    if (mpz_even_p(n)) { mpz_set_ui(d, 2); return 1; }

    /* The suboptions describe the INPUT (its base, its multiplier, its
     * coefficient): a cofactor left by the driver need not share that
     * structure, so when the configured search fails on one, the default
     * search (exact descent over the Stern-Brocot bases) is tried. */
    int cofactor = o->ncalls++ > 0;
    /* AllHits is exempt: its pieces are the answer, and re-running it on a
     * piece (which happens anyway) is what recovers nested algebraic
     * divisors; a generic search on huge cyclotomic leftovers is wasted. */
    int configured = o->variant != FI_RBD_ALLHITS &&
                     (o->nbases || o->mult_set || o->variant != FI_RBD_EXACT ||
                      o->height || o->max_bases || o->order || o->cmax_set);
    if (o->variant == FI_RBD_COPPERSMITH) {
        if (coppersmith(d, n, o)) return 1;
        return cofactor && default_search(d, n);
    }

    int found = 0;
    fi_rbd_hits_t hits;
    fi_rbd_hits_init(&hits);
    int variant = o->variant;
    int pass1 = variant == FI_RBD_TWOPASS ? FI_RBD_EXACT : variant;

    Setup S;
    setup_init(&S, n, o, pass1);
    Source src;
    source_init(&src, &S, 0);
    found = search_bases(d, &S, &src, &hits);
    source_clear(&src);

    if (found && variant == FI_RBD_ALLHITS) {
        fi_rbd_hits_t one;
        fi_rbd_hits_init(&one);
        if (hits.n == 0) fi_rbd_hits_add(&one, d, 0, 0);
        found = cache_hits(d, n, o, hits.n ? &hits : &one);
        fi_rbd_hits_clear(&one);
    }

    if (!found && variant == FI_RBD_TWOPASS) {
        /* P2: the lattice detector over the coarse grid (or the bases). */
        S.variant = FI_RBD_LATTICE;
        source_init(&src, &S, 1);
        found = search_bases(d, &S, &src, NULL);
        source_clear(&src);
    }

    setup_clear(&S);
    fi_rbd_hits_clear(&hits);
    if (!found && cofactor && configured) found = default_search(d, n);
    return found;
}
