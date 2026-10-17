/*
 * fi_smooth.c — threaded smooth-relation collection for CFRAC and Dixon.
 *
 * Both methods spend nearly all their time trial dividing candidates over
 * the factor base (fi_fb_factor), and every candidate is independent of the
 * others once it exists.  A region therefore hands the indices [0, count)
 * of a batch out in chunks (fi_claim_next), each worker trial-divides with
 * private scratch, and the hits go into a per-worker relation log tagged by
 * index.  The main thread then replays the logs merged by index, so the
 * relation store sees exactly the serial order: which relations are found,
 * when the matrix is solved and which divisor comes out do not depend on the
 * thread count or on scheduling.
 *
 * Chunks are claimed dynamically, not split statically, because candidate
 * costs are uneven (early abort drops most of them after a quarter of the
 * base) and because worker 0 may first run a `side` task — CFRAC generates
 * its next batch there, overlapping the sequential continued-fraction
 * recurrence with the trial division of the current batch.
 *
 * With one thread the same code runs on the calling thread, so the serial
 * and parallel paths cannot drift apart.
 */
#include "fi_internal.h"

#include <stdlib.h>

typedef struct {
    FiSmooth* S;
    int       id;
    mpz_t     X, v, tmp;
    uint32_t* col;
    uint32_t* exps;
} SmoothWorker;

struct FiSmooth {
    int          nt;
    FiClaim      claim;
    FiRelLog     log[FI_MAXTHREADS];
    SmoothWorker w[FI_MAXTHREADS];
    /* The current region. */
    size_t          chunk;
    FiSmoothPrep    prep;
    void          (*side)(void*);
    void*           user;
    const uint32_t* fb;
    size_t          nfb;
    uint64_t        lp_max;
    size_t          abort_bits;
};

FiSmooth* fi_smooth_new(int nthreads, size_t max_nfb) {
    FiSmooth* S = calloc(1, sizeof(FiSmooth));
    if (nthreads < 1) nthreads = 1;
    if (nthreads > FI_MAXTHREADS) nthreads = FI_MAXTHREADS;
    S->nt = nthreads;
    fi_claim_init(&S->claim);
    for (int t = 0; t < nthreads; t++) {
        SmoothWorker* w = &S->w[t];
        w->S = S;
        w->id = t;
        mpz_inits(w->X, w->v, w->tmp, NULL);
        w->col = malloc((max_nfb + 2) * sizeof(uint32_t));
        w->exps = malloc((max_nfb + 2) * sizeof(uint32_t));
    }
    /* Workers read the shared small-prime table: build it here, not lazily
     * inside the region. */
    size_t np;
    (void)fi_small_primes(&np);
    return S;
}

void fi_smooth_free(FiSmooth* S) {
    if (!S) return;
    for (int t = 0; t < S->nt; t++) {
        SmoothWorker* w = &S->w[t];
        mpz_clears(w->X, w->v, w->tmp, NULL);
        free(w->col);
        free(w->exps);
        fi_rellog_free(&S->log[t]);
    }
    fi_claim_destroy(&S->claim);
    free(S);
}

static void smooth_job(void* arg) {
    SmoothWorker* w = (SmoothWorker*)arg;
    FiSmooth* S = w->S;
    FiRelLog* lg = &S->log[w->id];
    if (w->id == 0 && S->side) S->side(S->user);
    size_t lo, hi;
    while (fi_claim_next(&S->claim, S->chunk, &lo, &hi)) {
        for (size_t i = lo; i < hi; i++) {
            size_t cnt = 0;
            uint64_t L;
            if (!S->prep(S->user, i, w->X, w->v, w->tmp, w->col, w->exps, &cnt)) continue;
            if (fi_fb_factor(w->v, S->fb, S->nfb, w->col, w->exps, &cnt, &L, S->lp_max, S->abort_bits))
                fi_rellog_push(lg, i, w->X, w->col, w->exps, cnt, L, 1);
        }
    }
}

void fi_smooth_run(FiSmooth* S, size_t count, size_t chunk, FiSmoothPrep prep,
                   void (*side)(void*), void* user, const uint32_t* fb, size_t nfb,
                   uint64_t lp_max, size_t abort_bits) {
    S->chunk = chunk ? chunk : 1;
    S->prep = prep;
    S->side = side;
    S->user = user;
    S->fb = fb;
    S->nfb = nfb;
    S->lp_max = lp_max;
    S->abort_bits = abort_bits;
    fi_claim_reset(&S->claim, count);
    /* Fewer candidates than threads: spare the thread starts (but keep
     * worker 0 when there is a side task). */
    int nt = S->nt;
    if ((size_t)nt > count / S->chunk + 1) nt = (int)(count / S->chunk + 1);
    if (nt < 1) nt = 1;
    fi_run_parallel(smooth_job, S->w, sizeof(SmoothWorker), nt);
}

int fi_smooth_replay(FiSmooth* S,
                     int (*visit)(void* user, const FiLogRec* r, const uint32_t* col,
                                  const uint32_t* exps),
                     void* user) {
    return fi_rellog_merge(S->log, S->nt, visit, user);
}
