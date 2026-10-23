/*
 * sparse_dot_kern.inc -- the typed Dot kernels, instantiated once per value
 * type by sparse_dot.c. Not a translation unit of its own.
 *
 * Parameters (defined by the includer, undefined at the end):
 *
 *   KT              element type: int64_t, double or SpCplx
 *   KS(name)        name ## _<type suffix>
 *   K_ZERO          the additive identity of KT
 *   K_FMA(acc,a,b)  acc += a * b; the int64 instance sets the local `ovf`
 *                   on overflow (every kernel declares `bool ovf = false`)
 *
 * Every matrix operand is CSR with 32-bit column indices (sp_dims_idx32 holds
 * for every rank-2 array whose dense partner fits in memory); dense operands
 * are row-major KT buffers. Each kernel returns false when the int64 instance
 * overflowed or an allocation failed, and the caller falls back to the exact
 * boxed path. Summation within one output element runs in storage order, so
 * an element's value does not depend on the thread count -- except for v . s,
 * whose per-thread accumulators are combined at the end.
 */

/* ---- s . v : y[i] = sum_k a[i,k] x[k] ------------------------------------ */

typedef struct {
    const int64_t* rp; const int32_t* ci; const KT* av;
    const KT* x; KT* y;
} KS(SpmvCtx);

static bool KS(spmv_rows)(void* c, int64_t r0, int64_t r1) {
    const KS(SpmvCtx)* t = c;
    bool ovf = false;
    const int64_t* rp = t->rp; const int32_t* ci = t->ci;
    const KT* av = t->av; const KT* x = t->x;
    for (int64_t i = r0; i < r1; i++) {
        KT acc = K_ZERO;
        for (int64_t k = rp[i]; k < rp[i + 1]; k++) K_FMA(acc, av[k], x[ci[k]]);
        t->y[i] = acc;
    }
    return !ovf;
}

static bool KS(spmv)(int64_t n, const int64_t* rp, const int32_t* ci, const KT* av,
                     const KT* x, KT* y) {
    KS(SpmvCtx) c = { rp, ci, av, x, y };
    return par_rows(n, rp, 1, KS(spmv_rows), &c);
}

/* ---- s . M : Y[i,:] = sum_k a[i,k] M[k,:]  (M is K x p) --------------------- */

typedef struct {
    const int64_t* rp; const int32_t* ci; const KT* av;
    const KT* m; KT* y; int64_t p;
} KS(SpmmCtx);

static bool KS(spmm_rows)(void* c, int64_t r0, int64_t r1) {
    const KS(SpmmCtx)* t = c;
    bool ovf = false;
    int64_t p = t->p;
    if (p <= 16) {
        /* a narrow dense operand: the output row in a local accumulator, and
         * the next gathered row requested ahead of its use */
        KT acc[16];
        for (int64_t i = r0; i < r1; i++) {
            for (int64_t j = 0; j < p; j++) acc[j] = K_ZERO;
            int64_t k1 = t->rp[i + 1];
            for (int64_t k = t->rp[i]; k < k1; k++) {
                KT a = t->av[k];
                const KT* mr = t->m + (int64_t)t->ci[k] * p;
#if defined(__GNUC__)
                if (k + 4 < k1) __builtin_prefetch(t->m + (int64_t)t->ci[k + 4] * p);
#endif
                for (int64_t j = 0; j < p; j++) K_FMA(acc[j], a, mr[j]);
            }
            KT* yr = t->y + i * p;
            for (int64_t j = 0; j < p; j++) yr[j] = acc[j];
        }
        return !ovf;
    }
    for (int64_t i = r0; i < r1; i++) {
        KT* yr = t->y + i * p;
        for (int64_t j = 0; j < p; j++) yr[j] = K_ZERO;
        for (int64_t k = t->rp[i]; k < t->rp[i + 1]; k++) {
            KT a = t->av[k];
            const KT* mr = t->m + (int64_t)t->ci[k] * p;
            for (int64_t j = 0; j < p; j++) K_FMA(yr[j], a, mr[j]);
        }
    }
    return !ovf;
}

static bool KS(spmm)(int64_t n, int64_t p, const int64_t* rp, const int32_t* ci,
                     const KT* av, const KT* m, KT* y) {
    KS(SpmmCtx) c = { rp, ci, av, m, y, p };
    /* a row costs its entries times the row length each one streams */
    return par_rows_w(n, rp, p, p, KS(spmm_rows), &c);
}

/* ---- v . s : y[j] = sum_i x[i] a[i,j] ------------------------------------- */
/* A scatter. The rows are cut into blocks of about equal entry count; each
 * block scatters into a private accumulator that spans only the columns the
 * block touches (a banded matrix -- a stencil -- gives every block a narrow
 * span), and the accumulators are then summed column-parallel. When the spans
 * are the full width (a random matrix) the block count is cut so the
 * accumulators cost no more than about half the scatter itself. */

typedef struct {
    int64_t r0, r1;          /* rows of the block */
    int64_t c0, c1;          /* column span [c0, c1) */
    KT* acc;
} KS(VBlock);

typedef struct {
    const int64_t* rp; const int32_t* ci; const KT* av;
    const KT* x;
    KS(VBlock)* blk;
    KT* region;              /* nb * p accumulators, one width-p slot per block */
    int64_t p;
} KS(VspmCtx);

static void KS(vblock_span)(const int64_t* rp, const int32_t* ci, KS(VBlock)* b) {
    int64_t lo = INT64_MAX, hi = -1;
    for (int64_t k = rp[b->r0]; k < rp[b->r1]; k++) {
        if (ci[k] < lo) lo = ci[k];
        if (ci[k] > hi) hi = ci[k];
    }
    b->c0 = hi < 0 ? 0 : lo;
    b->c1 = hi < 0 ? 0 : hi + 1;
}

static bool KS(vspm_blocks)(void* c, int64_t b0, int64_t b1) {
    KS(VspmCtx)* t = c;
    bool ovf = false;
    for (int64_t b = b0; b < b1; b++) {
        KS(VBlock)* bl = &t->blk[b];
        if (bl->c1 < 0) KS(vblock_span)(t->rp, t->ci, bl);
        int64_t w = bl->c1 - bl->c0;
        bl->acc = t->region + b * t->p;
        for (int64_t j = 0; j < w; j++) bl->acc[j] = K_ZERO;
        KT* y = bl->acc - bl->c0;
        for (int64_t i = bl->r0; i < bl->r1; i++) {
            KT xi = t->x[i];
            for (int64_t k = t->rp[i]; k < t->rp[i + 1]; k++) K_FMA(y[t->ci[k]], xi, t->av[k]);
        }
    }
    return !ovf;
}

typedef struct {
    const KS(VBlock)* blk; int64_t nb; KT* y;
} KS(VspmSum);

/* y[j0..j1) = the sum of the block accumulators over that range */
static bool KS(vspm_sum)(void* c, int64_t j0, int64_t j1) {
    const KS(VspmSum)* t = c;
    bool ovf = false;
    for (int64_t j = j0; j < j1; j++) t->y[j] = K_ZERO;
    for (int64_t b = 0; b < t->nb; b++) {
        const KS(VBlock)* bl = &t->blk[b];
        int64_t lo = bl->c0 > j0 ? bl->c0 : j0;
        int64_t hi = bl->c1 < j1 ? bl->c1 : j1;
        const KT* a = bl->acc - bl->c0;
        for (int64_t j = lo; j < hi; j++) K_FMA(t->y[j], K_ONE, a[j]);
    }
    return !ovf;
}

static bool KS(vspm)(int64_t nrows, int64_t p, const int64_t* rp, const int32_t* ci,
                     const KT* av, const KT* x, KT* y) {
    int64_t nnz = rp[nrows];
    int64_t nb = (nnz + nrows) / (int64_t)NDARRAY_THREAD_THRESHOLD;
    if (nb > NDARRAY_MAX_THREADS) nb = NDARRAY_MAX_THREADS;
    /* the accumulator region is nb width-p slots: keep it under 256 MB */
    int64_t nb_mem = ((int64_t)256 << 20) / ((int64_t)sizeof(KT) * (p ? p : 1));
    if (nb > nb_mem) nb = nb_mem;
    if (nb < 1) nb = 1;
    KS(VBlock) blk[NDARRAY_MAX_THREADS];
    /* blocks of equal entry count (keys rp[r] + r, as par_rows deals them) */
    RowPar key = { nrows, rp, 1, NULL, NULL, 1 };
    int64_t total = nnz + nrows;
    for (int64_t b = 0; b < nb; b++) {
        blk[b].r0 = b == 0 ? 0 : rp_lower(&key, total / nb * b);
        blk[b].r1 = b == nb - 1 ? nrows : rp_lower(&key, total / nb * (b + 1));
        blk[b].c1 = -1;
        blk[b].acc = NULL;
    }
    if (nb > 1) {
        /* a full-width span costs its width in zeroing and summing: keep the
         * block count where that stays under half the scatter */
        KS(vblock_span)(rp, ci, &blk[0]);
        int64_t w = blk[0].c1 - blk[0].c0;
        if (w > 0 && w * nb > nnz / 2) {
            int64_t nb2 = nnz / (2 * w);
            if (nb2 < 1) nb2 = 1;
            if (nb2 < nb) {
                nb = nb2;
                for (int64_t b = 0; b < nb; b++) {
                    blk[b].r0 = b == 0 ? 0 : rp_lower(&key, total / nb * b);
                    blk[b].r1 = b == nb - 1 ? nrows : rp_lower(&key, total / nb * (b + 1));
                    blk[b].c1 = -1;
                }
            }
        }
    }
    bool ok;
    if (nb == 1) {
        /* serial: scatter straight into y */
        bool ovf = false;
        for (int64_t j = 0; j < p; j++) y[j] = K_ZERO;
        for (int64_t i = 0; i < nrows; i++) {
            KT xi = x[i];
            for (int64_t k = rp[i]; k < rp[i + 1]; k++) K_FMA(y[ci[k]], xi, av[k]);
        }
        return !ovf;
    }
    /* the accumulators live in the reusable scratch region: a fresh
     * allocation would fault in every page of every block on every call */
    bool cached;
    KT* region = scratch_get(sizeof(KT) * (size_t)nb * (size_t)p, &cached);
    if (!region) return false;
    KS(VspmCtx) c = { rp, ci, av, x, blk, region, p };
    /* one block per thread: weight each block past the threading threshold */
    ok = par_rows(nb, NULL, (int64_t)NDARRAY_THREAD_THRESHOLD, KS(vspm_blocks), &c);
    if (ok) {
        KS(VspmSum) sc = { blk, nb, y };
        ok = par_rows(p, NULL, nb, KS(vspm_sum), &sc);
    }
    scratch_put(region, cached);
    return ok;
}

/* ---- M . s : Y[r,:] = sum_k M[r,k] a[k,:]  (M is n x K, s is K x p) -------- */

typedef struct {
    const int64_t* rp; const int32_t* ci; const KT* av;
    const KT* m; KT* y; int64_t K, p;
} KS(MspmCtx);

static bool KS(mspm_rows)(void* c, int64_t r0, int64_t r1) {
    const KS(MspmCtx)* t = c;
    bool ovf = false;
    for (int64_t r = r0; r < r1; r++) {
        KT* yr = t->y + r * t->p;
        const KT* mr = t->m + r * t->K;
        for (int64_t j = 0; j < t->p; j++) yr[j] = K_ZERO;
        for (int64_t k = 0; k < t->K; k++) {
            KT mk = mr[k];
            for (int64_t q = t->rp[k]; q < t->rp[k + 1]; q++) K_FMA(yr[t->ci[q]], mk, t->av[q]);
        }
    }
    return !ovf;
}

static bool KS(mspm)(int64_t n, int64_t K, int64_t p, const int64_t* rp, const int32_t* ci,
                     const KT* av, const KT* m, KT* y) {
    KS(MspmCtx) c = { rp, ci, av, m, y, K, p };
    int64_t w = rp[K] + K;
    return par_rows(n, NULL, w > 0 ? w : 1, KS(mspm_rows), &c);
}

/* ---- M . sv : y[r] = sum over the vector's entries M[r,k] x_k --------------- */

typedef struct {
    int64_t nnz; const int32_t* ci; const KT* xv;
    const KT* m; KT* y; int64_t K;
} KS(MsvCtx);

static bool KS(msv_rows)(void* c, int64_t r0, int64_t r1) {
    const KS(MsvCtx)* t = c;
    bool ovf = false;
    for (int64_t r = r0; r < r1; r++) {
        const KT* mr = t->m + r * t->K;
        KT acc = K_ZERO;
        for (int64_t q = 0; q < t->nnz; q++) K_FMA(acc, mr[t->ci[q]], t->xv[q]);
        t->y[r] = acc;
    }
    return !ovf;
}

static bool KS(msv)(int64_t n, int64_t K, int64_t nnz, const int32_t* ci, const KT* xv,
                    const KT* m, KT* y) {
    KS(MsvCtx) c = { nnz, ci, xv, m, y, K };
    return par_rows(n, NULL, nnz > 0 ? nnz : 1, KS(msv_rows), &c);
}

/* ---- sv . x and sv . sv : scalars ------------------------------------------ */

static bool KS(svdot_dense)(int64_t nnz, const int32_t* ci, const KT* xv, const KT* d, KT* out) {
    bool ovf = false;
    KT acc = K_ZERO;
    for (int64_t q = 0; q < nnz; q++) K_FMA(acc, xv[q], d[ci[q]]);
    *out = acc;
    return !ovf;
}

static bool KS(svdot_sparse)(int64_t na, const int32_t* ca, const KT* va,
                             int64_t nb, const int32_t* cb, const KT* vb, KT* out) {
    bool ovf = false;
    KT acc = K_ZERO;
    int64_t i = 0, j = 0;
    while (i < na && j < nb) {
        if (ca[i] == cb[j]) { K_FMA(acc, va[i], vb[j]); i++; j++; }
        else if (ca[i] < cb[j]) i++;
        else j++;
    }
    *out = acc;
    return !ovf;
}

/* ---- s . sv : sparse result; row i is explicit when it meets an entry ------- */

typedef struct {
    const int64_t* rp; const int32_t* ci; const KT* av;
    const KT* xd; const unsigned char* xm;      /* dense values + presence mask */
    KT* yv; unsigned char* ym;
} KS(SpsvCtx);

static bool KS(spsv_rows)(void* c, int64_t r0, int64_t r1) {
    const KS(SpsvCtx)* t = c;
    bool ovf = false;
    for (int64_t i = r0; i < r1; i++) {
        KT acc = K_ZERO;
        unsigned char hit = 0;
        for (int64_t k = t->rp[i]; k < t->rp[i + 1]; k++) {
            int32_t j = t->ci[k];
            if (t->xm[j]) { K_FMA(acc, t->av[k], t->xd[j]); hit = 1; }
        }
        t->yv[i] = acc;
        t->ym[i] = hit;
    }
    return !ovf;
}

static bool KS(spsv)(int64_t n, const int64_t* rp, const int32_t* ci, const KT* av,
                     const KT* xd, const unsigned char* xm, KT* yv, unsigned char* ym) {
    KS(SpsvCtx) c = { rp, ci, av, xd, xm, yv, ym };
    return par_rows(n, rp, 1, KS(spsv_rows), &c);
}

/* ---- s1 . s2 : Gustavson SpGEMM ------------------------------------------- */
/* Two passes over row chunks balanced by multiply count. The first counts each
 * output row's distinct columns with a marker array; the exact output is then
 * allocated once and the second pass fills each row in place. Building the
 * output in per-chunk buffers and concatenating doubled the fresh pages
 * touched and cost more than the arithmetic.
 *
 * The fill works row by row in one of two ways:
 *   - a row with few products (the common case) is expanded into a small
 *     buffer of (column, value) pairs, sorted by column, and compressed --
 *     everything stays in L1;
 *   - a row with many products uses a dense accumulator over the output
 *     columns with the marker, the touched columns sorted afterwards.
 * A dense accumulator per thread is width-sized; with every thread streaming
 * random updates into its own, the working set overflows the shared cache and
 * each product pays a memory round trip -- which is why it is kept for the
 * heavy rows only. Either way every row comes out in canonical order.
 *
 * Markers are int32 stamps relative to the chunk, from calloc, so a banded
 * product touches only the pages of its band. Both passes cut the rows
 * identically (same keys and total), so a chunk's marker is reused by the
 * second pass under higher stamps without clearing. */

#define KGEMM_ESC_MAX 512        /* most products in a row for the pair path */

typedef struct {
    int32_t c;
    KT v;
} KS(GPair);

typedef struct {
    int32_t* mark;           /* p stamps, by column */
    int64_t base;            /* the first pass's highest stamp */
} KS(GemmScratch);

typedef struct {
    int64_t n, p;
    const int64_t* arp; const int32_t* aci; const KT* av;
    const int64_t* brp; const int32_t* bci; const KT* bv;
    const int64_t* fl;           /* multiply-count prefix, n + 1 (the keys) */
    int64_t* rowptr;             /* n + 1: counts after pass 1, offsets after */
    KS(GemmScratch)* scr;        /* by first row of the chunk */
    int32_t* oc;                 /* pass 2 output */
    KT* ov;
} KS(GemmCtx);

static int32_t* KS(gemm_mark)(KS(GemmCtx)* t, int64_t r0) {
    if (!t->scr[r0].mark) t->scr[r0].mark = calloc((size_t)(t->p ? t->p : 1), sizeof(int32_t));
    return t->scr[r0].mark;
}

static bool KS(gemm_count)(void* c, int64_t r0, int64_t r1) {
    KS(GemmCtx)* t = c;
    int32_t* mark = KS(gemm_mark)(t, r0);
    if (!mark) return false;
    t->scr[r0].base = r1 - r0;
    for (int64_t i = r0; i < r1; i++) {
        int32_t stamp = (int32_t)(i - r0 + 1);
        int64_t cnt = 0;
        for (int64_t ka = t->arp[i]; ka < t->arp[i + 1]; ka++) {
            int32_t k = t->aci[ka];
            for (int64_t kb = t->brp[k]; kb < t->brp[k + 1]; kb++) {
                int32_t j = t->bci[kb];
                if (mark[j] != stamp) { mark[j] = stamp; cnt++; }
            }
        }
        t->rowptr[i + 1] = cnt;
    }
    return true;
}

/* In-place ascending sort of pairs by column. */
static void KS(gpair_sort)(KS(GPair)* a, int64_t n) {
    while (n > 24) {
        int32_t x = a[0].c, y = a[n / 2].c, z = a[n - 1].c;
        int32_t piv = x < y ? (y < z ? y : (x < z ? z : x)) : (x < z ? x : (y < z ? z : y));
        int64_t i = -1, j = n;
        for (;;) {
            do i++; while (a[i].c < piv);
            do j--; while (a[j].c > piv);
            if (i >= j) break;
            KS(GPair) t = a[i]; a[i] = a[j]; a[j] = t;
        }
        int64_t left = j + 1;
        if (left < n - left) { KS(gpair_sort)(a, left); a += left; n -= left; }
        else { KS(gpair_sort)(a + left, n - left); n = left; }
    }
    for (int64_t i = 1; i < n; i++) {
        KS(GPair) v = a[i];
        int64_t j = i - 1;
        while (j >= 0 && a[j].c > v.c) { a[j + 1] = a[j]; j--; }
        a[j + 1] = v;
    }
}

static bool KS(gemm_fill)(void* c, int64_t r0, int64_t r1) {
    KS(GemmCtx)* t = c;
    bool ovf = false;
    int64_t p = t->p;
    int32_t* mark = NULL;
    KT* acc = NULL;
    KS(GPair)* pr = malloc(sizeof(KS(GPair)) * KGEMM_ESC_MAX);
    if (!pr) return false;
    int32_t* oc = t->oc;
    KT* ov = t->ov;
    for (int64_t i = r0; i < r1; i++) {
        int64_t start = t->rowptr[i], cnt = start;
        int64_t prods = t->fl[i + 1] - t->fl[i] - 1;
        if (prods <= KGEMM_ESC_MAX) {
            /* expand, sort, compress */
            int64_t m = 0;
            for (int64_t ka = t->arp[i]; ka < t->arp[i + 1]; ka++) {
                int32_t k = t->aci[ka];
                KT a = t->av[ka];
                for (int64_t kb = t->brp[k]; kb < t->brp[k + 1]; kb++) {
                    pr[m].c = t->bci[kb];
                    pr[m].v = K_ZERO;
                    K_FMA(pr[m].v, a, t->bv[kb]);
                    m++;
                }
            }
            if (m > 1) KS(gpair_sort)(pr, m);
            for (int64_t q = 0; q < m; q++) {
                if (cnt > start && oc[cnt - 1] == pr[q].c) {
                    K_FMA(ov[cnt - 1], K_ONE, pr[q].v);
                } else {
                    oc[cnt] = pr[q].c;
                    ov[cnt] = pr[q].v;
                    cnt++;
                }
            }
            continue;
        }
        /* a heavy row: dense accumulator */
        if (!acc) {
            mark = KS(gemm_mark)(t, r0);
            acc = malloc(sizeof(KT) * (size_t)(p ? p : 1));
            if (!mark || !acc) { free(acc); free(pr); return false; }
        }
        int32_t stamp = (int32_t)(t->scr[r0].base + (i - r0) + 1);
        for (int64_t ka = t->arp[i]; ka < t->arp[i + 1]; ka++) {
            int32_t k = t->aci[ka];
            KT a = t->av[ka];
            for (int64_t kb = t->brp[k]; kb < t->brp[k + 1]; kb++) {
                int32_t j = t->bci[kb];
                if (mark[j] != stamp) {
                    mark[j] = stamp;
                    oc[cnt++] = j;
                    acc[j] = K_ZERO;
                }
                K_FMA(acc[j], a, t->bv[kb]);
            }
        }
        int64_t len = cnt - start;
        if (len > 1) {
            if (len * 16 > p) {
                /* a dense-ish row: the marker scan is already in order */
                int64_t q = start;
                for (int64_t j = 0; j < p; j++) if (mark[j] == stamp) oc[q++] = (int32_t)j;
            } else {
                sort_i32(oc + start, len);
            }
        }
        for (int64_t q = start; q < cnt; q++) ov[q] = acc[oc[q]];
    }
    free(acc);
    free(pr);
    return !ovf;
}

typedef struct {
    const int64_t* arp; const int32_t* aci; const int64_t* brp;
    int64_t* fl;
} KS(GemmFlops);

static bool KS(gemm_flops)(void* c, int64_t r0, int64_t r1) {
    const KS(GemmFlops)* t = c;
    for (int64_t i = r0; i < r1; i++) {
        int64_t f = 1;
        for (int64_t ka = t->arp[i]; ka < t->arp[i + 1]; ka++) {
            int32_t k = t->aci[ka];
            f += t->brp[k + 1] - t->brp[k];
        }
        t->fl[i + 1] = f;
    }
    return true;
}

/* C = A . B, A n x K and B K x p (CSR, int32 columns). On success *out_pat is a
 * fresh pattern (n rows, width 1, idx32) and *out_vals its KT values. */
static bool KS(spgemm)(int64_t n, int64_t p,
                       const int64_t* arp, const int32_t* aci, const KT* av,
                       const int64_t* brp, const int32_t* bci, const KT* bv,
                       SparsePattern** out_pat, KT** out_vals) {
    *out_pat = NULL; *out_vals = NULL;
    /* the int32 stamps count rows within a chunk */
    if (n >= ((int64_t)1 << 30)) return false;
    int64_t* fl = malloc(sizeof(int64_t) * (size_t)(n + 1));
    KS(GemmScratch)* scr = calloc((size_t)(n ? n : 1), sizeof(KS(GemmScratch)));
    SparsePattern* pat = NULL;
    KT* vals = NULL;
    bool ok = fl && scr;
    if (ok) {
        /* multiply counts: row i costs sum over its entries of |row k of B|,
         * +1 so an empty row still has a key of its own */
        KS(GemmFlops) fc = { arp, aci, brp, fl };
        par_rows(n, arp, 1, KS(gemm_flops), &fc);
        fl[0] = 0;
        for (int64_t i = 0; i < n; i++) fl[i + 1] += fl[i];
        pat = sp_pattern_new(n, 0, 1, true);        /* rowptr now, colidx below */
        ok = pat != NULL;
    }
    KS(GemmCtx) c = { n, p, arp, aci, av, brp, bci, bv, fl, pat ? pat->rowptr : NULL,
                      scr, NULL, NULL };
    if (ok) ok = par_rows(n, fl, 0, KS(gemm_count), &c);
    if (ok) {
        int64_t* rp = pat->rowptr;
        rp[0] = 0;
        for (int64_t i = 0; i < n; i++) rp[i + 1] += rp[i];
        int64_t total = rp[n];
        pat->nnz = total;
        pat->colidx = malloc(sizeof(int32_t) * (size_t)(total ? total : 1));
        vals = malloc(sizeof(KT) * (size_t)(total ? total : 1));
        ok = pat->colidx && vals;
        if (ok) {
            c.oc = pat->colidx;
            c.ov = vals;
            ok = par_rows(n, fl, 0, KS(gemm_fill), &c);
        }
        if (total == 0) { free(pat->colidx); pat->colidx = NULL; }
    }
    if (scr) for (int64_t r = 0; r < n; r++) free(scr[r].mark);
    free(scr);
    free(fl);
    if (!ok) { sp_pattern_release(pat); free(vals); return false; }
    *out_pat = pat;
    *out_vals = vals;
    return true;
}

#undef KGEMM_ESC_MAX

#undef KT
#undef KS
#undef K_ZERO
#undef K_ONE
#undef K_FMA
