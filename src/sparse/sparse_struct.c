/*
 * sparse_struct.c -- native structural operations on SparseArrays.
 *
 * Every head here rearranges, selects, concatenates or assembles the stored
 * entries of its SparseArray arguments and answers with a SparseArray built
 * straight from the index structure -- no dense intermediate, so the cost is
 * O(rows + nnz) of the operands (O(rows + nnz in range) for a slice), not
 * O(prod(dims)):
 *
 *   Transpose / ConjugateTranspose   counting-sort CSR -> CSR of the transpose
 *                                    (threaded); a general level permutation
 *                                    through the COO builder
 *   Flatten / ArrayReshape           re-key the row-major flat offsets: the
 *                                    entry order is unchanged, so the values
 *                                    are copied as one block
 *   Partition                        blocks of level-1 slices
 *   Join (levels 1 and 2)            concatenate CSR pieces
 *   Take / Drop / Rest / Most /      strided range selection per level (a
 *   First / Last / Reverse           row slice touches only its rows)
 *   Append / Prepend                 Join with a one-row block
 *   Diagonal / DiagonalMatrix        one binary search per row / one entry
 *                                    per row
 *   IdentityMatrix[n, SparseArray]   (called from linalg/construct.c)
 *   KroneckerProduct                 block-structured CSR, threaded by row
 *   ArrayFlatten                     block assembly
 *   Upper/LowerTriangularize         a cut per row
 *   ArrayPad / PadRight / PadLeft    shift (and crop) the coordinates; an
 *                                    explicit pad value is emitted in order
 *   Map over a vector                f over the stored values and the
 *                                    background (sparse_thread.c's tiers)
 *
 * Each is registered with symtab_set_sparse_native: it borrows the evaluated
 * call and returns the answer, or NULL to decline -- the head's gate mode
 * (sparse_gate.c: densify, evaluate, re-sparsify) then answers exactly as
 * before, and it is also what emits the head's messages for an invalid call.
 * So a kernel declines everything it has not verified: a spec it does not
 * implement, an out-of-range part, an empty result (which is a dense {}),
 * mismatched shapes or backgrounds.
 *
 * KIND OF THE ANSWER. A structural answer is sparse over the operand's
 * background. Where the answer is DENSE by the semantics -- a triangular part
 * or a Kronecker product over a non-zero background, an ArrayFlatten whose
 * blocks are not all sparse over one background, a Map whose values are
 * lists -- the kernel answers with the dense result (computed on the Normal
 * of the operands) rather than letting the gate re-sparsify it.
 */

#include "sparse_internal.h"
#include "ndarray_internal.h"
#include "sym_names.h"
#include "symtab.h"
#include "eval.h"
#include "checked_int.h"
#include "list/flatten.h"
#include <stdlib.h>
#include <string.h>

/* Largest dense element count a kernel materialises to give a dense answer
 * (the gate's SP_GATE_MAX). */
#define SPS_DENSE_MAX ((int64_t)1 << 27)
/* Largest number of explicit pad entries emitted, and largest row count a
 * reshaped CSR result may have (its rowptr is allocated in full). */
#define SPS_EXPLICIT_MAX ((int64_t)1 << 27)
#define SPS_ROWS_MAX ((int64_t)1 << 31)
/* Work (entries) below which a kernel runs serially. */
#define SPS_PAR_MIN ((int64_t)1 << 17)
/* Row blocks of the threaded transpose. */
#define SPS_BLOCKS 16

/* ======================================================================== */
/* Small helpers                                                            */
/* ======================================================================== */

#define NARGS(c) ((c)->data.function.arg_count)
#define ARG(c, i) ((c)->data.function.args[i])

static bool is_sym(const Expr* e, const char* s) {
    return e->type == EXPR_SYMBOL && e->data.symbol.name == s;
}

static bool head_is(const Expr* e, const char* s) {
    return e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL &&
           e->data.function.head->data.symbol.name == s;
}

/* An array-valued argument (List, packed or visible NDArray, SparseArray). */
static bool is_arraylike(const Expr* e) {
    return head_is(e, SYM_List) || e->type == EXPR_NDARRAY || e->type == EXPR_SPARSE;
}

static bool get_i64(const Expr* e, int64_t* v) {
    if (e->type != EXPR_INTEGER) return false;
    *v = e->data.integer;
    return true;
}

/* A List of machine Integers (or a packed int64 vector) of at most `max`
 * elements into out[]. Returns the count, or -1. */
static int int_list(const Expr* e, int64_t* out, int max) {
    if (e->type == EXPR_NDARRAY && e->data.ndarray.rank == 1 &&
        e->data.ndarray.dtype == NDT_INT64) {
        int64_t n = e->data.ndarray.dims[0];
        if (n > max) return -1;
        for (int64_t i = 0; i < n; i++) out[i] = ((const int64_t*)e->data.ndarray.data)[i];
        return (int)n;
    }
    if (!head_is(e, SYM_List) || NARGS(e) > (size_t)max) return -1;
    for (size_t i = 0; i < NARGS(e); i++)
        if (!get_i64(ARG(e, i), &out[i])) return -1;
    return (int)NARGS(e);
}

static Expr* mk_call1(const char* head, Expr* a) {
    Expr* args[1] = { a };
    return expr_new_function(expr_new_symbol(head), args, 1);
}

/* Evaluate head[a] (a borrowed). */
static Expr* eval1(const char* head, const Expr* a) {
    Expr* c = mk_call1(head, expr_copy((Expr*)a));
    Expr* r = evaluate(c);
    expr_free(c);
    return r;
}

static inline void pat_set(SparsePattern* p, int64_t k, int j, int64_t v) {
    if (p->idx32) ((int32_t*)p->colidx)[k * p->width + j] = (int32_t)v;
    else          ((int64_t*)p->colidx)[k * p->width + j] = v;
}

/* First entry in [lo, hi) of a row whose tuple coordinate 0 is >= v. */
static int64_t lower_bound0(const SparsePattern* p, int64_t lo, int64_t hi, int64_t v) {
    while (lo < hi) {
        int64_t mid = lo + (hi - lo) / 2;
        if (sp_col(p, mid, 0) < v) lo = mid + 1; else hi = mid;
    }
    return lo;
}

static int64_t row_lo(const SparseArray* sa, int64_t r) { return sa->pat->rowptr[r]; }
static int64_t row_hi(const SparseArray* sa, int64_t r) { return sa->pat->rowptr[r + 1]; }

/* ---- value buffers ------------------------------------------------------- */

/* The kind a single value would have in a typed buffer. */
static SpValKind expr_kind(const Expr* e) {
    if (e->type == EXPR_INTEGER) return SPV_INT64;
    if (e->type == EXPR_REAL) return SPV_F64;
    if (head_is(e, SYM_Complex) && NARGS(e) == 2 && ARG(e, 0)->type == EXPR_REAL &&
        ARG(e, 1)->type == EXPR_REAL) return SPV_C128;
    return SPV_EXPR;
}

static SpValKind kind_join(SpValKind a, SpValKind b) { return a == b ? a : SPV_EXPR; }

/* An output value buffer of one kind. Copying a value of another kind boxes
 * it, so the buffer's kind must be SPV_EXPR whenever the sources differ (the
 * callers decide that up front with kind_join); sp_new re-types a boxed
 * buffer that turns out uniform. */
typedef struct {
    SpValKind kind;
    void*     v;
    int64_t   cap;
} VBuf;

static bool vb_init(VBuf* b, SpValKind kind, int64_t n) {
    b->kind = kind;
    b->cap = n;
    size_t sz = kind == SPV_EXPR ? sizeof(Expr*) : sp_kind_size(kind);
    b->v = calloc((size_t)(n > 0 ? n : 1), sz);   /* zeroed: a partial EXPR buffer frees cleanly */
    return b->v != NULL;
}

static void vb_free(VBuf* b) {
    if (b->kind == SPV_EXPR && b->v)
        for (int64_t i = 0; i < b->cap; i++) if (((Expr**)b->v)[i]) expr_free(((Expr**)b->v)[i]);
    free(b->v);
    b->v = NULL;
}

/* Slot d <- value k of sa. */
static inline void vb_from(VBuf* b, int64_t d, const SparseArray* sa, int64_t k) {
    if (b->kind == sa->kind) {
        switch (b->kind) {
            case SPV_INT64: ((int64_t*)b->v)[d] = ((const int64_t*)sa->vals)[k]; return;
            case SPV_F64:   ((double*)b->v)[d] = ((const double*)sa->vals)[k]; return;
            case SPV_C128:
                ((double*)b->v)[2 * d] = ((const double*)sa->vals)[2 * k];
                ((double*)b->v)[2 * d + 1] = ((const double*)sa->vals)[2 * k + 1];
                return;
            default: ((Expr**)b->v)[d] = expr_copy(((Expr* const*)sa->vals)[k]); return;
        }
    }
    ((Expr**)b->v)[d] = sp_value_expr(sa, k);
}

/* Slots d .. d+n-1 <- values k .. k+n-1 of sa. */
/* memcpy, split across threads for a large block (a fresh destination's page
 * faults are then taken in parallel too). */
typedef struct { char* dst; const char* src; } CpyCtx;
static bool cpy_chunk(void* c, size_t lo, size_t hi) {
    CpyCtx* t = c;
    memcpy(t->dst + lo, t->src + lo, hi - lo);
    return true;
}
#define SPS_PAR_BYTES ((size_t)1 << 21)
static void par_memcpy(void* dst, const void* src, size_t n) {
    if (n < SPS_PAR_BYTES) { memcpy(dst, src, n); return; }
    CpyCtx t = { dst, src };
    nd_parallel_for(n, cpy_chunk, &t);
}

static void vb_range(VBuf* b, int64_t d, const SparseArray* sa, int64_t k, int64_t n) {
    if (n <= 0) return;
    if (b->kind == sa->kind && b->kind != SPV_EXPR) {
        size_t sz = sp_kind_size(b->kind);
        par_memcpy((char*)b->v + (size_t)d * sz, (const char*)sa->vals + (size_t)k * sz, (size_t)n * sz);
        return;
    }
    for (int64_t i = 0; i < n; i++) vb_from(b, d + i, sa, k + i);
}

/* Slot d <- a copy of e (whose expr_kind matches a typed buffer). */
static void vb_put(VBuf* b, int64_t d, const Expr* e) {
    switch (b->kind) {
        case SPV_INT64: ((int64_t*)b->v)[d] = e->data.integer; return;
        case SPV_F64:   ((double*)b->v)[d] = e->data.real; return;
        case SPV_C128:
            ((double*)b->v)[2 * d] = ARG(e, 0)->data.real;
            ((double*)b->v)[2 * d + 1] = ARG(e, 1)->data.real;
            return;
        default: ((Expr**)b->v)[d] = expr_copy((Expr*)e); return;
    }
}

/* Assemble the answer. Takes `bg`, the pattern reference and the buffer. */
static Expr* finish(int rank, const int64_t* dims, Expr* bg, SparsePattern* pat, VBuf* b) {
    void* v = b->v;
    b->v = NULL;
    return sp_to_expr(sp_new(rank, dims, bg, pat, b->kind, v));
}

/* ---- threading ----------------------------------------------------------- */

/* Run fn(ctx, b) for b in [0, nb), each block exactly once, in parallel when
 * threads are available. nd_parallel_for splits a flat range it chooses
 * itself; giving it nb virtual stretches of NDARRAY_THREAD_THRESHOLD slots and
 * running block b in whichever chunk holds slot b * THRESHOLD makes every
 * block land in exactly one chunk, with the chunks' block sets fixed by the
 * caller rather than by the split. */
typedef void (*BlockFn)(void* ctx, int b);
typedef struct { BlockFn fn; void* ctx; int nb; } BlockJob;

static bool block_chunk(void* c, size_t lo, size_t hi) {
    BlockJob* j = c;
    for (int b = 0; b < j->nb; b++) {
        size_t v = (size_t)b * NDARRAY_THREAD_THRESHOLD;
        if (v >= lo && v < hi) j->fn(j->ctx, b);
    }
    return true;
}

static void run_blocks(int nb, BlockFn fn, void* ctx) {
    if (nb <= 1) { if (nb == 1) fn(ctx, 0); return; }
    BlockJob j = { fn, ctx, nb };
    nd_parallel_for((size_t)nb * NDARRAY_THREAD_THRESHOLD, block_chunk, &j);
}

/* fn over [0, n) in SPS_BLOCKS contiguous pieces run in parallel when `work`
 * (the entries the pieces touch) is large -- for loops whose length is a row
 * count too small for nd_parallel_for's own threshold. fn must not fail. */
typedef struct { nd_chunk_fn fn; void* ctx; size_t n; int nb; } RangeJob;
static void range_block(void* c, int b) {
    RangeJob* j = c;
    size_t lo = j->n * (size_t)b / (size_t)j->nb, hi = j->n * (size_t)(b + 1) / (size_t)j->nb;
    if (lo < hi) j->fn(j->ctx, lo, hi);
}
static void par_range(size_t n, int64_t work, nd_chunk_fn fn, void* ctx) {
    if (work < SPS_PAR_MIN || n < 2) { fn(ctx, 0, n); return; }
    RangeJob j = { fn, ctx, n, n < SPS_BLOCKS ? (int)n : SPS_BLOCKS };
    run_blocks(j.nb, range_block, &j);
}

/* ---- dense answers ------------------------------------------------------- */

/* `e` with every SparseArray (looking `depth` List levels down) replaced by
 * its Normal; NULL with *fail set when one is too large. */
static Expr* densify(const Expr* e, int depth, bool* fail) {
    if (e->type == EXPR_SPARSE) {
        const SparseArray* sa = sp_of(e);
        int64_t t = sp_total(sa->rank, sa->dims);
        Expr* d = (t >= 0 && t <= SPS_DENSE_MAX) ? sp_normal(sa) : NULL;
        if (!d) *fail = true;
        return d;
    }
    if (depth > 0 && head_is(e, SYM_List)) {
        size_t n = NARGS(e);
        Expr** it = malloc(sizeof(Expr*) * (n ? n : 1));
        if (!it) { *fail = true; return NULL; }
        for (size_t i = 0; i < n; i++) {
            it[i] = densify(ARG(e, i), depth - 1, fail);
            if (*fail) {
                for (size_t q = 0; q < i; q++) expr_free(it[q]);
                free(it);
                return NULL;
            }
        }
        Expr* out = expr_new_function(expr_copy(e->data.function.head), it, n);
        free(it);
        return out;
    }
    return expr_copy((Expr*)e);
}

/* The call evaluated on the Normal of its SparseArray arguments, with the
 * dense answer kept dense -- for the cases whose answer is dense by the
 * semantics. NULL when an operand is too large or the head declined (the
 * gate then answers, and reports). */
static Expr* dense_answer(const Expr* call, int depth) {
    bool fail = false;
    size_t n = NARGS(call);
    Expr** a = malloc(sizeof(Expr*) * (n ? n : 1));
    if (!a) return NULL;
    for (size_t i = 0; i < n; i++) {
        a[i] = densify(ARG(call, i), depth, &fail);
        if (fail) {
            for (size_t q = 0; q < i; q++) expr_free(a[q]);
            free(a);
            return NULL;
        }
    }
    Expr* d = expr_new_function(expr_copy(call->data.function.head), a, n);
    free(a);
    Expr* r = evaluate(d);
    expr_free(d);
    if (r->type == EXPR_FUNCTION && expr_eq(r->data.function.head, call->data.function.head)) {
        expr_free(r);
        return NULL;
    }
    return r;
}

static bool has_machine(const SparseArray* sa);
static bool all_numbers(const SparseArray* sa);

/* ======================================================================== */
/* Transpose / ConjugateTranspose                                           */
/* ======================================================================== */

typedef struct {
    const SparseArray* sa;
    SparsePattern*     out;
    void*              ovals;
    int64_t*           cnt;     /* nb x ncols: per block column counts, then cursors */
    int64_t            m, n;    /* source rows, columns */
    int                nb;
} TrCtx;

static void tr_count(void* c, int b) {
    TrCtx* t = c;
    const SparsePattern* p = t->sa->pat;
    int64_t r0 = t->m * b / t->nb, r1 = t->m * (b + 1) / t->nb;
    int64_t* cnt = t->cnt + (int64_t)b * t->n;
    int64_t k0 = p->rowptr[r0], k1 = p->rowptr[r1];
    if (p->idx32) { const int32_t* ci = p->colidx; for (int64_t k = k0; k < k1; k++) cnt[ci[k]]++; }
    else          { const int64_t* ci = p->colidx; for (int64_t k = k0; k < k1; k++) cnt[ci[k]]++; }
}

/* Scatter block b's entries to their transposed slots. The slot of an entry
 * is its column's cursor, so within a new row the entries arrive in source
 * row order: the result is canonical without sorting. */
#define TR_SCATTER(IN_T, OUT_T, VAL_STMT)                                       \
    do {                                                                        \
        const IN_T* ci = p->colidx;                                             \
        OUT_T* co = t->out->colidx;                                             \
        for (int64_t r = r0; r < r1; r++)                                       \
            for (int64_t k = p->rowptr[r]; k < p->rowptr[r + 1]; k++) {         \
                int64_t pos = cnt[ci[k]]++;                                     \
                co[pos] = (OUT_T)r;                                             \
                VAL_STMT;                                                       \
            }                                                                   \
    } while (0)

#define TR_SCATTER_KIND(IN_T, OUT_T)                                                        \
    do {                                                                                    \
        switch (sa->kind) {                                                                 \
            case SPV_INT64: { const int64_t* sv = sa->vals; int64_t* dv = t->ovals;         \
                TR_SCATTER(IN_T, OUT_T, dv[pos] = sv[k]); break; }                          \
            case SPV_F64: { const double* sv = sa->vals; double* dv = t->ovals;             \
                TR_SCATTER(IN_T, OUT_T, dv[pos] = sv[k]); break; }                          \
            case SPV_C128: { const double* sv = sa->vals; double* dv = t->ovals;            \
                TR_SCATTER(IN_T, OUT_T, (dv[2 * pos] = sv[2 * k], dv[2 * pos + 1] = sv[2 * k + 1])); \
                break; }                                                                    \
            default: { Expr* const* sv = sa->vals; Expr** dv = t->ovals;                    \
                TR_SCATTER(IN_T, OUT_T, dv[pos] = expr_copy(sv[k])); break; }               \
        }                                                                                   \
    } while (0)

static void tr_scatter(void* c, int b) {
    TrCtx* t = c;
    const SparseArray* sa = t->sa;
    const SparsePattern* p = sa->pat;
    int64_t r0 = t->m * b / t->nb, r1 = t->m * (b + 1) / t->nb;
    int64_t* cnt = t->cnt + (int64_t)b * t->n;
    if (p->idx32 && t->out->idx32)       TR_SCATTER_KIND(int32_t, int32_t);
    else if (p->idx32)                   TR_SCATTER_KIND(int32_t, int64_t);
    else if (t->out->idx32)              TR_SCATTER_KIND(int64_t, int32_t);
    else                                 TR_SCATTER_KIND(int64_t, int64_t);
}

/* The transpose of a matrix: a counting sort of the entries by column. */
static SparseArray* transpose2(const SparseArray* sa) {
    int64_t m = sa->dims[0], n = sa->dims[1], nnz = sa->pat->nnz;
    int64_t odims[2] = { n, m };
    if (n > SPS_ROWS_MAX) return NULL;
    /* Boxed values are reference-counted without atomics: copy them on one
     * thread. */
    int nb = (nnz >= SPS_PAR_MIN && sa->kind != SPV_EXPR) ? SPS_BLOCKS : 1;
    if (nb > m) nb = m > 0 ? (int)m : 1;
    /* per-block column counts cost nb * n; fewer blocks when n dwarfs nnz */
    while (nb > 1 && (int64_t)nb * n > 4 * nnz) nb /= 2;
    SparsePattern* out = sp_pattern_new(n, nnz, 1, sp_dims_idx32(2, odims));
    int64_t* cnt = calloc((size_t)nb * (size_t)(n > 0 ? n : 1), sizeof(int64_t));
    VBuf vb;
    if (!out || !cnt || !vb_init(&vb, sa->kind, nnz)) {
        sp_pattern_release(out); free(cnt);
        if (out && cnt) vb_free(&vb);
        return NULL;
    }
    TrCtx t = { sa, out, vb.v, cnt, m, n, nb };
    run_blocks(nb, tr_count, &t);
    /* rowptr from the column totals, then each block's cursor per column */
    out->rowptr[0] = 0;
    for (int64_t c = 0; c < n; c++) {
        int64_t s = 0;
        for (int b = 0; b < nb; b++) s += cnt[(int64_t)b * n + c];
        out->rowptr[c + 1] = out->rowptr[c] + s;
    }
    if (nb == 1) {
        memcpy(cnt, out->rowptr, sizeof(int64_t) * (size_t)n);
    } else {
        int64_t* cur = malloc(sizeof(int64_t) * (size_t)(n > 0 ? n : 1));
        if (!cur) { sp_pattern_release(out); free(cnt); vb_free(&vb); return NULL; }
        memcpy(cur, out->rowptr, sizeof(int64_t) * (size_t)n);
        for (int b = 0; b < nb; b++) {
            int64_t* cb = cnt + (int64_t)b * n;
            for (int64_t c = 0; c < n; c++) { int64_t x = cb[c]; cb[c] = cur[c]; cur[c] += x; }
        }
        free(cur);
    }
    run_blocks(nb, tr_scatter, &t);
    free(cnt);
    void* v = vb.v;
    return sp_new(2, odims, expr_copy(sa->bg), out, sa->kind, v);
}

/* A general level permutation: old level l becomes new level perm[l]. */
static SparseArray* transpose_perm(const SparseArray* sa, const int* perm) {
    int rank = sa->rank;
    int64_t nd[SP_MAX_RANK];
    for (int l = 0; l < rank; l++) nd[perm[l]] = sa->dims[l];
    SpCoo c;
    sp_coo_init(&c, rank, nd, sa->kind);
    const SparsePattern* p = sa->pat;
    if (!sp_coo_reserve(&c, p->nnz)) { sp_coo_free(&c); return NULL; }
    int64_t old[SP_MAX_RANK], nw[SP_MAX_RANK];
    bool ok = true;
    for (int64_t r = 0; r < p->nrows && ok; r++) {
        old[0] = r;
        for (int64_t k = p->rowptr[r]; k < p->rowptr[r + 1] && ok; k++) {
            for (int j = 0; j < p->width; j++) old[1 + j] = sp_col(p, k, j);
            for (int l = 0; l < rank; l++) nw[perm[l]] = old[l];
            switch (sa->kind) {
                case SPV_INT64: ok = sp_coo_push_i64(&c, nw, ((const int64_t*)sa->vals)[k]); break;
                case SPV_F64:   ok = sp_coo_push_f64(&c, nw, ((const double*)sa->vals)[k]); break;
                case SPV_C128:  ok = sp_coo_push_c128(&c, nw, ((const double*)sa->vals)[2 * k],
                                                      ((const double*)sa->vals)[2 * k + 1]); break;
                default:        ok = sp_coo_push_expr(&c, nw, sp_value_expr(sa, k)); break;
            }
        }
    }
    if (!ok) { sp_coo_free(&c); return NULL; }
    return sp_coo_finish(&c, expr_copy(sa->bg), false);
}

/* Conjugate the values and the background of a freshly built array (the
 * caller's only reference), giving the replacement. */
static SparseArray* conjugate_payload(SparseArray* sa) {
    int64_t nnz = sa->pat->nnz;
    if (sa->kind == SPV_C128) {
        double* v = sa->vals;
        for (int64_t k = 0; k < nnz; k++) v[2 * k + 1] = -v[2 * k + 1];
    } else if (sa->kind == SPV_EXPR) {
        Expr** v = sa->vals;
        for (int64_t k = 0; k < nnz; k++) {
            Expr* c = eval1(SYM_Conjugate, v[k]);
            expr_free(v[k]);
            v[k] = c;
        }
    }
    Expr* bg = eval1(SYM_Conjugate, sa->bg);
    /* Move the values into a fresh payload so that sp_new re-types them. */
    SparseArray* r = sp_new(sa->rank, sa->dims, bg, sp_pattern_retain(sa->pat), sa->kind, sa->vals);
    sa->vals = NULL;
    sp_release(sa);
    return r;
}

/* The permutation of Transpose[s, spec] (old level -> new level, 0-based), or
 * false for a spec this path does not take (repeated levels collapse to a
 * diagonal; the builtin handles those and the messages). */
static bool transpose_spec(const Expr* call, int rank, int* perm) {
    for (int l = 0; l < rank; l++) perm[l] = l;
    if (NARGS(call) == 1) {
        if (rank >= 2) { perm[0] = 1; perm[1] = 0; }
        return true;
    }
    if (NARGS(call) != 2) return false;
    int64_t sp[SP_MAX_RANK];
    int k = int_list(ARG(call, 1), sp, SP_MAX_RANK);
    if (k < 1 || k > rank) return false;
    bool seen[SP_MAX_RANK] = { false };
    for (int l = 0; l < k; l++) {
        if (sp[l] < 1 || sp[l] > k || seen[sp[l] - 1]) return false;
        seen[sp[l] - 1] = true;
        perm[l] = (int)sp[l] - 1;
    }
    return true;
}

static Expr* transpose_common(const Expr* call, bool conj) {
    const SparseArray* sa = sp_of(ARG(call, 0));
    if (!sa || NARGS(call) < 1) return NULL;
    for (size_t i = 1; i < NARGS(call); i++) if (ARG(call, i)->type == EXPR_SPARSE) return NULL;
    int perm[SP_MAX_RANK];
    if (!transpose_spec(call, sa->rank, perm)) return NULL;
    bool ident = true;
    for (int l = 0; l < sa->rank; l++) if (perm[l] != l) ident = false;
    SparseArray* r;
    if (ident) r = sp_retain((SparseArray*)sa);
    else if (sa->rank == 2) r = transpose2(sa);
    else r = transpose_perm(sa, perm);
    if (!r) return NULL;
    if (conj) {
        if (r == sa) {
            /* a copy we own, sharing the pattern */
            VBuf vb;
            if (!vb_init(&vb, sa->kind, sa->pat->nnz)) { sp_release(r); return NULL; }
            vb_range(&vb, 0, sa, 0, sa->pat->nnz);
            sp_release(r);
            r = sp_new(sa->rank, sa->dims, expr_copy(sa->bg), sp_pattern_retain(sa->pat), vb.kind, vb.v);
            if (!r) return NULL;
        }
        r = conjugate_payload(r);
    }
    return sp_to_expr(r);
}

static Expr* nat_transpose(const Expr* call) { return transpose_common(call, false); }
static Expr* nat_conjugate_transpose(const Expr* call) { return transpose_common(call, true); }

/* ======================================================================== */
/* Flatten / ArrayReshape: re-keying the row-major offsets                  */
/* ======================================================================== */

typedef struct {
    const SparseArray* sa;
    SparsePattern*     out;
    int64_t*           nrow;          /* new row of each entry (new rank >= 2) */
    int64_t            keep, oinner, ninner;
    const int64_t*     nd;
    int                nrank, nw, ofirst;
} RkCtx;

/* Re-key the surviving entries of source rows [lo, hi). */
static bool rk_rows(void* c, size_t lo, size_t hi) {
    RkCtx* t = c;
    const SparseArray* sa = t->sa;
    const SparsePattern* p = sa->pat;
    int64_t sub[SP_MAX_RANK];
    for (size_t r = lo; r < hi; r++) {
        int64_t base = sa->rank >= 2 ? (int64_t)r * t->oinner : 0;
        int64_t k1 = p->rowptr[r + 1] < t->keep ? p->rowptr[r + 1] : t->keep;
        for (int64_t k = p->rowptr[r]; k < k1; k++) {
            int64_t off;
            if (p->width == 1) off = base + sp_col(p, k, 0);
            else {
                off = 0;
                for (int j = 0; j < p->width; j++) off = off * sa->dims[t->ofirst + j] + sp_col(p, k, j);
                off += base;
            }
            if (t->nrank == 1) { pat_set(t->out, k, 0, off); continue; }
            int64_t nr = off / t->ninner, rem = off - nr * t->ninner;
            if (t->nw == 1) pat_set(t->out, k, 0, rem);
            else {
                for (int j = t->nw - 1; j >= 0; j--) { sub[j] = rem % t->nd[1 + j]; rem /= t->nd[1 + j]; }
                for (int j = 0; j < t->nw; j++) pat_set(t->out, k, j, sub[j]);
            }
            t->nrow[k] = nr;
        }
    }
    return true;
}

/* sa re-keyed to (nrank, nd), keeping the row-major order: the entries whose
 * flat offset is below the new total survive, and when the new total is
 * larger the positions past the old total hold `pad` (explicit unless it is
 * the background). */
static Expr* rekey(const SparseArray* sa, int nrank, const int64_t* nd, const Expr* pad) {
    int64_t ot = sp_total(sa->rank, sa->dims), nt = sp_total(nrank, nd);
    if (ot < 0 || nt <= 0) return NULL;
    int64_t nrows = nrank >= 2 ? nd[0] : 1;
    if (nrows > SPS_ROWS_MAX) return NULL;
    const SparsePattern* p = sa->pat;
    /* offsets of the old inner tuple, and the new inner stride */
    int64_t oinner = sa->rank >= 2 ? ot / sa->dims[0] : 0;
    int ofirst = sa->rank >= 2 ? 1 : 0;
    int nw = nrank >= 2 ? nrank - 1 : 1;
    int64_t ninner = nrank >= 2 ? nt / nd[0] : 1;

    /* how many entries survive (offsets increase with k) */
    int64_t keep = p->nnz;
    if (nt < ot) {
        keep = 0;
        for (int64_t r = 0; r < p->nrows; r++) {
            int64_t base = sa->rank >= 2 ? r * oinner : 0;
            if (base >= nt) break;
            int64_t k = p->rowptr[r];
            for (; k < p->rowptr[r + 1]; k++) {
                int64_t off = 0;
                for (int j = 0; j < p->width; j++) off = off * sa->dims[ofirst + j] + sp_col(p, k, j);
                if (base + off >= nt) break;
            }
            keep = k;
            if (k < p->rowptr[r + 1]) break;
        }
    }
    int64_t npad = 0;
    if (nt > ot && pad && !sp_value_is_background(pad, sa->bg)) npad = nt - ot;
    if (npad > SPS_EXPLICIT_MAX) return NULL;
    int64_t nnz = keep + npad;
    SpValKind kind = npad ? kind_join(sa->kind, expr_kind(pad)) : sa->kind;
    SparsePattern* out = sp_pattern_new(nrows, nnz, nw, sp_dims_idx32(nrank, nd));
    VBuf vb;
    if (!out || !vb_init(&vb, kind, nnz)) { sp_pattern_release(out); return NULL; }
    memset(out->rowptr, 0, sizeof(int64_t) * (size_t)(nrows + 1));

    /* the surviving entries keep their index k: re-key them by source row
     * (threaded), recording each one's new row, then count the rows */
    int64_t* nrow = NULL;
    if (nrank >= 2 && keep > 0) {
        nrow = malloc(sizeof(int64_t) * (size_t)keep);
        if (!nrow) { sp_pattern_release(out); vb_free(&vb); return NULL; }
    }
    RkCtx t = { sa, out, nrow, keep, oinner, ninner, nd, nrank, nw, ofirst };
    par_range((size_t)p->nrows, keep, rk_rows, &t);
    if (nrow) for (int64_t k = 0; k < keep; k++) out->rowptr[nrow[k] + 1]++;
    free(nrow);
    int64_t q = keep;
    int64_t sub[SP_MAX_RANK];
    vb_range(&vb, 0, sa, 0, keep);
    for (int64_t off = ot; q < nnz; off++, q++) {
        if (nrank == 1) pat_set(out, q, 0, off);
        else {
            int64_t nr = off / ninner, rem = off - nr * ninner;
            for (int j = nw - 1; j >= 0; j--) { sub[j] = rem % nd[1 + j]; rem /= nd[1 + j]; }
            for (int j = 0; j < nw; j++) pat_set(out, q, j, sub[j]);
            out->rowptr[nr + 1]++;
        }
        vb_put(&vb, q, pad);
    }
    if (nrank == 1) { out->rowptr[0] = 0; out->rowptr[1] = nnz; }
    else for (int64_t r = 0; r < nrows; r++) out->rowptr[r + 1] += out->rowptr[r];
    return finish(nrank, nd, expr_copy(sa->bg), out, &vb);
}

/* Flatten[s, {{l11, l12, ...}, {l21, ...}, ...}]: output level g combines the
 * listed input levels (row-major in the listed order); levels not listed
 * follow, one each, in order. A grouping that keeps the levels in order is a
 * re-keying; any other one goes through the COO builder. */
static Expr* flatten_groups(const SparseArray* sa, const Expr* spec) {
    int rank = sa->rank;
    int ng = 0, nl = 0;
    int order[SP_MAX_RANK];                    /* levels in output order */
    int gsize[SP_MAX_RANK];
    bool seen[SP_MAX_RANK] = { false };
    if (NARGS(spec) < 1 || NARGS(spec) > (size_t)rank) return NULL;
    for (size_t g = 0; g < NARGS(spec); g++) {
        int64_t lv[SP_MAX_RANK];
        int k = int_list(ARG(spec, g), lv, SP_MAX_RANK);
        if (k < 1) return NULL;
        for (int i = 0; i < k; i++) {
            if (lv[i] < 1 || lv[i] > rank || seen[lv[i] - 1]) return NULL;
            seen[lv[i] - 1] = true;
            order[nl++] = (int)lv[i] - 1;
        }
        gsize[ng++] = k;
    }
    for (int l = 0; l < rank; l++)
        if (!seen[l]) { gsize[ng++] = 1; order[nl++] = l; }
    int64_t nd[SP_MAX_RANK];
    for (int g = 0, i = 0; g < ng; g++) {
        nd[g] = 1;
        for (int q = 0; q < gsize[g]; q++, i++)
            if (ci_mul_i64(nd[g], sa->dims[order[i]], &nd[g])) return NULL;
    }
    bool in_order = true;
    for (int i = 0; i < rank; i++) if (order[i] != i) in_order = false;
    if (in_order) return rekey(sa, ng, nd, NULL);
    SpCoo c;
    sp_coo_init(&c, ng, nd, sa->kind);
    const SparsePattern* p = sa->pat;
    if (!sp_coo_reserve(&c, p->nnz)) { sp_coo_free(&c); return NULL; }
    int64_t old[SP_MAX_RANK], nw[SP_MAX_RANK];
    int first = rank >= 2 ? 1 : 0;
    bool ok = true;
    for (int64_t r = 0; r < p->nrows && ok; r++) {
        old[0] = r;
        for (int64_t k = p->rowptr[r]; k < p->rowptr[r + 1] && ok; k++) {
            for (int j = 0; j < p->width; j++) old[first + j] = sp_col(p, k, j);
            for (int g = 0, i = 0; g < ng; g++) {
                nw[g] = 0;
                for (int q = 0; q < gsize[g]; q++, i++) nw[g] = nw[g] * sa->dims[order[i]] + old[order[i]];
            }
            switch (sa->kind) {
                case SPV_INT64: ok = sp_coo_push_i64(&c, nw, ((const int64_t*)sa->vals)[k]); break;
                case SPV_F64:   ok = sp_coo_push_f64(&c, nw, ((const double*)sa->vals)[k]); break;
                case SPV_C128:  ok = sp_coo_push_c128(&c, nw, ((const double*)sa->vals)[2 * k],
                                                      ((const double*)sa->vals)[2 * k + 1]); break;
                default:        ok = sp_coo_push_expr(&c, nw, sp_value_expr(sa, k)); break;
            }
        }
    }
    if (!ok) { sp_coo_free(&c); return NULL; }
    return sp_to_expr(sp_coo_finish(&c, expr_copy(sa->bg), false));
}

/* Flatten[s], Flatten[s, n] and Flatten[s, {{levels}, ...}]. */
static Expr* nat_flatten(const Expr* call) {
    const SparseArray* sa = sp_of(ARG(call, 0));
    if (!sa || NARGS(call) > 2) return NULL;
    int64_t n = SP_MAX_RANK;
    if (NARGS(call) == 2) {
        const Expr* a = ARG(call, 1);
        if (head_is(a, SYM_List)) {
            for (size_t g = 0; g < NARGS(a); g++) if (!head_is(ARG(a, g), SYM_List)) return NULL;
            return flatten_groups(sa, a);
        }
        if (is_sym(a, SYM_Infinity) || (head_is(a, SYM_DirectedInfinity) && NARGS(a) == 1 &&
                                        ARG(a, 0)->type == EXPR_INTEGER && ARG(a, 0)->data.integer == 1))
            n = SP_MAX_RANK;
        else if (a->type == EXPR_INTEGER && a->data.integer < 0) {
            /* Flatten::flev, reported against the call as written; the call
             * rests (the gate must not re-ask on the dense form). */
            flatten_flev_message(call);
            return expr_copy((Expr*)call);
        }
        else if (!get_i64(a, &n) || n < 0) return NULL;
    }
    if (n == 0 || sa->rank == 1) return expr_copy(ARG(call, 0));
    if (n > sa->rank - 1) n = sa->rank - 1;
    int64_t nd[SP_MAX_RANK];
    int nrank = sa->rank - (int)n;
    int64_t lead = sp_total((int)n + 1, sa->dims);
    if (lead <= 0) return NULL;
    nd[0] = lead;
    for (int l = 1; l < nrank; l++) nd[l] = sa->dims[n + l];
    return rekey(sa, nrank, nd, NULL);
}

/* ArrayReshape[s, dims] and ArrayReshape[s, dims, pad]. */
static Expr* nat_arrayreshape(const Expr* call) {
    const SparseArray* sa = sp_of(ARG(call, 0));
    if (!sa || NARGS(call) < 2 || NARGS(call) > 3) return NULL;
    int64_t nd[SP_MAX_RANK];
    int nrank = int_list(ARG(call, 1), nd, SP_MAX_RANK);
    if (nrank < 1) return NULL;
    for (int l = 0; l < nrank; l++) if (nd[l] < 1) return NULL;
    Expr* zero = NULL;
    const Expr* pad;
    if (NARGS(call) == 3) {
        pad = ARG(call, 2);
        if (is_arraylike(pad)) return NULL;
    } else {
        pad = zero = expr_new_integer(0);
    }
    Expr* r = rekey(sa, nrank, nd, pad);
    if (zero) expr_free(zero);
    return r;
}

/* ======================================================================== */
/* Join (levels 1 and 2), Append, Prepend                                   */
/* ======================================================================== */

/* Join of equal-shaped arrays at level 1 over the background `bg` (taken). */
static Expr* join1(const SparseArray* const* ops, int m, Expr* bg) {
    const SparseArray* a0 = ops[0];
    int rank = a0->rank;
    int64_t nd[SP_MAX_RANK];
    memcpy(nd, a0->dims, sizeof(int64_t) * (size_t)rank);
    int64_t nnz = 0;
    SpValKind kind = a0->kind;
    nd[0] = 0;
    for (int j = 0; j < m; j++) {
        if (ci_add_i64(nd[0], ops[j]->dims[0], &nd[0])) { expr_free(bg); return NULL; }
        nnz += ops[j]->pat->nnz;
        if (ops[j]->pat->nnz) kind = kind_join(kind, ops[j]->kind);
    }
    if (a0->pat->nnz == 0) {
        kind = SPV_INT64;
        bool any = false;
        for (int j = 0; j < m; j++)
            if (ops[j]->pat->nnz) { kind = any ? kind_join(kind, ops[j]->kind) : ops[j]->kind; any = true; }
    }
    if (rank >= 2 && nd[0] > SPS_ROWS_MAX) { expr_free(bg); return NULL; }
    int width = rank >= 2 ? rank - 1 : 1;
    SparsePattern* out = sp_pattern_new(rank >= 2 ? nd[0] : 1, nnz, width, sp_dims_idx32(rank, nd));
    VBuf vb;
    if (!out || !vb_init(&vb, kind, nnz)) { sp_pattern_release(out); expr_free(bg); return NULL; }
    int64_t q = 0, row = 0, shift = 0;
    out->rowptr[0] = 0;
    for (int j = 0; j < m; j++) {
        const SparseArray* a = ops[j];
        const SparsePattern* p = a->pat;
        if (rank >= 2) {
            for (int64_t r = 0; r < p->nrows; r++) out->rowptr[row + r + 1] = q + p->rowptr[r + 1];
            row += p->nrows;
            if (p->idx32 == out->idx32 && p->nnz) {
                size_t isz = p->idx32 ? sizeof(int32_t) : sizeof(int64_t);
                par_memcpy((char*)out->colidx + (size_t)q * (size_t)width * isz, p->colidx,
                       (size_t)p->nnz * (size_t)width * isz);
            } else {
                for (int64_t k = 0; k < p->nnz; k++)
                    for (int w = 0; w < width; w++) pat_set(out, q + k, w, sp_col(p, k, w));
            }
        } else {
            for (int64_t k = 0; k < p->nnz; k++) pat_set(out, q + k, 0, sp_col(p, k, 0) + shift);
            shift += a->dims[0];
        }
        vb_range(&vb, q, a, 0, p->nnz);
        q += p->nnz;
    }
    if (rank == 1) out->rowptr[1] = nnz;
    return finish(rank, nd, bg, out, &vb);
}

/* Join at level 2 (rank >= 2): row r is the operands' rows r side by side. */
static Expr* join2(const SparseArray* const* ops, int m, Expr* bg) {
    const SparseArray* a0 = ops[0];
    int rank = a0->rank;
    int64_t nd[SP_MAX_RANK];
    memcpy(nd, a0->dims, sizeof(int64_t) * (size_t)rank);
    int64_t nnz = 0, off[64];
    SpValKind kind = SPV_INT64;
    bool any = false;
    nd[1] = 0;
    for (int j = 0; j < m; j++) {
        off[j] = nd[1];
        if (ci_add_i64(nd[1], ops[j]->dims[1], &nd[1])) { expr_free(bg); return NULL; }
        nnz += ops[j]->pat->nnz;
        if (ops[j]->pat->nnz) { kind = any ? kind_join(kind, ops[j]->kind) : ops[j]->kind; any = true; }
    }
    int width = rank - 1;
    SparsePattern* out = sp_pattern_new(nd[0], nnz, width, sp_dims_idx32(rank, nd));
    VBuf vb;
    if (!out || !vb_init(&vb, kind, nnz)) { sp_pattern_release(out); expr_free(bg); return NULL; }
    int64_t q = 0;
    out->rowptr[0] = 0;
    for (int64_t r = 0; r < nd[0]; r++) {
        for (int j = 0; j < m; j++) {
            const SparseArray* a = ops[j];
            int64_t lo = row_lo(a, r), hi = row_hi(a, r);
            for (int64_t k = lo; k < hi; k++, q++) {
                pat_set(out, q, 0, sp_col(a->pat, k, 0) + off[j]);
                for (int w = 1; w < width; w++) pat_set(out, q, w, sp_col(a->pat, k, w));
            }
            vb_range(&vb, q - (hi - lo), a, lo, hi - lo);
        }
        out->rowptr[r + 1] = q;
    }
    return finish(rank, nd, bg, out, &vb);
}

/* Do the backgrounds agree (SameQ, or numerically equal)? The answer's
 * background is the first operand's. */
static bool bg_compatible(const SparseArray* const* ops, int m) {
    for (int j = 1; j < m; j++)
        if (!expr_eq(ops[j]->bg, ops[0]->bg) && !sp_value_is_background(ops[j]->bg, ops[0]->bg))
            return false;
    return true;
}

static Expr* nat_join(const Expr* call) {
    size_t n = NARGS(call);
    int64_t level = 1;
    size_t na = n;
    if (n >= 2 && ARG(call, n - 1)->type == EXPR_INTEGER) { level = ARG(call, n - 1)->data.integer; na = n - 1; }
    if (na < 1 || na > 64 || (level != 1 && level != 2)) return NULL;
    const SparseArray* ops[64];
    for (size_t i = 0; i < na; i++) {
        ops[i] = sp_of(ARG(call, i));
        if (!ops[i]) return NULL;                 /* a dense array makes the answer dense */
    }
    int m = (int)na;
    int rank = ops[0]->rank;
    if (level == 2 && rank < 2) return NULL;
    for (int j = 0; j < m; j++) {
        /* arrays of other shapes join to a ragged, dense List */
        if (ops[j]->rank != rank) return dense_answer(call, 0);
        for (int l = 0; l < rank; l++) {
            if (l == level - 1) continue;
            if (ops[j]->dims[l] != ops[0]->dims[l]) return dense_answer(call, 0);
        }
    }
    if (!bg_compatible(ops, m)) return NULL;
    if (m == 1) return expr_copy(ARG(call, 0));
    Expr* bg = expr_copy(ops[0]->bg);
    return level == 1 ? join1(ops, m, bg) : join2(ops, m, bg);
}

/* `row` (rank r - 1) as a one-row block of rank r. */
static SparseArray* as_block(const SparseArray* row) {
    int rank = row->rank + 1;
    int64_t nd[SP_MAX_RANK];
    nd[0] = 1;
    memcpy(nd + 1, row->dims, sizeof(int64_t) * (size_t)row->rank);
    const SparsePattern* p = row->pat;
    int width = rank - 1;
    SparsePattern* out = sp_pattern_new(1, p->nnz, width, sp_dims_idx32(rank, nd));
    VBuf vb;
    if (!out || !vb_init(&vb, row->kind, p->nnz)) { sp_pattern_release(out); return NULL; }
    out->rowptr[0] = 0;
    out->rowptr[1] = p->nnz;
    int rf = row->rank >= 2 ? 1 : 0;
    for (int64_t r = 0; r < p->nrows; r++)
        for (int64_t k = p->rowptr[r]; k < p->rowptr[r + 1]; k++) {
            if (rf) pat_set(out, k, 0, r);
            for (int j = 0; j < p->width; j++) pat_set(out, k, rf + j, sp_col(p, k, j));
        }
    vb_range(&vb, 0, row, 0, p->nnz);
    void* v = vb.v;
    return sp_new(rank, nd, expr_copy(row->bg), out, vb.kind, v);
}

/* The element to append, as a SparseArray of rank sa->rank - 1 over sa's
 * background (a fresh reference), or NULL when it is not one this path takes
 * (a row of the wrong shape makes a ragged, dense answer). */
static SparseArray* row_operand(const SparseArray* sa, const Expr* e) {
    int rr = sa->rank - 1;
    const SparseArray* re = sp_of(e);
    if (re) {
        if (re->rank != rr) return NULL;
        for (int l = 0; l < rr; l++) if (re->dims[l] != sa->dims[l + 1]) return NULL;
        if (expr_eq(re->bg, sa->bg) || sp_value_is_background(re->bg, sa->bg))
            return sp_retain((SparseArray*)re);
    } else {
        if (!is_arraylike(e)) return NULL;
        /* a row of another shape makes a ragged (dense) answer */
        Expr* de = eval1("Dimensions", e);
        int64_t dd[SP_MAX_RANK];
        int nd = int_list(de, dd, SP_MAX_RANK);
        expr_free(de);
        if (nd != rr) return NULL;
        for (int l = 0; l < rr; l++) if (dd[l] != sa->dims[l + 1]) return NULL;
    }
    /* SparseArray[dense row, dims, bg]: the stored entries are those that
     * differ from sa's background */
    int64_t t = sp_total(rr, sa->dims + 1);
    if (t < 0 || t > SPS_DENSE_MAX) return NULL;
    Expr* dense = re ? sp_normal(re) : expr_copy((Expr*)e);
    if (!dense) return NULL;
    Expr** dl = malloc(sizeof(Expr*) * (size_t)rr);
    if (!dl) { expr_free(dense); return NULL; }
    for (int l = 0; l < rr; l++) dl[l] = expr_new_integer(sa->dims[l + 1]);
    Expr* dims = expr_new_function(expr_new_symbol(SYM_List), dl, (size_t)rr);
    free(dl);
    Expr* args[3] = { dense, dims, expr_copy(sa->bg) };
    Expr* c = expr_new_function(expr_new_symbol(SYM_SparseArray), args, 3);
    Expr* r = evaluate(c);
    expr_free(c);
    SparseArray* out = NULL;
    if (r->type == EXPR_SPARSE && r->data.sparse->rank == rr) {
        bool ok = true;
        for (int l = 0; l < rr; l++) if (r->data.sparse->dims[l] != sa->dims[l + 1]) ok = false;
        if (ok) out = sp_retain(r->data.sparse);
    }
    expr_free(r);
    return out;
}

static Expr* append_common(const Expr* call, bool at_end) {
    if (NARGS(call) != 2) return NULL;
    const SparseArray* sa = sp_of(ARG(call, 0));
    if (!sa) return NULL;
    const Expr* e = ARG(call, 1);
    SparseArray* blk;
    if (sa->rank == 1) {
        if (is_arraylike(e)) return NULL;
        int64_t one = 1;
        bool implicit = sp_value_is_background(e, sa->bg);
        SparsePattern* p = sp_pattern_new(1, implicit ? 0 : 1, 1, true);
        VBuf vb;
        if (!p || !vb_init(&vb, expr_kind(e), implicit ? 0 : 1)) { sp_pattern_release(p); return NULL; }
        p->rowptr[0] = 0;
        p->rowptr[1] = implicit ? 0 : 1;
        if (!implicit) { pat_set(p, 0, 0, 0); vb_put(&vb, 0, e); }
        void* v = vb.v;
        blk = sp_new(1, &one, expr_copy(sa->bg), p, vb.kind, v);
    } else {
        SparseArray* row = row_operand(sa, e);
        /* a row of another shape: the answer is a ragged, dense List */
        if (!row) return dense_answer(call, 0);
        blk = as_block(row);
        sp_release(row);
    }
    if (!blk) return NULL;
    const SparseArray* ops[2];
    ops[0] = at_end ? sa : blk;
    ops[1] = at_end ? blk : sa;
    Expr* r = join1(ops, 2, expr_copy(sa->bg));
    sp_release(blk);
    return r;
}

static Expr* nat_append(const Expr* call) { return append_common(call, true); }
static Expr* nat_prepend(const Expr* call) { return append_common(call, false); }

/* ======================================================================== */
/* Slicing: Take, Drop, Rest, Most, First, Last, Reverse                     */
/* ======================================================================== */

/* One level's selection: a strided range of `count` positions from `start`
 * (0-based), one fixed position, or an explicit index list. */
enum { SEL_RANGE, SEL_FIXED, SEL_LIST };
typedef struct {
    int      mode;
    int64_t  start, step, count;
    int64_t* list;               /* SEL_LIST: 0-based, `count` of them (owned) */
} Sel;

static void sel_all(Sel* s, int64_t d) { s->mode = SEL_RANGE; s->start = 0; s->step = 1; s->count = d; s->list = NULL; }

/* Position c of a level, mapped through a range selection: the new index, or
 * -1 when c is not selected. */
static inline int64_t sel_map(const Sel* s, int64_t c) {
    int64_t d = c - s->start;
    if (s->step == 1) return (d >= 0 && d < s->count) ? d : -1;
    if (s->step > 0) {
        if (d < 0 || d % s->step) return -1;
        d /= s->step;
    } else {
        if (d > 0 || (-d) % (-s->step)) return -1;
        d /= s->step;
    }
    return d < s->count ? d : -1;
}

/* The selection through the general Part path (sparse_part.c). */
static Expr* slice_part(const Expr* s, const Sel* sel, int n) {
    Expr* idx[SP_MAX_RANK] = { NULL };
    for (int l = 0; l < n; l++) {
        const Sel* q = &sel[l];
        if (q->mode == SEL_FIXED) idx[l] = expr_new_integer(q->start + 1);
        else if (q->mode == SEL_LIST) {
            Expr** it = malloc(sizeof(Expr*) * (size_t)(q->count ? q->count : 1));
            for (int64_t i = 0; i < q->count; i++) it[i] = expr_new_integer(q->list[i] + 1);
            idx[l] = expr_new_function(expr_new_symbol(SYM_List), it, (size_t)q->count);
            free(it);
        } else {
            Expr* a[3] = { expr_new_integer(q->start + 1),
                           expr_new_integer(q->start + (q->count - 1) * q->step + 1),
                           expr_new_integer(q->step) };
            idx[l] = expr_new_function(expr_new_symbol(SYM_Span), a, 3);
        }
    }
    bool err = false;
    Expr* r = sp_part(s, idx, (size_t)n, &err);
    for (int l = 0; l < n; l++) expr_free(idx[l]);
    return r;
}

typedef struct {
    const SparseArray* sa;
    const Sel*         s0;
    SparsePattern*     out;
    VBuf*              vb;
} RowsCtx;

/* Copy output rows [lo, hi) of a row selection (out->rowptr is filled). */
static bool rows_copy(void* c, size_t lo, size_t hi) {
    RowsCtx* t = c;
    const SparsePattern* p = t->sa->pat;
    size_t isz = (p->idx32 ? sizeof(int32_t) : sizeof(int64_t)) * (size_t)p->width;
    for (size_t i = lo; i < hi; i++) {
        int64_t R = t->s0->start + (int64_t)i * t->s0->step;
        int64_t src = p->rowptr[R], len = p->rowptr[R + 1] - src, q = t->out->rowptr[i];
        if (!len) continue;
        memcpy((char*)t->out->colidx + (size_t)q * isz, (const char*)p->colidx + (size_t)src * isz,
               (size_t)len * isz);
        vb_range(t->vb, q, t->sa, src, len);
    }
    return true;
}

/* Rows of a matrix-or-higher array (every deeper level whole): each selected
 * row's index and value segments are copied as blocks, rows in parallel.
 * Take / Drop / Rest / Most / Reverse of rows. */
static Expr* slice_rows(const SparseArray* sa, const Sel* s0) {
    const SparsePattern* p = sa->pat;
    int64_t nr = s0->count;
    int64_t od[SP_MAX_RANK];
    memcpy(od, sa->dims, sizeof(int64_t) * (size_t)sa->rank);
    od[0] = nr;
    int64_t nnz = 0;
    for (int64_t i = 0; i < nr; i++) {
        int64_t R = s0->start + i * s0->step;
        nnz += p->rowptr[R + 1] - p->rowptr[R];
    }
    SparsePattern* out = sp_pattern_new(nr, nnz, p->width, p->idx32);
    VBuf vb;
    if (!out || !vb_init(&vb, sa->kind, nnz)) { sp_pattern_release(out); return NULL; }
    out->rowptr[0] = 0;
    for (int64_t i = 0; i < nr; i++) {
        int64_t R = s0->start + i * s0->step;
        out->rowptr[i + 1] = out->rowptr[i] + p->rowptr[R + 1] - p->rowptr[R];
    }
    RowsCtx t = { sa, s0, out, &vb };
    /* boxed values are reference-counted without atomics: one thread */
    par_range((size_t)nr, sa->kind != SPV_EXPR ? nnz : 0, rows_copy, &t);
    return finish(sa->rank, od, expr_copy(sa->bg), out, &vb);
}

/* A vector range of step +1 or -1: a binary search for its ends, then the
 * entries in order (or reversed). */
static Expr* slice_vec_unit(const SparseArray* sa, const Sel* s0) {
    const SparsePattern* p = sa->pat;
    int64_t a = s0->step > 0 ? s0->start : s0->start - s0->count + 1;   /* lowest selected */
    int64_t k0 = lower_bound0(p, 0, p->nnz, a);
    int64_t k1 = lower_bound0(p, k0, p->nnz, a + s0->count);
    int64_t nnz = k1 - k0, len = s0->count;
    SparsePattern* out = sp_pattern_new(1, nnz, 1, sp_dims_idx32(1, &len));
    VBuf vb;
    if (!out || !vb_init(&vb, sa->kind, nnz)) { sp_pattern_release(out); return NULL; }
    out->rowptr[0] = 0;
    out->rowptr[1] = nnz;
    if (s0->step > 0) {
        for (int64_t k = k0; k < k1; k++) pat_set(out, k - k0, 0, sp_col(p, k, 0) - a);
        vb_range(&vb, 0, sa, k0, nnz);
    } else {
        for (int64_t k = k1 - 1, q = 0; k >= k0; k--, q++) {
            pat_set(out, q, 0, s0->start - sp_col(p, k, 0));
            vb_from(&vb, q, sa, k);
        }
    }
    return finish(1, &len, expr_copy(sa->bg), out, &vb);
}

/* The direct path: level 0 a range or a fixed index; level 1 (the first
 * tuple coordinate) a range of either direction; deeper levels ascending
 * ranges. Entries are emitted in canonical order, so no sort is needed and a
 * row slice reads only its rows. */
static Expr* slice_fast(const SparseArray* sa, const Sel* sel) {
    int rank = sa->rank;
    const SparsePattern* p = sa->pat;
    bool fixed0 = rank >= 2 && sel[0].mode == SEL_FIXED;
    /* selections of the tuple coordinates, j = 0 .. width-1 */
    const Sel* ts = rank >= 2 ? sel + 1 : sel;
    int width = p->width;
    int orank = fixed0 ? rank - 1 : rank;
    int64_t od[SP_MAX_RANK];
    for (int l = 0, o = 0; l < rank; l++) if (!(l == 0 && fixed0)) od[o++] = sel[l].count;
    int64_t onrows = orank >= 2 ? od[0] : 1;
    int owidth = orank >= 2 ? orank - 1 : 1;
    int64_t nout_rows = rank >= 2 ? (fixed0 ? 1 : sel[0].count) : 1;
    if (onrows > SPS_ROWS_MAX) return NULL;

    /* tuple coordinate 0 bounds, for the binary search */
    const Sel* t0 = &ts[0];
    int64_t lo0 = t0->step > 0 ? t0->start : t0->start + (t0->count - 1) * t0->step;
    int64_t hi0 = t0->step > 0 ? t0->start + (t0->count - 1) * t0->step : t0->start;

    /* pass 1: an upper bound on the entries */
    int64_t* rlo = malloc(sizeof(int64_t) * (size_t)nout_rows);
    int64_t* rhi = malloc(sizeof(int64_t) * (size_t)nout_rows);
    if (!rlo || !rhi) { free(rlo); free(rhi); return NULL; }
    int64_t cap = 0;
    for (int64_t i = 0; i < nout_rows; i++) {
        int64_t R = rank >= 2 ? (fixed0 ? sel[0].start : sel[0].start + i * sel[0].step) : 0;
        rlo[i] = lower_bound0(p, p->rowptr[R], p->rowptr[R + 1], lo0);
        rhi[i] = lower_bound0(p, rlo[i], p->rowptr[R + 1], hi0 + 1);
        cap += rhi[i] - rlo[i];
    }
    SparsePattern* out = sp_pattern_new(onrows, cap, owidth, sp_dims_idx32(orank, od));
    int64_t* src = malloc(sizeof(int64_t) * (size_t)(cap ? cap : 1));
    if (!out || !src) { sp_pattern_release(out); free(src); free(rlo); free(rhi); return NULL; }
    memset(out->rowptr, 0, sizeof(int64_t) * (size_t)(onrows + 1));

    int64_t q = 0;
    int64_t nc[SP_MAX_RANK];
    for (int64_t i = 0; i < nout_rows; i++) {
        int64_t a = rlo[i], b = rhi[i];
        /* walk the row's groups of equal tuple coordinate 0 in the selection's
         * direction, each group forwards */
        int64_t k = t0->step > 0 ? a : b;
        while (t0->step > 0 ? k < b : k > a) {
            int64_t g0, g1;
            if (t0->step > 0) {
                g0 = k; g1 = k + 1;
                if (width > 1) while (g1 < b && sp_col(p, g1, 0) == sp_col(p, g0, 0)) g1++;
                k = g1;
            } else {
                g1 = k; g0 = k - 1;
                if (width > 1) while (g0 > a && sp_col(p, g0 - 1, 0) == sp_col(p, g1 - 1, 0)) g0--;
                k = g0;
            }
            int64_t m0 = sel_map(t0, sp_col(p, g0, 0));
            if (m0 < 0) continue;
            for (int64_t e = g0; e < g1; e++) {
                bool keep = true;
                for (int j = 1; j < width && keep; j++) {
                    nc[j] = sel_map(&ts[j], sp_col(p, e, j));
                    keep = nc[j] >= 0;
                }
                if (!keep) continue;
                nc[0] = m0;
                /* output coordinates: (row, tuple) */
                if (fixed0 && orank >= 2) {
                    out->rowptr[nc[0] + 1]++;
                    for (int j = 1; j < width; j++) pat_set(out, q, j - 1, nc[j]);
                } else {
                    if (orank >= 2) out->rowptr[i + 1]++;
                    for (int j = 0; j < width; j++) pat_set(out, q, j, nc[j]);
                }
                src[q++] = e;
            }
        }
    }
    free(rlo); free(rhi);
    out->nnz = q;
    if (orank >= 2) for (int64_t r = 0; r < onrows; r++) out->rowptr[r + 1] += out->rowptr[r];
    else { out->rowptr[0] = 0; out->rowptr[1] = q; }
    VBuf vb;
    if (!vb_init(&vb, sa->kind, q)) { sp_pattern_release(out); free(src); return NULL; }
    if (sa->kind == SPV_EXPR) for (int64_t e = 0; e < q; e++) vb_from(&vb, e, sa, src[e]);
    else {
        /* runs of consecutive sources copy as blocks */
        int64_t e = 0;
        while (e < q) {
            int64_t f = e + 1;
            while (f < q && src[f] == src[f - 1] + 1) f++;
            vb_range(&vb, e, sa, src[e], f - e);
            e = f;
        }
    }
    free(src);
    return finish(orank, od, expr_copy(sa->bg), out, &vb);
}

/* Select sel[0..n) of s (levels past n: All). Takes the LIST buffers. An
 * empty selection is a dense empty List: the call is answered on the Normal. */
static Expr* slice(const Expr* call, const Expr* s, Sel* sel, int n) {
    const SparseArray* sa = sp_of(s);
    int rank = sa->rank;
    Expr* r = NULL;
    Sel full[SP_MAX_RANK];
    for (int l = 0; l < rank; l++) {
        if (l < n) full[l] = sel[l];
        else sel_all(&full[l], sa->dims[l]);
    }
    bool fast = true, empty = false, allfixed = true;
    for (int l = 0; l < rank; l++) {
        if (full[l].mode != SEL_FIXED) { allfixed = false; if (full[l].count <= 0) empty = true; }
        if (full[l].mode == SEL_LIST) fast = false;
        if (full[l].mode == SEL_FIXED && !(l == 0 && rank >= 2)) fast = false;
        if (full[l].mode == SEL_RANGE && l >= 2 && full[l].step < 0) fast = false;
    }
    if (empty) {
        r = dense_answer(call, 0);
    } else if (allfixed) {
        int64_t c[SP_MAX_RANK];
        for (int l = 0; l < rank; l++) c[l] = full[l].start;
        r = sp_get(sa, c);
    } else if (fast) {
        bool whole_inner = true;
        for (int l = 1; l < rank; l++)
            if (!(full[l].mode == SEL_RANGE && full[l].start == 0 && full[l].step == 1 &&
                  full[l].count == sa->dims[l])) whole_inner = false;
        if (rank >= 2 && whole_inner && full[0].mode == SEL_RANGE) r = slice_rows(sa, &full[0]);
        else if (rank == 1 && (full[0].step == 1 || full[0].step == -1)) r = slice_vec_unit(sa, &full[0]);
        else r = slice_fast(sa, full);
    } else {
        r = slice_part(s, full, rank);
    }
    for (int l = 0; l < n; l++) free(sel[l].list);
    return r;
}

/* Normalise a 1-based, possibly negative position against a length. */
static bool pos_norm(int64_t i, int64_t d, int64_t* out) {
    if (i < 0) i = d + i + 1;
    if (i < 1 || i > d) return false;
    *out = i;
    return true;
}

/* Take[..., spec] for one level of length d. */
static bool take_spec(const Expr* e, int64_t d, Sel* s) {
    int64_t n, sp[3];
    s->list = NULL;
    if (get_i64(e, &n)) {
        if (n > d || n < -d) return false;
        s->mode = SEL_RANGE; s->step = 1;
        if (n >= 0) { s->start = 0; s->count = n; }
        else { s->start = d + n; s->count = -n; }
        return true;
    }
    if (is_sym(e, SYM_All)) { sel_all(s, d); return true; }
    if (is_sym(e, SYM_None)) { sel_all(s, d); s->count = 0; return true; }
    int k = int_list(e, sp, 3);
    if (k < 1) return false;
    int64_t a, b, st = k == 3 ? sp[2] : 1;
    if (!pos_norm(sp[0], d, &a) || !pos_norm(k >= 2 ? sp[1] : sp[0], d, &b) || st == 0) return false;
    s->mode = SEL_RANGE; s->start = a - 1; s->step = st;
    if ((st > 0 && a > b) || (st < 0 && a < b)) { s->count = 0; return true; }
    s->count = (st > 0 ? b - a : a - b) / (st > 0 ? st : -st) + 1;
    return true;
}

/* The complement of positions keep[] (0-based flags) as a selection. */
static bool sel_from_flags(const bool* drop, int64_t d, Sel* s) {
    int64_t cnt = 0;
    for (int64_t i = 0; i < d; i++) if (!drop[i]) cnt++;
    if (cnt == 0) { sel_all(s, d); s->count = 0; return true; }
    /* one contiguous run is a range */
    int64_t f = 0;
    while (drop[f]) f++;
    bool run = true;
    for (int64_t i = f; i < f + cnt; i++) if (drop[i]) run = false;
    if (run) { s->mode = SEL_RANGE; s->start = f; s->step = 1; s->count = cnt; s->list = NULL; return true; }
    s->mode = SEL_LIST; s->count = cnt; s->step = 1; s->start = 0;
    s->list = malloc(sizeof(int64_t) * (size_t)cnt);
    if (!s->list) return false;
    for (int64_t i = 0, q = 0; i < d; i++) if (!drop[i]) s->list[q++] = i;
    return true;
}

/* Drop[..., spec] for one level of length d. */
static bool drop_spec(const Expr* e, int64_t d, Sel* s) {
    int64_t n, sp[3];
    s->list = NULL;
    if (get_i64(e, &n)) {
        if (n > d || n < -d) return false;
        s->mode = SEL_RANGE; s->step = 1;
        if (n >= 0) { s->start = n; s->count = d - n; }
        else { s->start = 0; s->count = d + n; }
        return true;
    }
    if (is_sym(e, SYM_None)) { sel_all(s, d); return true; }
    if (is_sym(e, SYM_All)) { sel_all(s, d); s->count = 0; return true; }
    int k = int_list(e, sp, 3);
    if (k < 1) return false;
    int64_t a, b, st = k == 3 ? sp[2] : 1;
    if (!pos_norm(sp[0], d, &a) || !pos_norm(k >= 2 ? sp[1] : sp[0], d, &b) || st <= 0 || a > b)
        return false;
    if (d > SPS_EXPLICIT_MAX) return false;
    bool* drop = calloc((size_t)d, sizeof(bool));
    if (!drop) return false;
    for (int64_t i = a; i <= b; i += st) drop[i - 1] = true;
    bool ok = sel_from_flags(drop, d, s);
    free(drop);
    return ok;
}

static Expr* take_drop(const Expr* call, bool is_take) {
    const SparseArray* sa = sp_of(ARG(call, 0));
    size_t n = NARGS(call);
    if (!sa || n < 2 || n - 1 > (size_t)sa->rank) return NULL;
    Sel sel[SP_MAX_RANK];
    for (size_t l = 0; l + 1 < n; l++) {
        const Expr* e = ARG(call, l + 1);
        bool ok = e->type != EXPR_SPARSE &&
                  (is_take ? take_spec(e, sa->dims[l], &sel[l]) : drop_spec(e, sa->dims[l], &sel[l]));
        if (!ok) {
            for (size_t q = 0; q < l; q++) free(sel[q].list);
            return NULL;
        }
    }
    return slice(call, ARG(call, 0), sel, (int)n - 1);
}

static Expr* nat_take(const Expr* call) { return take_drop(call, true); }
static Expr* nat_drop(const Expr* call) { return take_drop(call, false); }

/* Rest / Most / First / Last: one level-1 selection. */
static Expr* level1(const Expr* call, int which) {
    const SparseArray* sa = sp_of(ARG(call, 0));
    if (!sa || NARGS(call) != 1) return NULL;
    int64_t d = sa->dims[0];
    Sel s;
    s.list = NULL;
    switch (which) {
        case 0: if (d < 1) return NULL; s.mode = SEL_RANGE; s.start = 1; s.step = 1; s.count = d - 1; break;
        case 1: if (d < 1) return NULL; s.mode = SEL_RANGE; s.start = 0; s.step = 1; s.count = d - 1; break;
        case 2: if (d < 1) return NULL; s.mode = SEL_FIXED; s.start = 0; s.step = 1; s.count = 1; break;
        default: if (d < 1) return NULL; s.mode = SEL_FIXED; s.start = d - 1; s.step = 1; s.count = 1; break;
    }
    return slice(call, ARG(call, 0), &s, 1);
}

static Expr* nat_rest(const Expr* call) { return level1(call, 0); }
static Expr* nat_most(const Expr* call) { return level1(call, 1); }
static Expr* nat_first(const Expr* call) { return level1(call, 2); }
static Expr* nat_last(const Expr* call) { return level1(call, 3); }

/* Reverse[s], Reverse[s, n], Reverse[s, {n1, n2, ...}]. */
static Expr* nat_reverse(const Expr* call) {
    const SparseArray* sa = sp_of(ARG(call, 0));
    if (!sa || NARGS(call) > 2) return NULL;
    int64_t lv[SP_MAX_RANK];
    int nl = 1;
    lv[0] = 1;
    if (NARGS(call) == 2) {
        if (get_i64(ARG(call, 1), &lv[0])) nl = 1;
        else if ((nl = int_list(ARG(call, 1), lv, SP_MAX_RANK)) < 0) return NULL;
    }
    Sel sel[SP_MAX_RANK];
    for (int l = 0; l < sa->rank; l++) sel_all(&sel[l], sa->dims[l]);
    bool any = false;
    for (int i = 0; i < nl; i++) {
        if (lv[i] < 1) return NULL;
        if (lv[i] > sa->rank) continue;              /* deeper than the array: no change */
        Sel* s = &sel[lv[i] - 1];
        s->start = s->count - 1;
        s->step = -1;
        any = true;
    }
    if (!any) return expr_copy(ARG(call, 0));
    return slice(call, ARG(call, 0), sel, sa->rank);
}

/* ======================================================================== */
/* Partition                                                                */
/* ======================================================================== */

/* Partition[s, n] and Partition[s, n, d]: blocks of n level-1 parts at
 * offset d, as one more level. */
static Expr* nat_partition(const Expr* call) {
    const SparseArray* sa = sp_of(ARG(call, 0));
    if (!sa || NARGS(call) < 2 || NARGS(call) > 3) return NULL;
    int64_t n, d;
    /* a block per level gives a dense nest of Lists */
    if (head_is(ARG(call, 1), SYM_List)) return dense_answer(call, 0);
    if (!get_i64(ARG(call, 1), &n) || n < 1) return NULL;
    if (NARGS(call) == 3) { if (!get_i64(ARG(call, 2), &d) || d < 1) return NULL; }
    else d = n;
    if (sa->rank + 1 > SP_MAX_RANK) return NULL;
    int64_t len = sa->dims[0];
    if (len < n) return dense_answer(call, 0);           /* {} */
    int64_t nb = (len - n) / d + 1;
    int rank = sa->rank + 1;
    int64_t nd[SP_MAX_RANK];
    nd[0] = nb; nd[1] = n;
    for (int l = 1; l < sa->rank; l++) nd[l + 1] = sa->dims[l];
    const SparsePattern* p = sa->pat;
    int width = rank - 1;
    /* entries of level-1 position i, as an index range */
    bool vec = sa->rank == 1;
    /* pass 1: count */
    int64_t nnz = 0;
    for (int64_t b = 0; b < nb; b++) {
        int64_t i0 = b * d, i1 = i0 + n;
        if (vec) nnz += lower_bound0(p, 0, p->nnz, i1) - lower_bound0(p, 0, p->nnz, i0);
        else nnz += p->rowptr[i1] - p->rowptr[i0];
    }
    SparsePattern* out = sp_pattern_new(nb, nnz, width, sp_dims_idx32(rank, nd));
    VBuf vb;
    if (!out || !vb_init(&vb, sa->kind, nnz)) { sp_pattern_release(out); return NULL; }
    int64_t q = 0;
    out->rowptr[0] = 0;
    for (int64_t b = 0; b < nb; b++) {
        int64_t i0 = b * d, i1 = i0 + n;
        if (vec) {
            int64_t k0 = lower_bound0(p, 0, p->nnz, i0), k1 = lower_bound0(p, k0, p->nnz, i1);
            for (int64_t k = k0; k < k1; k++) pat_set(out, q + k - k0, 0, sp_col(p, k, 0) - i0);
            vb_range(&vb, q, sa, k0, k1 - k0);
            q += k1 - k0;
        } else {
            for (int64_t r = i0; r < i1; r++) {
                int64_t k0 = p->rowptr[r], k1 = p->rowptr[r + 1];
                for (int64_t k = k0; k < k1; k++) {
                    pat_set(out, q + k - k0, 0, r - i0);
                    for (int j = 0; j < p->width; j++) pat_set(out, q + k - k0, j + 1, sp_col(p, k, j));
                }
                vb_range(&vb, q, sa, k0, k1 - k0);
                q += k1 - k0;
            }
        }
        out->rowptr[b + 1] = q;
    }
    return finish(rank, nd, expr_copy(sa->bg), out, &vb);
}

/* ======================================================================== */
/* Diagonal, DiagonalMatrix, IdentityMatrix                                 */
/* ======================================================================== */

typedef struct { const SparsePattern* p; int64_t* src; int64_t i0, k; } DiagCtx;

/* src[i] = the entry at (i0 + i, i0 + i + k), or -1, for i in [lo, hi). */
static bool diag_find(void* c, size_t lo, size_t hi) {
    DiagCtx* t = c;
    const SparsePattern* p = t->p;
    for (size_t i = lo; i < hi; i++) {
        int64_t r = t->i0 + (int64_t)i, col = r + t->k;
        int64_t e = lower_bound0(p, p->rowptr[r], p->rowptr[r + 1], col);
        t->src[i] = (e < p->rowptr[r + 1] && sp_col(p, e, 0) == col) ? e : -1;
    }
    return true;
}

static Expr* nat_diagonal(const Expr* call) {
    const SparseArray* sa = sp_of(ARG(call, 0));
    if (!sa || sa->rank != 2 || NARGS(call) > 2) return NULL;
    int64_t k = 0;
    if (NARGS(call) == 2 && !get_i64(ARG(call, 1), &k)) return NULL;
    if (k > SPS_ROWS_MAX || k < -SPS_ROWS_MAX) return NULL;
    int64_t m = sa->dims[0], n = sa->dims[1];
    int64_t len = k >= 0 ? (m < n - k ? m : n - k) : (m + k < n ? m + k : n);
    if (len <= 0) return dense_answer(call, 0);          /* {} */
    int64_t i0 = k >= 0 ? 0 : -k;
    const SparsePattern* p = sa->pat;
    /* one binary search per row (threaded), then compact */
    int64_t* src = malloc(sizeof(int64_t) * (size_t)len);
    if (!src) return NULL;
    DiagCtx t = { p, src, i0, k };
    par_range((size_t)len, 8 * len, diag_find, &t);
    int64_t q = 0;
    for (int64_t i = 0; i < len; i++) if (src[i] >= 0) q++;
    SparsePattern* out = sp_pattern_new(1, q, 1, sp_dims_idx32(1, &len));
    if (!out) { free(src); return NULL; }
    q = 0;
    for (int64_t i = 0; i < len; i++)
        if (src[i] >= 0) { pat_set(out, q, 0, i); src[q++] = src[i]; }
    out->nnz = q;
    out->rowptr[0] = 0;
    out->rowptr[1] = q;
    VBuf vb;
    if (!vb_init(&vb, sa->kind, q)) { sp_pattern_release(out); free(src); return NULL; }
    for (int64_t e = 0; e < q; e++) vb_from(&vb, e, sa, src[e]);
    free(src);
    return finish(1, &len, expr_copy(sa->bg), out, &vb);
}

/* DiagonalMatrix[v], [v, k], [v, k, n], [v, k, {m, n}] for a sparse vector:
 * the answer's background is 0, so where v's own background is not zero its
 * implicit elements become explicit diagonal entries. */
static Expr* nat_diagonalmatrix(const Expr* call) {
    const SparseArray* sa = sp_of(ARG(call, 0));
    size_t na = NARGS(call);
    if (!sa || sa->rank != 1 || na > 3) return NULL;
    int64_t k = 0, len = sa->dims[0];
    if (na >= 2 && !get_i64(ARG(call, 1), &k)) return NULL;
    int64_t sz = len + (k >= 0 ? k : -k);
    int64_t dims[2] = { sz, sz };
    if (na == 3) {
        int64_t mn[2];
        if (get_i64(ARG(call, 2), &mn[0])) { dims[0] = dims[1] = mn[0]; }
        else if (int_list(ARG(call, 2), mn, 2) == 2) { dims[0] = mn[0]; dims[1] = mn[1]; }
        else return NULL;
    }
    if (dims[0] < 1 || dims[1] < 1 || dims[0] > SPS_ROWS_MAX) return NULL;
    Expr* zero = expr_new_integer(0);
    bool fill = !sp_value_is_background(sa->bg, zero);
    if (fill && len > SPS_EXPLICIT_MAX) { expr_free(zero); return NULL; }
    const SparsePattern* p = sa->pat;
    int64_t cap = fill ? len : p->nnz;
    SpValKind kind = fill ? (p->nnz ? kind_join(sa->kind, expr_kind(sa->bg)) : expr_kind(sa->bg)) : sa->kind;
    SparsePattern* out = sp_pattern_new(dims[0], cap, 1, sp_dims_idx32(2, dims));
    VBuf vb;
    if (!out || !vb_init(&vb, kind, cap)) { sp_pattern_release(out); expr_free(zero); return NULL; }
    memset(out->rowptr, 0, sizeof(int64_t) * (size_t)(dims[0] + 1));
    int64_t q = 0, e = 0;
    for (int64_t i = 0; i < len; i++) {
        bool stored = e < p->nnz && sp_col(p, e, 0) == i;
        int64_t r = k >= 0 ? i : i - k, c = k >= 0 ? i + k : i;
        bool inside = r < dims[0] && c < dims[1];
        if ((stored || fill) && inside) {
            pat_set(out, q, 0, c);
            out->rowptr[r + 1]++;
            if (stored) vb_from(&vb, q, sa, e);
            else vb_put(&vb, q, sa->bg);
            q++;
        }
        if (stored) e++;
    }
    out->nnz = q;
    for (int64_t r = 0; r < dims[0]; r++) out->rowptr[r + 1] += out->rowptr[r];
    if (has_machine(sa)) {
        /* a machine vector gives a machine matrix: background 0., and every
         * value numericalised unless a symbolic value keeps them exact */
        expr_free(zero);
        zero = expr_new_real(0.0);
        if (!all_numbers(sa)) {
            /* values stay as they are */
        } else if (vb.kind == SPV_INT64) {
            double* d = malloc(sizeof(double) * (size_t)(q ? q : 1));
            if (!d) { sp_pattern_release(out); vb_free(&vb); expr_free(zero); return NULL; }
            for (int64_t i = 0; i < q; i++) d[i] = (double)((int64_t*)vb.v)[i];
            free(vb.v);
            vb.v = d;
            vb.kind = SPV_F64;
        } else if (vb.kind == SPV_EXPR) {
            Expr** ev = vb.v;
            for (int64_t i = 0; i < q; i++) {
                Expr* nv = eval1(SYM_N, ev[i]);
                expr_free(ev[i]);
                ev[i] = nv;
            }
        }
    }
    return finish(2, dims, zero, out, &vb);
}

typedef struct { SparsePattern* p; int64_t* v; int64_t d; } IdCtx;

/* Rows [lo, hi) of the identity: row i holds column i when i < d. */
static bool id_rows(void* c, size_t lo, size_t hi) {
    IdCtx* t = c;
    int64_t d = t->d;
    for (size_t i = lo; i < hi; i++) t->p->rowptr[i] = (int64_t)i < d ? (int64_t)i : d;
    int64_t e = (int64_t)hi < d ? (int64_t)hi : d;
    if (t->p->idx32) { int32_t* cc = t->p->colidx; for (int64_t i = (int64_t)lo; i < e; i++) cc[i] = (int32_t)i; }
    else             { int64_t* cc = t->p->colidx; for (int64_t i = (int64_t)lo; i < e; i++) cc[i] = i; }
    for (int64_t i = (int64_t)lo; i < e; i++) t->v[i] = 1;
    return true;
}

Expr* sp_identity_matrix(int64_t m, int64_t n) {
    if (m < 1 || n < 1 || m > SPS_ROWS_MAX) return NULL;
    int64_t dims[2] = { m, n };
    int64_t d = m < n ? m : n;
    SparsePattern* p = sp_pattern_new(m, d, 1, sp_dims_idx32(2, dims));
    int64_t* v = malloc(sizeof(int64_t) * (size_t)(d ? d : 1));
    if (!p || !v) { sp_pattern_release(p); free(v); return NULL; }
    IdCtx t = { p, v, d };
    nd_parallel_for((size_t)m, id_rows, &t);
    p->rowptr[m] = d;
    return sp_to_expr(sp_new(2, dims, expr_new_integer(0), p, SPV_INT64, v));
}

/* ======================================================================== */
/* KroneckerProduct                                                         */
/* ======================================================================== */

/* Does a boxed buffer hold a machine number? */
static bool has_machine(const SparseArray* sa) {
    if (sa->bg->type == EXPR_REAL || expr_kind(sa->bg) == SPV_C128) return true;
    if (sa->kind != SPV_EXPR) return sa->kind != SPV_INT64;
    for (int64_t k = 0; k < sa->pat->nnz; k++) {
        const Expr* v = ((Expr* const*)sa->vals)[k];
        if (v->type == EXPR_REAL || expr_kind(v) == SPV_C128) return true;
    }
    return false;
}

/* An exact or machine number (not a symbolic constant such as Pi). */
static bool is_number(const Expr* e) {
    if (e->type == EXPR_INTEGER || e->type == EXPR_BIGINT || e->type == EXPR_REAL) return true;
    if ((head_is(e, SYM_Rational) || head_is(e, SYM_Complex)) && NARGS(e) == 2)
        return is_number(ARG(e, 0)) && is_number(ARG(e, 1));
    return false;
}

/* Are the background and every value numbers? Such an array with a machine
 * value is a machine array as a whole: the operations that combine it
 * numericalise every value (a symbolic value keeps them exact). */
static bool all_numbers(const SparseArray* sa) {
    if (!is_number(sa->bg)) return false;
    if (sa->kind != SPV_EXPR) return true;
    for (int64_t k = 0; k < sa->pat->nnz; k++)
        if (!is_number(((Expr* const*)sa->vals)[k])) return false;
    return true;
}

typedef struct {
    const SparseArray *a, *b;
    SparsePattern*     out;
    void*              vals;
    SpValKind          kind;
    int64_t            mb, nb;
} KpCtx;

static inline void kp_val(const SparseArray* s, int64_t k, SpValKind kind, double* re, double* im) {
    switch (s->kind) {
        case SPV_INT64: *re = (double)((const int64_t*)s->vals)[k]; *im = 0; break;
        case SPV_F64:   *re = ((const double*)s->vals)[k]; *im = 0; break;
        default:        *re = ((const double*)s->vals)[2 * k]; *im = ((const double*)s->vals)[2 * k + 1]; break;
    }
    (void)kind;
}

static bool kp_rows(void* c, size_t lo, size_t hi) {
    KpCtx* t = c;
    const SparsePattern* pa = t->a->pat;
    const SparsePattern* pb = t->b->pat;
    for (size_t row = lo; row < hi; row++) {
        int64_t ia = (int64_t)row / t->mb, ib = (int64_t)row % t->mb;
        int64_t q = t->out->rowptr[row];
        for (int64_t ka = pa->rowptr[ia]; ka < pa->rowptr[ia + 1]; ka++) {
            int64_t cbase = sp_col(pa, ka, 0) * t->nb;
            double are = 0, aim = 0;
            int64_t ai = 0;
            if (t->kind == SPV_INT64) ai = ((const int64_t*)t->a->vals)[ka];
            else kp_val(t->a, ka, t->kind, &are, &aim);
            for (int64_t kb = pb->rowptr[ib]; kb < pb->rowptr[ib + 1]; kb++, q++) {
                pat_set(t->out, q, 0, cbase + sp_col(pb, kb, 0));
                if (t->kind == SPV_INT64) {
                    int64_t r;
                    if (ci_mul_i64(ai, ((const int64_t*)t->b->vals)[kb], &r)) return false;
                    ((int64_t*)t->vals)[q] = r;
                } else {
                    double bre, bim;
                    kp_val(t->b, kb, t->kind, &bre, &bim);
                    if (t->kind == SPV_F64) ((double*)t->vals)[q] = are * bre;
                    else {
                        ((double*)t->vals)[2 * q] = are * bre - aim * bim;
                        ((double*)t->vals)[2 * q + 1] = are * bim + aim * bre;
                    }
                }
            }
        }
    }
    return true;
}

/* A (x) B for two matrices over zero backgrounds. */
static SparseArray* kp2(const SparseArray* a, const SparseArray* b) {
    int64_t dims[2];
    if (ci_mul_i64(a->dims[0], b->dims[0], &dims[0]) || ci_mul_i64(a->dims[1], b->dims[1], &dims[1]))
        return NULL;
    int64_t nnz;
    if (ci_mul_i64(a->pat->nnz, b->pat->nnz, &nnz) || nnz > ((int64_t)1 << 31) ||
        dims[0] > SPS_ROWS_MAX) return NULL;
    /* A machine operand makes the whole product machine: every value is
     * numericalised and the background is 0. (as for the dense product of a
     * machine matrix). */
    bool boxed = a->kind == SPV_EXPR || b->kind == SPV_EXPR;
    bool boxed_machine = boxed && (has_machine(a) || has_machine(b)) && all_numbers(a) && all_numbers(b);
    SpValKind kind;
    Expr* bg;
    if (boxed) {
        kind = SPV_EXPR;
        bg = boxed_machine ? expr_new_real(0.0) : expr_new_integer(0);
    } else {
        bool machine = a->bg->type == EXPR_REAL || b->bg->type == EXPR_REAL ||
                       (a->kind != SPV_INT64 && a->pat->nnz) || (b->kind != SPV_INT64 && b->pat->nnz);
        bool cplx = (a->kind == SPV_C128 && a->pat->nnz) || (b->kind == SPV_C128 && b->pat->nnz);
        kind = cplx ? SPV_C128 : machine ? SPV_F64 : SPV_INT64;
        bg = machine ? expr_new_real(0.0) : expr_new_integer(0);
    }
    SparsePattern* out = sp_pattern_new(dims[0], nnz, 1, sp_dims_idx32(2, dims));
    VBuf vb;
    if (!out || !vb_init(&vb, kind, nnz)) { sp_pattern_release(out); expr_free(bg); return NULL; }
    int64_t mb = b->dims[0];
    out->rowptr[0] = 0;
    for (int64_t ia = 0; ia < a->dims[0]; ia++) {
        int64_t la = a->pat->rowptr[ia + 1] - a->pat->rowptr[ia];
        for (int64_t ib = 0; ib < mb; ib++) {
            int64_t row = ia * mb + ib;
            out->rowptr[row + 1] = out->rowptr[row] + la * (b->pat->rowptr[ib + 1] - b->pat->rowptr[ib]);
        }
    }
    bool ok;
    if (boxed) {
        ok = true;
        int64_t q = 0;
        for (int64_t row = 0; row < dims[0]; row++) {
            int64_t ia = row / mb, ib = row % mb;
            for (int64_t ka = a->pat->rowptr[ia]; ka < a->pat->rowptr[ia + 1]; ka++)
                for (int64_t kb = b->pat->rowptr[ib]; kb < b->pat->rowptr[ib + 1]; kb++, q++) {
                    pat_set(out, q, 0, sp_col(a->pat, ka, 0) * b->dims[1] + sp_col(b->pat, kb, 0));
                    Expr* ta[2] = { sp_value_expr(a, ka), sp_value_expr(b, kb) };
                    Expr* tc = expr_new_function(expr_new_symbol(SYM_Times), ta, 2);
                    if (boxed_machine) tc = mk_call1(SYM_N, tc);
                    ((Expr**)vb.v)[q] = evaluate(tc);
                    expr_free(tc);
                }
        }
    } else {
        KpCtx t = { a, b, out, vb.v, kind, mb, b->dims[1] };
        ok = nnz >= SPS_PAR_MIN ? nd_parallel_for((size_t)dims[0], kp_rows, &t)
                                : kp_rows(&t, 0, (size_t)dims[0]);
    }
    if (!ok) { sp_pattern_release(out); vb_free(&vb); expr_free(bg); return NULL; }
    sp_c128_canonicalise(&vb.kind, &vb.v, nnz);
    void* v = vb.v;
    return sp_new(2, dims, bg, out, vb.kind, v);
}

/* The dimensions of a dense array argument (List or NDArray), or -1 when it
 * is not a full array. */
static int dense_dims(const Expr* e, int64_t* d) {
    if (!head_is(e, SYM_List) && e->type != EXPR_NDARRAY) return -1;
    Expr* de = eval1("Dimensions", e);
    int n = int_list(de, d, SP_MAX_RANK);
    expr_free(de);
    return n;
}

/* sa over the background 0: the same array with its background elements
 * stored explicitly (a fresh reference), NULL when that is too large. */
static SparseArray* with_zero_bg(const SparseArray* sa) {
    if (sp_is_zero(sa->bg)) return sp_retain((SparseArray*)sa);
    int64_t t = sp_total(sa->rank, sa->dims);
    if (t < 0 || t > SPS_DENSE_MAX) return NULL;
    Expr* d = sp_normal(sa);
    if (!d) return NULL;
    Expr* r = eval1(SYM_SparseArray, d);
    expr_free(d);
    SparseArray* out = r->type == EXPR_SPARSE ? sp_retain(r->data.sparse) : NULL;
    expr_free(r);
    return out;
}

/* A matrix operand over the background 0: a SparseArray, or a dense matrix
 * made sparse. */
static SparseArray* matrix_operand(const Expr* e, bool first) {
    SparseArray* m = NULL;
    if (e->type == EXPR_SPARSE) {
        m = with_zero_bg(sp_of(e));
    } else {
        int64_t d[SP_MAX_RANK];
        int r = dense_dims(e, d);
        if (r != 1 && r != 2) return NULL;
        Expr* x = eval1(SYM_SparseArray, e);
        if (x->type == EXPR_SPARSE) m = sp_retain(x->data.sparse);
        expr_free(x);
    }
    if (!m || m->rank == 2) return m;
    if (m->rank != 1) { sp_release(m); return NULL; }
    /* a vector is a column when it leads the product, a row after that */
    int64_t d2[2] = { first ? m->dims[0] : 1, first ? 1 : m->dims[0] };
    Expr* mx = rekey(m, 2, d2, NULL);
    sp_release(m);
    if (!mx) return NULL;
    SparseArray* out = sp_retain(mx->data.sparse);
    expr_free(mx);
    return out;
}

static Expr* nat_kronecker(const Expr* call) {
    size_t n = NARGS(call);
    if (n < 2) return NULL;
    /* A non-zero background makes the product dense (every block is a
     * multiple of the other operand): computed over the explicit form and
     * answered as a dense matrix. */
    bool dense = false;
    int64_t total = 1;
    for (size_t i = 0; i < n; i++) {
        const SparseArray* s = sp_of(ARG(call, i));
        int64_t d[SP_MAX_RANK];
        int r;
        if (s) {
            r = s->rank;
            if (!sp_is_zero(s->bg)) dense = true;
            memcpy(d, s->dims, sizeof(int64_t) * (size_t)r);
        } else {
            r = dense_dims(ARG(call, i), d);
        }
        if (r == 1) d[1] = 1;
        else if (r != 2) return NULL;
        if (total >= 0 && (ci_mul_i64(total, d[0], &total) || ci_mul_i64(total, d[1], &total))) total = -1;
    }
    if (dense && (total < 0 || total > SPS_DENSE_MAX)) return NULL;
    SparseArray* acc = matrix_operand(ARG(call, 0), true);
    for (size_t i = 1; i < n && acc; i++) {
        SparseArray* b = matrix_operand(ARG(call, i), false);
        SparseArray* nx = b ? kp2(acc, b) : NULL;
        if (b) sp_release(b);
        sp_release(acc);
        acc = nx;
    }
    if (!acc) return NULL;
    if (dense) {
        Expr* r = sp_normal(acc);
        sp_release(acc);
        return r;
    }
    return sp_to_expr(acc);
}

/* ======================================================================== */
/* ArrayFlatten                                                             */
/* ======================================================================== */

static Expr* nat_arrayflatten(const Expr* call) {
    if (NARGS(call) != 1) return NULL;
    const Expr* M = ARG(call, 0);
    if (!head_is(M, SYM_List) || NARGS(M) == 0) return NULL;
    size_t br = NARGS(M), bc = 0;
    for (size_t i = 0; i < br; i++) {
        const Expr* row = ARG(M, i);
        if (!head_is(row, SYM_List) || NARGS(row) == 0) return NULL;
        if (i == 0) bc = NARGS(row);
        else if (NARGS(row) != bc) return NULL;
    }
    /* every block a matrix SparseArray over one background, or a scalar equal
     * to that background; anything else is a dense answer */
    const SparseArray* first = NULL;
    bool dense = false;
    for (size_t i = 0; i < br && !dense; i++)
        for (size_t j = 0; j < bc && !dense; j++) {
            const Expr* b = ARG(ARG(M, i), j);
            const SparseArray* s = sp_of(b);
            if (s) {
                if (s->rank != 2) return NULL;
                if (!first) first = s;
                else if (!expr_eq(s->bg, first->bg)) dense = true;
            } else if (is_arraylike(b)) {
                dense = true;
            }
        }
    if (!first) return NULL;
    for (size_t i = 0; i < br && !dense; i++)
        for (size_t j = 0; j < bc && !dense; j++) {
            const Expr* b = ARG(ARG(M, i), j);
            if (b->type != EXPR_SPARSE && !expr_eq(b, first->bg)) dense = true;
        }
    if (dense) return dense_answer(call, 2);
    if (br > SPS_EXPLICIT_MAX || bc > SPS_EXPLICIT_MAX) return NULL;
    /* block heights and widths, from the sparse blocks */
    int64_t* h = calloc(br, sizeof(int64_t));
    int64_t* w = calloc(bc, sizeof(int64_t));
    int64_t* roff = malloc(sizeof(int64_t) * (br + 1));
    int64_t* coff = malloc(sizeof(int64_t) * (bc + 1));
    Expr* result = NULL;
    SparsePattern* out = NULL;
    if (!h || !w || !roff || !coff) goto done;
    int64_t nnz = 0;
    SpValKind kind = SPV_INT64;
    bool any = false;
    for (size_t i = 0; i < br; i++)
        for (size_t j = 0; j < bc; j++) {
            const SparseArray* s = sp_of(ARG(ARG(M, i), j));
            if (!s) continue;
            if ((h[i] && h[i] != s->dims[0]) || (w[j] && w[j] != s->dims[1])) goto done;
            h[i] = s->dims[0];
            w[j] = s->dims[1];
            nnz += s->pat->nnz;
            if (s->pat->nnz) { kind = any ? kind_join(kind, s->kind) : s->kind; any = true; }
        }
    roff[0] = coff[0] = 0;
    for (size_t i = 0; i < br; i++) { if (!h[i]) goto done; roff[i + 1] = roff[i] + h[i]; }
    for (size_t j = 0; j < bc; j++) { if (!w[j]) goto done; coff[j + 1] = coff[j] + w[j]; }
    int64_t dims[2] = { roff[br], coff[bc] };
    if (dims[0] > SPS_ROWS_MAX) goto done;
    out = sp_pattern_new(dims[0], nnz, 1, sp_dims_idx32(2, dims));
    VBuf vb;
    if (!out || !vb_init(&vb, kind, nnz)) goto done;
    int64_t q = 0;
    out->rowptr[0] = 0;
    for (size_t i = 0; i < br; i++)
        for (int64_t r = 0; r < h[i]; r++) {
            for (size_t j = 0; j < bc; j++) {
                const SparseArray* s = sp_of(ARG(ARG(M, i), j));
                if (!s) continue;
                int64_t lo = row_lo(s, r), hi = row_hi(s, r);
                const SparsePattern* p = s->pat;
                if (out->idx32 && p->idx32) {
                    const int32_t* ci = p->colidx;
                    int32_t* co = out->colidx;
                    int32_t o = (int32_t)coff[j];
                    for (int64_t k = lo; k < hi; k++) co[q + k - lo] = ci[k] + o;
                } else {
                    for (int64_t k = lo; k < hi; k++) pat_set(out, q + k - lo, 0, sp_col(p, k, 0) + coff[j]);
                }
                vb_range(&vb, q, s, lo, hi - lo);
                q += hi - lo;
            }
            out->rowptr[roff[i] + r + 1] = q;
        }
    result = finish(2, dims, expr_copy(first->bg), out, &vb);
    out = NULL;
done:
    if (out) sp_pattern_release(out);
    free(h); free(w); free(roff); free(coff);
    return result;
}

/* ======================================================================== */
/* Upper/LowerTriangularize                                                 */
/* ======================================================================== */

/* The triangular part of a matrix over a zero background: a cut per row. */
static Expr* tri_core(const SparseArray* sa, int64_t k, bool upper) {
    const SparsePattern* p = sa->pat;
    int64_t m = sa->dims[0];
    int64_t* lo = malloc(sizeof(int64_t) * (size_t)(m ? m : 1));
    int64_t* hi = malloc(sizeof(int64_t) * (size_t)(m ? m : 1));
    if (!lo || !hi) { free(lo); free(hi); return NULL; }
    int64_t nnz = 0;
    for (int64_t r = 0; r < m; r++) {
        /* upper keeps c >= r + k, lower keeps c <= r + k */
        int64_t cut = lower_bound0(p, p->rowptr[r], p->rowptr[r + 1], upper ? r + k : r + k + 1);
        lo[r] = upper ? cut : p->rowptr[r];
        hi[r] = upper ? p->rowptr[r + 1] : cut;
        nnz += hi[r] - lo[r];
    }
    SparsePattern* out = sp_pattern_new(m, nnz, 1, p->idx32);
    VBuf vb;
    if (!out || !vb_init(&vb, sa->kind, nnz)) { sp_pattern_release(out); free(lo); free(hi); return NULL; }
    size_t isz = p->idx32 ? sizeof(int32_t) : sizeof(int64_t);
    int64_t q = 0;
    out->rowptr[0] = 0;
    for (int64_t r = 0; r < m; r++) {
        int64_t n = hi[r] - lo[r];
        if (n) memcpy((char*)out->colidx + (size_t)q * isz, (const char*)p->colidx + (size_t)lo[r] * isz,
                      (size_t)n * isz);
        vb_range(&vb, q, sa, lo[r], n);
        q += n;
        out->rowptr[r + 1] = q;
    }
    free(lo); free(hi);
    return finish(2, sa->dims, expr_copy(sa->bg), out, &vb);
}

static Expr* triangularize(const Expr* call, bool upper) {
    const SparseArray* sa = sp_of(ARG(call, 0));
    if (!sa || sa->rank != 2 || NARGS(call) > 2) return NULL;
    int64_t k = 0;
    if (NARGS(call) == 2 && !get_i64(ARG(call, 1), &k)) return NULL;
    if (k > SPS_ROWS_MAX || k < -SPS_ROWS_MAX) return NULL;   /* beyond any matrix: no overflow below */
    if (sp_is_zero(sa->bg)) return tri_core(sa, k, upper);
    /* the removed part becomes 0, so over a non-zero background the answer
     * is a dense matrix: cut the explicit form, then densify */
    SparseArray* z = with_zero_bg(sa);
    if (!z) return NULL;
    Expr* r = tri_core(z, k, upper);
    sp_release(z);
    if (!r) return NULL;
    Expr* d = sp_normal(sp_of(r));
    expr_free(r);
    return d;
}

static Expr* nat_upper(const Expr* call) { return triangularize(call, true); }
static Expr* nat_lower(const Expr* call) { return triangularize(call, false); }

/* ======================================================================== */
/* ArrayPad, PadRight, PadLeft                                              */
/* ======================================================================== */

/* sa with left[l] / right[l] positions added at level l (negative: removed),
 * the new positions holding `pad`. */
static Expr* pad_array(const SparseArray* sa, const int64_t* left, const int64_t* right, const Expr* pad) {
    int rank = sa->rank;
    int64_t nd[SP_MAX_RANK];
    for (int l = 0; l < rank; l++) {
        if (left[l] > SPS_ROWS_MAX || left[l] < -SPS_ROWS_MAX || right[l] > SPS_ROWS_MAX || right[l] < -SPS_ROWS_MAX) return NULL;
        nd[l] = sa->dims[l] + left[l] + right[l];
        if (nd[l] < 1) return NULL;
    }
    if (rank >= 2 && nd[0] > SPS_ROWS_MAX) return NULL;
    bool explicit_pad = !sp_value_is_background(pad, sa->bg);
    bool grows = false;
    for (int l = 0; l < rank; l++) if (left[l] > 0 || right[l] > 0) grows = true;
    if (!grows) explicit_pad = false;
    if (explicit_pad && rank > 2) return NULL;
    int64_t total = sp_total(rank, nd);
    if (explicit_pad && (total < 0 || total > SPS_EXPLICIT_MAX + sa->pat->nnz)) return NULL;
    const SparsePattern* p = sa->pat;
    int width = p->width;
    int64_t onrows = rank >= 2 ? nd[0] : 1;
    /* the source rows that land inside, for rank >= 2 */
    SpValKind kind = explicit_pad ? (p->nnz ? kind_join(sa->kind, expr_kind(pad)) : expr_kind(pad)) : sa->kind;
    /* capacity: kept entries plus the pad positions */
    int64_t cap = p->nnz + (explicit_pad ? total : 0);
    SparsePattern* out = sp_pattern_new(onrows, cap, width, sp_dims_idx32(rank, nd));
    VBuf vb;
    if (!out || !vb_init(&vb, kind, cap)) { sp_pattern_release(out); return NULL; }
    int64_t q = 0;
    out->rowptr[0] = 0;
    const int first = rank >= 2 ? 1 : 0;
    /* the inner coordinate range [ilo, ihi) of the new array the old one covers
     * at level `first` (the only inner level that may carry explicit pads) */
    int64_t ilo = left[first] > 0 ? left[first] : 0;
    int64_t ihi = sa->dims[first] + left[first];
    if (ihi > nd[first]) ihi = nd[first];
    for (int64_t nr = 0; nr < onrows; nr++) {
        int64_t r = rank >= 2 ? nr - left[0] : 0;        /* the source row */
        bool src_row = rank < 2 || (r >= 0 && r < sa->dims[0]);
        if (explicit_pad && !src_row) {
            for (int64_t c = 0; c < nd[1]; c++, q++) { pat_set(out, q, 0, c); vb_put(&vb, q, pad); }
            out->rowptr[nr + 1] = q;
            continue;
        }
        if (explicit_pad) for (int64_t c = 0; c < ilo; c++, q++) { pat_set(out, q, 0, c); vb_put(&vb, q, pad); }
        if (src_row) {
            for (int64_t k = p->rowptr[r]; k < p->rowptr[r + 1]; k++) {
                bool in = true;
                for (int j = 0; j < width && in; j++) {
                    int64_t c = sp_col(p, k, j) + left[first + j];
                    in = c >= 0 && c < nd[first + j];
                }
                if (!in) continue;
                for (int j = 0; j < width; j++) pat_set(out, q, j, sp_col(p, k, j) + left[first + j]);
                vb_from(&vb, q, sa, k);
                q++;
            }
        }
        if (explicit_pad) for (int64_t c = ihi; c < nd[first]; c++, q++) { pat_set(out, q, 0, c); vb_put(&vb, q, pad); }
        if (rank >= 2) out->rowptr[nr + 1] = q;
    }
    out->nnz = q;
    if (rank == 1) out->rowptr[1] = q;
    return finish(rank, nd, expr_copy(sa->bg), out, &vb);
}

/* ArrayPad[s, m], [s, {m, n}], [s, {{m1, n1}, ...}], each with an optional
 * scalar pad value (default 0). */
static Expr* nat_arraypad(const Expr* call) {
    const SparseArray* sa = sp_of(ARG(call, 0));
    if (!sa || NARGS(call) < 2 || NARGS(call) > 3) return NULL;
    int rank = sa->rank;
    int64_t left[SP_MAX_RANK] = { 0 }, right[SP_MAX_RANK] = { 0 };
    const Expr* spec = ARG(call, 1);
    int64_t m, two[2];
    if (get_i64(spec, &m)) {
        for (int l = 0; l < rank; l++) left[l] = right[l] = m;
    } else if (int_list(spec, two, 2) == 2) {
        for (int l = 0; l < rank; l++) { left[l] = two[0]; right[l] = two[1]; }
    } else if (head_is(spec, SYM_List) && NARGS(spec) >= 1 && NARGS(spec) <= (size_t)rank) {
        for (size_t l = 0; l < NARGS(spec); l++) {
            if (int_list(ARG(spec, l), two, 2) != 2) return NULL;
            left[l] = two[0]; right[l] = two[1];
        }
    } else {
        return NULL;
    }
    Expr* zero = NULL;
    const Expr* pad;
    if (NARGS(call) == 3) {
        pad = ARG(call, 2);
        if (is_arraylike(pad) || pad->type == EXPR_STRING) return NULL;
    } else {
        pad = zero = expr_new_integer(0);
    }
    Expr* r = pad_array(sa, left, right, pad);
    if (zero) expr_free(zero);
    return r;
}

/* PadRight[s, n], [s, {n1, ...}], each with an optional scalar pad. */
static Expr* pad_side(const Expr* call, bool at_right) {
    const SparseArray* sa = sp_of(ARG(call, 0));
    if (!sa || NARGS(call) < 2 || NARGS(call) > 3) return NULL;
    int rank = sa->rank;
    int64_t want[SP_MAX_RANK];
    int nw;
    /* a scalar length pads or cuts level 1 with scalars, not rows: on a
     * matrix that is a dense (ragged) List */
    if (get_i64(ARG(call, 1), &want[0])) {
        if (rank >= 2) return dense_answer(call, 0);
        nw = 1;
    }
    else if ((nw = int_list(ARG(call, 1), want, SP_MAX_RANK)) < 1 || nw > rank) return NULL;
    int64_t left[SP_MAX_RANK] = { 0 }, right[SP_MAX_RANK] = { 0 };
    for (int l = 0; l < nw; l++) {
        if (want[l] < 1) return NULL;
        int64_t delta = want[l] - sa->dims[l];
        if (at_right) right[l] = delta; else left[l] = delta;
    }
    Expr* zero = NULL;
    const Expr* pad;
    if (NARGS(call) == 3) {
        pad = ARG(call, 2);
        if (is_arraylike(pad)) return NULL;
    } else {
        pad = zero = expr_new_integer(0);
    }
    Expr* r = pad_array(sa, left, right, pad);
    if (zero) expr_free(zero);
    return r;
}

static Expr* nat_padright(const Expr* call) { return pad_side(call, true); }
static Expr* nat_padleft(const Expr* call) { return pad_side(call, false); }

/* ======================================================================== */
/* Map over a vector                                                        */
/* ======================================================================== */

/* Map[f, v] for a sparse vector: f over the stored values and once over the
 * background, by the Listable-threading tiers (a head's ND kernel over a
 * typed buffer, else f per value). A matrix is the gate's: Map walks its
 * rows. */
static Expr* nat_map(const Expr* call) {
    if (NARGS(call) != 2) return NULL;
    const SparseArray* sa = sp_of(ARG(call, 1));
    if (!sa || sa->rank != 1 || ARG(call, 0)->type == EXPR_SPARSE) return NULL;
    Expr* a[1] = { expr_copy(ARG(call, 1)) };
    Expr* fcall = expr_new_function(expr_copy(ARG(call, 0)), a, 1);
    Expr* r = sp_thread_listable(fcall);
    expr_free(fcall);
    if (r) return r;
    /* f giving a List on the background: the answer is a dense List of
     * Lists, not an array with list-valued elements */
    Expr* fa[1] = { expr_copy(sa->bg) };
    Expr* fb = expr_new_function(expr_copy(ARG(call, 0)), fa, 1);
    Expr* v = evaluate(fb);
    expr_free(fb);
    bool listy = is_arraylike(v);
    expr_free(v);
    return listy ? dense_answer(call, 0) : NULL;
}

/* ======================================================================== */
/* Registration                                                             */
/* ======================================================================== */

void sparse_struct_init(void) {
    symtab_set_sparse_native("Transpose", (const void*)nat_transpose);
    symtab_set_sparse_native("ConjugateTranspose", (const void*)nat_conjugate_transpose);
    symtab_set_sparse_native("Flatten", (const void*)nat_flatten);
    symtab_set_sparse_native("ArrayReshape", (const void*)nat_arrayreshape);
    symtab_set_sparse_native("Join", (const void*)nat_join);
    symtab_set_sparse_native("Append", (const void*)nat_append);
    symtab_set_sparse_native("Prepend", (const void*)nat_prepend);
    symtab_set_sparse_native("Take", (const void*)nat_take);
    symtab_set_sparse_native("Drop", (const void*)nat_drop);
    symtab_set_sparse_native("Rest", (const void*)nat_rest);
    symtab_set_sparse_native("Most", (const void*)nat_most);
    symtab_set_sparse_native("First", (const void*)nat_first);
    symtab_set_sparse_native("Last", (const void*)nat_last);
    symtab_set_sparse_native("Reverse", (const void*)nat_reverse);
    symtab_set_sparse_native("Partition", (const void*)nat_partition);
    symtab_set_sparse_native("Diagonal", (const void*)nat_diagonal);
    symtab_set_sparse_native("DiagonalMatrix", (const void*)nat_diagonalmatrix);
    symtab_set_sparse_native("KroneckerProduct", (const void*)nat_kronecker);
    symtab_set_sparse_native("ArrayFlatten", (const void*)nat_arrayflatten);
    symtab_set_sparse_native("UpperTriangularize", (const void*)nat_upper);
    symtab_set_sparse_native("LowerTriangularize", (const void*)nat_lower);
    symtab_set_sparse_native("ArrayPad", (const void*)nat_arraypad);
    symtab_set_sparse_native("PadRight", (const void*)nat_padright);
    symtab_set_sparse_native("PadLeft", (const void*)nat_padleft);
    symtab_set_sparse_native("Map", (const void*)nat_map);
}
