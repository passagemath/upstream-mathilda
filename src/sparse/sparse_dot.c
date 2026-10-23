/*
 * sparse_dot.c -- native Dot, MatrixPower and Inner[Times, _, _, Plus] over
 * SparseArray operands.
 *
 * Shapes (rank-1 and rank-2 operands, at least one of them sparse):
 *
 *   s . v, s . m, v . s, m . s   dense answer (a packed List when machine)
 *   s1 . s2                      SparseArray (Gustavson SpGEMM)
 *   s . sv, sv . s               SparseArray vector
 *   sv1 . sv2, sv . v            scalar
 *   Dot[a, b, c, ...]            pairwise, left to right; a chain that ends in
 *                                a dense vector runs right to left instead, so
 *                                every intermediate is a vector
 *
 * Explicit entries of a sparse result are STRUCTURAL: an output position is
 * explicit when some product a[i,k] b[k,j] of two explicit entries lands on
 * it, even if those products cancel to zero. That keeps the result's pattern
 * a function of the operands' patterns alone (and is what the two-pass sizes
 * below rely on).
 *
 * Value types. The operands' values join on Integer < Real < Complex (machine
 * values), and the kernels run in that type: int64 with overflow checks, then
 * double, then (re, im) pairs. A Rational, a big integer, an exact complex or
 * a symbolic value anywhere makes the product BOXED: per output element the
 * products are formed as Times[a, b] and summed by one evaluated Plus, which is
 * also where an int64 overflow lands. Only zero backgrounds take these paths;
 * a nonzero background makes every element of the product depend on whole
 * rows and columns, the answer is dense, and the transparency gate computes
 * it from the dense forms.
 *
 * Backgrounds of a sparse result:
 *   - typed result with explicit entries: the zero of the result type (0 or
 *     0.; a machine complex zero is the Real 0. here);
 *   - typed result without explicit entries, and any boxed result: exact 0;
 *   - an operand with no explicit entries at all short-circuits the product to
 *     an empty SparseArray (even against a dense operand), whose background is
 *     0. when either operand carries machine numbers and 0 otherwise.
 *
 * A shape mismatch reports Dot::dotsh here, with the operands as the user
 * wrote them, and leaves the call unevaluated -- it never densifies a large
 * array just to find out the shapes disagree.
 *
 * Threading: every kernel partitions rows (by stored entries, or by multiply
 * count for SpGEMM) through nd_parallel_for; the boxed path is serial.
 */

#include "sparse_internal.h"
#include "ndarray_internal.h"
#include "sym_names.h"
#include "symtab.h"
#include "eval.h"
#include "message.h"
#include "print.h"
#include "pack.h"
#include "checked_int.h"
#include "arithmetic.h"
#include <stdlib.h>
#include <string.h>

/* Largest dense answer (elements) these kernels will build. */
#define SPD_MAX_DENSE ((int64_t)1 << 30)

typedef struct { double re, im; } SpCplx;

enum { TY_INT = 0, TY_REAL = 1, TY_CPLX = 2, TY_BOX = 3 };

/* ======================================================================== */
/* Row-parallel driver                                                      */
/* ======================================================================== */

/*
 * Row r carries the key  mul * (key ? key[r] : 0) + r * scale  (mul 0 means
 * 1), strictly increasing
 * in r, with key(n) the total weight. nd_parallel_for splits [0, total) into
 * contiguous ranges and a range [lo, hi) owns the rows whose key falls in it,
 * so the rows are dealt out in contiguous blocks of about equal weight (stored
 * entries plus one per row, multiply count, ...). Below the threading
 * threshold the whole range runs serially in one call.
 */
typedef bool (*RowFn)(void* ctx, int64_t r0, int64_t r1);

typedef struct {
    int64_t n;
    const int64_t* key;
    int64_t scale;
    RowFn fn;
    void* ctx;
    int64_t mul;
} RowPar;

static inline int64_t rp_key(const RowPar* p, int64_t r) {
    return (p->key ? (p->mul > 1 ? p->mul : 1) * p->key[r] : 0) + r * p->scale;
}

/* First row whose key is >= x (n when none). */
static int64_t rp_lower(const RowPar* p, int64_t x) {
    int64_t lo = 0, hi = p->n;
    while (lo < hi) {
        int64_t mid = lo + (hi - lo) / 2;
        if (rp_key(p, mid) < x) lo = mid + 1;
        else hi = mid;
    }
    return lo;
}

static bool rowpar_chunk(void* c, size_t lo, size_t hi) {
    const RowPar* p = c;
    int64_t r0 = rp_lower(p, (int64_t)lo);
    int64_t r1 = (hi >= (size_t)rp_key(p, p->n)) ? p->n : rp_lower(p, (int64_t)hi);
    if (r0 >= r1) return true;
    return p->fn(p->ctx, r0, r1);
}

static bool par_rows_w(int64_t n, const int64_t* key, int64_t mul, int64_t scale,
                       RowFn fn, void* ctx) {
    if (n <= 0) return true;
    RowPar p = { n, key, scale, fn, ctx, mul };
    int64_t total = rp_key(&p, n);
    if (total <= 0) return fn(ctx, 0, n);
    return nd_parallel_for((size_t)total, rowpar_chunk, &p);
}

static bool par_rows(int64_t n, const int64_t* key, int64_t scale, RowFn fn, void* ctx) {
    return par_rows_w(n, key, 1, scale, fn, ctx);
}

/*
 * A scratch region reused across calls (the v . s accumulators). Evaluation is
 * single-threaded, so one region with a busy flag is enough; a nested or
 * oversized request gets a plain allocation. Reusing it is the point: the
 * accumulators are rewritten every call, and fresh pages cost a fault each.
 */
#define SPD_SCRATCH_MAX ((size_t)32 << 20)
static void*  g_scratch = NULL;
static size_t g_scratch_size = 0;
static bool   g_scratch_busy = false;

static void scratch_release(void) {
    free(g_scratch);
    g_scratch = NULL;
    g_scratch_size = 0;
}

static void* scratch_get(size_t size, bool* cached) {
    *cached = false;
    if (g_scratch_busy || size > SPD_SCRATCH_MAX) return malloc(size ? size : 1);
    if (!g_scratch_size && !g_scratch) atexit(scratch_release);
    if (g_scratch_size < size) {
        free(g_scratch);
        g_scratch = malloc(size ? size : 1);
        g_scratch_size = g_scratch ? size : 0;
        if (!g_scratch) return NULL;
    }
    g_scratch_busy = true;
    *cached = true;
    return g_scratch;
}

static void scratch_put(void* p, bool cached) {
    if (cached) g_scratch_busy = false;
    else free(p);
}

/* In-place ascending sort of a short int32 run (a row's output columns). */
static void sort_i32(int32_t* a, int64_t n) {
    while (n > 24) {
        /* median of three, Hoare partition */
        int32_t x = a[0], y = a[n / 2], z = a[n - 1];
        int32_t piv = x < y ? (y < z ? y : (x < z ? z : x)) : (x < z ? x : (y < z ? z : y));
        int64_t i = -1, j = n;
        for (;;) {
            do i++; while (a[i] < piv);
            do j--; while (a[j] > piv);
            if (i >= j) break;
            int32_t t = a[i]; a[i] = a[j]; a[j] = t;
        }
        /* recurse into the smaller half, loop on the larger */
        int64_t left = j + 1;
        if (left < n - left) { sort_i32(a, left); a += left; n -= left; }
        else { sort_i32(a + left, n - left); n = left; }
    }
    for (int64_t i = 1; i < n; i++) {
        int32_t v = a[i];
        int64_t j = i - 1;
        while (j >= 0 && a[j] > v) { a[j + 1] = a[j]; j--; }
        a[j + 1] = v;
    }
}

/* ======================================================================== */
/* Typed kernels (one instance per value type)                              */
/* ======================================================================== */

#define KT int64_t
#define KS(name) name##_i64
#define K_ZERO ((int64_t)0)
#define K_ONE ((int64_t)1)
#define K_FMA(acc, a, b) do { int64_t t_; \
        if (ci_mul_i64((a), (b), &t_) || ci_add_i64((acc), t_, &(acc))) ovf = true; } while (0)
#include "sparse_dot_kern.inc"

#define KT double
#define KS(name) name##_f64
#define K_ZERO 0.0
#define K_ONE 1.0
#define K_FMA(acc, a, b) ((acc) += (a) * (b))
#include "sparse_dot_kern.inc"

static const SpCplx SPC_ZERO = { 0.0, 0.0 };
static const SpCplx SPC_ONE = { 1.0, 0.0 };
#define KT SpCplx
#define KS(name) name##_c128
#define K_ZERO SPC_ZERO
#define K_ONE SPC_ONE
#define K_FMA(acc, a, b) do { SpCplx a_ = (a), b_ = (b); \
        (acc).re += a_.re * b_.re - a_.im * b_.im; \
        (acc).im += a_.re * b_.im + a_.im * b_.re; } while (0)
#include "sparse_dot_kern.inc"

/* ======================================================================== */
/* Operands                                                                 */
/* ======================================================================== */

typedef struct {
    const Expr* e;
    bool sparse;
    int rank;                /* 1 or 2 */
    int64_t dims[2];
    bool machine;            /* carries a machine number (value or background) */
    bool allnum;             /* every value a number (exact or machine) */
    bool exactnum;           /* ... some of them exact non-Integers */
    bool cplx;               /* ... some of them complex */
    const SparseArray* sa;   /* sparse */
    const Expr* nd;          /* dense: an EXPR_NDARRAY */
    Expr** leaves;           /* dense plain List: borrowed leaves, row-major */
} Opd;

static bool is_list(const Expr* e) {
    return e && e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL &&
           e->data.function.head->data.symbol.name == SYM_List;
}

static bool is_machine_complex(const Expr* e) {
    return e->type == EXPR_FUNCTION && e->data.function.arg_count == 2 &&
           e->data.function.head->type == EXPR_SYMBOL &&
           e->data.function.head->data.symbol.name == SYM_Complex &&
           e->data.function.args[0]->type == EXPR_REAL &&
           e->data.function.args[1]->type == EXPR_REAL;
}

/* Fold one value into the operand's flags. */
static void opd_note(Opd* o, const Expr* e) {
    double re, im;
    if (e->type == EXPR_INTEGER) return;
    if (e->type == EXPR_REAL) { o->machine = true; return; }
    if (is_machine_complex(e)) { o->machine = true; o->cplx = true; return; }
    if (power_scalar_components(e, &re, &im)) {
        o->exactnum = true;
        if (e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL &&
            e->data.function.head->data.symbol.name == SYM_Complex)
            o->cplx = true;
        return;
    }
    o->allnum = false;
}

/* The machine type a pair computes in, or TY_BOX. When both operands are
 * all-numeric, a machine number anywhere in the pair makes the product
 * machine (exact Rationals, big integers and exact complexes are converted);
 * without one, only Integers stay on a typed path (as exact int64). A
 * symbolic value anywhere makes the product boxed, and boxed products keep
 * every exact value exact. */
static int pair_ty(const Opd* a, const Opd* b) {
    if (!a->allnum || !b->allnum) return TY_BOX;
    bool mach = a->machine || b->machine;
    if (!mach && (a->exactnum || b->exactnum)) return TY_BOX;
    if (a->cplx || b->cplx) return TY_CPLX;
    return mach ? TY_REAL : TY_INT;
}

static void opd_free(Opd* o) {
    free(o->leaves);
    o->leaves = NULL;
}

/* Classify one Dot operand. False when it is not an array these kernels take
 * (wrong rank, ragged, a nonzero background, a float32 / bool buffer, ...). */
static bool opd_init(Opd* o, const Expr* e) {
    memset(o, 0, sizeof *o);
    o->e = e;
    if (e->type == EXPR_SPARSE) {
        const SparseArray* sa = sp_of(e);
        if (sa->rank < 1 || sa->rank > 2 || !sp_is_zero(sa->bg) || !sa->pat->idx32) return false;
        o->sparse = true;
        o->sa = sa;
        o->rank = sa->rank;
        for (int k = 0; k < sa->rank; k++) {
            o->dims[k] = sa->dims[k];
            if (sa->dims[k] <= 0 || sa->dims[k] > INT32_MAX) return false;
        }
        o->allnum = true;
        /* a zero background: 0 adds nothing, 0. makes the operand machine */
        if (sa->bg->type == EXPR_REAL || is_machine_complex(sa->bg)) o->machine = true;
        switch (sa->kind) {
            case SPV_INT64: break;
            case SPV_F64:  o->machine = true; break;
            case SPV_C128: o->machine = true; o->cplx = true; break;
            default: {
                Expr* const* v = sa->vals;
                for (int64_t k = 0; k < sa->pat->nnz; k++) opd_note(o, v[k]);
            }
        }
        return true;
    }
    if (e->type == EXPR_NDARRAY) {
        const NDArrayData* d = &e->data.ndarray;
        if (d->rank < 1 || d->rank > 2) return false;
        if (d->dtype == NDT_FLOAT64) o->machine = true;
        else if (d->dtype == NDT_COMPLEX64) o->machine = o->cplx = true;
        else if (d->dtype != NDT_INT64) return false;
        o->allnum = true;
        o->rank = d->rank;
        for (int k = 0; k < d->rank; k++) {
            o->dims[k] = d->dims[k];
            if (d->dims[k] <= 0) return false;
        }
        o->nd = e;
        return true;
    }
    if (!is_list(e)) return false;
    size_t n = e->data.function.arg_count;
    if (n == 0) return false;
    const Expr* first = e->data.function.args[0];
    if (is_list(first)) {
        size_t m = first->data.function.arg_count;
        if (m == 0) return false;
        o->rank = 2; o->dims[0] = (int64_t)n; o->dims[1] = (int64_t)m;
        if ((int64_t)n > SPD_MAX_DENSE / (int64_t)m) return false;
        o->leaves = malloc(sizeof(Expr*) * n * m);
        if (!o->leaves) return false;
        for (size_t i = 0; i < n; i++) {
            const Expr* row = e->data.function.args[i];
            if (!is_list(row) || row->data.function.arg_count != m) { opd_free(o); return false; }
            for (size_t j = 0; j < m; j++) {
                Expr* x = row->data.function.args[j];
                if (is_list(x) || x->type == EXPR_NDARRAY || x->type == EXPR_SPARSE) {
                    opd_free(o); return false;
                }
                o->leaves[i * m + j] = x;
            }
        }
    } else {
        o->rank = 1; o->dims[0] = (int64_t)n;
        o->leaves = malloc(sizeof(Expr*) * n);
        if (!o->leaves) return false;
        for (size_t i = 0; i < n; i++) {
            Expr* x = e->data.function.args[i];
            if (is_list(x) || x->type == EXPR_NDARRAY || x->type == EXPR_SPARSE) {
                opd_free(o); return false;
            }
            o->leaves[i] = x;
        }
    }
    int64_t total = o->rank == 2 ? o->dims[0] * o->dims[1] : o->dims[0];
    o->allnum = true;
    for (int64_t k = 0; k < total; k++) opd_note(o, o->leaves[k]);
    return true;
}

static int64_t opd_count(const Opd* o) {
    if (o->sparse) return o->sa->pat->nnz;
    return o->rank == 2 ? o->dims[0] * o->dims[1] : o->dims[0];
}

/* The operand's values (sparse: the explicit values in storage order; dense:
 * every element, row-major) as a buffer of type `ty`. Returns a borrowed
 * pointer when the storage already has that type, else a fresh buffer that is
 * also stored in *owned for the caller to free. NULL on allocation failure. */
static void* opd_typed(const Opd* o, int ty, void** owned) {
    *owned = NULL;
    int64_t n = opd_count(o);
    size_t esz = ty == TY_INT ? sizeof(int64_t) : ty == TY_REAL ? sizeof(double) : sizeof(SpCplx);
    if (o->sparse) {
        SpValKind k = o->sa->kind;
        if ((k == SPV_INT64 && ty == TY_INT) || (k == SPV_F64 && ty == TY_REAL) ||
            (k == SPV_C128 && ty == TY_CPLX))
            return o->sa->vals;
    } else if (o->nd) {
        NDType dt = o->nd->data.ndarray.dtype;
        if ((dt == NDT_INT64 && ty == TY_INT) || (dt == NDT_FLOAT64 && ty == TY_REAL) ||
            (dt == NDT_COMPLEX64 && ty == TY_CPLX))
            return o->nd->data.ndarray.data;
    }
    void* buf = malloc(esz * (size_t)(n ? n : 1));
    if (!buf) return NULL;
    for (int64_t k = 0; k < n; k++) {
        int64_t iv = 0; double re = 0, im = 0;
        if (o->sparse && o->sa->kind != SPV_EXPR) {
            switch (o->sa->kind) {
                case SPV_INT64: iv = ((const int64_t*)o->sa->vals)[k]; re = (double)iv; break;
                case SPV_F64:   re = ((const double*)o->sa->vals)[k]; break;
                default:        re = ((const double*)o->sa->vals)[2 * k];
                                im = ((const double*)o->sa->vals)[2 * k + 1]; break;
            }
        } else if (o->nd) {
            NDType dt = o->nd->data.ndarray.dtype;
            const void* d = o->nd->data.ndarray.data;
            if (dt == NDT_INT64) { iv = ((const int64_t*)d)[k]; re = (double)iv; }
            else if (dt == NDT_FLOAT64) { re = ((const double*)d)[k]; }
            else { re = ((const double*)d)[2 * k]; im = ((const double*)d)[2 * k + 1]; }
        } else {
            const Expr* x = o->sparse ? ((Expr* const*)o->sa->vals)[k] : o->leaves[k];
            if (x->type == EXPR_INTEGER) { iv = x->data.integer; re = (double)iv; }
            else if (is_machine_complex(x)) {
                re = x->data.function.args[0]->data.real;
                im = x->data.function.args[1]->data.real;
            } else if (!power_scalar_components(x, &re, &im)) {
                free(buf);
                return NULL;
            }
        }
        if (ty == TY_INT) ((int64_t*)buf)[k] = iv;
        else if (ty == TY_REAL) ((double*)buf)[k] = re;
        else { ((SpCplx*)buf)[k].re = re; ((SpCplx*)buf)[k].im = im; }
    }
    *owned = buf;
    return buf;
}

/* ======================================================================== */
/* Results                                                                  */
/* ======================================================================== */

static NDType ty_dtype(int ty) {
    return ty == TY_INT ? NDT_INT64 : ty == TY_REAL ? NDT_FLOAT64 : NDT_COMPLEX64;
}

static size_t ty_size(int ty) {
    return ty == TY_INT ? sizeof(int64_t) : ty == TY_REAL ? sizeof(double) : sizeof(SpCplx);
}

static Expr* scalar_expr(int ty, const void* v) {
    if (ty == TY_INT) return expr_new_integer(*(const int64_t*)v);
    if (ty == TY_REAL) return expr_new_real(*(const double*)v);
    const SpCplx* c = v;
    if (c->im == 0.0) return expr_new_real(c->re);
    Expr* args[2] = { expr_new_real(c->re), expr_new_real(c->im) };
    return expr_new_function(expr_new_symbol(SYM_Complex), args, 2);
}

/* A dense answer from a filled machine buffer (adopted): it inherits a visible
 * NDArray operand's presentation, and is a packed List otherwise. */
static Expr* dense_result(int rank, const int64_t* dims, int ty, void* buf, const Expr* like) {
    if (like && like->type == EXPR_NDARRAY)
        return expr_new_ndarray_like(like, rank, dims, buf, ty_dtype(ty));
    Expr* e = expr_new_ndarray_raw(rank, dims, buf, ty_dtype(ty));
    if (e) {
        e->data.ndarray.present_as = NDA_HEAD_LIST;
        pack_g_any_created = true;      /* arm the evaluator's packed gate */
    }
    return e;
}

/* The dense operand that decides a dense answer's presentation. */
static const Expr* like_of(const Opd* a, const Opd* b) {
    if (a->nd && a->nd->data.ndarray.present_as == NDA_HEAD_NDARRAY) return a->nd;
    if (b->nd && b->nd->data.ndarray.present_as == NDA_HEAD_NDARRAY) return b->nd;
    if (a->nd) return a->nd;
    if (b->nd) return b->nd;
    return NULL;
}

static Expr* zero_of(bool machine) {
    return machine ? expr_new_real(0.0) : expr_new_integer(0);
}

/* An empty SparseArray of the product's shape (rank 0: the scalar zero). */
static Expr* empty_product(const Opd* a, const Opd* b) {
    bool machine = a->machine || b->machine;
    int64_t dims[2];
    int rank = 0;
    if (a->rank == 2) dims[rank++] = a->dims[0];
    if (b->rank == 2) dims[rank++] = b->dims[1];
    if (rank == 0) return zero_of(machine);
    SparsePattern* pat = sp_pattern_new(rank == 2 ? dims[0] : 1, 0, 1, true);
    if (!pat) return NULL;
    for (int64_t r = 0; r <= pat->nrows; r++) pat->rowptr[r] = 0;
    return sp_to_expr(sp_new(rank, dims, zero_of(machine), pat, SPV_INT64, NULL));
}

/* A sparse answer over `pat` (adopted) with typed `vals` (adopted). */
static Expr* sparse_result(int rank, const int64_t* dims, SparsePattern* pat, int ty, void* vals) {
    SpValKind kind = ty == TY_INT ? SPV_INT64 : ty == TY_REAL ? SPV_F64 : SPV_C128;
    bool any = pat->nnz > 0;
    if (kind == SPV_C128) sp_c128_canonicalise(&kind, &vals, pat->nnz);
    Expr* bg = (any && ty != TY_INT) ? expr_new_real(0.0) : expr_new_integer(0);
    SparseArray* sa = sp_new(rank, dims, bg, pat, kind, vals);
    return sa ? sp_to_expr(sa) : NULL;
}

/* ======================================================================== */
/* Typed products                                                           */
/* ======================================================================== */

/* Dispatch one kernel over the three machine types. */
#define TY_CALL(ty, fn, ...) \
    ((ty) == TY_INT ? fn##_i64(__VA_ARGS__) : (ty) == TY_REAL ? fn##_f64(__VA_ARGS__) : fn##_c128(__VA_ARGS__))

/* a . b for machine values of type ty. NULL when the int64 instance overflowed
 * (or an allocation failed): the caller then takes the boxed path. */
static Expr* typed_product(const Opd* a, const Opd* b, int ty) {
    void* oa = NULL;
    void* ob = NULL;
    void* va = opd_typed(a, ty, &oa);
    void* vb = va ? opd_typed(b, ty, &ob) : NULL;
    Expr* out = NULL;
    size_t esz = ty_size(ty);
    if (!va || !vb) goto done;

    if (a->sparse && b->sparse) {
        const SparsePattern* pa = a->sa->pat;
        const SparsePattern* pb = b->sa->pat;
        if (a->rank == 1 && b->rank == 1) {
            unsigned char acc[sizeof(SpCplx)];
            if (TY_CALL(ty, svdot_sparse, pa->nnz, pa->colidx, va, pb->nnz, pb->colidx, vb,
                        (void*)acc))
                out = scalar_expr(ty, acc);
            goto done;
        }
        if (b->rank == 2) {
            /* rank-1 left: a vector is a one-row CSR matrix */
            SparsePattern* pat;
            void* vals;
            bool ok = TY_CALL(ty, spgemm, pa->nrows, b->dims[1], pa->rowptr, pa->colidx, va,
                              pb->rowptr, pb->colidx, vb, &pat, (void*)&vals);
            if (!ok) goto done;
            int64_t dims[2] = { a->rank == 2 ? a->dims[0] : b->dims[1], b->dims[1] };
            out = a->rank == 2 ? sparse_result(2, dims, pat, ty, vals)
                               : sparse_result(1, dims + 1, pat, ty, vals);
            goto done;
        }
        /* s . sv: a mask over the vector, then one pass over the matrix */
        int64_t n = a->dims[0], K = a->dims[1];
        void* xd = calloc((size_t)K, esz);
        unsigned char* xm = calloc((size_t)K, 1);
        void* yv = malloc(esz * (size_t)n);
        unsigned char* ym = malloc((size_t)n);
        bool ok = xd && xm && yv && ym;
        if (ok) {
            for (int64_t q = 0; q < pb->nnz; q++) {
                int32_t k = ((const int32_t*)pb->colidx)[q];
                memcpy((char*)xd + (size_t)k * esz, (const char*)vb + (size_t)q * esz, esz);
                xm[k] = 1;
            }
            ok = TY_CALL(ty, spsv, n, pa->rowptr, pa->colidx, va, xd, xm, yv, ym);
        }
        if (ok) {
            int64_t cnt = 0;
            for (int64_t i = 0; i < n; i++) cnt += ym[i];
            SparsePattern* pat = sp_pattern_new(1, cnt, 1, true);
            void* vals = malloc(esz * (size_t)(cnt ? cnt : 1));
            if (pat && vals) {
                pat->rowptr[0] = 0; pat->rowptr[1] = cnt;
                int64_t q = 0;
                for (int64_t i = 0; i < n; i++) {
                    if (!ym[i]) continue;
                    ((int32_t*)pat->colidx)[q] = (int32_t)i;
                    memcpy((char*)vals + (size_t)q * esz, (char*)yv + (size_t)i * esz, esz);
                    q++;
                }
                int64_t dims[1] = { n };
                out = sparse_result(1, dims, pat, ty, vals);
            } else {
                sp_pattern_release(pat);
                free(vals);
            }
        }
        free(xd); free(xm); free(yv); free(ym);
        goto done;
    }

    /* one sparse and one dense operand: a dense answer */
    const Opd* s = a->sparse ? a : b;
    const void* sv = a->sparse ? va : vb;
    const void* dv = a->sparse ? vb : va;
    const SparsePattern* ps = s->sa->pat;
    int64_t dims[2];
    int rank = 0;
    if (a->rank == 2) dims[rank++] = a->dims[0];
    if (b->rank == 2) dims[rank++] = b->dims[1];
    if (rank == 0) {
        unsigned char acc[sizeof(SpCplx)];
        if (TY_CALL(ty, svdot_dense, ps->nnz, ps->colidx, sv, dv, (void*)acc))
            out = scalar_expr(ty, acc);
        goto done;
    }
    int64_t total = rank == 2 ? dims[0] * dims[1] : dims[0];
    if (rank == 2 && dims[0] > SPD_MAX_DENSE / dims[1]) goto done;
    if (total > SPD_MAX_DENSE) goto done;
    void* y = malloc(esz * (size_t)total);
    if (!y) goto done;
    bool ok;
    if (a->sparse) {
        if (b->rank == 1)        /* s . v  (a is a matrix here: rank 0 handled above) */
            ok = TY_CALL(ty, spmv, ps->nrows, ps->rowptr, ps->colidx, sv, dv, y);
        else                     /* s . m, sv . m (a one-row CSR) */
            ok = TY_CALL(ty, spmm, ps->nrows, b->dims[1], ps->rowptr, ps->colidx, sv, dv, y);
    } else if (b->rank == 2) {
        if (a->rank == 1)        /* v . s */
            ok = TY_CALL(ty, vspm, ps->nrows, b->dims[1], ps->rowptr, ps->colidx, sv, dv, y);
        else                     /* m . s */
            ok = TY_CALL(ty, mspm, a->dims[0], a->dims[1], b->dims[1], ps->rowptr, ps->colidx,
                         sv, dv, y);
    } else {                     /* m . sv */
        ok = TY_CALL(ty, msv, a->dims[0], a->dims[1], ps->nnz, ps->colidx, sv, dv, y);
    }
    if (!ok) { free(y); goto done; }
    out = dense_result(rank, dims, ty, y, like_of(a, b));

done:
    free(oa);
    free(ob);
    return out;
}

/* ======================================================================== */
/* Boxed products                                                           */
/* ======================================================================== */

/* A generic CSR matrix with owned Expr values: the left operand as n x K
 * (a vector as one row), the right as K x p (a vector as one column). */
typedef struct {
    int64_t nr, nc, nnz;
    int64_t* rp;
    int64_t* ci;
    Expr** v;
} GCsr;

static void gcsr_free(GCsr* g) {
    if (g->v) for (int64_t k = 0; k < g->nnz; k++) expr_free(g->v[k]);
    free(g->rp); free(g->ci); free(g->v);
    memset(g, 0, sizeof *g);
}

static Expr* dense_elem(const Opd* o, int64_t k) {
    if (o->nd) return ndarray_buffer_element_to_expr(o->nd->data.ndarray.data, (size_t)k,
                                                     o->nd->data.ndarray.dtype);
    return expr_copy(o->leaves[k]);
}

static bool gcsr_build(GCsr* g, const Opd* o, bool left) {
    memset(g, 0, sizeof *g);
    bool col = !left && o->rank == 1;          /* a right-hand vector: one column */
    g->nr = o->rank == 2 ? o->dims[0] : (left ? 1 : o->dims[0]);
    g->nc = o->rank == 2 ? o->dims[1] : (left ? o->dims[0] : 1);
    g->nnz = opd_count(o);
    g->rp = calloc((size_t)g->nr + 1, sizeof(int64_t));
    g->ci = malloc(sizeof(int64_t) * (size_t)(g->nnz ? g->nnz : 1));
    g->v = calloc((size_t)(g->nnz ? g->nnz : 1), sizeof(Expr*));
    if (!g->rp || !g->ci || !g->v) { gcsr_free(g); return false; }
    if (o->sparse) {
        const SparsePattern* p = o->sa->pat;
        for (int64_t k = 0; k < p->nnz; k++) g->v[k] = sp_value_expr(o->sa, k);
        if (col) {
            /* vector entries k -> (row k, column 0) */
            for (int64_t q = 0; q < p->nnz; q++) {
                int64_t k = sp_col(p, q, 0);
                g->rp[k + 1] = 1;
                g->ci[q] = 0;
            }
            for (int64_t r = 0; r < g->nr; r++) g->rp[r + 1] += g->rp[r];
        } else {
            for (int64_t r = 0; r <= g->nr; r++) g->rp[r] = p->rowptr[r];
            for (int64_t q = 0; q < p->nnz; q++) g->ci[q] = sp_col(p, q, 0);
        }
    } else {
        for (int64_t k = 0; k < g->nnz; k++) g->v[k] = dense_elem(o, k);
        for (int64_t r = 0; r <= g->nr; r++) g->rp[r] = r * g->nc;
        for (int64_t k = 0; k < g->nnz; k++) g->ci[k] = k % g->nc;
    }
    for (int64_t k = 0; k < g->nnz; k++) if (!g->v[k]) { gcsr_free(g); return false; }
    return true;
}

static int cmp_i64(const void* x, const void* y) {
    int64_t a = *(const int64_t*)x, b = *(const int64_t*)y;
    return a < b ? -1 : (a > b ? 1 : 0);
}

/* Sum of the terms Times[a, b] as one evaluated Plus (adopts the terms). */
static Expr* boxed_sum(Expr** terms, size_t n) {
    Expr* e = n == 1 ? terms[0] : expr_new_function(expr_new_symbol(SYM_Plus), terms, n);
    Expr* r = evaluate(e);
    expr_free(e);
    return r;
}

/* a . b over exact or symbolic values: a structural product, each output
 * element evaluated once as Plus of its Times terms. */
static Expr* boxed_product(const Opd* a, const Opd* b) {
    GCsr L, R;
    if (!gcsr_build(&L, a, true)) return NULL;
    if (!gcsr_build(&R, b, false)) { gcsr_free(&L); return NULL; }
    int64_t n = L.nr, p = R.nc;
    bool sparse_out = a->sparse && b->sparse;
    int64_t dims[2];
    int rank = 0;
    if (a->rank == 2) dims[rank++] = a->dims[0];
    if (b->rank == 2) dims[rank++] = b->dims[1];
    if (rank == 2 && n > SPD_MAX_DENSE / p && !sparse_out) { gcsr_free(&L); gcsr_free(&R); return NULL; }

    /* per output column: a chain of terms (head/next over a growing pool) */
    int64_t* head = malloc(sizeof(int64_t) * (size_t)p);
    int64_t* mark = malloc(sizeof(int64_t) * (size_t)p);
    int64_t* touched = malloc(sizeof(int64_t) * (size_t)p);
    int64_t* oc = NULL;
    int64_t* orow = NULL;
    Expr** ov = NULL;
    int64_t ocnt = 0, ocap = 0;
    Expr** pool = NULL;
    int64_t* next = NULL;
    int64_t pcap = 0;
    Expr** terms = NULL;
    int64_t tcap = 0;
    bool ok = head && mark && touched;
    for (int64_t j = 0; ok && j < p; j++) mark[j] = -1;
    for (int64_t i = 0; ok && i < n; i++) {
        int64_t nt = 0, npool = 0;
        for (int64_t ka = L.rp[i]; ok && ka < L.rp[i + 1]; ka++) {
            int64_t k = L.ci[ka];
            for (int64_t kb = R.rp[k]; kb < R.rp[k + 1]; kb++) {
                int64_t j = R.ci[kb];
                if (npool == pcap) {
                    pcap = pcap ? 2 * pcap : 64;
                    Expr** np = realloc(pool, sizeof(Expr*) * (size_t)pcap);
                    int64_t* nn = np ? realloc(next, sizeof(int64_t) * (size_t)pcap) : NULL;
                    if (np) pool = np;
                    if (nn) next = nn;
                    if (!np || !nn) { ok = false; break; }
                }
                Expr* targs[2] = { expr_copy(L.v[ka]), expr_copy(R.v[kb]) };
                pool[npool] = expr_new_function(expr_new_symbol(SYM_Times), targs, 2);
                if (mark[j] != i) { mark[j] = i; touched[nt++] = j; head[j] = -1; }
                next[npool] = head[j];
                head[j] = npool;
                npool++;
            }
        }
        if (!ok) { for (int64_t q = 0; q < npool; q++) expr_free(pool[q]); break; }
        /* output columns in order */
        int64_t* tt = touched;
        if (nt > 32) qsort(tt, (size_t)nt, sizeof(int64_t), cmp_i64);
        else for (int64_t x = 1; x < nt; x++) {
            int64_t v = tt[x], y = x - 1;
            while (y >= 0 && tt[y] > v) { tt[y + 1] = tt[y]; y--; }
            tt[y + 1] = v;
        }
        for (int64_t x = 0; x < nt; x++) {
            int64_t j = tt[x];
            int64_t m = 0;
            for (int64_t q = head[j]; q >= 0; q = next[q]) m++;
            if (m > tcap) {
                Expr** nt2 = realloc(terms, sizeof(Expr*) * (size_t)m);
                if (!nt2) { ok = false; break; }
                terms = nt2; tcap = m;
            }
            /* the chain is newest-first: restore storage order */
            int64_t w = m;
            for (int64_t q = head[j]; q >= 0; q = next[q]) terms[--w] = pool[q];
            for (int64_t q = head[j]; q >= 0; q = next[q]) pool[q] = NULL;
            Expr* val = boxed_sum(terms, (size_t)m);
            if (ocnt == ocap) {
                ocap = ocap ? 2 * ocap : 64;
                int64_t* a1 = realloc(oc, sizeof(int64_t) * (size_t)ocap);
                int64_t* a2 = a1 ? realloc(orow, sizeof(int64_t) * (size_t)ocap) : NULL;
                Expr** a3 = a2 ? realloc(ov, sizeof(Expr*) * (size_t)ocap) : NULL;
                if (a1) oc = a1;
                if (a2) orow = a2;
                if (a3) ov = a3;
                if (!a1 || !a2 || !a3) { expr_free(val); ok = false; break; }
            }
            oc[ocnt] = j; orow[ocnt] = i; ov[ocnt] = val; ocnt++;
        }
        for (int64_t q = 0; q < npool; q++) if (pool[q]) expr_free(pool[q]);
    }
    free(head); free(mark); free(touched); free(pool); free(next); free(terms);
    gcsr_free(&L);
    gcsr_free(&R);
    Expr* out = NULL;
    if (ok) {
        if (rank == 0) {
            out = ocnt ? expr_copy(ov[0]) : expr_new_integer(0);
        } else if (sparse_out) {
            /* rank 2: rows i, columns j; a vector result is one row over the
             * free index (j for sv . s, i for s . sv) */
            bool vec_rows = rank == 1 && a->rank == 2;
            SparsePattern* pat = sp_pattern_new(rank == 2 ? n : 1, ocnt, 1, true);
            Expr** vals = malloc(sizeof(Expr*) * (size_t)(ocnt ? ocnt : 1));
            if (pat && vals) {
                for (int64_t r = 0; r <= pat->nrows; r++) pat->rowptr[r] = 0;
                for (int64_t q = 0; q < ocnt; q++) {
                    ((int32_t*)pat->colidx)[q] = (int32_t)(vec_rows ? orow[q] : oc[q]);
                    vals[q] = expr_copy(ov[q]);
                    if (rank == 2) pat->rowptr[orow[q] + 1]++;
                }
                if (rank == 2) for (int64_t r = 0; r < n; r++) pat->rowptr[r + 1] += pat->rowptr[r];
                else pat->rowptr[1] = ocnt;
                SparseArray* sa = sp_new(rank, dims, expr_new_integer(0), pat, SPV_EXPR, vals);
                out = sa ? sp_to_expr(sa) : NULL;
            } else {
                sp_pattern_release(pat);
                free(vals);
            }
        } else {
            int64_t total = n * p;
            Expr** grid = malloc(sizeof(Expr*) * (size_t)total);
            if (grid) {
                for (int64_t k = 0; k < total; k++) grid[k] = NULL;
                for (int64_t q = 0; q < ocnt; q++) grid[orow[q] * p + oc[q]] = expr_copy(ov[q]);
                for (int64_t k = 0; k < total; k++) if (!grid[k]) grid[k] = expr_new_integer(0);
                if (rank == 1) {
                    out = expr_new_function(expr_new_symbol(SYM_List), grid, (size_t)total);
                } else {
                    Expr** rows = malloc(sizeof(Expr*) * (size_t)n);
                    if (rows) {
                        for (int64_t i = 0; i < n; i++)
                            rows[i] = expr_new_function(expr_new_symbol(SYM_List), grid + i * p, (size_t)p);
                        out = expr_new_function(expr_new_symbol(SYM_List), rows, (size_t)n);
                        free(rows);
                    } else {
                        for (int64_t k = 0; k < total; k++) expr_free(grid[k]);
                    }
                }
                free(grid);
            }
        }
    }
    for (int64_t q = 0; q < ocnt; q++) expr_free(ov[q]);
    free(oc); free(orow); free(ov);
    return out;
}

/* ======================================================================== */
/* Dot                                                                      */
/* ======================================================================== */

enum { DOT_OK, DOT_DECLINE, DOT_SHAPE };

/* a . b for one pair with at least one sparse operand. */
static Expr* dot_pair(const Expr* ea, const Expr* eb, int* status) {
    Opd a, b;
    *status = DOT_DECLINE;
    if (!opd_init(&a, ea)) return NULL;
    if (!opd_init(&b, eb)) { opd_free(&a); return NULL; }
    Expr* out = NULL;
    if (!a.sparse && !b.sparse) goto done;
    if (a.dims[a.rank - 1] != b.dims[0]) { *status = DOT_SHAPE; goto done; }
    if ((a.sparse && a.sa->pat->nnz == 0) || (b.sparse && b.sa->pat->nnz == 0)) {
        out = empty_product(&a, &b);
    } else {
        int ty = pair_ty(&a, &b);
        if (ty != TY_BOX) out = typed_product(&a, &b, ty);
        if (!out) out = boxed_product(&a, &b);   /* exact values, or int64 overflow */
    }
    if (out) *status = DOT_OK;
done:
    opd_free(&a);
    opd_free(&b);
    return out;
}

static void dotsh(const Expr* a, const Expr* b) {
    char* as = expr_to_string((Expr*)a);
    char* bs = expr_to_string((Expr*)b);
    mth_message("Dot", "dotsh", "Tensors %s and %s have incompatible shapes.",
                as ? as : "", bs ? bs : "");
    free(as);
    free(bs);
}

static bool any_sparse(const Expr* x, const Expr* y) {
    return x->type == EXPR_SPARSE || y->type == EXPR_SPARSE;
}

/* One step of a chain: a native pair when a SparseArray is involved, the
 * ordinary Dot otherwise. */
static Expr* chain_step(const Expr* x, const Expr* y) {
    if (any_sparse(x, y)) {
        int st;
        return dot_pair(x, y, &st);
    }
    Expr* args[2] = { expr_copy((Expr*)x), expr_copy((Expr*)y) };
    Expr* call = expr_new_function(expr_new_symbol(SYM_Dot), args, 2);
    Expr* r = evaluate(call);
    expr_free(call);
    if (r->type == EXPR_FUNCTION && r->data.function.head->type == EXPR_SYMBOL &&
        r->data.function.head->data.symbol.name == SYM_Dot) {
        expr_free(r);
        return NULL;
    }
    return r;
}

/* The rank and dims of a Dot operand (rank 0: not an array these take). */
static int arr_shape(const Expr* e, int64_t* d0, int64_t* dl) {
    Opd o;
    if (!opd_init(&o, e)) return 0;
    *d0 = o.dims[0];
    *dl = o.dims[o.rank - 1];
    int r = o.rank;
    bool empty = o.sparse && o.sa->pat->nnz == 0;
    opd_free(&o);
    return empty ? -r : r;
}

static Expr* sparse_dot_native(const Expr* call) {
    size_t n = call->data.function.arg_count;
    Expr* const* args = call->data.function.args;
    if (n == 1) return expr_copy(args[0]);
    if (n == 2) {
        int st;
        Expr* r = dot_pair(args[0], args[1], &st);
        if (r) return r;
        if (st == DOT_SHAPE) {
            /* reported here, with the arrays as written; the call stays */
            dotsh(args[0], args[1]);
            return expr_copy((Expr*)call);
        }
        return NULL;
    }
    /* a chain: every operand a vector or matrix, the inner ones matrices
     * (a vector inside the chain would make a scalar), adjacent shapes
     * agreeing; otherwise the gate answers */
    bool right_to_left = true;
    int64_t prev_last = -1;
    for (size_t i = 0; i < n; i++) {
        int64_t d0, dl;
        int r = arr_shape(args[i], &d0, &dl);
        if (r == 0) return NULL;
        if (r < 0) { right_to_left = false; r = -r; }   /* an empty factor: keep the order */
        if (i > 0 && i + 1 < n && r != 2) return NULL;
        if (i > 0 && prev_last != d0) return NULL;
        prev_last = dl;
        if (i + 1 == n && (r != 1 || args[i]->type == EXPR_SPARSE)) right_to_left = false;
    }
    Expr* acc;
    if (right_to_left) {
        acc = expr_copy(args[n - 1]);
        for (size_t i = n - 1; i-- > 0;) {
            Expr* r = chain_step(args[i], acc);
            expr_free(acc);
            if (!r) return NULL;
            acc = r;
        }
    } else {
        acc = expr_copy(args[0]);
        for (size_t i = 1; i < n; i++) {
            Expr* r = chain_step(acc, args[i]);
            expr_free(acc);
            if (!r) return NULL;
            acc = r;
        }
    }
    return acc;
}

/* ======================================================================== */
/* MatrixPower[s, k] and MatrixPower[s, k, v]                               */
/* ======================================================================== */

/* MatrixPower on the dense form, for the powers whose answer is dense: a
 * nonzero background, a negative, non-integer or symbolic exponent. A call the
 * dense routine leaves alone stays as written (its message already given). */
static Expr* matpow_dense(const Expr* call, const SparseArray* sa) {
    int64_t total = sp_total(2, sa->dims);
    if (total < 0 || total > ((int64_t)1 << 24)) return NULL;
    Expr* d = sp_normal(sa);
    if (!d) return NULL;
    if (d->type == EXPR_NDARRAY && call->data.function.args[1]->type != EXPR_INTEGER) {
        /* a non-integer power is declined with a message; a packed operand
         * would be materialised and the call asked (and the message given)
         * a second time */
        Expr* u = pack_unpack(d);
        expr_free(d);
        d = u;
        if (!d) return NULL;
    }
    size_t n = call->data.function.arg_count;
    Expr* cargs[3];
    cargs[0] = d;
    for (size_t i = 1; i < n; i++) {
        const Expr* x = call->data.function.args[i];
        cargs[i] = x->type == EXPR_SPARSE ? sp_normal(sp_of(x)) : expr_copy((Expr*)x);
        if (!cargs[i]) {
            for (size_t j = 0; j < i; j++) expr_free(cargs[j]);
            return NULL;
        }
    }
    Expr* c = expr_new_function(expr_copy(call->data.function.head), cargs, n);
    Expr* r = evaluate(c);
    expr_free(c);
    if (r->type == EXPR_FUNCTION && expr_eq(r->data.function.head, call->data.function.head)) {
        expr_free(r);
        return expr_copy((Expr*)call);
    }
    return r;
}

static Expr* sparse_matrixpower_native(const Expr* call) {
    size_t n = call->data.function.arg_count;
    if (n < 2 || n > 3) return NULL;
    const Expr* s = call->data.function.args[0];
    const Expr* k = call->data.function.args[1];
    if (s->type != EXPR_SPARSE) return NULL;
    const SparseArray* sa = sp_of(s);
    if (sa->rank != 2 || sa->dims[0] != sa->dims[1] || sa->dims[0] <= 0) {
        char* shown = expr_to_string((Expr*)s);
        mth_message("MatrixPower", "matsq",
                    "Argument %s at position 1 is not a non-empty square matrix.",
                    shown ? shown : "");
        free(shown);
        return expr_copy((Expr*)call);
    }
    if (k->type != EXPR_INTEGER || k->data.integer < 0 || !sp_is_zero(sa->bg))
        return matpow_dense(call, sa);
    int64_t e = k->data.integer;
    if (n == 3) {
        /* s . (s . ( ... v)): k sparse-matrix-vector products; the answer is
         * dense, even for a sparse v */
        const Expr* v = call->data.function.args[2];
        Expr* acc = v->type == EXPR_SPARSE ? sp_normal(sp_of(v)) : expr_copy((Expr*)v);
        if (!acc) return NULL;
        for (int64_t i = 0; i < e; i++) {
            int st;
            Expr* r = dot_pair(s, acc, &st);
            if (!r && st == DOT_SHAPE) {
                mth_message("MatrixPower", "lslc",
                            "Coefficient matrix and target vector(s) or matrix do not have the same dimensions.");
                expr_free(acc);
                return expr_copy((Expr*)call);
            }
            expr_free(acc);
            if (!r) return NULL;
            acc = r;
        }
        return acc;
    }
    if (e == 0) return NULL;            /* the identity: the gate builds it */
    if (e == 1) return expr_copy((Expr*)s);
    /* binary powering: s^e from repeated squares */
    Expr* result = NULL;
    Expr* base = expr_copy((Expr*)s);
    bool ok = true;
    for (;;) {
        int st;
        if (e & 1) {
            if (!result) result = expr_copy(base);
            else {
                Expr* r = dot_pair(result, base, &st);
                expr_free(result);
                result = r;
                ok = r != NULL;
            }
        }
        e >>= 1;
        if (!ok || e == 0) break;
        Expr* sq = dot_pair(base, base, &st);
        expr_free(base);
        base = sq;
        if (!sq) { ok = false; break; }
    }
    if (base) expr_free(base);
    if (!ok) { if (result) expr_free(result); return NULL; }
    return result;
}

/* ======================================================================== */
/* Inner[Times, a, b, Plus]                                                 */
/* ======================================================================== */

static bool is_sym(const Expr* e, const char* name) {
    return e->type == EXPR_SYMBOL && e->data.symbol.name == name;
}

/* The generalised inner product with the ordinary operations is Dot -- and
 * its answer is always dense. */
static Expr* sparse_inner_native(const Expr* call) {
    if (call->data.function.arg_count != 4) return NULL;
    Expr* const* a = call->data.function.args;
    if (!is_sym(a[0], SYM_Times) || !is_sym(a[3], SYM_Plus)) return NULL;
    Opd x, y;
    if (!opd_init(&x, a[1])) return NULL;
    if (!opd_init(&y, a[2])) { opd_free(&x); return NULL; }
    bool empty = (x.sparse && x.sa->pat->nnz == 0) || (y.sparse && y.sa->pat->nnz == 0);
    bool fit = x.dims[x.rank - 1] == y.dims[0];
    opd_free(&x);
    opd_free(&y);
    if (empty || !fit) return NULL;
    int st;
    Expr* r = dot_pair(a[1], a[2], &st);
    if (r && r->type == EXPR_SPARSE) {
        Expr* d = sp_normal(sp_of(r));
        expr_free(r);
        r = d;
    }
    return r;
}

void sparse_dot_init(void) {
    symtab_set_sparse_native("Dot", (const void*)sparse_dot_native);
    symtab_set_sparse_native("MatrixPower", (const void*)sparse_matrixpower_native);
    symtab_set_sparse_native("Inner", (const void*)sparse_inner_native);
}
