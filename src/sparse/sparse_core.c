/*
 * sparse_core.c -- the SparseArray payload: lifecycle, canonical value kinds,
 * element access, identity (SameQ / hash / order), the canonical
 * SparseArray[Automatic, ...] form, Normal, and the COO -> CSR builder every
 * constructor and operation funnels through.
 *
 * See sparse.h for the representation and the invariants it guarantees.
 */

#include "sparse.h"
#include "sym_names.h"
#include "arithmetic.h"
#include "ndarray.h"
#include "pack.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

/* Largest dense array Normal builds: boxed (an Expr* per element) and typed (a
 * machine buffer). Beyond these Normal is left unevaluated rather than
 * attempting a multi-gigabyte allocation. */
#define SP_NORMAL_MAX_BOXED ((int64_t)1 << 27)
#define SP_NORMAL_MAX_TYPED ((int64_t)1 << 30)

/* ======================================================================== */
/* Patterns                                                                 */
/* ======================================================================== */

SparsePattern* sp_pattern_new(int64_t nrows, int64_t nnz, int width, bool idx32) {
    SparsePattern* p = malloc(sizeof *p);
    if (!p) return NULL;
    p->refcount = 1;
    p->nrows = nrows;
    p->nnz = nnz;
    p->width = width;
    p->idx32 = idx32;
    p->rowptr = malloc(sizeof(int64_t) * (size_t)(nrows + 1));
    size_t isz = idx32 ? sizeof(int32_t) : sizeof(int64_t);
    p->colidx = nnz > 0 ? malloc(isz * (size_t)nnz * (size_t)width) : NULL;
    if (!p->rowptr || (nnz > 0 && !p->colidx)) {
        free(p->rowptr);
        free(p->colidx);
        free(p);
        return NULL;
    }
    return p;
}

SparsePattern* sp_pattern_retain(SparsePattern* p) {
    if (p) p->refcount++;
    return p;
}

void sp_pattern_release(SparsePattern* p) {
    if (!p) return;
    if (--p->refcount > 0) return;
    free(p->rowptr);
    free(p->colidx);
    free(p);
}

bool sp_pattern_same(const SparsePattern* a, const SparsePattern* b) {
    if (a == b) return true;
    if (a->nrows != b->nrows || a->nnz != b->nnz || a->width != b->width ||
        a->idx32 != b->idx32) return false;
    if (memcmp(a->rowptr, b->rowptr, sizeof(int64_t) * (size_t)(a->nrows + 1)) != 0)
        return false;
    size_t isz = a->idx32 ? sizeof(int32_t) : sizeof(int64_t);
    return a->nnz == 0 ||
           memcmp(a->colidx, b->colidx, isz * (size_t)a->nnz * (size_t)a->width) == 0;
}

bool sp_dims_idx32(int rank, const int64_t* dims) {
    int first = rank >= 2 ? 1 : 0;
    for (int k = first; k < rank; k++)
        if (dims[k] > INT32_MAX) return false;
    return true;
}

/* ======================================================================== */
/* Values                                                                   */
/* ======================================================================== */

size_t sp_kind_size(SpValKind kind) {
    switch (kind) {
        case SPV_INT64: return sizeof(int64_t);
        case SPV_F64:   return sizeof(double);
        case SPV_C128:  return 2 * sizeof(double);
        default:        return 0;
    }
}

/* Is `e` a machine complex number Complex[Real, Real]? */
static bool is_machine_complex(const Expr* e, double* re, double* im) {
    if (!e || e->type != EXPR_FUNCTION || e->data.function.arg_count != 2 ||
        e->data.function.head->type != EXPR_SYMBOL ||
        e->data.function.head->data.symbol.name != SYM_Complex) return false;
    const Expr* a = e->data.function.args[0];
    const Expr* b = e->data.function.args[1];
    if (a->type != EXPR_REAL || b->type != EXPR_REAL) return false;
    if (re) *re = a->data.real;
    if (im) *im = b->data.real;
    return true;
}

static Expr* box_c128(double re, double im) {
    Expr* args[2] = { expr_new_real(re), expr_new_real(im) };
    return expr_new_function(expr_new_symbol(SYM_Complex), args, 2);
}

Expr* sp_value_expr(const SparseArray* sa, int64_t k) {
    switch (sa->kind) {
        case SPV_INT64: return expr_new_integer(((const int64_t*)sa->vals)[k]);
        case SPV_F64:   return expr_new_real(((const double*)sa->vals)[k]);
        case SPV_C128: {
            const double* v = (const double*)sa->vals;
            return box_c128(v[2 * k], v[2 * k + 1]);
        }
        default:        return expr_copy(((Expr* const*)sa->vals)[k]);
    }
}

/* Convert a boxed value buffer to the canonical typed kind when every value
 * shares one machine type. Frees the boxed buffer on conversion. */
static void normalise_kind(SpValKind* kind, void** vals, int64_t nnz) {
    if (nnz == 0) {
        /* The empty array has one canonical kind (and no values to free). */
        free(*vals);
        *vals = NULL;
        *kind = SPV_INT64;
        return;
    }
    if (*kind != SPV_EXPR) return;
    Expr** ev = (Expr**)*vals;
    ExprType t0 = ev[0]->type;
    bool all_int = t0 == EXPR_INTEGER, all_real = t0 == EXPR_REAL;
    bool all_cplx = !all_int && !all_real && is_machine_complex(ev[0], NULL, NULL);
    if (!all_int && !all_real && !all_cplx) return;
    for (int64_t k = 1; k < nnz; k++) {
        const Expr* v = ev[k];
        if (all_int && v->type != EXPR_INTEGER) return;
        if (all_real && v->type != EXPR_REAL) return;
        if (all_cplx && !is_machine_complex(v, NULL, NULL)) return;
    }
    void* tv;
    if (all_int) {
        int64_t* b = malloc(sizeof(int64_t) * (size_t)nnz);
        if (!b) return;
        for (int64_t k = 0; k < nnz; k++) b[k] = ev[k]->data.integer;
        tv = b; *kind = SPV_INT64;
    } else if (all_real) {
        double* b = malloc(sizeof(double) * (size_t)nnz);
        if (!b) return;
        for (int64_t k = 0; k < nnz; k++) b[k] = ev[k]->data.real;
        tv = b; *kind = SPV_F64;
    } else {
        double* b = malloc(2 * sizeof(double) * (size_t)nnz);
        if (!b) return;
        for (int64_t k = 0; k < nnz; k++) is_machine_complex(ev[k], &b[2 * k], &b[2 * k + 1]);
        tv = b; *kind = SPV_C128;
    }
    for (int64_t k = 0; k < nnz; k++) expr_free(ev[k]);
    free(ev);
    *vals = tv;
}

static void free_vals(SpValKind kind, void* vals, int64_t nnz) {
    if (kind == SPV_EXPR && vals) {
        Expr** ev = (Expr**)vals;
        for (int64_t k = 0; k < nnz; k++) if (ev[k]) expr_free(ev[k]);
    }
    free(vals);
}

/* ======================================================================== */
/* Arrays                                                                   */
/* ======================================================================== */

bool sp_g_any_created = false;

SparseArray* sp_new(int rank, const int64_t* dims, Expr* bg, SparsePattern* pat,
                    SpValKind kind, void* vals) {
    sp_g_any_created = true;
    SparseArray* sa = malloc(sizeof *sa);
    int64_t* d = malloc(sizeof(int64_t) * (size_t)(rank > 0 ? rank : 1));
    if (!sa || !d) {
        free(sa); free(d);
        free_vals(kind, vals, pat ? pat->nnz : 0);
        sp_pattern_release(pat);
        if (bg) expr_free(bg);
        return NULL;
    }
    memcpy(d, dims, sizeof(int64_t) * (size_t)rank);
    normalise_kind(&kind, &vals, pat->nnz);
    sa->refcount = 1;
    sa->rank = rank;
    sa->dims = d;
    sa->bg = bg;
    sa->pat = pat;
    sa->kind = kind;
    sa->vals = vals;
    return sa;
}

SparseArray* sp_retain(SparseArray* sa) {
    if (sa) sa->refcount++;
    return sa;
}

void sp_release(SparseArray* sa) {
    if (!sa) return;
    if (--sa->refcount > 0) return;
    free_vals(sa->kind, sa->vals, sa->pat->nnz);
    sp_pattern_release(sa->pat);
    if (sa->bg) expr_free(sa->bg);
    free(sa->dims);
    free(sa);
}

Expr* sp_to_expr(SparseArray* sa) {
    if (!sa) return NULL;
    Expr* e = expr_new_sparse(sa);
    if (!e) sp_release(sa);
    return e;
}

int64_t sp_total(int rank, const int64_t* dims) {
    int64_t t = 1;
    for (int k = 0; k < rank; k++) {
        if (dims[k] < 0) return -1;
        if (dims[k] == 0) return 0;
        if (t > INT64_MAX / dims[k]) return -1;
        t *= dims[k];
    }
    return t;
}

/* Lexicographic comparison of entry k's tuple with the 0-based tuple `c`. */
static int tuple_cmp(const SparsePattern* p, int64_t k, const int64_t* c) {
    for (int j = 0; j < p->width; j++) {
        int64_t v = sp_col(p, k, j);
        if (v != c[j]) return v < c[j] ? -1 : 1;
    }
    return 0;
}

/* Index of the stored entry at 0-based coordinates `c`, or -1. */
static int64_t sp_find(const SparseArray* sa, const int64_t* c) {
    const SparsePattern* p = sa->pat;
    int64_t row = sa->rank >= 2 ? c[0] : 0;
    const int64_t* tup = sa->rank >= 2 ? c + 1 : c;
    int64_t lo = p->rowptr[row], hi = p->rowptr[row + 1] - 1;
    while (lo <= hi) {
        int64_t mid = lo + (hi - lo) / 2;
        int r = tuple_cmp(p, mid, tup);
        if (r == 0) return mid;
        if (r < 0) lo = mid + 1; else hi = mid - 1;
    }
    return -1;
}

Expr* sp_get(const SparseArray* sa, const int64_t* c) {
    int64_t k = sp_find(sa, c);
    return k >= 0 ? sp_value_expr(sa, k) : expr_copy(sa->bg);
}

/* ======================================================================== */
/* Identity                                                                 */
/* ======================================================================== */

static bool real_same(double a, double b) {
    /* expr_eq's rule for two Reals: equal values, and NaN matches NaN. */
    return a == b || (isnan(a) && isnan(b));
}

bool sp_eq(const SparseArray* a, const SparseArray* b) {
    if (a == b) return true;
    if (a->rank != b->rank) return false;
    for (int k = 0; k < a->rank; k++) if (a->dims[k] != b->dims[k]) return false;
    if (!sp_pattern_same(a->pat, b->pat)) return false;
    if (a->kind != b->kind) return false;   /* the kind is canonical */
    if (!expr_eq(a->bg, b->bg)) return false;
    int64_t n = a->pat->nnz;
    switch (a->kind) {
        case SPV_INT64:
            return n == 0 || memcmp(a->vals, b->vals, sizeof(int64_t) * (size_t)n) == 0;
        case SPV_F64: {
            const double* x = a->vals; const double* y = b->vals;
            for (int64_t k = 0; k < n; k++) if (!real_same(x[k], y[k])) return false;
            return true;
        }
        case SPV_C128: {
            const double* x = a->vals; const double* y = b->vals;
            for (int64_t k = 0; k < 2 * n; k++) if (!real_same(x[k], y[k])) return false;
            return true;
        }
        default: {
            Expr* const* x = a->vals; Expr* const* y = b->vals;
            for (int64_t k = 0; k < n; k++) if (!expr_eq(x[k], y[k])) return false;
            return true;
        }
    }
}

#define SP_FNV_PRIME 1099511628211ULL
#define SP_MIX(h, v) do { (h) ^= (uint64_t)(v); (h) *= SP_FNV_PRIME; } while (0)

static uint64_t double_bits(double d) {
    uint64_t u;
    if (d == 0.0) d = 0.0;            /* -0.0 == 0.0, so they must hash alike */
    if (isnan(d)) d = NAN;            /* every NaN is SameQ to every other */
    memcpy(&u, &d, sizeof u);
    return u;
}

uint64_t sp_hash(const SparseArray* sa) {
    uint64_t h = 14695981039346656037ULL;
    SP_MIX(h, 0x5A5A0001u);           /* type tag: never collides with a List */
    SP_MIX(h, sa->rank);
    for (int k = 0; k < sa->rank; k++) SP_MIX(h, sa->dims[k]);
    SP_MIX(h, expr_hash(sa->bg));
    const SparsePattern* p = sa->pat;
    SP_MIX(h, p->nnz);
    for (int64_t r = 0; r <= p->nrows; r++) SP_MIX(h, p->rowptr[r]);
    for (int64_t k = 0; k < p->nnz; k++)
        for (int j = 0; j < p->width; j++) SP_MIX(h, sp_col(p, k, j));
    switch (sa->kind) {
        case SPV_INT64:
            for (int64_t k = 0; k < p->nnz; k++) SP_MIX(h, ((const int64_t*)sa->vals)[k]);
            break;
        case SPV_F64:
            for (int64_t k = 0; k < p->nnz; k++) SP_MIX(h, double_bits(((const double*)sa->vals)[k]));
            break;
        case SPV_C128:
            for (int64_t k = 0; k < 2 * p->nnz; k++) SP_MIX(h, double_bits(((const double*)sa->vals)[k]));
            break;
        default:
            for (int64_t k = 0; k < p->nnz; k++) SP_MIX(h, expr_hash(((Expr* const*)sa->vals)[k]));
            break;
    }
    return h;
}

int sp_compare(const SparseArray* a, const SparseArray* b) {
    if (a == b) return 0;
    if (a->rank != b->rank) return a->rank < b->rank ? -1 : 1;
    for (int k = 0; k < a->rank; k++)
        if (a->dims[k] != b->dims[k]) return a->dims[k] < b->dims[k] ? -1 : 1;
    int64_t na = a->pat->nnz, nb = b->pat->nnz;
    if (na != nb) return na < nb ? -1 : 1;
    int c = expr_compare(a->bg, b->bg);
    if (c) return c;
    /* Same size: positions first, then values, in canonical entry order. */
    const SparsePattern* pa = a->pat; const SparsePattern* pb = b->pat;
    for (int64_t r = 0; r <= pa->nrows; r++)
        if (pa->rowptr[r] != pb->rowptr[r]) return pa->rowptr[r] > pb->rowptr[r] ? -1 : 1;
    for (int64_t k = 0; k < na; k++)
        for (int j = 0; j < pa->width; j++) {
            int64_t x = sp_col(pa, k, j), y = sp_col(pb, k, j);
            if (x != y) return x < y ? -1 : 1;
        }
    for (int64_t k = 0; k < na; k++) {
        Expr* x = sp_value_expr(a, k);
        Expr* y = sp_value_expr(b, k);
        c = expr_compare(x, y);
        expr_free(x); expr_free(y);
        if (c) return c;
    }
    return 0;
}

/* ======================================================================== */
/* Canonical form and Normal                                                */
/* ======================================================================== */

static Expr* mk_list_adopt(Expr** items, size_t n) {
    Expr* out = expr_new_function(expr_new_symbol(SYM_List), items, n);
    free(items);
    return out;
}

static Expr** alloc_items(size_t n) {
    return malloc(sizeof(Expr*) * (n ? n : 1));
}

static Expr* dims_list(int rank, const int64_t* dims) {
    Expr** it = alloc_items((size_t)rank);
    for (int k = 0; k < rank; k++) it[k] = expr_new_integer(dims[k]);
    return mk_list_adopt(it, (size_t)rank);
}

Expr* sp_canonical_expr(const SparseArray* sa) {
    const SparsePattern* p = sa->pat;
    int64_t n = p->nnz;

    Expr** rp = alloc_items((size_t)(p->nrows + 1));
    for (int64_t r = 0; r <= p->nrows; r++) rp[r] = expr_new_integer(p->rowptr[r]);
    Expr* rowptr = mk_list_adopt(rp, (size_t)(p->nrows + 1));

    Expr** ci = alloc_items((size_t)n);
    for (int64_t k = 0; k < n; k++) {
        Expr** t = alloc_items((size_t)p->width);
        for (int j = 0; j < p->width; j++) t[j] = expr_new_integer(sp_col(p, k, j) + 1);
        ci[k] = mk_list_adopt(t, (size_t)p->width);
    }
    Expr* colidx = mk_list_adopt(ci, (size_t)n);

    Expr** vs = alloc_items((size_t)n);
    for (int64_t k = 0; k < n; k++) vs[k] = sp_value_expr(sa, k);
    Expr* vals = mk_list_adopt(vs, (size_t)n);

    Expr* idx_args[2] = { rowptr, colidx };
    Expr* idx = expr_new_function(expr_new_symbol(SYM_List), idx_args, 2);
    Expr* data_args[3] = { expr_new_integer(1), idx, vals };
    Expr* data = expr_new_function(expr_new_symbol(SYM_List), data_args, 3);
    Expr* args[4] = { expr_new_symbol(SYM_Automatic), dims_list(sa->rank, sa->dims),
                      expr_copy(sa->bg), data };
    return expr_new_function(expr_new_symbol(SYM_SparseArray), args, 4);
}

/* Row-major 0-based offset of entry k (row r) in the dense array. */
static int64_t entry_offset(const SparseArray* sa, int64_t r, int64_t k) {
    const SparsePattern* p = sa->pat;
    int64_t off = sa->rank >= 2 ? r : 0;
    int first = sa->rank >= 2 ? 1 : 0;
    for (int j = 0; j < p->width; j++)
        off = off * sa->dims[first + j] + sp_col(p, k, j);
    return off;
}

/* Fold a flat row-major slot buffer into nested Lists, moving the slots. */
static Expr* fold_slots(Expr** slots, int rank, const int64_t* dims, int level, int64_t* off) {
    int64_t n = dims[level];
    Expr** items = alloc_items((size_t)n);
    for (int64_t i = 0; i < n; i++) {
        if (level + 1 == rank) items[i] = slots[(*off)++];
        else items[i] = fold_slots(slots, rank, dims, level + 1, off);
    }
    return mk_list_adopt(items, (size_t)n);
}

/* Try to build Normal[s] straight into a packed buffer. NULL when the kinds do
 * not share one machine type, or packing declines. */
static Expr* normal_packed(const SparseArray* sa, int64_t total) {
    const Expr* bg = sa->bg;
    double bre = 0, bim = 0;
    NDType dt;
    bool empty = sa->pat->nnz == 0;
    if (bg->type == EXPR_REAL && (empty || sa->kind == SPV_F64)) {
        dt = NDT_FLOAT64; bre = bg->data.real;
    } else if (bg->type == EXPR_INTEGER && (empty || sa->kind == SPV_INT64)) {
        dt = NDT_INT64;
    } else if (is_machine_complex(bg, &bre, &bim) && (empty || sa->kind == SPV_C128)) {
        dt = NDT_COMPLEX64;
    } else {
        return NULL;
    }
    if (total > SP_NORMAL_MAX_TYPED) return NULL;
    /* A zero background (+0.0, 0, 0.+0.I, all-zero bits): start from zeroed
     * pages and write only the explicit entries. */
    bool zero_bg = (dt == NDT_INT64) ? bg->data.integer == 0
                                     : (bre == 0.0 && !signbit(bre) && bim == 0.0 && !signbit(bim));
    void* buf = NULL;
    Expr* node = zero_bg ? ndbuild_open_zeroed(sa->rank, sa->dims, dt, &buf)
                         : ndbuild_open(sa->rank, sa->dims, dt, &buf);
    if (!node) return NULL;
    const SparsePattern* p = sa->pat;
    if (dt == NDT_FLOAT64) {
        double* b = buf;
        if (!zero_bg) for (int64_t i = 0; i < total; i++) b[i] = bre;
        const double* v = sa->vals;
        for (int64_t r = 0; r < p->nrows; r++)
            for (int64_t k = p->rowptr[r]; k < p->rowptr[r + 1]; k++)
                b[entry_offset(sa, r, k)] = v[k];
    } else if (dt == NDT_INT64) {
        int64_t* b = buf;
        int64_t bv = bg->data.integer;
        if (!zero_bg) for (int64_t i = 0; i < total; i++) b[i] = bv;
        const int64_t* v = sa->vals;
        for (int64_t r = 0; r < p->nrows; r++)
            for (int64_t k = p->rowptr[r]; k < p->rowptr[r + 1]; k++)
                b[entry_offset(sa, r, k)] = v[k];
    } else {
        double* b = buf;
        if (!zero_bg) for (int64_t i = 0; i < total; i++) { b[2 * i] = bre; b[2 * i + 1] = bim; }
        const double* v = sa->vals;
        for (int64_t r = 0; r < p->nrows; r++)
            for (int64_t k = p->rowptr[r]; k < p->rowptr[r + 1]; k++) {
                int64_t o = entry_offset(sa, r, k);
                b[2 * o] = v[2 * k]; b[2 * o + 1] = v[2 * k + 1];
            }
    }
    return node;
}

Expr* sp_normal(const SparseArray* sa) {
    int64_t total = sp_total(sa->rank, sa->dims);
    if (total < 0) return NULL;
    if (total > 0) {
        Expr* packed = normal_packed(sa, total);
        if (packed) return packed;
    }
    if (total > SP_NORMAL_MAX_BOXED) return NULL;
    /* A zero dimension: the array is a nest of empty Lists down to it, and
     * the fold never reaches a leaf. */
    if (total == 0) {
        int z = 0;
        while (sa->dims[z] != 0) z++;
        int64_t off = 0;
        return fold_slots(NULL, z + 1, sa->dims, 0, &off);
    }
    Expr** slots = malloc(sizeof(Expr*) * (size_t)total);
    if (!slots) return NULL;
    for (int64_t i = 0; i < total; i++) slots[i] = expr_copy(sa->bg);
    const SparsePattern* p = sa->pat;
    for (int64_t r = 0; r < p->nrows; r++)
        for (int64_t k = p->rowptr[r]; k < p->rowptr[r + 1]; k++) {
            int64_t o = entry_offset(sa, r, k);
            expr_free(slots[o]);
            slots[o] = sp_value_expr(sa, k);
        }
    int64_t off = 0;
    Expr* out = fold_slots(slots, sa->rank, sa->dims, 0, &off);
    free(slots);
    return out;
}

/* ======================================================================== */
/* Background equality                                                      */
/* ======================================================================== */

/* A machine-or-exact number's value as a double pair, when it has one. */
static bool number_value(const Expr* e, double* re, double* im) {
    *im = 0;
    switch (e->type) {
        case EXPR_INTEGER: *re = (double)e->data.integer; return true;
        case EXPR_REAL:    *re = e->data.real; return true;
        case EXPR_FUNCTION:
            if (is_machine_complex(e, re, im)) return true;
            return false;
        default: return false;
    }
}

bool sp_value_is_background(const Expr* v, const Expr* bg) {
    if (v == bg) return true;
    /* An inexact number against a number: compare values, so 0. is 0. */
    bool vi = v->type == EXPR_REAL || is_machine_complex(v, NULL, NULL);
    bool bi = bg->type == EXPR_REAL || is_machine_complex(bg, NULL, NULL);
    if (vi || bi) {
        double vr, vm, br, bm;
        if (number_value(v, &vr, &vm) && number_value(bg, &br, &bm))
            return vr == br && vm == bm;
    }
    return expr_eq(v, bg);
}

/* ======================================================================== */
/* COO builder                                                              */
/* ======================================================================== */

void sp_coo_init(SpCoo* c, int rank, const int64_t* dims, SpValKind kind) {
    memset(c, 0, sizeof *c);
    c->rank = rank;
    memcpy(c->dims, dims, sizeof(int64_t) * (size_t)rank);
    c->kind = kind;
}

bool sp_coo_reserve(SpCoo* c, int64_t extra) {
    if (c->n + extra <= c->cap) return true;
    int64_t cap = c->cap ? c->cap : 16;
    while (cap < c->n + extra) cap *= 2;
    size_t vsz = c->kind == SPV_EXPR ? sizeof(Expr*) : sp_kind_size(c->kind);
    int64_t* co = realloc(c->coords, sizeof(int64_t) * (size_t)cap * (size_t)(c->rank ? c->rank : 1));
    if (!co) return false;
    c->coords = co;
    void* va = realloc(c->vals, vsz * (size_t)cap);
    if (!va) return false;
    c->vals = va;
    c->cap = cap;
    return true;
}

static bool coo_slot(SpCoo* c, const int64_t* coords) {
    if (!sp_coo_reserve(c, 1)) return false;
    memcpy(c->coords + c->n * c->rank, coords, sizeof(int64_t) * (size_t)c->rank);
    return true;
}

bool sp_coo_push_expr(SpCoo* c, const int64_t* coords, Expr* v) {
    if (!coo_slot(c, coords)) { expr_free(v); return false; }
    ((Expr**)c->vals)[c->n++] = v;
    return true;
}
bool sp_coo_push_f64(SpCoo* c, const int64_t* coords, double v) {
    if (!coo_slot(c, coords)) return false;
    ((double*)c->vals)[c->n++] = v;
    return true;
}
bool sp_coo_push_i64(SpCoo* c, const int64_t* coords, int64_t v) {
    if (!coo_slot(c, coords)) return false;
    ((int64_t*)c->vals)[c->n++] = v;
    return true;
}
bool sp_coo_push_c128(SpCoo* c, const int64_t* coords, double re, double im) {
    if (!coo_slot(c, coords)) return false;
    double* b = c->vals;
    b[2 * c->n] = re; b[2 * c->n + 1] = im;
    c->n++;
    return true;
}

void sp_coo_free(SpCoo* c) {
    if (c->kind == SPV_EXPR && c->vals)
        for (int64_t i = 0; i < c->n; i++) if (((Expr**)c->vals)[i]) expr_free(((Expr**)c->vals)[i]);
    free(c->vals);
    free(c->coords);
    free(c->ck);
    memset(c, 0, sizeof *c);
}

/* ---- sorting ------------------------------------------------------------- */

/* Stable LSD radix sort of (key, idx) pairs by key, 11 bits per pass and only
 * as many passes as the largest key needs. `tk`/`ti` are scratch of length n. */
static void radix_sort_pairs(uint64_t* key, int64_t* idx, int64_t n,
                             uint64_t* tk, int64_t* ti) {
    uint64_t maxk = 0;
    for (int64_t i = 0; i < n; i++) if (key[i] > maxk) maxk = key[i];
    enum { BITS = 11, BUCKETS = 1 << 11 };
    int64_t count[BUCKETS];
    uint64_t* sk = key; int64_t* si = idx;
    uint64_t* dk = tk;  int64_t* di = ti;
    int shift = 0;
    while (shift < 64 && (shift == 0 || (maxk >> shift) != 0)) {
        memset(count, 0, sizeof count);
        for (int64_t i = 0; i < n; i++) count[(sk[i] >> shift) & (BUCKETS - 1)]++;
        int64_t sum = 0;
        for (int b = 0; b < BUCKETS; b++) { int64_t t = count[b]; count[b] = sum; sum += t; }
        for (int64_t i = 0; i < n; i++) {
            int64_t pos = count[(sk[i] >> shift) & (BUCKETS - 1)]++;
            dk[pos] = sk[i]; di[pos] = si[i];
        }
        uint64_t* t1 = sk; sk = dk; dk = t1;
        int64_t* t2 = si; si = di; di = t2;
        shift += BITS;
    }
    if (sk != key) {
        memcpy(key, sk, sizeof(uint64_t) * (size_t)n);
        memcpy(idx, si, sizeof(int64_t) * (size_t)n);
    }
}

/* Tuple mode (inner index space too large for a 64-bit key): compare entries'
 * tuples, ties broken by input order so the sort is stable. */
typedef struct { const int64_t* coords; int rank; int first; } TupleCtx;

static int tuple_order(const TupleCtx* t, int64_t a, int64_t b) {
    const int64_t* x = t->coords + a * t->rank;
    const int64_t* y = t->coords + b * t->rank;
    for (int j = t->first; j < t->rank; j++)
        if (x[j] != y[j]) return x[j] < y[j] ? -1 : 1;
    return a < b ? -1 : (a > b ? 1 : 0);
}

static void merge_sort_idx(int64_t* idx, int64_t* tmp, int64_t n, const TupleCtx* t) {
    if (n < 2) return;
    if (n <= 16) {
        for (int64_t i = 1; i < n; i++) {
            int64_t v = idx[i], j = i - 1;
            while (j >= 0 && tuple_order(t, idx[j], v) > 0) { idx[j + 1] = idx[j]; j--; }
            idx[j + 1] = v;
        }
        return;
    }
    int64_t h = n / 2;
    merge_sort_idx(idx, tmp, h, t);
    merge_sort_idx(idx + h, tmp, n - h, t);
    int64_t i = 0, j = h, o = 0;
    while (i < h && j < n) tmp[o++] = tuple_order(t, idx[i], idx[j]) <= 0 ? idx[i++] : idx[j++];
    while (i < h) tmp[o++] = idx[i++];
    while (j < n) tmp[o++] = idx[j++];
    memcpy(idx, tmp, sizeof(int64_t) * (size_t)n);
}

/* ---- finish -------------------------------------------------------------- */

/* Is value i of a typed or boxed buffer equal to the background? The typed
 * arms are sp_value_is_background specialised: an inexact value equals a
 * numeric background by value; an exact Integer equals an Integer exactly, or
 * an inexact background by value. */
typedef struct {
    const Expr* bg;
    bool num;               /* the background has a machine value */
    double re, im;
} BgTest;

static inline bool val_is_bg(SpValKind kind, const void* vals, int64_t i, const BgTest* t) {
    switch (kind) {
        case SPV_F64:
            return t->num && t->im == 0 && ((const double*)vals)[i] == t->re;
        case SPV_INT64: {
            int64_t v = ((const int64_t*)vals)[i];
            return t->bg->type == EXPR_INTEGER ? v == t->bg->data.integer
                                               : t->num && t->im == 0 && (double)v == t->re;
        }
        case SPV_C128: {
            const double* b = vals;
            return t->num && b[2 * i] == t->re && b[2 * i + 1] == t->im;
        }
        default:
            return sp_value_is_background(((Expr* const*)vals)[i], t->bg);
    }
}

/* Move value `from` of buffer `src` to slot `to` of `dst` (element size sz). */
static inline void val_move(void* dst, int64_t to, const void* src, int64_t from, size_t sz) {
    if (sz == 8) ((uint64_t*)dst)[to] = ((const uint64_t*)src)[from];
    else memcpy((char*)dst + (size_t)to * sz, (const char*)src + (size_t)from * sz, sz);
}

/* The general (tuple-compared) path: an inner index space too large for a
 * 64-bit key. Rare; a comparison sort over a permutation. */
static SparseArray* coo_finish_tuples(SpCoo* c, Expr* bg, bool drop_bg) {
    int rank = c->rank;
    int64_t n = c->n;
    int first = rank >= 2 ? 1 : 0;
    int width = rank >= 2 ? rank - 1 : 1;
    int64_t nrows = rank >= 2 ? c->dims[0] : 1;
    size_t vsz = c->kind == SPV_EXPR ? sizeof(Expr*) : sp_kind_size(c->kind);
    int64_t* perm = malloc(sizeof(int64_t) * (size_t)(n ? n : 1));
    int64_t* rp = calloc((size_t)(nrows + 1), sizeof(int64_t));
    int64_t* fill = malloc(sizeof(int64_t) * (size_t)(nrows + 1));
    int64_t* tmp = malloc(sizeof(int64_t) * (size_t)(n ? n : 1));
    if (!perm || !rp || !fill || !tmp) {
        free(perm); free(rp); free(fill); free(tmp);
        sp_coo_free(c); expr_free(bg);
        return NULL;
    }
    for (int64_t i = 0; i < n; i++) rp[(first ? c->coords[i * rank] : 0) + 1]++;
    for (int64_t r = 0; r < nrows; r++) rp[r + 1] += rp[r];
    memcpy(fill, rp, sizeof(int64_t) * (size_t)(nrows + 1));
    for (int64_t i = 0; i < n; i++) perm[fill[first ? c->coords[i * rank] : 0]++] = i;
    TupleCtx t = { c->coords, rank, first };
    for (int64_t r = 0; r < nrows; r++) merge_sort_idx(perm + rp[r], tmp, rp[r + 1] - rp[r], &t);
    free(tmp);

    BgTest bt = { bg, false, 0, 0 };
    bt.num = number_value(bg, &bt.re, &bt.im);
    SparsePattern* pat = sp_pattern_new(nrows, n, width, sp_dims_idx32(rank, c->dims));
    void* vals = malloc(vsz * (size_t)(n ? n : 1));
    if (!pat || !vals) {
        sp_pattern_release(pat); free(vals); free(perm); free(rp); free(fill);
        sp_coo_free(c); expr_free(bg);
        return NULL;
    }
    int64_t m = 0;
    pat->rowptr[0] = 0;
    for (int64_t r = 0; r < nrows; r++) {
        int64_t prev = -1;
        for (int64_t q = rp[r]; q < rp[r + 1]; q++) {
            int64_t i = perm[q];
            if (prev >= 0) {
                const int64_t* x = c->coords + prev * rank;
                const int64_t* y = c->coords + i * rank;
                bool same = true;
                for (int j = first; j < rank && same; j++) same = x[j] == y[j];
                if (same) continue;          /* a later duplicate: the first won */
            }
            prev = i;
            if (drop_bg && val_is_bg(c->kind, c->vals, i, &bt)) continue;
            const int64_t* x = c->coords + i * rank + first;
            for (int j = 0; j < width; j++) {
                if (pat->idx32) ((int32_t*)pat->colidx)[m * width + j] = (int32_t)x[j];
                else            ((int64_t*)pat->colidx)[m * width + j] = x[j];
            }
            val_move(vals, m, c->vals, i, vsz);
            if (c->kind == SPV_EXPR) ((Expr**)c->vals)[i] = NULL;   /* moved */
            m++;
        }
        pat->rowptr[r + 1] = m;
    }
    pat->nnz = m;
    free(perm); free(rp); free(fill);
    SpValKind kind = c->kind;
    int64_t dims[SP_MAX_RANK];
    memcpy(dims, c->dims, sizeof(int64_t) * (size_t)rank);
    sp_coo_free(c);
    return sp_new(rank, dims, bg, pat, kind, vals);
}

/* Number of bits needed to hold values 0 .. x. */
static int bits_for(uint64_t x) {
    int b = 0;
    while (x) { b++; x >>= 1; }
    return b;
}

/*
 * Composite path (the usual one): row and inner key fit together in 63 bits,
 * so every entry is one integer  ck = (row << kbits) | key  whose numeric
 * order IS canonical order. A stable LSD radix sort on ck (as few passes as
 * the bit count needs, at most 2048 buckets each, so every pass streams) puts
 * the entries in order with the specification order preserved among equal
 * positions; one sequential pass then keeps the first of each run, drops the
 * background, and splits ck back into row pointers and column tuples. A pass
 * whose digit is the same for every entry is skipped, and input already in
 * canonical order with no repeats skips the sort outright and has its value
 * buffer adopted rather than copied.
 */
static SparseArray* coo_finish_composite(SpCoo* c, Expr* bg, bool drop_bg, const int64_t* idims,
                                         int kbits, int total_bits) {
    int rank = c->rank;
    int64_t n = c->n;
    int first = rank >= 2 ? 1 : 0;
    int width = rank >= 2 ? rank - 1 : 1;
    int64_t nrows = rank >= 2 ? c->dims[0] : 1;
    size_t vsz = c->kind == SPV_EXPR ? sizeof(Expr*) : sp_kind_size(c->kind);
    uint64_t kmask = kbits >= 64 ? ~UINT64_C(0) : ((UINT64_C(1) << kbits) - 1);

    uint64_t* ck = c->ck;
    bool ordered = true;
    if (ck) {
        c->ck = NULL;                     /* the producer computed them */
        for (int64_t i = 1; i < n && ordered; i++) ordered = ck[i] > ck[i - 1];
    } else {
        ck = malloc(sizeof(uint64_t) * (size_t)(n ? n : 1));
        if (!ck) { sp_coo_free(c); expr_free(bg); return NULL; }
        for (int64_t i = 0; i < n; i++) {
            const int64_t* x = c->coords + i * rank;
            uint64_t k = 0;
            for (int j = 0; j < width; j++) k = k * (uint64_t)idims[j] + (uint64_t)x[first + j];
            uint64_t v = first ? (((uint64_t)x[0] << kbits) | k) : k;
            ck[i] = v;
            if (i > 0 && v <= ck[i - 1]) ordered = false;
        }
    }

    void* sv = c->vals;                 /* the buffer holding the sorted values */
    int ibits = bits_for((uint64_t)(n > 0 ? n - 1 : 0));
    if (!ordered && !c->ordered && n > 1 && total_bits + ibits <= 64) {
        /* Key and input index share one word, (ck << ibits) | i, so each
         * radix pass streams a single 8-byte array; the values follow in one
         * gather at the end. Stable: equal keys keep their input order. */
        uint64_t* a = malloc(sizeof(uint64_t) * (size_t)n);
        uint64_t* t = malloc(sizeof(uint64_t) * (size_t)n);
        void* nv = malloc(vsz * (size_t)n);
        if (!a || !t || !nv) { free(a); free(t); free(nv); free(ck); sp_coo_free(c); expr_free(bg); return NULL; }
        for (int64_t i = 0; i < n; i++) a[i] = (ck[i] << ibits) | (uint64_t)i;
        /* 12-bit digits: fewest passes while the 4096 bucket heads stay in
         * L1 (measured: three 12-bit passes beat four 9-bit ones by ~30% on
         * 10^6 keys). */
        int passes = (total_bits + 11) / 12;
        if (passes < 1) passes = 1;
        int dw = (total_bits + passes - 1) / passes;
        uint64_t dmask = (UINT64_C(1) << dw) - 1;
        int64_t count[4096];
        uint64_t* src = a; uint64_t* dst = t;
        for (int p = 0; p < passes; p++) {
            int shift = ibits + p * dw;
            memset(count, 0, sizeof(int64_t) * ((size_t)dmask + 1));
            for (int64_t i = 0; i < n; i++) count[(src[i] >> shift) & dmask]++;
            if (count[(src[0] >> shift) & dmask] == n) continue;   /* one digit everywhere */
            int64_t sum = 0;
            for (uint64_t d = 0; d <= dmask; d++) { int64_t q = count[d]; count[d] = sum; sum += q; }
            for (int64_t i = 0; i < n; i++) dst[count[(src[i] >> shift) & dmask]++] = src[i];
            uint64_t* sw = src; src = dst; dst = sw;
        }
        uint64_t imask = ibits >= 64 ? ~UINT64_C(0) : ((UINT64_C(1) << ibits) - 1);
        if (vsz == 8) {
            const uint64_t* ov = c->vals; uint64_t* w = nv;
            for (int64_t q = 0; q < n; q++) { ck[q] = src[q] >> ibits; w[q] = ov[src[q] & imask]; }
        } else {
            const uint64_t* ov = c->vals; uint64_t* w = nv;
            for (int64_t q = 0; q < n; q++) {
                int64_t i = (int64_t)(src[q] & imask);
                ck[q] = src[q] >> ibits;
                w[2 * q] = ov[2 * i]; w[2 * q + 1] = ov[2 * i + 1];
            }
        }
        free(a); free(t);
        free(c->vals);                  /* raw: its values were moved, not copied */
        c->vals = nv;
        sv = nv;
    } else if (!ordered && !c->ordered && n > 1) {
        uint64_t* tk = malloc(sizeof(uint64_t) * (size_t)n);
        void* tv = malloc(vsz * (size_t)n);
        if (!tk || !tv) { free(tk); free(tv); free(ck); sp_coo_free(c); expr_free(bg); return NULL; }
        int passes = (total_bits + 10) / 11;
        if (passes < 1) passes = 1;
        int dw = (total_bits + passes - 1) / passes;
        uint64_t dmask = (UINT64_C(1) << dw) - 1;
        int64_t count[2048];
        uint64_t* sk = ck; uint64_t* dk = tk;
        void* svv = c->vals; void* dvv = tv;
        for (int p = 0; p < passes; p++) {
            int shift = p * dw;
            memset(count, 0, sizeof(int64_t) * ((size_t)dmask + 1));
            for (int64_t i = 0; i < n; i++) count[(sk[i] >> shift) & dmask]++;
            if (count[(sk[0] >> shift) & dmask] == n) continue;   /* one digit everywhere */
            int64_t sum = 0;
            for (uint64_t d = 0; d <= dmask; d++) { int64_t t = count[d]; count[d] = sum; sum += t; }
            if (vsz == 8) {
                const uint64_t* a = svv; uint64_t* b = dvv;
                for (int64_t i = 0; i < n; i++) {
                    int64_t q = count[(sk[i] >> shift) & dmask]++;
                    dk[q] = sk[i]; b[q] = a[i];
                }
            } else {
                const uint64_t* a = svv; uint64_t* b = dvv;     /* 16-byte values */
                for (int64_t i = 0; i < n; i++) {
                    int64_t q = count[(sk[i] >> shift) & dmask]++;
                    dk[q] = sk[i]; b[2 * q] = a[2 * i]; b[2 * q + 1] = a[2 * i + 1];
                }
            }
            uint64_t* t1 = sk; sk = dk; dk = t1;
            void* t2 = svv; svv = dvv; dvv = t2;
        }
        /* keep whichever buffers hold the result; the others are scratch */
        uint64_t* scratch_k = sk == ck ? tk : ck;
        void* scratch_v = svv == c->vals ? tv : c->vals;   /* raw: its pointers were moved */
        ck = sk;
        c->vals = svv;
        free(scratch_k);
        free(scratch_v);
        sv = c->vals;
    }

    /* dedup (first wins), drop the background, decode */
    BgTest bt = { bg, false, 0, 0 };
    bt.num = number_value(bg, &bt.re, &bt.im);
    bool idx32 = sp_dims_idx32(rank, c->dims);
    SparsePattern* pat = sp_pattern_new(nrows, n, width, idx32);
    if (!pat) { free(ck); sp_coo_free(c); expr_free(bg); return NULL; }
    memset(pat->rowptr, 0, sizeof(int64_t) * (size_t)(nrows + 1));
    int64_t m = 0;
    for (int64_t q = 0; q < n; q++) {
        uint64_t v = ck[q];
        if ((q > 0 && v == ck[q - 1]) || (drop_bg && val_is_bg(c->kind, sv, q, &bt))) {
            if (c->kind == SPV_EXPR) { expr_free(((Expr**)sv)[q]); ((Expr**)sv)[q] = NULL; }
            continue;
        }
        uint64_t k = v & kmask;
        if (first) pat->rowptr[(v >> kbits) + 1]++;
        if (width == 1) {
            if (idx32) ((int32_t*)pat->colidx)[m] = (int32_t)k;
            else       ((int64_t*)pat->colidx)[m] = (int64_t)k;
        } else {
            for (int j = width - 1; j >= 0; j--) {
                int64_t cj = (int64_t)(k % (uint64_t)idims[j]);
                k /= (uint64_t)idims[j];
                if (idx32) ((int32_t*)pat->colidx)[m * width + j] = (int32_t)cj;
                else       ((int64_t*)pat->colidx)[m * width + j] = cj;
            }
        }
        if (m != q) val_move(sv, m, sv, q, vsz);
        m++;
    }
    if (first) for (int64_t r = 0; r < nrows; r++) pat->rowptr[r + 1] += pat->rowptr[r];
    else pat->rowptr[1] = m;
    pat->nnz = m;
    free(ck);
    if (m < n) {
        size_t isz = idx32 ? sizeof(int32_t) : sizeof(int64_t);
        if (m == 0) { free(pat->colidx); pat->colidx = NULL; }
        else {
            void* ci = realloc(pat->colidx, isz * (size_t)m * (size_t)width);
            if (ci) pat->colidx = ci;
        }
    }
    void* vals = c->vals;
    c->vals = NULL;                    /* adopted */
    c->n = 0;
    if (m == 0) { free(vals); vals = NULL; }
    else if (m < n) { void* vv = realloc(vals, vsz * (size_t)m); if (vv) vals = vv; }
    SpValKind kind = c->kind;
    int64_t dims[SP_MAX_RANK];
    memcpy(dims, c->dims, sizeof(int64_t) * (size_t)rank);
    sp_coo_free(c);
    return sp_new(rank, dims, bg, pat, kind, vals);
}

bool sp_coo_composite(int rank, const int64_t* dims, int* kbits) {
    int first = rank >= 2 ? 1 : 0;
    uint64_t span = 1;
    for (int j = first; j < rank; j++) {
        uint64_t d = (uint64_t)(dims[j] > 0 ? dims[j] : 1);
        if (span > (UINT64_C(1) << 62) / d) return false;
        span *= d;
    }
    int kb = bits_for(span - 1);
    int rb = first ? bits_for((uint64_t)(dims[0] > 0 ? dims[0] - 1 : 0)) : 0;
    if (kb + rb > 63) return false;
    *kbits = kb;
    return true;
}

/*
 * Keyed path (the normal one). Each entry's inner coordinates are folded into
 * one 64-bit row-major key, then:
 *   1. one pass scatters (key, value) into row order -- stable, so within a
 *      row the entries keep their specification order;
 *   2. each row is sorted (insertion for short rows, LSD radix for long ones,
 *      both stable) on contiguous, cache-resident memory;
 *   3. one sequential pass keeps the first of each run of equal keys (the
 *      first rule wins), drops background values, and decodes the key back
 *      into the column tuple.
 * Input already in canonical order with no repeats -- every dense conversion,
 * most generated rule lists -- skips 1 and 2, and its value buffer is adopted
 * rather than copied.
 */
SparseArray* sp_coo_finish(SpCoo* c, Expr* bg, bool drop_bg) {
    int rank = c->rank;
    int64_t n = c->n;
    int first = rank >= 2 ? 1 : 0;
    int width = rank >= 2 ? rank - 1 : 1;
    int64_t nrows = rank >= 2 ? c->dims[0] : 1;

    /* inner dims, and whether their product fits a 62-bit key */
    int64_t idims[SP_MAX_RANK];
    uint64_t span = 1;
    for (int j = 0; j < width; j++) {
        idims[j] = c->dims[first + j] > 0 ? c->dims[first + j] : 1;
        if (span > (UINT64_C(1) << 62) / (uint64_t)idims[j]) return coo_finish_tuples(c, bg, drop_bg);
        span *= (uint64_t)idims[j];
    }
    /* the composite (row, key) integer fits: the radix path */
    {
        int kbits = bits_for(span - 1);
        int rbits = first ? bits_for((uint64_t)(nrows > 0 ? nrows - 1 : 0)) : 0;
        if (kbits + rbits <= 63) return coo_finish_composite(c, bg, drop_bg, idims, kbits, kbits + rbits);
    }
    if (c->ck) { sp_coo_free(c); expr_free(bg); return NULL; }   /* producer contract broken */
    size_t vsz = c->kind == SPV_EXPR ? sizeof(Expr*) : sp_kind_size(c->kind);

    uint64_t* key = malloc(sizeof(uint64_t) * (size_t)(n ? n : 1));
    int64_t* rp = calloc((size_t)(nrows + 1), sizeof(int64_t));
    if (!key || !rp) { free(key); free(rp); sp_coo_free(c); expr_free(bg); return NULL; }

    /* keys and row counts, sequentially; and is the input already canonical? */
    bool ordered = true;
    int64_t prev_row = -1;
    uint64_t prev_key = 0;
    for (int64_t i = 0; i < n; i++) {
        const int64_t* x = c->coords + i * rank;
        int64_t row = first ? x[0] : 0;
        uint64_t k = 0;
        for (int j = 0; j < width; j++) k = k * (uint64_t)idims[j] + (uint64_t)x[first + j];
        key[i] = k;
        rp[row + 1]++;
        if (ordered && i > 0 && (row < prev_row || (row == prev_row && k <= prev_key))) ordered = false;
        prev_row = row;
        prev_key = k;
    }
    for (int64_t r = 0; r < nrows; r++) rp[r + 1] += rp[r];

    uint64_t* skey = key;
    void* svals = c->vals;
    if (!ordered && !c->ordered) {
        /* 1. scatter into row order */
        skey = malloc(sizeof(uint64_t) * (size_t)(n ? n : 1));
        svals = malloc(vsz * (size_t)(n ? n : 1));
        int64_t* fill = malloc(sizeof(int64_t) * (size_t)(nrows + 1));
        if (!skey || !svals || !fill) {
            free(skey); free(svals); free(fill); free(key); free(rp);
            sp_coo_free(c); expr_free(bg);
            return NULL;
        }
        memcpy(fill, rp, sizeof(int64_t) * (size_t)(nrows + 1));
        for (int64_t i = 0; i < n; i++) {
            int64_t p = fill[first ? c->coords[i * rank] : 0]++;
            skey[p] = key[i];
            val_move(svals, p, c->vals, i, vsz);
        }
        free(fill);
        free(key);
        key = NULL;
        if (c->kind == SPV_EXPR) memset(c->vals, 0, sizeof(Expr*) * (size_t)n);   /* moved */

        /* 2. sort each row (stable); a value permutation follows the keys */
        int64_t maxlen = 0;
        for (int64_t r = 0; r < nrows; r++)
            if (rp[r + 1] - rp[r] > maxlen) maxlen = rp[r + 1] - rp[r];
        int64_t* idx = NULL; uint64_t* tk = NULL; int64_t* ti = NULL; char* tv = NULL;
        for (int64_t r = 0; r < nrows; r++) {
            int64_t lo = rp[r], len = rp[r + 1] - lo;
            if (len < 2) continue;
            uint64_t* rk = skey + lo;
            bool sorted = true;
            for (int64_t i = 1; i < len && sorted; i++) sorted = rk[i - 1] <= rk[i];
            if (sorted) continue;
            if (len <= 32) {
                /* insertion sort carrying the values */
                char* rv = (char*)svals + (size_t)lo * vsz;
                char tmpv[16];
                for (int64_t i = 1; i < len; i++) {
                    uint64_t kk = rk[i];
                    memcpy(tmpv, rv + (size_t)i * vsz, vsz);
                    int64_t j = i - 1;
                    while (j >= 0 && rk[j] > kk) {
                        rk[j + 1] = rk[j];
                        memcpy(rv + (size_t)(j + 1) * vsz, rv + (size_t)j * vsz, vsz);
                        j--;
                    }
                    rk[j + 1] = kk;
                    memcpy(rv + (size_t)(j + 1) * vsz, tmpv, vsz);
                }
                continue;
            }
            if (!idx) {
                idx = malloc(sizeof(int64_t) * (size_t)maxlen);
                tk = malloc(sizeof(uint64_t) * (size_t)maxlen);
                ti = malloc(sizeof(int64_t) * (size_t)maxlen);
                tv = malloc(vsz * (size_t)maxlen);
                if (!idx || !tk || !ti || !tv) {
                    free(idx); free(tk); free(ti); free(tv);
                    free(skey); free(rp);
                    if (c->kind == SPV_EXPR) for (int64_t q = 0; q < n; q++) expr_free(((Expr**)svals)[q]);
                    free(svals);
                    sp_coo_free(c); expr_free(bg);
                    return NULL;
                }
            }
            for (int64_t i = 0; i < len; i++) idx[i] = i;
            radix_sort_pairs(rk, idx, len, tk, ti);
            char* rv = (char*)svals + (size_t)lo * vsz;
            for (int64_t i = 0; i < len; i++) memcpy(tv + (size_t)i * vsz, rv + (size_t)idx[i] * vsz, vsz);
            memcpy(rv, tv, (size_t)len * vsz);
        }
        free(idx); free(tk); free(ti); free(tv);
    }

    /* 3. dedup (first wins), drop the background, decode the keys */
    BgTest bt = { bg, false, 0, 0 };
    bt.num = number_value(bg, &bt.re, &bt.im);
    bool idx32 = sp_dims_idx32(rank, c->dims);
    SparsePattern* pat = sp_pattern_new(nrows, n, width, idx32);
    if (!pat) {
        if (skey != key) { free(skey);
            if (c->kind == SPV_EXPR) for (int64_t q = 0; q < n; q++) expr_free(((Expr**)svals)[q]);
            free(svals); }
        free(key); free(rp); sp_coo_free(c); expr_free(bg);
        return NULL;
    }
    int64_t m = 0;
    pat->rowptr[0] = 0;
    for (int64_t r = 0; r < nrows; r++) {
        int64_t lo = rp[r], hi = rp[r + 1];
        for (int64_t q = lo; q < hi; q++) {
            if (q > lo && skey[q] == skey[q - 1]) {
                /* a later duplicate: the first one won */
                if (c->kind == SPV_EXPR) { expr_free(((Expr**)svals)[q]); ((Expr**)svals)[q] = NULL; }
                continue;
            }
            if (drop_bg && val_is_bg(c->kind, svals, q, &bt)) {
                if (c->kind == SPV_EXPR) { expr_free(((Expr**)svals)[q]); ((Expr**)svals)[q] = NULL; }
                continue;
            }
            uint64_t k = skey[q];
            if (width == 1) {
                if (idx32) ((int32_t*)pat->colidx)[m] = (int32_t)k;
                else       ((int64_t*)pat->colidx)[m] = (int64_t)k;
            } else {
                for (int j = width - 1; j >= 0; j--) {
                    int64_t cj = (int64_t)(k % (uint64_t)idims[j]);
                    k /= (uint64_t)idims[j];
                    if (idx32) ((int32_t*)pat->colidx)[m * width + j] = (int32_t)cj;
                    else       ((int64_t*)pat->colidx)[m * width + j] = cj;
                }
            }
            if (m != q) val_move(svals, m, svals, q, vsz);
            m++;
        }
        pat->rowptr[r + 1] = m;
    }
    pat->nnz = m;
    if (m < n) {
        size_t isz = idx32 ? sizeof(int32_t) : sizeof(int64_t);
        if (m == 0) { free(pat->colidx); pat->colidx = NULL; }
        else {
            void* ci = realloc(pat->colidx, isz * (size_t)m * (size_t)width);
            if (ci) pat->colidx = ci;
        }
    }
    /* the value buffer: adopted from the builder when it was not scattered */
    void* vals = svals;
    if (svals == c->vals) c->vals = NULL;     /* adopted: the builder must not free it */
    if (m < n && m > 0) {
        void* vv = realloc(vals, vsz * (size_t)m);
        if (vv) vals = vv;
    }
    if (skey != key) free(skey);
    free(key);
    free(rp);
    SpValKind kind = c->kind;
    int64_t dims[SP_MAX_RANK];
    memcpy(dims, c->dims, sizeof(int64_t) * (size_t)rank);
    c->n = 0;                 /* boxed values were moved or freed above */
    sp_coo_free(c);
    if (m == 0) { free(vals); vals = NULL; }
    return sp_new(rank, dims, bg, pat, kind, vals);
}

/* ======================================================================== */
/* Value maps                                                               */
/* ======================================================================== */

Expr* sp_map_values(const Expr* s, SpValueFn fn, void* ctx, bool int_to_real) {
    const SparseArray* sa = sp_of(s);
    if (!sa) return NULL;
    int64_t n = sa->pat->nnz;
    Expr* bg = fn(sa->bg, ctx);
    if (!bg) return NULL;
    SpValKind kind = sa->kind;
    void* vals = NULL;
    if (int_to_real && kind != SPV_EXPR) {
        size_t sz = sp_kind_size(kind == SPV_INT64 ? SPV_F64 : kind) * (size_t)n;
        vals = n ? malloc(sz) : NULL;
        if (n && !vals) { expr_free(bg); return NULL; }
        if (kind == SPV_INT64) {
            for (int64_t k = 0; k < n; k++) ((double*)vals)[k] = (double)((const int64_t*)sa->vals)[k];
            kind = SPV_F64;
        } else if (n) {
            memcpy(vals, sa->vals, sz);
        }
    } else {
        Expr** ev = malloc(sizeof(Expr*) * (size_t)(n ? n : 1));
        if (!ev) { expr_free(bg); return NULL; }
        for (int64_t k = 0; k < n; k++) {
            Expr* v = sp_value_expr(sa, k);
            ev[k] = fn(v, ctx);
            expr_free(v);
            if (!ev[k]) {
                for (int64_t j = 0; j < k; j++) expr_free(ev[j]);
                free(ev); expr_free(bg);
                return NULL;
            }
        }
        vals = ev;
        kind = SPV_EXPR;
    }
    return sp_to_expr(sp_new(sa->rank, sa->dims, bg, sp_pattern_retain(sa->pat), kind, vals));
}
