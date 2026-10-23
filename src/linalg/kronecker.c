/* KroneckerProduct[a, b, ...] for dense arrays (List, packed List, NDArray).
 *
 * A SparseArray argument never reaches this builtin first: the head's native
 * sparse kernel (src/sparse/sparse_struct.c, nat_kronecker) is dispatched by
 * the evaluator before it, and only a call that kernel declines arrives here
 * (densified by the sparse gate).
 *
 * SHAPE RULE.  For two arrays of ranks ra and rb (R = ra + rb) the product is a
 * re-grouping of their tensor product T (whose R indices are a's followed by
 * b's):
 *
 *   - two vectors, or a vector beside an array of even rank >= 4: T itself
 *     (the outer product; a vector pair is a column times a row);
 *   - a matrix then a vector: the vector is a row, giving m x (n * len);
 *   - otherwise, with h = floor(R / 2): result axis j < h joins T-index j (the
 *     major part) with T-index h + j (the minor part), and when R is odd the
 *     last T-index is a final axis of its own.  For two matrices this is the
 *     usual block product (m p) x (n q); for equal ranks r it is the
 *     element-wise product of the dimensions.
 *
 * Three or more arguments fold from the left; every argument but the last must
 * then have rank <= 2 (otherwise the call stays unevaluated).
 *
 * VALUES.  Each entry is a_i * b_j.  When both operands consist of numbers and
 * either holds a machine number, the whole product is machine (an exact 2 next
 * to a 1. comes out as 2.): computed straight into a float64 / complex buffer.
 * Two exact-integer operands multiply in int64 (ci_mul_i64), abandoning to the
 * exact boxed path on the first overflow so big products stay exact.  Anything
 * else (symbols, Rationals beside exact numbers, ...) is Times[a_i, b_j],
 * evaluated.
 *
 * The buffer result inherits the presentation of an array operand (packed List
 * or visible NDArray), and is a fresh packed producer result for two plain
 * Lists.  Rank and layout are computed once as a stride per T-index, so the
 * fill is one multiply and one store per entry, whatever the ranks. */

#include "linalg.h"
#include "sym_names.h"
#include "common.h"
#include "eval.h"
#include "pack.h"
#include "ndarray.h"
#include "message.h"
#include "print.h"
#include "checked_int.h"
#include <stdlib.h>
#include <string.h>

#define KP_MAX_RANK 64

/* Leaf classes, as bits: what the product's value path needs to know. */
enum {
    KL_INT     = 1u << 0,   /* machine Integer                                */
    KL_REAL    = 1u << 1,   /* machine Real                                   */
    KL_CPLX_M  = 1u << 2,   /* Complex with a machine part                    */
    KL_CPLX_X  = 1u << 3,   /* exact Complex with machine-integer parts       */
    KL_NUM     = 1u << 4,   /* any other number (BigInt, Rational, MPFR, ...) */
    KL_OTHER   = 1u << 5    /* not a number                                   */
};

/* One operand: its shape and either a buffer or its leaves (row-major). */
typedef struct {
    int         rank;
    int64_t     dims[KP_MAX_RANK];
    int64_t     n;          /* product of dims                                */
    const Expr* nd;         /* the buffer operand (borrowed), or NULL         */
    Expr**      leaves;     /* owned, for a List operand (and on demand)      */
    unsigned    kinds;      /* union of KL_* over the leaves                  */
} KpOp;

typedef enum { KP_OK, KP_SILENT, KP_RECT } KpStatus;

static bool is_list(const Expr* e) { return head_is(e, SYM_List); }

static void kpop_free(KpOp* o) {
    if (o->leaves) {
        for (int64_t i = 0; i < o->n; i++) if (o->leaves[i]) expr_free(o->leaves[i]);
        free(o->leaves);
    }
    o->leaves = NULL;
}

static unsigned leaf_kind(const Expr* e) {
    switch (e->type) {
        case EXPR_INTEGER: return KL_INT;
        case EXPR_REAL:    return KL_REAL;
        case EXPR_BIGINT:  return KL_NUM;
#ifdef USE_MPFR
        case EXPR_MPFR:    return KL_NUM;
#endif
        default: break;
    }
    if (head_is(e, SYM_Rational) && e->data.function.arg_count == 2) return KL_NUM;
    if (head_is(e, SYM_Complex) && e->data.function.arg_count == 2) {
        const Expr* re = e->data.function.args[0];
        const Expr* im = e->data.function.args[1];
        if (re->type == EXPR_REAL || im->type == EXPR_REAL) {
            bool ok = (re->type == EXPR_REAL || re->type == EXPR_INTEGER) &&
                      (im->type == EXPR_REAL || im->type == EXPR_INTEGER);
            return ok ? KL_CPLX_M : KL_NUM;
        }
        if (re->type == EXPR_INTEGER && im->type == EXPR_INTEGER) return KL_CPLX_X;
        return KL_NUM;
    }
    return KL_OTHER;
}

static unsigned nd_kind(NDType dt) {
    switch (dt) {
        case NDT_INT64:     return KL_INT;
        case NDT_FLOAT64:
        case NDT_FLOAT32:   return KL_REAL;
        case NDT_COMPLEX64:
        case NDT_COMPLEX32: return KL_CPLX_M;
        default:            return KL_OTHER;   /* bool: True/False are not numbers */
    }
}

/* Walk a List level by level against the shape already read off its first
 * elements, collecting the leaves when `out` is non-NULL. A wrong length or a
 * non-List where a sub-array belongs is a ragged array (KP_RECT); a List (or a
 * nested buffer) where a leaf belongs is an array of a deeper ragged shape the
 * product does not take (KP_SILENT). */
static KpStatus walk(const Expr* e, const KpOp* o, int level, Expr** out, int64_t* pos,
                     unsigned* kinds) {
    if (level == o->rank) {
        if (is_list(e) || e->type == EXPR_NDARRAY) return KP_SILENT;
        *kinds |= leaf_kind(e);
        if (out) out[(*pos)++] = expr_copy((Expr*)e);
        return KP_OK;
    }
    if (e->type == EXPR_NDARRAY) return KP_SILENT;
    if (!is_list(e) || (int64_t)e->data.function.arg_count != o->dims[level]) return KP_RECT;
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        KpStatus s = walk(e->data.function.args[i], o, level + 1, out, pos, kinds);
        if (s != KP_OK) return s;
    }
    return KP_OK;
}

/* Read an operand's shape (and classify its leaves; collect them when
 * `collect`). A number is not an array (KP_RECT); any other non-List is a
 * symbolic operand the call keeps (KP_SILENT). */
static KpStatus kp_read(const Expr* e, KpOp* o, bool collect) {
    memset(o, 0, sizeof *o);
    if (e->type == EXPR_NDARRAY) {
        o->rank = e->data.ndarray.rank;
        if (o->rank > KP_MAX_RANK) return KP_SILENT;
        o->n = 1;
        for (int i = 0; i < o->rank; i++) {
            o->dims[i] = e->data.ndarray.dims[i];
            o->n *= o->dims[i];
        }
        o->nd = e;
        o->kinds = nd_kind(e->data.ndarray.dtype);
        return KP_OK;
    }
    if (!is_list(e)) return (leaf_kind(e) & KL_OTHER) ? KP_SILENT : KP_RECT;
    /* the shape along the first elements */
    const Expr* cur = e;
    o->n = 1;
    while (is_list(cur)) {
        if (o->rank == KP_MAX_RANK) return KP_SILENT;
        int64_t len = (int64_t)cur->data.function.arg_count;
        o->dims[o->rank++] = len;
        if (ci_mul_i64(o->n, len, &o->n)) return KP_SILENT;
        if (len == 0) break;
        cur = cur->data.function.args[0];
    }
    if (collect) {
        o->leaves = calloc((size_t)(o->n ? o->n : 1), sizeof(Expr*));
        if (!o->leaves) return KP_SILENT;
    }
    int64_t pos = 0;
    KpStatus s = walk(e, o, 0, o->leaves, &pos, &o->kinds);
    if (s != KP_OK) kpop_free(o);
    return s;
}

/* The leaves of a buffer operand, as boxed values (for the exact path). */
static bool kp_box(KpOp* o) {
    if (o->leaves || !o->nd) return true;
    o->leaves = calloc((size_t)(o->n ? o->n : 1), sizeof(Expr*));
    if (!o->leaves) return false;
    for (int64_t i = 0; i < o->n; i++) o->leaves[i] = ndarray_element_to_expr(o->nd, (size_t)i);
    return true;
}

/* ------------------------------------------------------------------ layout */

typedef struct {
    int     rank;
    int64_t dims[KP_MAX_RANK];
    int64_t total;
    int64_t sa[KP_MAX_RANK], sb[KP_MAX_RANK];   /* out stride per index of a, b */
} KpLayout;

static bool kp_layout(const KpOp* a, const KpOp* b, KpLayout* L) {
    int ra = a->rank, rb = b->rank;
    int64_t D[2 * KP_MAX_RANK];
    int64_t st[2 * KP_MAX_RANK];
    int R = 0;
    for (int i = 0; i < ra; i++) D[R++] = a->dims[i];
    int bstart = R;
    bool row = (ra == 2 && rb == 1);
    if (row) D[R++] = 1;                       /* the vector as a 1 x n row */
    for (int i = 0; i < rb; i++) D[R++] = b->dims[i];
    bool outer = (ra == 1 && rb == 1) || (ra == 1 && rb >= 4 && rb % 2 == 0) ||
                 (rb == 1 && ra >= 4 && ra % 2 == 0);
    if (outer) {
        if (R > KP_MAX_RANK) return false;
        L->rank = R;
        memcpy(L->dims, D, sizeof(int64_t) * (size_t)R);
        int64_t s = 1;
        for (int t = R - 1; t >= 0; t--) {
            st[t] = s;
            if (ci_mul_i64(s, D[t] ? D[t] : 1, &s)) return false;
        }
    } else {
        int h = R / 2, m = R - h;
        if (m > KP_MAX_RANK) return false;
        L->rank = m;
        for (int j = 0; j < h; j++)
            if (ci_mul_i64(D[j], D[h + j], &L->dims[j])) return false;
        if (R % 2) L->dims[h] = D[R - 1];
        int64_t S[KP_MAX_RANK];
        int64_t s = 1;
        for (int j = m - 1; j >= 0; j--) {
            S[j] = s;
            if (ci_mul_i64(s, L->dims[j] ? L->dims[j] : 1, &s)) return false;
        }
        for (int t = 0; t < h; t++) {
            if (ci_mul_i64(D[h + t], S[t], &st[t])) return false;   /* major part */
            st[h + t] = S[t];                                       /* minor part */
        }
        if (R % 2) st[R - 1] = S[h];
    }
    L->total = 1;
    for (int j = 0; j < L->rank; j++)
        if (ci_mul_i64(L->total, L->dims[j], &L->total)) return false;
    for (int i = 0; i < ra; i++) L->sa[i] = st[i];
    for (int i = 0; i < rb; i++) L->sb[i] = st[bstart + (row ? 1 : 0) + i];
    return true;
}

/* off[i] = the output offset contributed by element i (row-major) of an
 * operand with these dims and per-index strides. */
static int64_t* kp_offsets(int rank, const int64_t* dims, const int64_t* str, int64_t n) {
    int64_t* off = malloc(sizeof(int64_t) * (size_t)(n ? n : 1));
    if (!off || n == 0) return off;
    int64_t idx[KP_MAX_RANK] = { 0 };
    int64_t cur = 0;
    for (int64_t i = 0; i < n; i++) {
        off[i] = cur;
        for (int d = rank - 1; d >= 0; d--) {          /* odometer step */
            if (++idx[d] < dims[d]) { cur += str[d]; break; }
            cur -= str[d] * (dims[d] - 1);
            idx[d] = 0;
        }
    }
    return off;
}

/* ------------------------------------------------------------------ values */

static void leaf_c(const Expr* e, double* re, double* im) {
    if (e->type == EXPR_INTEGER) { *re = (double)e->data.integer; *im = 0; return; }
    if (e->type == EXPR_REAL)    { *re = e->data.real; *im = 0; return; }
    const Expr* p = e->data.function.args[0];
    const Expr* q = e->data.function.args[1];
    *re = p->type == EXPR_REAL ? p->data.real : (double)p->data.integer;
    *im = q->type == EXPR_REAL ? q->data.real : (double)q->data.integer;
}

static void op_c(const KpOp* o, int64_t i, double* re, double* im) {
    if (o->nd) ndt_get(o->nd->data.ndarray.data, (size_t)i, o->nd->data.ndarray.dtype, re, im);
    else leaf_c(o->leaves[i], re, im);
}

static int64_t op_i(const KpOp* o, int64_t i) {
    if (o->nd) return ((const int64_t*)o->nd->data.ndarray.data)[i];
    return o->leaves[i]->data.integer;
}

/* Fill a machine buffer. Returns false on an int64 overflow (int path only). */
static bool kp_fill(const KpOp* a, const KpOp* b, const int64_t* oa, const int64_t* ob,
                    void* buf, NDType dt) {
    if (dt == NDT_INT64) {
        int64_t* o = buf;
        for (int64_t i = 0; i < a->n; i++) {
            int64_t x = op_i(a, i);
            for (int64_t j = 0; j < b->n; j++)
                if (ci_mul_i64(x, op_i(b, j), &o[oa[i] + ob[j]])) return false;
        }
        return true;
    }
    bool fast = a->nd && b->nd && a->nd->data.ndarray.dtype == NDT_FLOAT64 &&
                b->nd->data.ndarray.dtype == NDT_FLOAT64 && dt == NDT_FLOAT64;
    if (fast) {
        const double* pa = a->nd->data.ndarray.data;
        const double* pb = b->nd->data.ndarray.data;
        double* o = buf;
        for (int64_t i = 0; i < a->n; i++) {
            double x = pa[i];
            double* oi = o + oa[i];
            for (int64_t j = 0; j < b->n; j++) oi[ob[j]] = x * pb[j];
        }
        return true;
    }
    double* bre = malloc(sizeof(double) * (size_t)(b->n ? b->n : 1));
    double* bim = malloc(sizeof(double) * (size_t)(b->n ? b->n : 1));
    if (!bre || !bim) { free(bre); free(bim); return false; }
    for (int64_t j = 0; j < b->n; j++) op_c(b, j, &bre[j], &bim[j]);
    bool cplx = ndt_is_complex(dt);
    for (int64_t i = 0; i < a->n; i++) {
        double xr, xi;
        op_c(a, i, &xr, &xi);
        for (int64_t j = 0; j < b->n; j++) {
            double re = xr * bre[j], im = 0;
            if (cplx) {
                re -= xi * bim[j];
                im = xr * bim[j] + xi * bre[j];
            }
            ndt_set(buf, (size_t)(oa[i] + ob[j]), dt, re, im);
        }
    }
    free(bre);
    free(bim);
    return true;
}

/* A nested List of the given dims over `flat` (consumed, row-major). */
static Expr* nest(const int64_t* dims, int rank, Expr** flat, int64_t* pos) {
    int64_t len = dims[0];
    Expr** args = malloc(sizeof(Expr*) * (size_t)(len ? len : 1));
    for (int64_t i = 0; i < len; i++)
        args[i] = rank == 1 ? flat[(*pos)++] : nest(dims + 1, rank - 1, flat, pos);
    Expr* r = expr_new_function(expr_new_symbol(SYM_List), args, (size_t)len);
    free(args);
    return r;
}

/* The exact (boxed) product. */
static Expr* kp_boxed(KpOp* a, KpOp* b, const KpLayout* L, const int64_t* oa, const int64_t* ob) {
    if (!kp_box(a) || !kp_box(b)) return NULL;
    unsigned num = KL_INT | KL_REAL | KL_CPLX_M | KL_CPLX_X | KL_NUM;
    unsigned mach = KL_REAL | KL_CPLX_M;
    bool numericize = !(a->kinds & ~num) && !(b->kinds & ~num) && ((a->kinds | b->kinds) & mach);
    Expr** out = calloc((size_t)(L->total ? L->total : 1), sizeof(Expr*));
    if (!out) return NULL;
    for (int64_t i = 0; i < a->n; i++) {
        Expr* x = a->leaves[i];
        for (int64_t j = 0; j < b->n; j++) {
            Expr* y = b->leaves[j];
            int64_t p;
            Expr* v;
            if (x->type == EXPR_INTEGER && y->type == EXPR_INTEGER && !numericize &&
                !ci_mul_i64(x->data.integer, y->data.integer, &p)) {
                v = expr_new_integer(p);
            } else {
                Expr* ta[2] = { expr_copy(x), expr_copy(y) };
                Expr* t = expr_new_function(expr_new_symbol(SYM_Times), ta, 2);
                if (numericize) {
                    Expr* na[1] = { t };
                    t = expr_new_function(expr_new_symbol(SYM_N), na, 1);
                }
                v = evaluate(t);
                expr_free(t);
            }
            out[oa[i] + ob[j]] = v;
        }
    }
    int64_t pos = 0;
    Expr* r = nest(L->dims, L->rank, out, &pos);
    free(out);
    return pack_offer(r);
}

/* KroneckerProduct of two (already validated) operands; NULL on failure. */
static Expr* kp_pair(const Expr* ea, const Expr* eb) {
    KpOp a, b;
    if (kp_read(ea, &a, true) != KP_OK) return NULL;
    if (kp_read(eb, &b, true) != KP_OK) { kpop_free(&a); return NULL; }
    KpLayout L;
    Expr* result = NULL;
    int64_t *oa = NULL, *ob = NULL;
    if (!kp_layout(&a, &b, &L)) goto done;
    oa = kp_offsets(a.rank, a.dims, L.sa, a.n);
    ob = kp_offsets(b.rank, b.dims, L.sb, b.n);
    if (!oa || !ob) goto done;

    /* The machine path: both operands numeric with only machine-convertible
     * leaves, and either both exact integers (int64) or one machine number. */
    unsigned conv = KL_INT | KL_REAL | KL_CPLX_M | KL_CPLX_X;
    unsigned u = a.kinds | b.kinds;
    bool machine = !(a.kinds & ~conv) && !(b.kinds & ~conv) &&
                   (u == KL_INT || (u & (KL_REAL | KL_CPLX_M))) && L.total > 0 &&
                   L.rank < NDARRAY_MAX_RANK;
    if (machine) {
        NDType dt;
        if (u == KL_INT) dt = NDT_INT64;
        else if (a.nd && b.nd) dt = ndt_promote(a.nd->data.ndarray.dtype, b.nd->data.ndarray.dtype);
        else dt = (u & (KL_CPLX_M | KL_CPLX_X)) ? NDT_COMPLEX64 : NDT_FLOAT64;
        if (dt != NDT_INT64 && (u & (KL_CPLX_M | KL_CPLX_X)) && !ndt_is_complex(dt))
            dt = NDT_COMPLEX64;
        const Expr* src = nd_present_src2(a.nd, b.nd);
        void* buf = NULL;
        Expr* nd = src ? ndbuild_open_like(src, L.rank, L.dims, dt, &buf)
                       : ndbuild_open(L.rank, L.dims, dt, &buf);
        void* temp = NULL;
        if (!nd) {
            temp = malloc(ndt_elem_size(dt) * (size_t)L.total);
            buf = temp;
        }
        if (buf && kp_fill(&a, &b, oa, ob, buf, dt)) {
            if (nd) result = nd;
            else {
                Expr* raw = expr_new_ndarray_raw(L.rank, L.dims, temp, dt);
                result = ndarray_to_nested_list(raw);
                expr_free(raw);
            }
            goto done;
        }
        if (nd) expr_free(nd);
        else free(temp);
        /* an int64 overflow: the exact path below */
    }
    result = kp_boxed(&a, &b, &L, oa, ob);
done:
    free(oa);
    free(ob);
    kpop_free(&a);
    kpop_free(&b);
    return result;
}

Expr* builtin_kroneckerproduct(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 2) {
        mth_message("KroneckerProduct", argc == 1 ? "argmu" : "argm",
                    "KroneckerProduct called with %zu argument%s; 2 or more arguments are expected.",
                    argc, argc == 1 ? "" : "s");
        return NULL;
    }
    /* Validate every operand before computing anything: a symbolic operand
     * keeps the call silently; a number or a ragged array is ::rect. */
    bool rect = false;
    for (size_t i = 0; i < argc; i++) {
        KpOp o;
        KpStatus s = kp_read(res->data.function.args[i], &o, false);
        if (s == KP_SILENT) return NULL;
        if (s == KP_RECT) { rect = true; continue; }
        if (argc > 2 && i + 1 < argc && o.rank > 2) return NULL;
    }
    if (rect) {
        /* reported on the List form, once: declining with a buffer in place
         * lets the evaluator's post-gate materialise it and ask again */
        for (size_t i = 0; i < argc; i++)
            if (res->data.function.args[i]->type == EXPR_NDARRAY) return ndarray_delist_and_reeval(res);
        mth_message("KroneckerProduct", "rect", "Nonrectangular array encountered.");
        return NULL;
    }
    Expr* acc = NULL;
    for (size_t i = 1; i < argc; i++) {
        const Expr* left = acc ? acc : res->data.function.args[0];
        Expr* next = kp_pair(left, res->data.function.args[i]);
        if (acc) expr_free(acc);
        acc = next;
        if (!acc) return NULL;
    }
    return acc;
}
