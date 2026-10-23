/*
 * sparse_arith.c -- Plus, Times and Power over machine-typed SparseArrays.
 *
 * Each routine reproduces exactly what the scalar arithmetic does element by
 * element (Integer op Integer exact and overflow-checked, anything with a Real
 * in it a Real, machine complex with a zero imaginary part collapsing to a
 * Real), so it is a faster spelling of sparse_thread.c's generic tier and
 * never a different answer. Whatever it cannot type exactly -- an overflow, an
 * exact zero factor, an exact non-integer background, a result that would mix
 * Integers with Reals -- it declines with NULL, and the generic tier answers.
 *
 * The background is evaluated through the ordinary Plus/Times AFTER the value
 * loop has succeeded, so a declined fast path never emits its messages twice.
 */

#include "sparse_internal.h"
#include "sym_names.h"
#include "eval.h"
#include "checked_int.h"
#include <stdlib.h>
#include <string.h>

/* ---- scalars ------------------------------------------------------------- */

enum { NS_INT, NS_REAL, NS_CPLX };

/* A machine scalar: Integer, Real or Complex[Real, Real]. */
static bool machine_scalar(const Expr* c, int* t, int64_t* iv, double* re, double* im) {
    *im = 0;
    if (c->type == EXPR_INTEGER) { *t = NS_INT; *iv = c->data.integer; *re = (double)*iv; return true; }
    if (c->type == EXPR_REAL) { *t = NS_REAL; *re = c->data.real; return true; }
    if (c->type == EXPR_FUNCTION && c->data.function.arg_count == 2 &&
        c->data.function.head->type == EXPR_SYMBOL &&
        c->data.function.head->data.symbol.name == SYM_Complex &&
        c->data.function.args[0]->type == EXPR_REAL && c->data.function.args[1]->type == EXPR_REAL) {
        *t = NS_CPLX;
        *re = c->data.function.args[0]->data.real;
        *im = c->data.function.args[1]->data.real;
        return true;
    }
    return false;
}

static Expr* eval2(const char* head, const Expr* a, const Expr* b) {
    Expr* args[2] = { expr_copy((Expr*)a), expr_copy((Expr*)b) };
    Expr* e = expr_new_function(expr_new_symbol(head), args, 2);
    Expr* r = evaluate(e);
    expr_free(e);
    return r;
}

static bool is_list_like(const Expr* v) {
    return v->type == EXPR_SPARSE || v->type == EXPR_NDARRAY ||
           (v->type == EXPR_FUNCTION && v->data.function.head->type == EXPR_SYMBOL &&
            v->data.function.head->data.symbol.name == SYM_List);
}

/* Assemble a result over `pat` (one reference adopted) with background
 * head[bga, bgb]. */
static Expr* finish(const SparseArray* like, const char* head, const Expr* bga, const Expr* bgb,
                    SparsePattern* pat, SpValKind kind, void* vals) {
    sp_c128_canonicalise(&kind, &vals, pat->nnz);
    Expr* bg = eval2(head, bga, bgb);
    if (!bg || is_list_like(bg)) {
        if (bg) expr_free(bg);
        if (kind == SPV_EXPR) for (int64_t k = 0; k < pat->nnz; k++) expr_free(((Expr**)vals)[k]);
        free(vals);
        sp_pattern_release(pat);
        return NULL;
    }
    return sp_to_expr(sp_new(like->rank, like->dims, bg, pat, kind, vals));
}

/* ======================================================================== */
/* s op c                                                                   */
/* ======================================================================== */

Expr* sp_arith_scalar(bool is_plus, const SparseArray* sa, const Expr* c) {
    int t; int64_t ci; double cr, cm;
    if (sa->kind == SPV_EXPR || !machine_scalar(c, &t, &ci, &cr, &cm)) return NULL;
    /* An exact zero factor makes every product exact (the scalar rule), which
     * no typed buffer can say; the generic tier answers. */
    if (!is_plus && t == NS_INT && ci == 0) return NULL;
    const char* head = is_plus ? SYM_Plus : SYM_Times;
    int64_t n = sa->pat->nnz;
    SparsePattern* pat = sp_pattern_retain(sa->pat);

    /* exact identities: the values are unchanged */
    if (t == NS_INT && ((is_plus && ci == 0) || (!is_plus && ci == 1))) {
        size_t sz = sp_kind_size(sa->kind) * (size_t)n;
        void* v = n ? malloc(sz) : NULL;
        if (n && !v) { sp_pattern_release(pat); return NULL; }
        if (n) memcpy(v, sa->vals, sz);
        return finish(sa, head, sa->bg, c, pat, sa->kind, v);
    }

    if (sa->kind == SPV_INT64 && t == NS_INT) {
        int64_t* out = malloc(sizeof(int64_t) * (size_t)(n ? n : 1));
        if (!out) { sp_pattern_release(pat); return NULL; }
        const int64_t* in = sa->vals;
        bool ok = true;
        if (is_plus) for (int64_t k = 0; k < n && ok; k++) ok = !ci_add_i64(in[k], ci, &out[k]);
        else         for (int64_t k = 0; k < n && ok; k++) ok = !ci_mul_i64(in[k], ci, &out[k]);
        if (!ok) { free(out); sp_pattern_release(pat); return NULL; }
        return finish(sa, head, sa->bg, c, pat, SPV_INT64, out);
    }

    bool cplx = t == NS_CPLX || sa->kind == SPV_C128;
    if (!cplx) {
        double* out = malloc(sizeof(double) * (size_t)(n ? n : 1));
        if (!out) { sp_pattern_release(pat); return NULL; }
        if (sa->kind == SPV_INT64) {
            const int64_t* in = sa->vals;
            if (is_plus) for (int64_t k = 0; k < n; k++) out[k] = (double)in[k] + cr;
            else         for (int64_t k = 0; k < n; k++) out[k] = (double)in[k] * cr;
        } else {
            const double* in = sa->vals;
            if (is_plus) for (int64_t k = 0; k < n; k++) out[k] = in[k] + cr;
            else         for (int64_t k = 0; k < n; k++) out[k] = in[k] * cr;
        }
        return finish(sa, head, sa->bg, c, pat, SPV_F64, out);
    }

    double* out = malloc(2 * sizeof(double) * (size_t)(n ? n : 1));
    if (!out) { sp_pattern_release(pat); return NULL; }
    for (int64_t k = 0; k < n; k++) {
        double a, b;
        if (sa->kind == SPV_C128) { a = ((const double*)sa->vals)[2 * k]; b = ((const double*)sa->vals)[2 * k + 1]; }
        else if (sa->kind == SPV_F64) { a = ((const double*)sa->vals)[k]; b = 0; }
        else { a = (double)((const int64_t*)sa->vals)[k]; b = 0; }
        if (is_plus) { out[2 * k] = a + cr; out[2 * k + 1] = b + cm; }
        else if (sa->kind != SPV_C128) { out[2 * k] = a * cr; out[2 * k + 1] = a * cm; }
        else if (t != NS_CPLX) { out[2 * k] = a * cr; out[2 * k + 1] = b * cr; }
        else { out[2 * k] = a * cr - b * cm; out[2 * k + 1] = a * cm + b * cr; }
    }
    return finish(sa, head, sa->bg, c, pat, SPV_C128, out);
}

/* ======================================================================== */
/* s1 op s2                                                                 */
/* ======================================================================== */

/* Operand values as the merge loops read them. */
typedef struct {
    const SparseArray* sa;
    bool  is_int;           /* SPV_INT64 */
    int64_t bgi;            /* background as an Integer (when is_int and exact) */
    double bgd;             /* background as a double */
} Opnd;

static inline double val_d(const Opnd* o, int64_t k) {
    return o->is_int ? (double)((const int64_t*)o->sa->vals)[k] : ((const double*)o->sa->vals)[k];
}

/* Two-pointer merge of one row's tuples: entries ka in [a0, a1), kb in
 * [b0, b1). Emits (ka or -1, kb or -1) in canonical order; with `inter` only
 * the positions both have. */
static inline int tuple_cmp2(const SparsePattern* pa, int64_t ka, const SparsePattern* pb, int64_t kb) {
    if (pa->width == 1) {
        int64_t x = sp_col(pa, ka, 0), y = sp_col(pb, kb, 0);
        return x < y ? -1 : (x > y ? 1 : 0);
    }
    for (int j = 0; j < pa->width; j++) {
        int64_t x = sp_col(pa, ka, j), y = sp_col(pb, kb, j);
        if (x != y) return x < y ? -1 : 1;
    }
    return 0;
}

static inline void copy_tuple(SparsePattern* out, int64_t ko, const SparsePattern* in, int64_t ki) {
    int w = in->width;
    if (out->idx32) {
        const int32_t* s = (const int32_t*)in->colidx + ki * w;
        int32_t* d = (int32_t*)out->colidx + ko * w;
        for (int j = 0; j < w; j++) d[j] = s[j];
    } else {
        const int64_t* s = (const int64_t*)in->colidx + ki * w;
        int64_t* d = (int64_t*)out->colidx + ko * w;
        for (int j = 0; j < w; j++) d[j] = s[j];
    }
}

Expr* sp_arith_binary(bool is_plus, const SparseArray* a, const SparseArray* b) {
    if (a->kind != SPV_INT64 && a->kind != SPV_F64) return NULL;
    if (b->kind != SPV_INT64 && b->kind != SPV_F64) return NULL;
    const Expr* ga = a->bg;
    const Expr* gb = b->bg;
    if ((ga->type != EXPR_INTEGER && ga->type != EXPR_REAL) ||
        (gb->type != EXPR_INTEGER && gb->type != EXPR_REAL)) return NULL;
    bool ai = a->kind == SPV_INT64, bi = b->kind == SPV_INT64;
    bool int_result = ai && bi && ga->type == EXPR_INTEGER && gb->type == EXPR_INTEGER;
    if (!int_result) {
        /* Every element must come out Real: an Integer entry may only meet the
         * other operand's REAL background or a Real entry. */
        bool ok = (!ai && !bi) || (ai && !bi && gb->type == EXPR_REAL) ||
                  (!ai && bi && ga->type == EXPR_REAL);
        if (!ok) return NULL;
    }
    Opnd A = { a, ai, ga->type == EXPR_INTEGER ? ga->data.integer : 0,
               ga->type == EXPR_INTEGER ? (double)ga->data.integer : ga->data.real };
    Opnd B = { b, bi, gb->type == EXPR_INTEGER ? gb->data.integer : 0,
               gb->type == EXPR_INTEGER ? (double)gb->data.integer : gb->data.real };
    const char* head = is_plus ? SYM_Plus : SYM_Times;
    /* Times over two zero backgrounds keeps only the common positions (the
     * product is the zero background everywhere else). */
    bool inter = !is_plus && sp_is_zero(ga) && sp_is_zero(gb);
    bool need_a = inter, need_b = inter;
    const SparsePattern* pa = a->pat;
    const SparsePattern* pb = b->pat;

    if (sp_pattern_same(pa, pb)) {
        int64_t n = pa->nnz;
        SparsePattern* pat = sp_pattern_retain((SparsePattern*)pa);
        if (int_result) {
            int64_t* out = malloc(sizeof(int64_t) * (size_t)(n ? n : 1));
            if (!out) { sp_pattern_release(pat); return NULL; }
            const int64_t* x = a->vals; const int64_t* y = b->vals;
            bool ok = true;
            if (is_plus) for (int64_t k = 0; k < n && ok; k++) ok = !ci_add_i64(x[k], y[k], &out[k]);
            else         for (int64_t k = 0; k < n && ok; k++) ok = !ci_mul_i64(x[k], y[k], &out[k]);
            if (!ok) { free(out); sp_pattern_release(pat); return NULL; }
            return finish(a, head, ga, gb, pat, SPV_INT64, out);
        }
        double* out = malloc(sizeof(double) * (size_t)(n ? n : 1));
        if (!out) { sp_pattern_release(pat); return NULL; }
        if (!ai && !bi) {
            const double* x = a->vals; const double* y = b->vals;
            if (is_plus) for (int64_t k = 0; k < n; k++) out[k] = x[k] + y[k];
            else         for (int64_t k = 0; k < n; k++) out[k] = x[k] * y[k];
        } else {
            if (is_plus) for (int64_t k = 0; k < n; k++) out[k] = val_d(&A, k) + val_d(&B, k);
            else         for (int64_t k = 0; k < n; k++) out[k] = val_d(&A, k) * val_d(&B, k);
        }
        return finish(a, head, ga, gb, pat, SPV_F64, out);
    }

    /* Different patterns: one fused merge pass per row. */
    int64_t cap;
    if (need_a && need_b) cap = pa->nnz < pb->nnz ? pa->nnz : pb->nnz;
    else if (need_a) cap = pa->nnz;
    else if (need_b) cap = pb->nnz;
    else cap = pa->nnz + pb->nnz;
    SparsePattern* pat = sp_pattern_new(pa->nrows, cap, pa->width, pa->idx32);
    void* out = malloc((int_result ? sizeof(int64_t) : sizeof(double)) * (size_t)(cap ? cap : 1));
    if (!pat || !out) { sp_pattern_release(pat); free(out); return NULL; }
    int64_t* oi = out;
    double* od = out;
    int64_t o = 0;
    bool ok = true;
    pat->rowptr[0] = 0;
    bool fast = pa->width == 1 && pa->idx32 && !int_result && !ai && !bi;
    if (fast) {
        /* The common case -- a matrix or vector of Reals with 32-bit column
         * indices -- as a branch-light merge on raw pointers. */
        const int32_t* ca = pa->colidx; const int32_t* cb = pb->colidx;
        int32_t* co = pat->colidx;
        const double* va = a->vals; const double* vb = b->vals;
        double ga_d = A.bgd, gb_d = B.bgd;
        for (int64_t r = 0; r < pa->nrows; r++) {
            int64_t ka = pa->rowptr[r], ea = pa->rowptr[r + 1];
            int64_t kb = pb->rowptr[r], eb = pb->rowptr[r + 1];
            if (inter) {
                while (ka < ea && kb < eb) {
                    int32_t x = ca[ka], y = cb[kb];
                    if (x == y) { co[o] = x; od[o++] = va[ka++] * vb[kb++]; }
                    else if (x < y) ka++;
                    else kb++;
                }
            } else if (is_plus) {
                while (ka < ea && kb < eb) {
                    int32_t x = ca[ka], y = cb[kb];
                    if (x < y)      { co[o] = x; od[o++] = va[ka++] + gb_d; }
                    else if (y < x) { co[o] = y; od[o++] = ga_d + vb[kb++]; }
                    else            { co[o] = x; od[o++] = va[ka++] + vb[kb++]; }
                }
                while (ka < ea) { co[o] = ca[ka]; od[o++] = va[ka++] + gb_d; }
                while (kb < eb) { co[o] = cb[kb]; od[o++] = ga_d + vb[kb++]; }
            } else {
                while (ka < ea && kb < eb) {
                    int32_t x = ca[ka], y = cb[kb];
                    if (x < y)      { co[o] = x; od[o++] = va[ka++] * gb_d; }
                    else if (y < x) { co[o] = y; od[o++] = ga_d * vb[kb++]; }
                    else            { co[o] = x; od[o++] = va[ka++] * vb[kb++]; }
                }
                while (ka < ea) { co[o] = ca[ka]; od[o++] = va[ka++] * gb_d; }
                while (kb < eb) { co[o] = cb[kb]; od[o++] = ga_d * vb[kb++]; }
            }
            pat->rowptr[r + 1] = o;
        }
    }
    for (int64_t r = 0; !fast && r < pa->nrows && ok; r++) {
        int64_t ka = pa->rowptr[r], ea = pa->rowptr[r + 1];
        int64_t kb = pb->rowptr[r], eb = pb->rowptr[r + 1];
        while (ka < ea || kb < eb) {
            int c;
            if (ka >= ea) c = 1;
            else if (kb >= eb) c = -1;
            else c = tuple_cmp2(pa, ka, pb, kb);
            int64_t xa = c <= 0 ? ka : -1, xb = c >= 0 ? kb : -1;
            if (c <= 0) ka++;
            if (c >= 0) kb++;
            if ((need_a && xa < 0) || (need_b && xb < 0)) continue;
            copy_tuple(pat, o, xa >= 0 ? pa : pb, xa >= 0 ? xa : xb);
            if (int_result) {
                int64_t x = xa >= 0 ? ((const int64_t*)a->vals)[xa] : A.bgi;
                int64_t y = xb >= 0 ? ((const int64_t*)b->vals)[xb] : B.bgi;
                ok = !(is_plus ? ci_add_i64(x, y, &oi[o]) : ci_mul_i64(x, y, &oi[o]));
                if (!ok) break;
            } else {
                double x = xa >= 0 ? val_d(&A, xa) : A.bgd;
                double y = xb >= 0 ? val_d(&B, xb) : B.bgd;
                od[o] = is_plus ? x + y : x * y;
            }
            o++;
        }
        pat->rowptr[r + 1] = o;
    }
    if (!ok) { sp_pattern_release(pat); free(out); return NULL; }
    pat->nnz = o;
    if (o < cap) {
        size_t isz = pat->idx32 ? sizeof(int32_t) : sizeof(int64_t);
        if (o == 0) { free(pat->colidx); pat->colidx = NULL; }
        else {
            void* ci = realloc(pat->colidx, isz * (size_t)o * (size_t)pat->width);
            if (ci) pat->colidx = ci;
            void* vv = realloc(out, (int_result ? sizeof(int64_t) : sizeof(double)) * (size_t)o);
            if (vv) out = vv;
        }
    }
    return finish(a, head, ga, gb, pat, int_result ? SPV_INT64 : SPV_F64, out);
}

/* ======================================================================== */
/* Power                                                                    */
/* ======================================================================== */

Expr* sp_power_scalar(const SparseArray* sa, const Expr* c, bool sparse_is_base) {
    if (sa->kind == SPV_EXPR) return NULL;
    if (c->type != EXPR_INTEGER && c->type != EXPR_REAL) return NULL;
    double cv = c->type == EXPR_INTEGER ? (double)c->data.integer : c->data.real;
    /* An Integer exponent past 2^53 does not survive the double. */
    if (c->type == EXPR_INTEGER && (c->data.integer > ((int64_t)1 << 53) ||
                                    c->data.integer < -((int64_t)1 << 53))) return NULL;
    Expr* view = sp_view(sa->kind, sa->vals, sa->pat->nnz);
    if (!view) return NULL;
    Expr* nd = sparse_is_base ? ndarray_scalar_power(view, cv, 0.0)
                              : ndarray_base_scalar_power(cv, 0.0, view);
    sp_view_free(view);
    SpValKind kind;
    void* vals;
    if (!sp_adopt_ndarray(nd, &kind, &vals)) return NULL;
    /* the caller replaces the background with Power[bg, c] (or Power[c, bg]) */
    return sp_to_expr(sp_new(sa->rank, sa->dims, expr_copy(sa->bg), sp_pattern_retain(sa->pat),
                             kind, vals));
}
