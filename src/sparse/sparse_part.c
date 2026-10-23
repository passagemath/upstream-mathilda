/*
 * sparse_part.c -- Part of a SparseArray, read and assignment.
 *
 * A read selects, per level, either one index (the level disappears) or a list
 * of indices -- All, a Span, or an explicit {i1, i2, ...}, which may repeat and
 * permute. Levels past the spec are All. All indices given: the element, by a
 * binary search in its row. Otherwise a SparseArray of the selection, built by
 * walking only the stored entries (and, when the first level is selected, only
 * the selected rows): a stored entry is kept, at every new position its
 * coordinates map to, and the background stays the background. A selection
 * with an empty dimension is the dense empty array, as for a List.
 *
 * Assignment of one element rebuilds the entry arrays with the entry inserted,
 * replaced, or -- when the new value is the background -- removed. Any other
 * assignment is left to the dense path in part.c.
 */

#include "sparse.h"
#include "sym_names.h"
#include "message.h"
#include "ndarray.h"
#include "part.h"
#include "print.h"
#include <stdlib.h>
#include <string.h>

/* One level's selector. */
typedef struct {
    bool     fixed;          /* a single index (the level is dropped) */
    int64_t  index;          /* 0-based, when fixed */
    bool     all;            /* every index, in order */
    int64_t  len;            /* number of selected positions (when !fixed) */
    int64_t* sel;            /* 0-based old index of each new position (!all) */
    /* inverse: new positions of old index i are inv[inv_ptr[i] .. inv_ptr[i+1]) */
    int64_t* inv_ptr;
    int64_t* inv;
} LevelSel;

static void levelsel_free(LevelSel* s) {
    free(s->sel); free(s->inv_ptr); free(s->inv);
    memset(s, 0, sizeof *s);
}

static bool head_is_sym(const Expr* e, const char* sym) {
    return e && e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL &&
           e->data.function.head->data.symbol.name == sym;
}

/* Parse one subscript against a dimension of length `dim`. Returns 1 ok, 0 not
 * a spec this path handles, -1 out of range. */
static int parse_level(const Expr* ix, int64_t dim, LevelSel* s) {
    memset(s, 0, sizeof *s);
    if (ix->type == EXPR_INTEGER) {
        int64_t k = ix->data.integer;
        if (k < 0) k = dim + k + 1;
        if (k < 1 || k > dim) return -1;
        s->fixed = true;
        s->index = k - 1;
        return 1;
    }
    if (ix->type == EXPR_SYMBOL && ix->data.symbol.name == SYM_All) {
        s->all = true;
        s->len = dim;
        return 1;
    }
    int64_t n;
    if (head_is_sym(ix, SYM_Span)) {
        int64_t start, step, count;
        if (!part_span_resolve(ix, dim, &start, &step, &count)) return 0;
        s->len = count;
        s->sel = malloc(sizeof(int64_t) * (size_t)(count ? count : 1));
        if (!s->sel) return 0;
        for (int64_t i = 0, k = start; i < count; i++, k += step) s->sel[i] = k - 1;
    } else if (head_is_sym(ix, SYM_List) || (is_packed_list(ix) && ix->data.ndarray.rank == 1 &&
                                             ix->data.ndarray.dtype == NDT_INT64)) {
        n = ix->type == EXPR_NDARRAY ? ix->data.ndarray.dims[0] : (int64_t)ix->data.function.arg_count;
        s->len = n;
        s->sel = malloc(sizeof(int64_t) * (size_t)(n ? n : 1));
        if (!s->sel) return 0;
        for (int64_t i = 0; i < n; i++) {
            int64_t k;
            if (ix->type == EXPR_NDARRAY) k = ndt_get_i(ix->data.ndarray.data, (size_t)i, NDT_INT64);
            else if (ix->data.function.args[i]->type == EXPR_INTEGER) k = ix->data.function.args[i]->data.integer;
            else { levelsel_free(s); return 0; }
            if (k < 0) k = dim + k + 1;
            if (k < 1 || k > dim) { levelsel_free(s); return -1; }
            s->sel[i] = k - 1;
        }
    } else {
        return 0;
    }
    /* inverse map, by counting */
    s->inv_ptr = calloc((size_t)(dim + 1), sizeof(int64_t));
    s->inv = malloc(sizeof(int64_t) * (size_t)(s->len ? s->len : 1));
    if (!s->inv_ptr || !s->inv) { levelsel_free(s); return 0; }
    for (int64_t i = 0; i < s->len; i++) s->inv_ptr[s->sel[i] + 1]++;
    for (int64_t i = 0; i < dim; i++) s->inv_ptr[i + 1] += s->inv_ptr[i];
    int64_t* fill = malloc(sizeof(int64_t) * (size_t)(dim ? dim : 1));
    if (!fill) { levelsel_free(s); return 0; }
    memcpy(fill, s->inv_ptr, sizeof(int64_t) * (size_t)dim);
    for (int64_t i = 0; i < s->len; i++) s->inv[fill[s->sel[i]]++] = i;
    free(fill);
    return 1;
}

static void part_msg(const char* tag, const Expr* s, Expr* const* idx, size_t n) {
    Expr** a = malloc(sizeof(Expr*) * (n + 1));
    a[0] = expr_copy((Expr*)s);
    for (size_t i = 0; i < n; i++) a[i + 1] = expr_copy(idx[i]);
    Expr* call = expr_new_function(expr_new_symbol(SYM_Part), a, n + 1);
    free(a);
    char* sc = expr_to_string_outputform(call);
    if (tag[4] == 'w') {
        /* Part::partw names the offending subscript */
        char* si = NULL;
        for (size_t i = 0; i < n && !si; i++) si = expr_to_string_outputform(idx[i]);
        mth_message("Part", "partw", "Part %s of %s does not exist.", si ? si : "", sc);
        free(si);
    } else {
        mth_message("Part", "partd", "Part specification %s is longer than depth of object.", sc);
    }
    free(sc);
    expr_free(call);
}

/* Recursively emit an entry at every new position its coordinates map to. */
static bool emit(SpCoo* c, const SparseArray* sa, int64_t k, const int64_t* coord,
                 const LevelSel* ls, int level, int out_level, int64_t* out) {
    if (level == sa->rank) {
        switch (sa->kind) {
            case SPV_INT64: return sp_coo_push_i64(c, out, ((const int64_t*)sa->vals)[k]);
            case SPV_F64:   return sp_coo_push_f64(c, out, ((const double*)sa->vals)[k]);
            case SPV_C128:  return sp_coo_push_c128(c, out, ((const double*)sa->vals)[2 * k],
                                                    ((const double*)sa->vals)[2 * k + 1]);
            default:        return sp_coo_push_expr(c, out, sp_value_expr(sa, k));
        }
    }
    const LevelSel* s = &ls[level];
    int64_t ci = coord[level];
    if (s->fixed) return emit(c, sa, k, coord, ls, level + 1, out_level, out);
    if (s->all) {
        out[out_level] = ci;
        return emit(c, sa, k, coord, ls, level + 1, out_level + 1, out);
    }
    for (int64_t q = s->inv_ptr[ci]; q < s->inv_ptr[ci + 1]; q++) {
        out[out_level] = s->inv[q];
        if (!emit(c, sa, k, coord, ls, level + 1, out_level + 1, out)) return false;
    }
    return true;
}

/* Does stored entry k (coords) survive the fixed / selected levels? */
static bool entry_selected(const int64_t* coord, const LevelSel* ls, int rank) {
    for (int l = 0; l < rank; l++) {
        const LevelSel* s = &ls[l];
        if (s->fixed) { if (coord[l] != s->index) return false; }
        else if (!s->all && s->inv_ptr[coord[l]] == s->inv_ptr[coord[l] + 1]) return false;
    }
    return true;
}

Expr* sp_part(const Expr* s, Expr* const* idx, size_t n, bool* err) {
    *err = false;
    const SparseArray* sa = sp_of(s);
    if (!sa || n == 0) return NULL;
    int rank = sa->rank;

    /* [[..., 0]]: the head of the part, which for an array is List. */
    for (size_t i = 0; i < n; i++) {
        if (idx[i]->type == EXPR_INTEGER && idx[i]->data.integer == 0) {
            for (size_t j = 0; j < i; j++) if (idx[j]->type != EXPR_INTEGER) return NULL;
            if (i > (size_t)rank) return NULL;
            if (i < (size_t)rank) return expr_new_symbol(SYM_List);
            /* every index given, then 0: the head of the element */
            int64_t c[SP_MAX_RANK];
            for (int l = 0; l < rank; l++) {
                int64_t k = idx[l]->data.integer;
                if (k < 0) k = sa->dims[l] + k + 1;
                if (k < 1 || k > sa->dims[l]) return NULL;
                c[l] = k - 1;
            }
            Expr* el = sp_get(sa, c);
            Expr* h = expr_head(el);
            expr_free(el);
            return h;
        }
    }
    if (n > (size_t)rank) { part_msg("partd", s, idx, n); *err = true; return NULL; }

    LevelSel ls[SP_MAX_RANK];
    memset(ls, 0, sizeof ls);
    int nfixed = 0;
    for (int l = 0; l < rank; l++) {
        int rc = (size_t)l < n ? parse_level(idx[l], sa->dims[l], &ls[l]) : 1;
        if ((size_t)l >= n) { ls[l].all = true; ls[l].len = sa->dims[l]; }
        if (rc <= 0) {
            for (int q = 0; q <= l; q++) levelsel_free(&ls[q]);
            if (rc < 0) { part_msg("partw", s, idx + l, n - (size_t)l); *err = true; }
            return NULL;
        }
        if (ls[l].fixed) nfixed++;
    }

    Expr* result = NULL;
    if (nfixed == rank) {
        int64_t c[SP_MAX_RANK];
        for (int l = 0; l < rank; l++) c[l] = ls[l].index;
        result = sp_get(sa, c);
        goto done;
    }

    int orank = rank - nfixed;
    int64_t odims[SP_MAX_RANK];
    for (int l = 0, o = 0; l < rank; l++) if (!ls[l].fixed) odims[o++] = ls[l].len;
    SpCoo c;
    sp_coo_init(&c, orank, odims, sa->kind);
    const SparsePattern* p = sa->pat;
    int64_t coord[SP_MAX_RANK], out[SP_MAX_RANK];
    int first = rank >= 2 ? 1 : 0;
    bool ok = true;
    for (int64_t r = 0; r < p->nrows && ok; r++) {
        if (first) {
            if (ls[0].fixed && r != ls[0].index) continue;
            if (!ls[0].fixed && !ls[0].all && ls[0].inv_ptr[r] == ls[0].inv_ptr[r + 1]) continue;
            coord[0] = r;
        }
        for (int64_t k = p->rowptr[r]; k < p->rowptr[r + 1] && ok; k++) {
            for (int j = 0; j < p->width; j++) coord[first + j] = sp_col(p, k, j);
            if (!entry_selected(coord, ls, rank)) continue;
            ok = emit(&c, sa, k, coord, ls, 0, 0, out);
        }
    }
    if (!ok) { sp_coo_free(&c); goto done; }
    if (sp_total(orank, odims) == 0) {
        /* an empty selection is the dense empty array */
        SparseArray* e = sp_coo_finish(&c, expr_copy(sa->bg), false);
        if (e) { result = sp_normal(e); sp_release(e); }
        goto done;
    }
    result = sp_to_expr(sp_coo_finish(&c, expr_copy(sa->bg), false));
done:
    for (int l = 0; l < rank; l++) levelsel_free(&ls[l]);
    return result;
}

/* ---- assignment ---------------------------------------------------------- */

static bool is_list_value(const Expr* v) {
    return head_is_sym(v, SYM_List) || v->type == EXPR_NDARRAY || v->type == EXPR_SPARSE;
}

Expr* sp_part_assign(const Expr* s, Expr* const* idx, size_t n, const Expr* rhs, bool* err) {
    *err = false;
    const SparseArray* sa = sp_of(s);
    if (!sa || n != (size_t)sa->rank || is_list_value(rhs)) return NULL;
    int rank = sa->rank;
    int64_t c[SP_MAX_RANK];
    for (int l = 0; l < rank; l++) {
        if (idx[l]->type != EXPR_INTEGER) return NULL;
        int64_t k = idx[l]->data.integer;
        if (k < 0) k = sa->dims[l] + k + 1;
        if (k < 1 || k > sa->dims[l]) { part_msg("partw", s, idx + l, n - (size_t)l); *err = true; return NULL; }
        c[l] = k - 1;
    }
    /* Rebuild through the COO builder: every stored entry, then the new one
     * placed FIRST so that it wins over the old value at that position. */
    const SparsePattern* p = sa->pat;
    bool is_bg = sp_value_is_background(rhs, sa->bg);
    /* Stay typed when the new value is (or removes to) the buffer's type. */
    SpValKind kind = SPV_EXPR;
    if (sa->kind == SPV_F64 && (is_bg || rhs->type == EXPR_REAL)) kind = SPV_F64;
    if (sa->kind == SPV_INT64 && (is_bg || rhs->type == EXPR_INTEGER)) kind = SPV_INT64;
    SpCoo coo;
    sp_coo_init(&coo, rank, sa->dims, kind);
    if (!sp_coo_reserve(&coo, p->nnz + 1)) { sp_coo_free(&coo); return NULL; }
    int first = rank >= 2 ? 1 : 0;
    int64_t coord[SP_MAX_RANK];
    bool ok = true;
    if (!is_bg) {
        if (kind == SPV_F64) ok = sp_coo_push_f64(&coo, c, rhs->data.real);
        else if (kind == SPV_INT64) ok = sp_coo_push_i64(&coo, c, rhs->data.integer);
        else ok = sp_coo_push_expr(&coo, c, expr_copy((Expr*)rhs));
    }
    for (int64_t r = 0; r < p->nrows && ok; r++) {
        if (first) coord[0] = r;
        for (int64_t k = p->rowptr[r]; k < p->rowptr[r + 1] && ok; k++) {
            for (int j = 0; j < p->width; j++) coord[first + j] = sp_col(p, k, j);
            if (is_bg && memcmp(coord, c, sizeof(int64_t) * (size_t)rank) == 0) continue;
            if (kind == SPV_F64) ok = sp_coo_push_f64(&coo, coord, ((const double*)sa->vals)[k]);
            else if (kind == SPV_INT64) ok = sp_coo_push_i64(&coo, coord, ((const int64_t*)sa->vals)[k]);
            else ok = sp_coo_push_expr(&coo, coord, sp_value_expr(sa, k));
        }
    }
    if (!ok) { sp_coo_free(&coo); return NULL; }
    /* Existing explicit entries keep their explicitness (no background drop);
     * only the assigned position is decided by the new value. */
    return sp_to_expr(sp_coo_finish(&coo, expr_copy(sa->bg), false));
}
