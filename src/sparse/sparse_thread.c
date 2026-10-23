/*
 * sparse_thread.c -- Listable functions over SparseArrays.
 *
 * f[s1, ..., c1, ...] with SparseArray arguments (and no List arguments, which
 * densify first -- a sparse plus a dense array is a dense array) threads over
 * the elements as it would over Normal[s]:
 *
 *   background  f[bg1, ..., c1, ...]   evaluated once
 *   entries     f[v1, ..., c1, ...]    at the union of the explicit positions,
 *                                      with an operand's background standing in
 *                                      where it has no entry
 *
 * Times keeps only the positions where every zero-background operand has an
 * entry: elsewhere the product is that zero, which is the result's background.
 * Nothing is re-minimised -- an entry that evaluates to the background stays
 * explicit -- so the result's pattern is a function of the operands' patterns
 * alone, and a single-operand map SHARES its operand's pattern.
 *
 * Three tiers, fastest first:
 *   1. Plus / Times / Power on typed (machine) values: C loops (sparse_arith.c)
 *   2. a head with a registered ND kernel on typed values: the SAME kernel the
 *      packed / NDArray surface uses, run over the value buffer as a vector
 *   3. everything else: evaluate f per entry (symbolic values, user functions)
 * A faster tier that declines (overflow, a domain escape it cannot represent)
 * hands the call to the next, so all three give the same answer.
 */

#include "sparse_internal.h"
#include "sym_names.h"
#include "symtab.h"
#include "eval.h"
#include "message.h"
#include "print.h"
#include "attr.h"
#include <stdlib.h>
#include <string.h>

/* ======================================================================== */
/* Value-buffer views                                                       */
/* ======================================================================== */

static bool kind_dtype(SpValKind kind, NDType* dt) {
    switch (kind) {
        case SPV_INT64: *dt = NDT_INT64; return true;
        case SPV_F64:   *dt = NDT_FLOAT64; return true;
        case SPV_C128:  *dt = NDT_COMPLEX64; return true;
        default:        return false;
    }
}

Expr* sp_view(SpValKind kind, void* vals, int64_t n) {
    NDType dt;
    if (!kind_dtype(kind, &dt)) return NULL;
    return expr_new_ndarray_raw(1, &n, vals, dt);
}

void sp_view_free(Expr* view) {
    if (!view) return;
    view->data.ndarray.data = NULL;   /* borrowed: not ours to free */
    expr_free(view);
}

bool sp_adopt_ndarray(Expr* nd, SpValKind* kind, void** vals) {
    if (!nd) return false;
    SpValKind k;
    switch (nd->data.ndarray.dtype) {
        case NDT_INT64:     k = SPV_INT64; break;
        case NDT_FLOAT64:   k = SPV_F64; break;
        case NDT_COMPLEX64: k = SPV_C128; break;
        default: expr_free(nd); return false;
    }
    if (nd->refcount != 1 || nd->data.ndarray.rank != 1) { expr_free(nd); return false; }
    int64_t n = nd->data.ndarray.dims[0];
    *kind = k;
    *vals = nd->data.ndarray.data;
    nd->data.ndarray.data = NULL;
    expr_free(nd);
    sp_c128_canonicalise(kind, vals, n);
    return true;
}

void sp_c128_canonicalise(SpValKind* kind, void** vals, int64_t n) {
    if (*kind != SPV_C128) return;
    const double* b = *vals;
    int64_t k = 0;
    while (k < n && b[2 * k + 1] != 0.0) k++;
    if (k == n) return;
    Expr** ev = malloc(sizeof(Expr*) * (size_t)n);
    if (!ev) return;
    for (int64_t i = 0; i < n; i++) {
        if (b[2 * i + 1] == 0.0) {
            ev[i] = expr_new_real(b[2 * i]);
        } else {
            Expr* args[2] = { expr_new_real(b[2 * i]), expr_new_real(b[2 * i + 1]) };
            ev[i] = expr_new_function(expr_new_symbol(SYM_Complex), args, 2);
        }
    }
    free(*vals);
    *vals = ev;
    *kind = SPV_EXPR;
}

bool sp_is_zero(const Expr* e) {
    if (e->type == EXPR_INTEGER) return e->data.integer == 0;
    if (e->type == EXPR_REAL) return e->data.real == 0.0;
    if (e->type == EXPR_FUNCTION && e->data.function.arg_count == 2 &&
        e->data.function.head->type == EXPR_SYMBOL &&
        e->data.function.head->data.symbol.name == SYM_Complex) {
        const Expr* a = e->data.function.args[0];
        const Expr* b = e->data.function.args[1];
        return a->type == EXPR_REAL && b->type == EXPR_REAL &&
               a->data.real == 0.0 && b->data.real == 0.0;
    }
    return false;
}

/* ======================================================================== */
/* Union of patterns                                                        */
/* ======================================================================== */

static int entry_cmp(const SparsePattern* pa, int64_t ka, const SparsePattern* pb, int64_t kb) {
    for (int j = 0; j < pa->width; j++) {
        int64_t x = sp_col(pa, ka, j), y = sp_col(pb, kb, j);
        if (x != y) return x < y ? -1 : 1;
    }
    return 0;
}

bool sp_union(const SparseArray* const* ops, int m, const bool* need_all, SpUnion* u) {
    memset(u, 0, sizeof *u);
    u->m = m;
    const SparsePattern* p0 = ops[0]->pat;
    bool same = true;
    for (int j = 1; j < m && same; j++) same = sp_pattern_same(p0, ops[j]->pat);
    if (same) {
        int64_t n = p0->nnz;
        u->pat = sp_pattern_retain((SparsePattern*)p0);
        u->src = malloc(sizeof(int64_t) * (size_t)(n * m > 0 ? n * m : 1));
        if (!u->src) { sp_pattern_release(u->pat); u->pat = NULL; return false; }
        for (int64_t k = 0; k < n; k++) for (int j = 0; j < m; j++) u->src[k * m + j] = k;
        return true;
    }
    int64_t cap = 0;
    bool any_need = false;
    for (int j = 0; j < m; j++) {
        if (need_all && need_all[j]) {
            if (!any_need || ops[j]->pat->nnz < cap) cap = ops[j]->pat->nnz;
            any_need = true;
        }
    }
    if (!any_need) for (int j = 0; j < m; j++) cap += ops[j]->pat->nnz;
    int64_t nrows = p0->nrows;
    int width = p0->width;
    SparsePattern* pat = sp_pattern_new(nrows, cap, width, p0->idx32);
    int64_t* src = malloc(sizeof(int64_t) * (size_t)(cap * m > 0 ? cap * m : 1));
    int64_t* cur = malloc(sizeof(int64_t) * (size_t)m);
    bool* here = malloc(sizeof(bool) * (size_t)m);
    if (!pat || !src || !cur || !here) {
        sp_pattern_release(pat); free(src); free(cur); free(here);
        return false;
    }
    int64_t out = 0;
    pat->rowptr[0] = 0;
    for (int64_t r = 0; r < nrows; r++) {
        for (int j = 0; j < m; j++) cur[j] = ops[j]->pat->rowptr[r];
        for (;;) {
            int lead = -1;
            for (int j = 0; j < m; j++) {
                if (cur[j] >= ops[j]->pat->rowptr[r + 1]) continue;
                if (lead < 0 || entry_cmp(ops[j]->pat, cur[j], ops[lead]->pat, cur[lead]) < 0) lead = j;
            }
            if (lead < 0) break;
            bool keep = true;
            for (int j = 0; j < m; j++) {
                here[j] = cur[j] < ops[j]->pat->rowptr[r + 1] &&
                          entry_cmp(ops[j]->pat, cur[j], ops[lead]->pat, cur[lead]) == 0;
                if (need_all && need_all[j] && !here[j]) keep = false;
            }
            if (keep) {
                for (int q = 0; q < width; q++) {
                    int64_t v = sp_col(ops[lead]->pat, cur[lead], q);
                    if (pat->idx32) ((int32_t*)pat->colidx)[out * width + q] = (int32_t)v;
                    else            ((int64_t*)pat->colidx)[out * width + q] = v;
                }
                for (int j = 0; j < m; j++) src[out * m + j] = here[j] ? cur[j] : -1;
                out++;
            }
            for (int j = m - 1; j >= 0; j--) if (here[j]) cur[j]++;
        }
        pat->rowptr[r + 1] = out;
    }
    free(cur); free(here);
    pat->nnz = out;
    if (out == 0) { free(pat->colidx); pat->colidx = NULL; }
    u->pat = pat;
    u->src = src;
    return true;
}

void sp_union_free(SpUnion* u) {
    sp_pattern_release(u->pat);
    free(u->src);
    memset(u, 0, sizeof *u);
}

/* ======================================================================== */
/* The generic (boxed) tier                                                 */
/* ======================================================================== */

/* Evaluate head[args] built from the call's arguments, with argument i of the
 * call replaced by repl[i] when that is non-NULL. */
static Expr* eval_with(const Expr* call, Expr* const* repl) {
    size_t n = call->data.function.arg_count;
    Expr** a = malloc(sizeof(Expr*) * (n ? n : 1));
    for (size_t i = 0; i < n; i++)
        a[i] = repl[i] ? repl[i] : expr_copy(call->data.function.args[i]);
    Expr* e = expr_new_function(expr_copy(call->data.function.head), a, n);
    free(a);
    Expr* r = evaluate(e);
    expr_free(e);
    return r;
}

static bool is_list_value(const Expr* v) {
    return (v->type == EXPR_FUNCTION && v->data.function.head->type == EXPR_SYMBOL &&
            v->data.function.head->data.symbol.name == SYM_List) || v->type == EXPR_NDARRAY ||
           v->type == EXPR_SPARSE;
}

static Expr* thread_generic(const Expr* call, const int* which, int m, bool is_times) {
    size_t n = call->data.function.arg_count;
    const SparseArray* ops[64];
    const SparseArray** opsv = m <= 64 ? ops : malloc(sizeof(SparseArray*) * (size_t)m);
    bool* need = malloc(sizeof(bool) * (size_t)m);
    /* Times over all-zero backgrounds: the product is the background wherever
     * any operand has none, so only the common positions are explicit. */
    bool all_zero = is_times;
    for (int j = 0; j < m; j++) {
        opsv[j] = call->data.function.args[which[j]]->data.sparse;
        if (!sp_is_zero(opsv[j]->bg)) all_zero = false;
    }
    for (int j = 0; j < m; j++) need[j] = all_zero;
    Expr* result = NULL;
    Expr** repl = calloc(n ? n : 1, sizeof(Expr*));
    SpUnion u;
    if (!repl || !sp_union(opsv, m, all_zero ? need : NULL, &u)) goto done;

    /* background */
    for (int j = 0; j < m; j++) repl[which[j]] = expr_copy(opsv[j]->bg);
    Expr* bg = eval_with(call, repl);
    for (int j = 0; j < m; j++) repl[which[j]] = NULL;
    if (!bg || is_list_value(bg)) { if (bg) expr_free(bg); sp_union_free(&u); goto done; }

    int64_t nnz = u.pat->nnz;
    Expr** vals = malloc(sizeof(Expr*) * (size_t)(nnz ? nnz : 1));
    if (!vals) { expr_free(bg); sp_union_free(&u); goto done; }
    for (int64_t k = 0; k < nnz; k++) {
        for (int j = 0; j < m; j++) {
            int64_t s = u.src[k * m + j];
            repl[which[j]] = s >= 0 ? sp_value_expr(opsv[j], s) : expr_copy(opsv[j]->bg);
        }
        vals[k] = eval_with(call, repl);
        for (int j = 0; j < m; j++) repl[which[j]] = NULL;
        if (!vals[k] || is_list_value(vals[k])) {
            /* An element that is itself a List cannot live in a SparseArray. */
            for (int64_t q = 0; q <= k; q++) if (vals[q]) expr_free(vals[q]);
            free(vals); expr_free(bg); sp_union_free(&u);
            goto done;
        }
    }
    const SparseArray* a0 = opsv[0];
    result = sp_to_expr(sp_new(a0->rank, a0->dims, bg, sp_pattern_retain(u.pat), SPV_EXPR, vals));
    sp_union_free(&u);
done:
    free(repl);
    free(need);
    if (opsv != ops) free(opsv);
    return result;
}

/* ======================================================================== */
/* The kernel tier                                                          */
/* ======================================================================== */

/* f[s] for a head with an NDUnaryKernel: run it over the value vector. */
static Expr* thread_unary_kernel(const Expr* call, const SparseArray* sa, const void* kernel) {
    if (sa->kind == SPV_EXPR) return NULL;
    int64_t n = sa->pat->nnz;
    Expr* view = sp_view(sa->kind, sa->vals, n);
    if (!view) return NULL;
    Expr* nd = ndarray_map_unary(view, (const NDUnaryKernel*)kernel);
    sp_view_free(view);
    SpValKind kind;
    void* vals;
    if (!sp_adopt_ndarray(nd, &kind, &vals)) return NULL;
    Expr* repl[1] = { expr_copy(sa->bg) };
    Expr* bg = eval_with(call, repl);
    if (!bg || is_list_value(bg)) {
        if (bg) expr_free(bg);
        free(vals);
        return NULL;
    }
    return sp_to_expr(sp_new(sa->rank, sa->dims, bg, sp_pattern_retain(sa->pat), kind, vals));
}

/* f[s, c] / f[c, s] for a head with an NDBinaryKernel and a machine scalar. */
static Expr* thread_binary_kernel(const Expr* call, int si, const void* kernel) {
    const SparseArray* sa = call->data.function.args[si]->data.sparse;
    const Expr* other = call->data.function.args[1 - si];
    if (sa->kind == SPV_EXPR) return NULL;
    if (other->type != EXPR_INTEGER && other->type != EXPR_REAL) return NULL;
    Expr* view = sp_view(sa->kind, sa->vals, sa->pat->nnz);
    if (!view) return NULL;
    Expr* nd = si == 0 ? ndarray_map_binary(view, other, (const NDBinaryKernel*)kernel)
                       : ndarray_map_binary(other, view, (const NDBinaryKernel*)kernel);
    sp_view_free(view);
    SpValKind kind;
    void* vals;
    if (!sp_adopt_ndarray(nd, &kind, &vals)) return NULL;
    Expr* repl[2] = { NULL, NULL };
    repl[si] = expr_copy(sa->bg);
    Expr* bg = eval_with(call, repl);
    if (!bg || is_list_value(bg)) {
        if (bg) expr_free(bg);
        free(vals);
        return NULL;
    }
    return sp_to_expr(sp_new(sa->rank, sa->dims, bg, sp_pattern_retain(sa->pat), kind, vals));
}

/* ======================================================================== */
/* Plus / Times folding                                                     */
/* ======================================================================== */

/* Plus/Times over typed operands: fold the sparse operands pairwise, then the
 * (pre-summed / pre-multiplied) scalars. NULL when any step declines. */
static Expr* thread_arith(const Expr* call, const int* which, int m, bool is_plus) {
    size_t n = call->data.function.arg_count;
    Expr* acc = expr_copy(call->data.function.args[which[0]]);
    for (int j = 1; j < m && acc; j++) {
        Expr* nx = sp_arith_binary(is_plus, acc->data.sparse,
                                   call->data.function.args[which[j]]->data.sparse);
        expr_free(acc);
        acc = nx;
    }
    if (!acc) return NULL;
    if ((size_t)m == n) return acc;
    /* the scalar operands, combined once */
    Expr** sc = malloc(sizeof(Expr*) * n);
    size_t ns = 0;
    for (size_t i = 0; i < n; i++)
        if (call->data.function.args[i]->type != EXPR_SPARSE)
            sc[ns++] = expr_copy(call->data.function.args[i]);
    Expr* c;
    if (ns == 1) { c = sc[0]; }
    else {
        Expr* e = expr_new_function(expr_new_symbol(is_plus ? SYM_Plus : SYM_Times), sc, ns);
        c = evaluate(e);
        expr_free(e);
    }
    free(sc);
    Expr* r = (c && c->type != EXPR_SPARSE) ? sp_arith_scalar(is_plus, acc->data.sparse, c) : NULL;
    if (c) expr_free(c);
    expr_free(acc);
    return r;
}

/* ======================================================================== */
/* Entry point                                                              */
/* ======================================================================== */

Expr* sp_thread_listable(const Expr* call) {
    if (call->type != EXPR_FUNCTION) return NULL;
    size_t n = call->data.function.arg_count;
    int which_buf[64];
    int* which = n <= 64 ? which_buf : malloc(sizeof(int) * n);
    int m = 0;
    for (size_t i = 0; i < n; i++)
        if (call->data.function.args[i]->type == EXPR_SPARSE) which[m++] = (int)i;
    Expr* result = NULL;
    if (m == 0) goto done;

    /* Equal dimensions, or Thread::tdlen. */
    const SparseArray* a0 = call->data.function.args[which[0]]->data.sparse;
    for (int j = 1; j < m; j++) {
        const SparseArray* aj = call->data.function.args[which[j]]->data.sparse;
        bool eq = aj->rank == a0->rank;
        for (int k = 0; eq && k < a0->rank; k++) eq = aj->dims[k] == a0->dims[k];
        if (!eq) {
            char* s = expr_to_string_outputform((Expr*)call);
            mth_message("Thread", "tdlen", "Objects of unequal length in %s cannot be combined.", s);
            free(s);
            goto done;
        }
    }

    const Expr* head = call->data.function.head;
    const char* hname = head->type == EXPR_SYMBOL ? head->data.symbol.name : NULL;
    bool is_plus = hname == SYM_Plus, is_times = hname == SYM_Times;
    /* N numericalises a SparseArray itself, converting a typed value buffer
     * in one pass (numeric.c); threading it per element would box them all. */
    if (hname == SYM_N) goto done;

    if (is_plus || is_times) {
        bool typed = true;
        for (int j = 0; j < m && typed; j++)
            typed = call->data.function.args[which[j]]->data.sparse->kind != SPV_EXPR;
        if (typed) result = thread_arith(call, which, m, is_plus);
    } else if (hname == SYM_Power && n == 2 && m == 1) {
        result = sp_power_scalar(call->data.function.args[which[0]]->data.sparse,
                                 call->data.function.args[1 - which[0]], which[0] == 0);
        if (result) {
            /* the background still goes through Power itself */
            Expr* repl[2] = { NULL, NULL };
            repl[which[0]] = expr_copy(call->data.function.args[which[0]]->data.sparse->bg);
            Expr* bg = eval_with(call, repl);
            SparseArray* r = result->data.sparse;
            if (!bg || is_list_value(bg)) {
                if (bg) expr_free(bg);
                expr_free(result);
                result = NULL;
            } else {
                expr_free(r->bg);
                r->bg = bg;
            }
        }
    } else if (hname) {
        SymbolDef* def = symtab_lookup(hname);
        if (def && n == 1 && def->ndarray_unary_kernel)
            result = thread_unary_kernel(call, a0, def->ndarray_unary_kernel);
        else if (def && n == 2 && m == 1 && def->ndarray_binary_kernel)
            result = thread_binary_kernel(call, which[0], def->ndarray_binary_kernel);
    }
    if (!result) result = thread_generic(call, which, m, is_times);
done:
    if (which != which_buf) free(which);
    return result;
}
