/*
 * sparse_func.c -- native SparseArray kernels for the functional heads.
 *
 * The gate's SPG_ROWS mode hands a head the List of an array's level-1 parts.
 * That is exactly right at level 1 (Scan[f, m], Apply[f, m], MapIndexed[f, m]
 * walk sparse rows), but the rows are SparseArrays -- atoms -- so a level
 * spec reaching below level 1 found nothing there: Map[f, m, {2}] and
 * Apply[f, m, {1}] gave m back, and Scan[f, m, {2}] visited nothing. These
 * kernels take the call first and pick the right representation per level
 * spec:
 *
 *   Map[f, s, lev]       leaf level only ({r}, or {-1} / Infinity over atomic
 *                        elements): f over the stored values and once over
 *                        the background, O(nnz), a SparseArray answer with
 *                        background f[bg]. Level 1 of a matrix: the gate's
 *                        rows. Any other level: the dense call.
 *   MapIndexed / Apply / Scan with an explicit level spec: the dense call
 *                        (every position is visited anyway; a background
 *                        position cannot share one f[...] value when its
 *                        index is an argument). Without one: the gate's rows.
 *   Fold / FoldList      the LIST argument becomes its level-1 parts, so a
 *                        SparseArray initial value stays sparse:
 *                        Fold[Plus, s0, {s1, s2}] and Fold[Plus, m] (a sum of
 *                        sparse rows) are SparseArrays.
 *   Outer                dense, except f = Times over a leading sparse array,
 *                        whose zeros survive: that answer is sparse.
 *   Clip / Chop          element-wise in the first argument: f over the
 *                        values and the background, O(nnz).
 *
 * "The dense call" is the gate's SPG_DENSE, under its implicit-densify cap.
 * Every kernel returns NULL to decline, and the head's gate mode is then the
 * fallback (sparse_gate.c).
 */

#include "sparse_internal.h"
#include "sym_names.h"
#include "symtab.h"
#include "eval.h"
#include <stdlib.h>
#include <string.h>

#define NARGS(c) ((c)->data.function.arg_count)
#define ARG(c, i) ((c)->data.function.args[i])

static bool head_is(const Expr* e, const char* sym) {
    return e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL &&
           e->data.function.head->data.symbol.name == sym;
}

static bool is_rule(const Expr* e) {
    return head_is(e, SYM_Rule) || head_is(e, SYM_RuleDelayed);
}

/* The call answered as on Normal of its sparse arguments (`deep`: also those
 * one and two List levels down). The call itself back (expr_copy) when the
 * head declined or densifying was refused: it then rests, without a second
 * attempt by the gate. NULL on allocation failure. */
static Expr* dense_call(const Expr* call, bool deep) {
    bool rest = false;
    Expr* r = sp_gate(call, SPG_DENSE | (deep ? SPG_DEEP : 0), &rest);
    if (r) return r;
    return rest ? expr_copy((Expr*)call) : NULL;
}

/* An element is atomic when Level[_, {-1}] stops at it: anything but a
 * compound expression, and the Rational / Complex numbers. */
static bool value_atomic(const Expr* v) {
    if (v->type != EXPR_FUNCTION) return v->type != EXPR_SPARSE && v->type != EXPR_NDARRAY;
    return head_is(v, SYM_Rational) || head_is(v, SYM_Complex);
}

static bool elements_atomic(const SparseArray* sa) {
    if (!value_atomic(sa->bg)) return false;
    if (sa->kind != SPV_EXPR) return true;
    Expr* const* vals = (Expr* const*)sa->vals;
    for (int64_t k = 0; k < sp_nnz(sa); k++)
        if (!value_atomic(vals[k])) return false;
    return true;
}

/* Is the level spec `ls` (NULL: the default {1}) the leaf level of `sa`, and
 * nothing else? Positive levels are exact in a rectangular array; a negative
 * level or an Infinity bound counts element depth, so it is the leaf level
 * only when every element is atomic. Anything unrecognised is "no". */
static bool is_leaf_level(const Expr* ls, const SparseArray* sa) {
    int64_t r = sa->rank, lo, hi;
    bool open_hi = false;
    if (!ls) {
        lo = hi = 1;
    } else if (ls->type == EXPR_INTEGER) {
        lo = 1; hi = ls->data.integer;
        if (hi >= r) open_hi = hi > r;
    } else if (ls->type == EXPR_SYMBOL && ls->data.symbol.name == SYM_Infinity) {
        lo = 1; hi = r; open_hi = true;
    } else if (head_is(ls, SYM_List) && (NARGS(ls) == 1 || NARGS(ls) == 2) &&
               ARG(ls, 0)->type == EXPR_INTEGER &&
               ARG(ls, NARGS(ls) - 1)->type == EXPR_INTEGER) {
        lo = ARG(ls, 0)->data.integer;
        hi = ARG(ls, NARGS(ls) - 1)->data.integer;
    } else {
        return false;
    }
    bool atomic = (lo < 0 || hi < 0 || open_hi) ? elements_atomic(sa) : true;
    if (!atomic) return false;
    /* over atomic elements, level -d is level r + 1 - d */
    if (lo < 0) lo = r + 1 + lo;
    if (hi < 0) hi = r + 1 + hi;
    return lo == r && hi == r;
}

/* f over every element: the stored values and once over the background (the
 * Listable-threading tiers: a head's ND kernel over a typed buffer, else f
 * per value). NULL when an element came out a List (the answer is then not
 * an array with those elements). */
static Expr* map_elements(const Expr* f, const Expr* s) {
    Expr* a[1] = { expr_copy((Expr*)s) };
    Expr* fcall = expr_new_function(expr_copy((Expr*)f), a, 1);
    Expr* r = sp_thread_listable(fcall);
    expr_free(fcall);
    return r;
}

static bool same_shape(const Expr* r, const SparseArray* sa, int depth) {
    if (r->type == EXPR_NDARRAY) {
        if (r->data.ndarray.rank != sa->rank - depth) return false;
        for (int k = 0; k < r->data.ndarray.rank; k++)
            if (r->data.ndarray.dims[k] != sa->dims[depth + k]) return false;
        return true;
    }
    bool lst = head_is(r, SYM_List);
    if (depth == sa->rank) return !lst && r->type != EXPR_SPARSE;
    if (!lst || (int64_t)NARGS(r) != sa->dims[depth]) return false;
    for (size_t i = 0; i < NARGS(r); i++)
        if (!same_shape(ARG(r, i), sa, depth + 1)) return false;
    return true;
}

bool sp_same_shape(const Expr* r, const SparseArray* sa) { return same_shape(r, sa, 0); }

/* ======================================================================== */
/* Map                                                                      */
/* ======================================================================== */

static Expr* nat_map(const Expr* call) {
    size_t n = NARGS(call);
    if (n < 2 || n > 4 || ARG(call, 0)->type == EXPR_SPARSE) return NULL;
    const SparseArray* sa = sp_of(ARG(call, 1));
    if (!sa) return NULL;
    /* Map[f, s, lev, Heads -> ...]: the heads of an array's parts are List,
     * which a SparseArray does not expose -- the option changes nothing. */
    const Expr* ls = NULL;
    for (size_t i = 2; i < n; i++) {
        if (is_rule(ARG(call, i))) continue;
        if (ls) return NULL;
        ls = ARG(call, i);
    }
    if (is_leaf_level(ls, sa)) {
        Expr* r = map_elements(ARG(call, 0), ARG(call, 1));
        if (r) return r;
        /* f gives a List on some element: an array of depth > rank */
        return dense_call(call, false);
    }
    /* Level 1 of a matrix: the gate maps over the sparse rows. */
    if (!ls || (ls->type == EXPR_INTEGER && ls->data.integer == 1) ||
        (head_is(ls, SYM_List) && NARGS(ls) == 1 && ARG(ls, 0)->type == EXPR_INTEGER &&
         ARG(ls, 0)->data.integer == 1))
        return NULL;
    /* Any other level. Mapping only inside the elements ({-1} over x^2)
     * keeps the array's shape, and then the answer is an array again. */
    Expr* r = dense_call(call, false);
    if (r && r != call && sp_same_shape(r, sa)) r = sp_gate_resparse(r, sa->bg);
    return r;
}

/* ======================================================================== */
/* MapIndexed / Apply / Scan                                                */
/* ======================================================================== */

/* With an explicit level spec (the third argument, not an option) the head
 * walks Normal[s]; without one -- or Apply's own level {0} -- the gate's
 * rows serve it. */
static Expr* nat_levelled(const Expr* call) {
    size_t n = NARGS(call);
    if (n < 3 || is_rule(ARG(call, 2)) || ARG(call, 0)->type == EXPR_SPARSE) return NULL;
    const Expr* ls = ARG(call, 2);
    const Expr* h = call->data.function.head;
    if (h->type == EXPR_SYMBOL && h->data.symbol.name == SYM_Apply &&
        head_is(ls, SYM_List) && NARGS(ls) == 1 && ARG(ls, 0)->type == EXPR_INTEGER &&
        ARG(ls, 0)->data.integer == 0)
        return NULL;
    return dense_call(call, false);
}

/* ======================================================================== */
/* Fold / FoldList                                                          */
/* ======================================================================== */

/* Fold[f, list] / Fold[f, x, list]: only the list argument is a sequence;
 * an initial value that is a SparseArray is a value, and stays one. */
static Expr* nat_fold(const Expr* call) {
    size_t n = NARGS(call);
    if (n != 2 && n != 3) return NULL;
    const Expr* lst = ARG(call, n - 1);
    if (lst->type != EXPR_SPARSE) {
        /* the sparse argument is f or the initial value: nothing to walk */
        return NULL;
    }
    bool refused = false;
    Expr* parts = sp_gate_rows(lst, &refused);
    if (!parts) return refused ? expr_copy((Expr*)call) : NULL;
    Expr** a = malloc(sizeof(Expr*) * n);
    if (!a) { expr_free(parts); return NULL; }
    for (size_t i = 0; i + 1 < n; i++) a[i] = expr_copy(ARG(call, i));
    a[n - 1] = parts;
    Expr* out = expr_new_function(expr_copy(call->data.function.head), a, n);
    free(a);
    return out;
}

/* ======================================================================== */
/* Outer                                                                    */
/* ======================================================================== */

/* Outer[f, a1, a2, ...]: the dense outer product -- except for f = Times with
 * a sparse a1, whose zero background survives every product: that answer is
 * sparse over a1's background (Outer[Times, s, v], Outer[Times, s, s]). */
static Expr* nat_outer(const Expr* call) {
    size_t n = NARGS(call);
    if (n < 2) return NULL;
    /* Outer[f, s]: f at every element, an element map that stays sparse */
    if (n == 2 && ARG(call, 1)->type == EXPR_SPARSE && ARG(call, 0)->type != EXPR_SPARSE) {
        Expr* r = map_elements(ARG(call, 0), ARG(call, 1));
        if (r) return r;
    }
    Expr* r = dense_call(call, true);
    if (!r || r == call) return r;
    const Expr* f = ARG(call, 0);
    const Expr* a1 = n >= 2 ? ARG(call, 1) : NULL;
    if (f->type == EXPR_SYMBOL && f->data.symbol.name == SYM_Times &&
        a1 && a1->type == EXPR_SPARSE && sp_is_zero(sp_of(a1)->bg))
        r = sp_gate_resparse(r, sp_of(a1)->bg);
    return r;
}

/* ======================================================================== */
/* Clip / Chop: element-wise in the first argument                          */
/* ======================================================================== */

static Expr* nat_elementwise1(const Expr* call) {
    size_t n = NARGS(call);
    if (n < 1 || ARG(call, 0)->type != EXPR_SPARSE) return NULL;
    for (size_t i = 1; i < n; i++)
        if (ARG(call, i)->type == EXPR_SPARSE) return NULL;
    Expr* r = sp_thread_listable(call);
    /* A spec the head does not accept (Chop[s, x]) leaves every element
     * unevaluated: the call is not element-wise then -- decline, so it stays
     * unevaluated as a whole, as on the dense form. */
    if (r && head_is(sp_of(r)->bg, call->data.function.head->data.symbol.name)) {
        expr_free(r);
        return NULL;
    }
    return r;
}

/* ======================================================================== */
/* Registration                                                             */
/* ======================================================================== */

void sparse_func_init(void) {
    symtab_set_sparse_native("Map", (const void*)nat_map);
    symtab_set_sparse_native("MapIndexed", (const void*)nat_levelled);
    symtab_set_sparse_native("Apply", (const void*)nat_levelled);
    symtab_set_sparse_native("Scan", (const void*)nat_levelled);
    symtab_set_sparse_native("Fold", (const void*)nat_fold);
    symtab_set_sparse_native("FoldList", (const void*)nat_fold);
    symtab_set_sparse_native("Outer", (const void*)nat_outer);
    symtab_set_sparse_native("Clip", (const void*)nat_elementwise1);
    symtab_set_sparse_native("Chop", (const void*)nat_elementwise1);
}
