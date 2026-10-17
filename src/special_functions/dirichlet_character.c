/* dirichlet_character.c -- DirichletCharacter[k, j, n].
 *
 * The Dirichlet character chi_{k,j}(n) in Mathematica's labelling (see
 * dirichlet.h). Values are exact: 0, a unit (+-1, +-I) or unit * E^(I Pi r).
 *
 * Attributes: Protected (as in Mathematica). DirichletCharacter is NOT Listable,
 * yet it threads over lists in every argument; that is done here.
 *
 * Packed / NDArray surface: an int64 buffer of n with scalar k, j is answered
 * straight from the character table. A real character (values 0, +-1) returns
 * an int64 array of the same shape and presentation; a complex character
 * returns a List of exact roots of unity -- those are not machine numbers, so
 * there is no buffer to return (recorded in the array-exactness EXEMPT list).
 * DirichletCharacter is exact-valued, so it has no Compile[] lowering
 * (compile-coverage EXEMPT).
 */
#include "dirichlet_character.h"
#include "dirichlet.h"
#include "symtab.h"
#include "attr.h"
#include "eval.h"
#include "common.h"        /* builtin_arg_error */
#include "message.h"
#include "print.h"
#include "sym_names.h"
#include "ndarray.h"
#include "arithmetic.h"
#include <gmp.h>
#include <stdlib.h>
#include <string.h>

static bool is_list(const Expr* e) {
    return e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == SYM_List;
}

/* n mod k for an exact integer n (int64 or BigInt, any sign). */
static bool int_mod_k(const Expr* n, uint64_t k, uint64_t* out) {
    if (n->type == EXPR_INTEGER) {
        int64_t v = n->data.integer;
        if (k == 1) { *out = 0; return true; }
        int64_t r = v % (int64_t)k;
        if (r < 0) r += (int64_t)k;
        *out = (uint64_t)r;
        return true;
    }
    if (n->type == EXPR_BIGINT) {
        mpz_t r; mpz_init(r);
        mpz_t K; mpz_init(K);
        mpz_import(K, 1, 1, sizeof(uint64_t), 0, 0, &k);
        mpz_fdiv_r(r, n->data.bigint, K);
        uint64_t v = 0;
        mpz_export(&v, NULL, 1, sizeof(uint64_t), 0, 0, r);
        mpz_clear(r); mpz_clear(K);
        *out = v;
        return true;
    }
    return false;
}

/* Real character over an int64 buffer -> int64 buffer of 0, +-1. */
static Expr* dc_ndarray_real(const DchTable* T, const Expr* nd) {
    size_t n = ndarray_size(nd);
    int64_t* out = (int64_t*)malloc((n ? n : 1) * sizeof(int64_t));
    if (!out) return NULL;
    uint64_t k = (uint64_t)T->k;
    for (size_t i = 0; i < n; i++) {
        int64_t v = ndt_get_i(nd->data.ndarray.data, i, nd->data.ndarray.dtype);
        int64_t r = (k == 1) ? 0 : v % (int64_t)k;
        if (r < 0) r += (int64_t)k;
        uint64_t e;
        if (!dch_exp(T, (uint64_t)r, &e)) out[i] = 0;
        else out[i] = (e == 0) ? 1 : -1;
    }
    return expr_new_ndarray_like(nd, nd->data.ndarray.rank, nd->data.ndarray.dims,
                                 out, NDT_INT64);
}

/* Complex character over an int64 buffer -> nested List of exact values.
 * Values are memoised per residue (an Expr copy is a refcount bump). */
static Expr* dc_build_list(const DchTable* T, const Expr* nd, int axis,
                           size_t* idx, Expr** memo, size_t memo_n) {
    const NDArrayData* A = &nd->data.ndarray;
    int64_t len = A->dims[axis];
    Expr** items = (Expr**)malloc((size_t)(len ? len : 1) * sizeof(Expr*));
    for (int64_t i = 0; i < len; i++) {
        if (axis + 1 < A->rank) {
            items[i] = dc_build_list(T, nd, axis + 1, idx, memo, memo_n);
            continue;
        }
        int64_t v = ndt_get_i(A->data, (*idx)++, A->dtype);
        uint64_t k = (uint64_t)T->k;
        int64_t r = (k == 1) ? 0 : v % (int64_t)k;
        if (r < 0) r += (int64_t)k;
        Expr* val;
        if ((uint64_t)r < memo_n) {
            if (!memo[r]) memo[r] = dch_value_expr_mpz_mod(T, (uint64_t)r);
            val = memo[r] ? expr_copy(memo[r]) : NULL;
        } else {
            val = dch_value_expr_mpz_mod(T, (uint64_t)r);
        }
        items[i] = val ? val : expr_new_integer(0);
    }
    Expr* out = expr_new_function(expr_new_symbol(SYM_List), items, (size_t)len);
    free(items);
    return out;
}

static Expr* dc_ndarray(Expr* res) {
    Expr** a = res->data.function.args;
    if (!is_ndarray(a[2]) || is_ndarray(a[0]) || is_ndarray(a[1])) return NULL;
    if (a[2]->data.ndarray.dtype != NDT_INT64) return NULL;
    if (a[0]->type != EXPR_INTEGER || a[1]->type != EXPR_INTEGER) return NULL;
    DchTable* owned = NULL;
    const DchTable* T = dch_validate_kj("DirichletCharacter", res, false, &owned);
    if (!T) return NULL;
    Expr* out;
    if (T->real) {
        out = dc_ndarray_real(T, a[2]);
    } else {
        size_t memo_n = (T->k <= 65536) ? (size_t)T->k : 0;
        Expr** memo = memo_n ? (Expr**)calloc(memo_n, sizeof(Expr*)) : NULL;
        size_t idx = 0;
        out = dc_build_list(T, a[2], 0, &idx, memo, memo ? memo_n : 0);
        if (memo) { for (size_t i = 0; i < memo_n; i++) if (memo[i]) expr_free(memo[i]); free(memo); }
    }
    dch_free(owned);
    return out;
}

/* Thread over List arguments: DirichletCharacter[{3,4}, 2, {1,2}] ->
 * {DirichletCharacter[3,2,1], DirichletCharacter[4,2,2]}. Lists must agree in
 * length (otherwise the call is left alone). */
static Expr* dc_thread(Expr* res) {
    Expr** a = res->data.function.args;
    int64_t len = -1;
    for (int i = 0; i < 3; i++) {
        if (!is_list(a[i])) continue;
        int64_t L = (int64_t)a[i]->data.function.arg_count;
        if (len >= 0 && L != len) return NULL;
        len = L;
    }
    if (len < 0) return NULL;
    Expr** items = (Expr**)malloc((size_t)(len ? len : 1) * sizeof(Expr*));
    for (int64_t t = 0; t < len; t++) {
        Expr* ca[3];
        for (int i = 0; i < 3; i++)
            ca[i] = expr_copy(is_list(a[i]) ? a[i]->data.function.args[t] : a[i]);
        items[t] = expr_new_function(expr_new_symbol(SYM_DirichletCharacter), ca, 3);
    }
    Expr* out = expr_new_function(expr_new_symbol(SYM_List), items, (size_t)len);
    free(items);
    return eval_and_free(out);
}

Expr* builtin_dirichletcharacter(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc != 3) return builtin_arg_error("DirichletCharacter", argc, 3, 3);
    Expr** a = res->data.function.args;

    if (is_ndarray(a[0]) || is_ndarray(a[1]) || is_ndarray(a[2])) {
        Expr* out = dc_ndarray(res);
        if (out) return out;
        return ndarray_delist_and_reeval(res);
    }
    if (is_list(a[0]) || is_list(a[1]) || is_list(a[2])) return dc_thread(res);

    DchTable* owned = NULL;
    const DchTable* T = dch_validate_kj("DirichletCharacter", res, false, &owned);
    if (!T) return NULL;

    Expr* out = NULL;
    uint64_t r;
    if (int_mod_k(a[2], (uint64_t)T->k, &r)) {
        out = dch_value_expr_mpz_mod(T, r);
    } else if (a[2]->type == EXPR_REAL || dch_is_mpfr(a[2])
               || is_rational_like(a[2])) {
        char* s = expr_to_string(res);
        mth_message("DirichletCharacter", "int", "Integer expected at position 3 in %s.",
                    s ? s : "?");
        free(s);
    } else {
        Expr *re, *im;
        if (is_complex(a[2], &re, &im)) {
            char* s = expr_to_string(res);
            mth_message("DirichletCharacter", "int", "Integer expected at position 3 in %s.",
                        s ? s : "?");
            free(s);
        }
    }
    dch_free(owned);
    return out;
}

void dirichletcharacter_init(void) {
    symtab_add_builtin("DirichletCharacter", builtin_dirichletcharacter);
    symtab_get_def("DirichletCharacter")->attributes |= ATTR_PROTECTED;
    /* Docstring lives in info.c (info_init). */
}
