/*
 * bitor.c
 *
 * BitOr[n1, n2, ...] -- the bitwise OR of the integers ni: the integer whose
 * binary representation has a one wherever any of the ni does. Negative
 * integers are read in two's complement (GMP's mpz_ior has exactly that
 * semantics, so bignums of either sign are exact).
 *
 * Attributes Flat, Listable, OneIdentity, Orderless, Protected. By the time the
 * builtin runs the evaluator has flattened nested BitOr, threaded over plain
 * Lists and sorted the arguments canonically, so identical arguments are
 * adjacent.
 *
 * Symbolic reduction (OR is idempotent, with identity 0 and absorbing -1) --
 * the exact dual of BitAnd:
 *   1. identical arguments collapse to one (x | x = x), numeric or not:
 *        BitOr[x, y, y, x] -> BitOr[x, y],   BitOr[2.5, 2.5] -> 2.5
 *   2. the integers fold into one constant c (BitOr[1, 2, x] -> BitOr[3, x]);
 *      c == -1 (all bits on) absorbs everything, symbols and non-integers
 *      included (BitOr[-1, 2, x] -> -1, BitOr[-1, 2.5] -> -1), and c == 0 is
 *      dropped (BitOr[0, x] -> x). Any other constant stays inline.
 *   3. OneIdentity rebuild: nothing left -> 0, one argument -> itself,
 *      otherwise BitOr[...].
 * A non-integer number (2.5, 1/2, I) is NOT an error: it is left in the call
 * exactly as a symbol is, silently (BitOr[2.5, 3] stays unevaluated with no
 * message). There are no absorption/complement rules: BitOr[x, BitAnd[x, y]]
 * and BitOr[x, BitNot[x]] are left alone.
 *
 * Surfaces:
 *   - packed / NDArray: the int64 binary kernel ndk_BitOr_ii (+ whole-buffer
 *     twins) in src/ndinteger.c handles the two-argument call before this
 *     builtin runs. The n-ary call (an array with several operands) is folded
 *     here, one buffer pass per operand, through the same kernel.
 *   - Compile[]: OP_BOR_I, folded left in src/compile/compile_emit_int.c.
 * All three compute the same int64 `a | b` as the machine-integer path below,
 * so every surface agrees bit for bit.
 */

#include "bitwise.h"
#include "symtab.h"
#include "sym_names.h"
#include "ndarray.h"

#include <stdlib.h>
#include <stdint.h>
#include <gmp.h>

/* OR accumulator over exact integers. Machine operands fold into the int64
 * `i` (starting from the identity 0; the OR of int64 values is an int64).
 * BigInt operands fold into `big`, the RESULT node itself: allocated on the
 * first BigInt (set to it, not ORed with 0), private (refcount 1, not yet
 * published), and ORed into in place. A bignum operand is read where it lives
 * and the answer is never copied out again. `i` is folded into `big` once, at
 * the end -- a negative `i` can shrink the result back into int64 range, which
 * expr_bigint_normalize then demotes. */
typedef struct {
    int64_t i;
    Expr*   big;      /* owned EXPR_BIGINT result node, or NULL */
    size_t  reserve;  /* bits to presize `big` with, 0 = let GMP grow it */
} OrAcc;

static void oacc_init(OrAcc* a) { a->i = 0; a->big = NULL; a->reserve = 0; }

/* Presize for operands of up to `bits` bits, so the in-place ORs never
 * reallocate. One spare limb covers the sign extension of a negative. */
static void oacc_reserve(OrAcc* a, size_t bits) { a->reserve = bits + GMP_NUMB_BITS; }

static void oacc_add(OrAcc* a, const Expr* n) {
    if (n->type == EXPR_INTEGER) { a->i |= n->data.integer; return; }
    if (!a->big) {
        a->big = expr_new_bigint_from_int64(0);
        if (a->reserve) mpz_realloc2(a->big->data.bigint, a->reserve);
        mpz_set(a->big->data.bigint, n->data.bigint);
        return;
    }
    mpz_ior(a->big->data.bigint, a->big->data.bigint, n->data.bigint);
}

/* Fold the machine part into the bignum part (no-op when either is trivial). */
static void oacc_fold(OrAcc* a) {
    if (!a->big || a->i == 0) return;
    Expr* m = expr_new_integer(a->i);
    mpz_t t;
    expr_to_mpz(m, t);        /* the codebase's int64 -> mpz, LLP64-safe */
    mpz_ior(a->big->data.bigint, a->big->data.bigint, t);
    mpz_clear(t);
    expr_free(m);
    a->i = 0;
}

/* `v` is only ever 0 or -1, both well inside `long` on every target. */
static bool oacc_is(OrAcc* a, int64_t v) {
    oacc_fold(a);
    return a->big ? (mpz_cmp_si(a->big->data.bigint, (long)v) == 0) : (a->i == v);
}

/* Consumes the accumulator. */
static Expr* oacc_take(OrAcc* a) {
    oacc_fold(a);
    if (!a->big) return expr_new_integer(a->i);
    Expr* r = expr_bigint_normalize(a->big);   /* back to int64 when it fits */
    a->big = NULL;
    return r;
}

static void oacc_clear(OrAcc* a) {
    if (a->big) expr_free(a->big);
    a->big = NULL;
    a->i = 0;
}

/* An array operand among the arguments (arity >= 3 reaches here; the
 * two-argument call is taken by the evaluator's kernel map first). Folds every
 * int64 array and machine-integer scalar through the registered binary kernel,
 * one buffer pass per array plus one for the combined scalar. Anything else --
 * a BigInt or symbolic operand, a non-int64 buffer, mismatched shapes -- takes
 * the List path, which answers exactly (and emits Thread::tdlen for ragged
 * shapes, as the List form would). */
static Expr* bitor_ndarray(Expr* res) {
    size_t argc = res->data.function.arg_count;
    Expr** args = res->data.function.args;
    SymbolDef* d = symtab_lookup(SYM_BitOr);
    const NDBinaryKernel* k = d ? (const NDBinaryKernel*)d->ndarray_binary_kernel : NULL;
    if (!k) return ndarray_delist_and_reeval(res);

    int64_t c = 0;
    for (size_t i = 0; i < argc; i++) {
        Expr* a = args[i];
        if (is_ndarray(a)) {
            if (a->data.ndarray.dtype != NDT_INT64) return ndarray_delist_and_reeval(res);
        } else if (a->type == EXPR_INTEGER) {
            c |= a->data.integer;
        } else {
            return ndarray_delist_and_reeval(res);
        }
    }

    Expr* acc = NULL;          /* owned, or NULL while it is still args[first] */
    const Expr* cur = NULL;    /* current value: acc, or a borrowed argument */
    for (size_t i = 0; i < argc; i++) {
        if (!is_ndarray(args[i])) continue;
        if (!cur) { cur = args[i]; continue; }
        Expr* next = ndarray_map_binary2(cur, args[i], k);
        if (acc) expr_free(acc);
        if (!next) return ndarray_delist_and_reeval(res);
        acc = next; cur = next;
    }
    /* A 0 constant (or none) is the identity -- skip the pass. */
    if (c != 0) {
        Expr* s = expr_new_integer(c);
        Expr* next = ndarray_map_binary(cur, s, k);
        expr_free(s);
        if (acc) expr_free(acc);
        if (!next) return ndarray_delist_and_reeval(res);
        acc = next;
    }
    return acc ? acc : expr_copy((Expr*)cur);
}

Expr* builtin_bitor(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    Expr** args = res->data.function.args;
    if (argc == 0) return expr_new_integer(0);
    if (argc == 1) return expr_copy(args[0]);          /* OneIdentity */

    for (size_t i = 0; i < argc; i++)
        if (is_ndarray(args[i])) return bitor_ndarray(res);

    /* Pure-integer call: the common numeric case, no allocation beyond the
     * result unless a BigInt is involved. */
    bool all_int = true;
    for (size_t i = 0; i < argc && all_int; i++)
        all_int = expr_is_integer_like(args[i]);
    if (all_int) {
        OrAcc a; oacc_init(&a);
        size_t bits = 0;
        for (size_t i = 0; i < argc; i++)
            if (args[i]->type == EXPR_BIGINT) {
                size_t b = mpz_size(args[i]->data.bigint) * GMP_NUMB_BITS;
                if (b > bits) bits = b;
            }
        if (bits) oacc_reserve(&a, bits);
        for (size_t i = 0; i < argc; i++) oacc_add(&a, args[i]);
        return oacc_take(&a);
    }

    /* Mixed / symbolic. One scan over the sorted arguments: integers go to the
     * accumulator, each run of identical non-integers is kept once. */
    OrAcc a; oacc_init(&a);
    size_t n_int = 0;
    bool int_trivial = false;      /* a 0 or -1 among the integers */
    bool collapsed = false;        /* a run of identical non-integers */
    Expr** keep = malloc(sizeof(Expr*) * argc);        /* borrowed */
    size_t nk = 0;
    for (size_t i = 0; i < argc; ) {
        Expr* e = args[i];
        if (expr_is_integer_like(e)) {
            n_int++;
            if (e->type == EXPR_INTEGER && (e->data.integer == 0 || e->data.integer == -1))
                int_trivial = true;
            oacc_add(&a, e);
            i++;
            continue;
        }
        size_t j = i + 1;
        while (j < argc && expr_eq(args[j], e)) j++;
        keep[nk++] = e;
        if (j - i > 1) collapsed = true;
        i = j;
    }

    /* -1 (all bits on) absorbs everything (a folded -1 too: BitOr[-2, 1, x]). */
    if (n_int > 0 && oacc_is(&a, -1)) {
        oacc_clear(&a);
        free(keep);
        return expr_new_integer(-1);
    }

    /* Nothing to simplify: no collapse, at most one integer and that one
     * neither 0 nor -1. Returning NULL keeps the fixed-point loop from
     * rebuilding an identical call. */
    if (!collapsed && n_int <= 1 && !int_trivial) {
        oacc_clear(&a);
        free(keep);
        return NULL;
    }

    /* nk >= 1 here: not every argument is an integer. */
    bool have_c = n_int > 0 && !oacc_is(&a, 0);
    size_t m = nk + (have_c ? 1 : 0);
    Expr* out;
    if (m == 1) {
        out = expr_copy(keep[0]);
        oacc_clear(&a);
    } else {
        Expr** ca = malloc(sizeof(Expr*) * m);
        size_t p = 0;
        if (have_c) ca[p++] = oacc_take(&a);   /* Orderless re-sorts it */
        else oacc_clear(&a);
        for (size_t t = 0; t < nk; t++) ca[p++] = expr_copy(keep[t]);
        out = expr_new_function(expr_new_symbol(SYM_BitOr), ca, m);
        free(ca);
    }
    free(keep);
    return out;
}
