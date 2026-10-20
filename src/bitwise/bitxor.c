/*
 * bitxor.c
 *
 * BitXor[n1, n2, ...] -- the bitwise XOR of the integers ni: the integer whose
 * binary representation has a one wherever an odd number of the ni do.
 * Negative integers are read in two's complement (GMP's mpz_xor has exactly
 * that semantics, so bignums of either sign are exact).
 *
 * Attributes Flat, Listable, OneIdentity, Orderless, Protected. By the time the
 * builtin runs the evaluator has flattened nested BitXor, threaded over plain
 * Lists and sorted the arguments canonically, so identical arguments are
 * adjacent and integers lead.
 *
 * Symbolic reduction, in this order (the order is observable -- see step 2):
 *   1. identical arguments cancel in pairs (x^x = 0), numeric or not:
 *        BitXor[x, x, y] -> y,   BitXor[x, 2.5, 2.5] -> x
 *   2. an odd number of -1 arguments is pulled out as BitNot (x^-1 = ~x),
 *      BEFORE the integer constants are combined:
 *        BitXor[1, -1, x] -> BitNot[BitXor[1, x]]
 *   3. the remaining integers fold into one constant c; c == 0 is dropped and
 *      c == -1 is pulled out as BitNot too (BitXor[-2, 1, x] -> BitNot[x]);
 *      any other negative constant stays inline (BitXor[-5, x] is left alone).
 *   4. OneIdentity rebuild: nothing left -> 0 (or -1 under BitNot), one
 *      argument -> itself, otherwise BitXor[...].
 * A concrete non-integer number (2.5, 1/2, I) that survives pair cancellation
 * is an error: BitXor::int names its position and the call stays unevaluated.
 * Symbolic arguments are silent -- they may yet become integers, and they are
 * what the cancellation rules act on. (A cancelled pair is never reported:
 * BitXor[x, 2.5, 2.5] is x, because 2.5 ^ 2.5 is 0 whatever 2.5 means.)
 *
 * Surfaces:
 *   - packed / NDArray: the int64 binary kernel ndk_BitXor_ii (+ whole-buffer
 *     twin) in src/ndinteger.c handles the two-argument call before this
 *     builtin runs. The n-ary call (an array with several operands) is folded
 *     here, one buffer pass per operand, through the same kernel.
 *   - Compile[]: OP_BXOR_I, folded left in src/compile/compile_emit_int.c.
 * All three compute the same int64 `a ^ b` as the machine-integer path below,
 * so every surface agrees bit for bit.
 */

#include "bitwise.h"
#include "symtab.h"
#include "sym_names.h"
#include "ndarray.h"

#include "print.h"        /* expr_to_string */
#include "message.h"      /* mth_message: Quiet/Check funnel */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <gmp.h>

/* `BitXor::int: Integer expected at position <pos> in <call>.`  Matches the
 * BitLength::int / IntegerLength::int surface diagnostic. */
static Expr* bitxor_emit_int(size_t pos, Expr* res) {
    char* call_str = expr_to_string(res);
    mth_message("BitXor", "int", "Integer expected at position %zu in %s.", pos, call_str ? call_str : "?");
    free(call_str);
    return NULL;
}

/* XOR accumulator over exact integers: stays in int64 until the first BigInt
 * operand arrives, then continues in GMP. XOR of two int64 values is always an
 * int64, so the machine path can never overflow. */
typedef struct {
    int64_t i;
    bool    big;    /* true once `z` is live and authoritative */
    mpz_t   z;
} XorAcc;

static void xacc_init(XorAcc* a) { a->i = 0; a->big = false; }

static void xacc_add(XorAcc* a, const Expr* n) {
    if (n->type == EXPR_INTEGER && !a->big) { a->i ^= n->data.integer; return; }
    if (!a->big) {
        /* Promote the machine accumulator. Via an Expr so the int64 -> mpz
         * conversion never narrows through `long` (32-bit on LLP64). */
        Expr* cur = expr_new_integer(a->i);
        expr_to_mpz(cur, a->z);
        expr_free(cur);
        a->big = true;
    }
    mpz_t t;
    expr_to_mpz(n, t);
    mpz_xor(a->z, a->z, t);
    mpz_clear(t);
}

static bool xacc_is(const XorAcc* a, int64_t v) {
    return a->big ? (mpz_cmp_si(a->z, (long)v) == 0) : (a->i == v);
}

/* Consumes the accumulator. */
static Expr* xacc_take(XorAcc* a) {
    if (!a->big) return expr_new_integer(a->i);
    Expr* r = expr_bigint_normalize(expr_new_bigint_from_mpz(a->z));
    mpz_clear(a->z);
    a->big = false;
    return r;
}

static void xacc_clear(XorAcc* a) {
    if (a->big) mpz_clear(a->z);
    a->big = false;
}

static bool is_minus_one(const Expr* e) {
    return e->type == EXPR_INTEGER && e->data.integer == -1;
}

/* An array operand among the arguments (arity >= 3 reaches here; the
 * two-argument call is taken by the evaluator's kernel map first). Folds every
 * int64 array and machine-integer scalar through the registered binary kernel,
 * one buffer pass per array plus one for the combined scalar. Anything else --
 * a BigInt or symbolic operand, a non-int64 buffer, mismatched shapes -- takes
 * the List path, which answers exactly (and emits Thread::tdlen for ragged
 * shapes, as the List form would). */
static Expr* bitxor_ndarray(Expr* res) {
    size_t argc = res->data.function.arg_count;
    Expr** args = res->data.function.args;
    SymbolDef* d = symtab_lookup(SYM_BitXor);
    const NDBinaryKernel* k = d ? (const NDBinaryKernel*)d->ndarray_binary_kernel : NULL;
    if (!k) return ndarray_delist_and_reeval(res);

    int64_t c = 0;
    bool any_scalar = false;
    for (size_t i = 0; i < argc; i++) {
        Expr* a = args[i];
        if (is_ndarray(a)) {
            if (a->data.ndarray.dtype != NDT_INT64) return ndarray_delist_and_reeval(res);
        } else if (a->type == EXPR_INTEGER) {
            c ^= a->data.integer;
            any_scalar = true;
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
    /* A zero constant is the identity -- skip the pass (BitXor[v, 1, 1]). */
    if (any_scalar && c != 0) {
        Expr* s = expr_new_integer(c);
        Expr* next = ndarray_map_binary(cur, s, k);
        expr_free(s);
        if (acc) expr_free(acc);
        if (!next) return ndarray_delist_and_reeval(res);
        acc = next;
    }
    return acc ? acc : expr_copy((Expr*)cur);
}

Expr* builtin_bitxor(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    Expr** args = res->data.function.args;
    if (argc == 0) return expr_new_integer(0);
    if (argc == 1) return expr_copy(args[0]);          /* OneIdentity */

    for (size_t i = 0; i < argc; i++)
        if (is_ndarray(args[i])) return bitxor_ndarray(res);

    /* Pure-integer call: the common numeric case, no allocation beyond the
     * result unless a BigInt is involved. */
    bool all_int = true;
    for (size_t i = 0; i < argc && all_int; i++)
        all_int = expr_is_integer_like(args[i]);
    if (all_int) {
        XorAcc a; xacc_init(&a);
        for (size_t i = 0; i < argc; i++) xacc_add(&a, args[i]);
        return xacc_take(&a);
    }

    /* Mixed / symbolic. One scan over the sorted arguments: integers go to the
     * accumulator (-1s counted separately, step 2), non-integers are kept when
     * their run of identical neighbours has odd length (step 1). */
    XorAcc a; xacc_init(&a);
    size_t n_int = 0, n_m1 = 0;
    bool int_trivial = false;      /* a 0 or -1 among the integers */
    bool cancelled = false;
    Expr** keep = malloc(sizeof(Expr*) * argc);        /* borrowed */
    size_t nk = 0;
    for (size_t i = 0; i < argc; ) {
        Expr* e = args[i];
        if (expr_is_integer_like(e)) {
            n_int++;
            if (is_minus_one(e)) { n_m1++; int_trivial = true; }
            else {
                if (e->type == EXPR_INTEGER && e->data.integer == 0) int_trivial = true;
                xacc_add(&a, e);
            }
            i++;
            continue;
        }
        size_t j = i + 1;
        while (j < argc && expr_eq(args[j], e)) j++;
        if ((j - i) & 1) {
            if (expr_is_numeric_like(e)) {
                xacc_clear(&a);
                free(keep);
                return bitxor_emit_int(i + 1, res);
            }
            keep[nk++] = e;
        }
        if (j - i > 1) cancelled = true;
        i = j;
    }

    /* Nothing to simplify: no cancellation, at most one integer and that one
     * neither 0 nor -1. Returning NULL keeps the fixed-point loop from
     * rebuilding an identical call. */
    if (!cancelled && n_int <= 1 && !int_trivial) {
        xacc_clear(&a);
        free(keep);
        return NULL;
    }

    bool bitnot = (n_m1 & 1) != 0;
    if (nk == 0) {
        /* Every non-integer cancelled: the answer is the integer constant. */
        free(keep);
        if (!bitnot) return xacc_take(&a);
        /* c ^ -1 = -c - 1 */
        Expr* c = xacc_take(&a);
        XorAcc b; xacc_init(&b);
        xacc_add(&b, c);
        Expr* m1 = expr_new_integer(-1);
        xacc_add(&b, m1);
        expr_free(m1);
        expr_free(c);
        return xacc_take(&b);
    }
    if (xacc_is(&a, -1)) { bitnot = !bitnot; xacc_clear(&a); xacc_init(&a); }
    bool have_c = !xacc_is(&a, 0);

    size_t m = nk + (have_c ? 1 : 0);
    Expr* core;
    if (m == 1) {
        core = expr_copy(keep[0]);
        xacc_clear(&a);
    } else {
        Expr** ca = malloc(sizeof(Expr*) * m);
        size_t p = 0;
        if (have_c) ca[p++] = xacc_take(&a);   /* integers sort first */
        else xacc_clear(&a);
        for (size_t t = 0; t < nk; t++) ca[p++] = expr_copy(keep[t]);
        core = expr_new_function(expr_new_symbol(SYM_BitXor), ca, m);
        free(ca);
    }
    free(keep);
    if (bitnot) {
        Expr* na[1] = { core };
        return expr_new_function(expr_new_symbol(SYM_BitNot), na, 1);
    }
    return core;
}
