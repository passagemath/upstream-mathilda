/*
 * bitset.c
 *
 * BitSet[n, k] -- the integer n with the bit k places from the right (the
 * coefficient of 2^k) set to 1. Integers are taken in two's complement with
 * an unlimited number of digits, so a negative n already has every high bit
 * on (BitSet[-1, 100] is -1), and a non-negative n grows to whatever size the
 * bit needs (BitSet[42, 9] is 554, BitSet[0, 64] is 2^64). GMP's mpz_setbit
 * has exactly that semantics for either sign.
 *
 * Attributes Listable, Protected. The evaluator threads plain Lists before the
 * builtin runs, so it only ever sees one (n, k) pair.
 *
 * Everything else is left unevaluated SILENTLY, with no message, including
 * the wrong argument count: k < 0, a bignum k, a non-integer n or k (2.5, x),
 * BitSet[5], BitSet[5, 1, 2]. There are no symbolic rules
 * (BitSet[BitSet[x, 1], 1] stays as it is).
 *
 * Surfaces:
 *   - packed / NDArray: the int64 binary kernel ndk_BitSet_ii and its
 *     whole-buffer twins in src/ndinteger.c answer before this builtin runs.
 *   - Compile[]: OP_BSET_I in src/compile/compile_emit_int.c.
 * All three, and the machine path below, call ndint_bitset_i64, so every
 * surface agrees bit for bit; where it declines (a bignum answer) each
 * surface falls back to this builtin's GMP path.
 */

#include "bitwise.h"
#include "ndinteger.h"   /* ndint_bitset_i64: the shared machine helper */

#include <limits.h>
#include <stdbool.h>
#include <stdint.h>
#include <gmp.h>

/* The largest k mpz_setbit can honour: an mpz holds at most INT_MAX limbs,
 * and mp_bitcnt_t is an unsigned long, which is only 32 bits on LLP64
 * Windows. A bit beyond either has no representation at all, so the call is
 * left unevaluated rather than letting GMP abort the process. */
static bool bitset_k_representable(int64_t k) {
    if ((uint64_t)k > (uint64_t)ULONG_MAX) return false;
    return k < (int64_t)INT_MAX * GMP_NUMB_BITS;
}

Expr* builtin_bitset(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2)
        return NULL;

    Expr* n = res->data.function.args[0];
    Expr* kx = res->data.function.args[1];
    if (kx->type != EXPR_INTEGER || kx->data.integer < 0) return NULL;
    int64_t k = kx->data.integer;

    if (n->type == EXPR_INTEGER) {
        int64_t r;
        if (ndint_bitset_i64(n->data.integer, k, &r)) return expr_new_integer(r);
    } else if (n->type != EXPR_BIGINT) {
        return NULL;
    }

    /* A bignum n, or an int64 n >= 0 whose bit k >= 63 makes a bignum. Set the
     * bit in a fresh node; the operand may be shared and is never written.
     * Setting a bit of a negative bignum moves it toward zero and can land it
     * back in int64 range (BitSet[-2^64, 63] is -2^63), which
     * expr_bigint_normalize demotes. */
    if (!bitset_k_representable(k)) return NULL;
    Expr* r;
    if (n->type == EXPR_INTEGER) {
        r = expr_new_bigint_from_int64(n->data.integer);
    } else {
        r = expr_new_bigint_from_int64(0);
        mpz_set(r->data.bigint, n->data.bigint);
    }
    mpz_setbit(r->data.bigint, (mp_bitcnt_t)k);
    return expr_bigint_normalize(r);
}
