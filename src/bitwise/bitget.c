/*
 * bitget.c
 *
 * BitGet[n, k] -- the bit of the integer n k places from the right (the
 * coefficient of 2^k), as 0 or 1. Integers are taken in two's complement with
 * an unlimited number of digits, so every bit of a negative n above its
 * magnitude is 1 (BitGet[-1, 100] is 1, BitGet[-42, 3] is 0). GMP's
 * mpz_tstbit has exactly that semantics for either sign.
 *
 * Attributes Listable, Protected. The evaluator threads plain Lists before the
 * builtin runs, so it only ever sees one (n, k) pair.
 *
 * Everything else is left unevaluated SILENTLY, with no message, including
 * the wrong argument count: k < 0, a bignum k, a non-integer n or k (2.5, x),
 * BitGet[5], BitGet[5, 1, 2]. There are no symbolic rules
 * (BitGet[BitSet[x, 1], 1] stays as it is).
 *
 * The answer is always 0 or 1, so unlike BitSet there is no bignum result and
 * no allocation beyond the returned integer.
 *
 * Surfaces:
 *   - packed / NDArray: the int64 binary kernel ndk_BitGet_ii and its
 *     whole-buffer twins in src/ndinteger.c answer before this builtin runs.
 *   - Compile[]: OP_BGET_I in src/compile/compile_emit_int.c.
 * All three, and the machine path below, call ndint_bitget_i64, so every
 * surface agrees bit for bit; it declines only for k < 0, where every surface
 * leaves BitGet unevaluated.
 */

#include "bitwise.h"
#include "ndinteger.h"   /* ndint_bitget_i64: the shared machine helper */

#include <limits.h>
#include <stdint.h>
#include <gmp.h>

Expr* builtin_bitget(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2)
        return NULL;

    Expr* n = res->data.function.args[0];
    Expr* kx = res->data.function.args[1];
    if (kx->type != EXPR_INTEGER || kx->data.integer < 0) return NULL;
    int64_t k = kx->data.integer;

    if (n->type == EXPR_INTEGER) {
        int64_t r;
        if (ndint_bitget_i64(n->data.integer, k, &r)) return expr_new_integer(r);
        return NULL;
    }
    if (n->type != EXPR_BIGINT) return NULL;

    /* A bignum n: read the bit in place (the operand may be shared and is
     * never written). mp_bitcnt_t is an unsigned long, only 32 bits on LLP64
     * Windows; a k past it lies above every limb an mpz can hold, where the
     * bit is the sign bit. */
    if ((uint64_t)k > (uint64_t)ULONG_MAX)
        return expr_new_integer(mpz_sgn(n->data.bigint) < 0 ? 1 : 0);
    return expr_new_integer(mpz_tstbit(n->data.bigint, (mp_bitcnt_t)k));
}
