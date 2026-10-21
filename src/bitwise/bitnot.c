/*
 * bitnot.c
 *
 * BitNot[n] -- the bitwise NOT of the integer n: ones become zeros and vice
 * versa. Integers are taken in two's complement with an unlimited number of
 * digits, so BitNot[n] is exactly -1 - n (BitNot[0] is -1, all bits on).
 * GMP's mpz_com has precisely that semantics, so bignums of either sign are
 * exact; the int64 complement ~n is total (BitNot[-2^63] is 2^63 - 1).
 *
 * Attributes Listable, Protected. The evaluator threads plain Lists before the
 * builtin runs, so it only ever sees one element.
 *
 * Symbolic reduction: BitNot is an involution, so BitNot[BitNot[x]] -> x.
 * Nothing else: BitNot[-1 - x] stays as it is. A non-integer number (2.5,
 * 1/2, I, Pi) is NOT an error: it is left unevaluated silently, like a symbol.
 * The wrong argument count emits BitNot::argx.
 *
 * Surfaces:
 *   - packed / NDArray: the int64 unary kernel ndk_BitNot_ii and its
 *     whole-buffer twin ndk_BitNot_whole in src/ndinteger.c answer before this
 *     builtin runs.
 *   - Compile[]: OP_BNOT_I in src/compile/compile_emit_int.c.
 * All three compute the same int64 `~n` as the machine-integer path below,
 * so every surface agrees bit for bit.
 */

#include "bitwise.h"
#include "sym_names.h"
#include "message.h"   /* mth_message: Quiet/Check funnel */

#include <stdint.h>
#include <gmp.h>

Expr* builtin_bitnot(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc != 1) {
        /* argx never fires at argc == 1, so the count is always plural. */
        mth_message("BitNot", "argx",
                    "BitNot called with %zu arguments; 1 argument is expected.",
                    argc);
        return NULL;
    }

    Expr* n = res->data.function.args[0];
    switch (n->type) {
    case EXPR_INTEGER:
        return expr_new_integer(~n->data.integer);
    case EXPR_BIGINT: {
        /* Complement into a fresh node; the operand may be shared and is
         * never written. -1 - n can land back in int64 range (n = 2^63 gives
         * -2^63), which expr_bigint_normalize demotes. */
        Expr* r = expr_new_bigint_from_int64(0);
        mpz_com(r->data.bigint, n->data.bigint);
        return expr_bigint_normalize(r);
    }
    case EXPR_FUNCTION:
        /* BitNot[BitNot[x]] -> x. The inner node may be shared, so it is not
         * stolen from; expr_copy is a refcount bump. */
        if (n->data.function.arg_count == 1
            && n->data.function.head->type == EXPR_SYMBOL
            && n->data.function.head->data.symbol.name == SYM_BitNot)
            return expr_copy(n->data.function.args[0]);
        return NULL;
    default:
        return NULL;
    }
}
