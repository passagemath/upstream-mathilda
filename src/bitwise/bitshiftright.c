/*
 * bitshiftright.c
 *
 * BitShiftRight[n, k] -- the bits of the integer n shifted k places to the
 * right, dropping the bits shifted past the units position: IntegerPart[n /
 * 2^k], which truncates toward zero (BitShiftRight[-5, 1] is -2, not the
 * arithmetic shift's -3). BitShiftRight[n] shifts by one place, and a
 * negative k shifts left. Works on integers of any size.
 *
 * BitShiftRight[n, k] is BitShiftLeft[n, -k] on every argument, so the whole
 * implementation -- the domain, messages and surfaces -- is bitshift_eval in
 * bitshiftleft.c, called here with the opposite direction. The machine
 * surfaces call ndint_bitshiftright_i64 (src/ndinteger.h), which is
 * ndint_bitshiftleft_i64 with k negated.
 *
 * Attributes Listable, Protected.
 */

#include "bitwise.h"

Expr* builtin_bitshiftright(Expr* res) {
    return bitshift_eval(res, "BitShiftRight", -1);
}
