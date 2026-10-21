/*
 * bitwise.h
 *
 * Bit-level integer builtins. Currently:
 *   - BitLength[n]          number of binary bits needed to represent n.
 *   - BitXor[n1, n2, ...]   bitwise XOR (two's complement for negatives).
 *   - BitAnd[n1, n2, ...]   bitwise AND (two's complement for negatives).
 *   - BitOr[n1, n2, ...]    bitwise OR (two's complement for negatives).
 *   - BitNot[n]             bitwise NOT, -1 - n (two's complement).
 *
 * The rest of the Bit* family (BitShiftLeft / BitShiftRight, ...) is still
 * to come. Each builtin lives in its own <name>.c
 * file and is registered by bitwise_init(); docstrings live in src/info.c.
 */

#ifndef BITWISE_H
#define BITWISE_H

#include "expr.h"

Expr* builtin_bitlength(Expr* res);
Expr* builtin_bitxor(Expr* res);
Expr* builtin_bitand(Expr* res);
Expr* builtin_bitor(Expr* res);
Expr* builtin_bitnot(Expr* res);

void bitwise_init(void);

#endif /* BITWISE_H */
