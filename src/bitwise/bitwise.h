/*
 * bitwise.h
 *
 * Bit-level integer builtins. Currently:
 *   - BitLength[n]          number of binary bits needed to represent n.
 *   - BitXor[n1, n2, ...]   bitwise XOR (two's complement for negatives).
 *
 * This is the seed of the Bit* family (BitAnd / BitOr / BitNot /
 * BitShiftLeft / BitShiftRight, ...). Each builtin lives in its own <name>.c
 * file and is registered by bitwise_init(); docstrings live in src/info.c.
 */

#ifndef BITWISE_H
#define BITWISE_H

#include "expr.h"

Expr* builtin_bitlength(Expr* res);
Expr* builtin_bitxor(Expr* res);

void bitwise_init(void);

#endif /* BITWISE_H */
