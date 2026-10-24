/*
 * permutations.h -- Cycles, PermutationGroup and the permutation builtins
 * (PermutationList, PermutationCycles, GroupElements, GroupOrder).
 */
#ifndef MATHILDA_PERMUTATIONS_H
#define MATHILDA_PERMUTATIONS_H

#include <stddef.h>
#include <stdint.h>
#include "expr.h"

/* Every element of PermutationGroup[{g1, ...}] (concrete Cycles generators)
 * as 0-based image arrays of length *deg (the largest moved point), *order of
 * them, sorted by image list -- GroupElements' order.  Caller frees.  NULL
 * when group is not such a group, or its elements would not fit in memory. */
uint32_t* perm_group_elements(const Expr* group, size_t* order, size_t* deg);

/* Registers the builtins.  Called from core_init(). */
void permutations_init(void);

#endif /* MATHILDA_PERMUTATIONS_H */
