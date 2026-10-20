/* Mathilda -- KroneckerDelta[n1, n2, ...]: 1 when every ni is equal, 0 when
 * any two are unequal.
 *
 *   KroneckerDelta[]            1
 *   KroneckerDelta[n]           KroneckerDelta[n, 0]
 *   KroneckerDelta[x, x]        1   (identical arguments, symbolic or not)
 *   KroneckerDelta[1, 2, x]     0   (one provably unequal pair settles it)
 *   KroneckerDelta[1/2, 0.5]    1   (numbers compare by value, as Equal does)
 *   KroneckerDelta[x, y]        stays unevaluated
 *
 * A pair is decided only when both sides are numbers or numeric quantities
 * (Pi, Sqrt[2]); Equal settles it, exactly where an exact zero test can. A
 * symbolic pair that is not identical is never guessed.
 *
 * Attributes: NumericFunction, Orderless, Protected. Not Listable: a List
 * argument is not a number, so KroneckerDelta[{1, 2}, 1] stays unevaluated.
 */
#include "kroneckerdelta.h"
#include "sym_names.h"

#include <stdbool.h>
#include <stdlib.h>

#include "attr.h"
#include "eval.h"
#include "expr.h"
#include "numeric.h"
#include "symtab.h"

/* 1 equal, 0 unequal, -1 undecided. */
static int kd_pair(Expr* a, Expr* b) {
    if (expr_eq(a, b)) return 1;
    if (!expr_is_numeric_quantity(a) || !expr_is_numeric_quantity(b)) return -1;
    Expr* args[2] = { expr_copy(a), expr_copy(b) };
    Expr* r = eval_and_free(expr_new_function(expr_new_symbol(SYM_Equal), args, 2));
    int v = -1;
    if (r->type == EXPR_SYMBOL) {
        if (r->data.symbol.name == SYM_True) v = 1;
        else if (r->data.symbol.name == SYM_False) v = 0;
    }
    expr_free(r);
    return v;
}

Expr* builtin_kroneckerdelta(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t n = res->data.function.arg_count;
    Expr** a = res->data.function.args;
    if (n == 0) return expr_new_integer(1);
    if (n == 1) {
        Expr* zero = expr_new_integer(0);
        int v = kd_pair(a[0], zero);
        expr_free(zero);
        return v < 0 ? NULL : expr_new_integer(v);
    }
    /* Any provably unequal pair gives 0; all pairs equal gives 1. */
    bool all_equal = true;
    for (size_t i = 0; i < n; i++) {
        for (size_t j = i + 1; j < n; j++) {
            int v = kd_pair(a[i], a[j]);
            if (v == 0) return expr_new_integer(0);
            if (v < 0) all_equal = false;
        }
    }
    return all_equal ? expr_new_integer(1) : NULL;
}

void kroneckerdelta_init(void) {
    symtab_add_builtin("KroneckerDelta", builtin_kroneckerdelta);
    symtab_get_def("KroneckerDelta")->attributes |=
        (ATTR_NUMERICFUNCTION | ATTR_ORDERLESS | ATTR_PROTECTED);
    /* Docstring lives in info.c (info_init). */
}
