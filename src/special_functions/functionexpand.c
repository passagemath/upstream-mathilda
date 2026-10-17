/*
 * functionexpand.c -- FunctionExpand for Mathilda (minimal).
 *
 *   FunctionExpand[expr]
 *   FunctionExpand[expr, assum]
 *
 * expands special functions in expr into simpler ones. In this first version
 * the only rule is MeijerG -> its Slater expansion (src/special_functions/
 * meijer_g.c, meijerg_slater_expand): a sum of HypergeometricPFQ terms in z
 * (p < q, or p == q with the shorter expansion) or in 1/z (p > q), each pFq
 * rewritten to a classical function where one exists. A MeijerG with no such
 * expansion -- the log case of integer-spaced parameters, or p == q with
 * m + n == p, where the two expansions are not continuations of each other --
 * is left as it is, as is everything else. The assumptions argument is
 * accepted for compatibility and not yet used.
 *
 * Memory contract: the builtin takes ownership of res but must not free it.
 */

#include "functionexpand.h"
#include "meijer_g.h"
#include "symtab.h"
#include "attr.h"
#include "eval.h"
#include "expr.h"
#include "sym_names.h"
#include <stdbool.h>
#include <stdlib.h>

/* A copy of e with every expandable MeijerG node replaced, bottom-up; sets
 * *changed when anything was rewritten. Returns an owned Expr. */
static Expr* fe_walk(const Expr* e, bool* changed) {
    if (!e || e->type != EXPR_FUNCTION) return expr_copy((Expr*)e);
    size_t n = e->data.function.arg_count;
    Expr** a = malloc(sizeof(Expr*) * (n ? n : 1));
    bool any = false;
    for (size_t i = 0; i < n; i++) {
        bool c = false;
        a[i] = fe_walk(e->data.function.args[i], &c);
        any |= c;
    }
    Expr* node = expr_new_function(expr_copy(e->data.function.head), a, n);
    free(a);
    if (any) {
        *changed = true;
        node = eval_and_free(node);
    }
    if (node->type == EXPR_FUNCTION && node->data.function.head->type == EXPR_SYMBOL &&
        node->data.function.head->data.symbol.name == SYM_MeijerG) {
        Expr* s = meijerg_slater_expand(node, true);
        if (s) { expr_free(node); node = s; *changed = true; }
    }
    return node;
}

Expr* builtin_function_expand(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t na = res->data.function.arg_count;
    if (na != 1 && na != 2) return NULL;
    bool changed = false;
    Expr* out = fe_walk(res->data.function.args[0], &changed);
    /* Unchanged: hand back the (already evaluated) argument itself. */
    return out;
}

void functionexpand_init(void) {
    symtab_add_builtin("FunctionExpand", builtin_function_expand);
    symtab_get_def("FunctionExpand")->attributes |= ATTR_PROTECTED;
}
