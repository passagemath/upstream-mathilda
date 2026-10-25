/*
 * minus.c -- Minus[x], the functional form of unary negation.
 *
 * The parser already reads `-x` as Times[-1, x], so Minus only has to exist as
 * a callable head: Minus[x] rewrites to Times[-1, x] and lets Times do the
 * arithmetic (numbers, Rationals, Complex, Plus distribution, packed arrays).
 * This is what makes SortBy[list, Minus] and KeySortBy[a, Minus] work.
 *
 * Attributes match Mathematica: Listable, NumericFunction, Protected.  Any
 * other argument count emits Minus::argx and stays unevaluated (WL prints
 * Minus[x, y] as x − y but does not evaluate it).
 */
#include "expr.h"
#include "symtab.h"
#include "attr.h"
#include "common.h"
#include "sym_names.h"

Expr* builtin_minus(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc != 1) return builtin_arg_error("Minus", argc, 1, 1);
    Expr* targs[2] = { expr_new_integer(-1), expr_copy(res->data.function.args[0]) };
    return expr_new_function(expr_new_symbol(SYM_Times), targs, 2);
}

/* Internal`SyntacticNegativeQ[e] -- does `e` carry a minus sign in the LEADING
 * position of its printed form? A *syntactic* question, not a sign test, read
 * as the reference reads it:
 *
 *   -3, -2.5, -1/2, -x, -2 x y, -I, 1 - I, -a - b, -1 + x    True
 *   3, x, x - y, 1 - 2 x, -1 + I, (-x)^2                      False
 *
 * Times and Plus defer to their first argument in canonical order (the
 * numeric coefficient of a product, the leading term of a sum); a Complex is
 * decided by its imaginary part alone (Complex[-1, 2] is False). */
static bool syntactic_negative(const Expr* e) {
    if (!e) return false;
    switch (e->type) {
        case EXPR_INTEGER: return e->data.integer < 0;
        case EXPR_REAL:    return e->data.real < 0.0;
        case EXPR_BIGINT:  return mpz_sgn(e->data.bigint) < 0;
        case EXPR_FUNCTION: break;
        default: return false;
    }
    if (e->data.function.head->type != EXPR_SYMBOL) return false;
    const char* h = e->data.function.head->data.symbol.name;
    size_t argc  = e->data.function.arg_count;
    if (h == SYM_Rational && argc == 2)
        return syntactic_negative(e->data.function.args[0]);   /* denominator > 0 */
    /* A complex number by the sign of its imaginary part alone, a sum by its
     * first term in canonical order -- the reference's reading:
     * Complex[-1, 2] is False and Complex[1, -2] True; -x - y and -1 + x are
     * True, x - y is False. (FeynCalc's SUNSimplify rewrites
     * colorFactor[a] /; SyntacticNegativeQ[a] to -colorFactor[ExpandAll[-a]].) */
    if (h == SYM_Complex && argc == 2)
        return syntactic_negative(e->data.function.args[1]);
    if ((h == SYM_Times || h == SYM_Plus) && argc >= 1)
        return syntactic_negative(e->data.function.args[0]);
    /* A literal Minus[x] head needs no clause: the parser reads -x as
     * Times[-1, x], and builtin_minus rewrites any explicit Minus[x] to the
     * same, so nothing reaches here still wearing that head. */
    return false;
}

static Expr* builtin_syntacticnegativeq(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc != 1) return builtin_arg_error("Internal`SyntacticNegativeQ", argc, 1, 1);
    return expr_new_symbol(syntactic_negative(res->data.function.args[0])
                           ? SYM_True : SYM_False);
}

void minus_init(void) {
    symtab_add_builtin("Minus", builtin_minus);
    symtab_get_def("Minus")->attributes |= ATTR_LISTABLE | ATTR_NUMERICFUNCTION | ATTR_PROTECTED;
    symtab_set_docstring("Minus",
        "Minus[x] is the arithmetic negation of x, equivalent to -x (Times[-1, x]).\n"
        "\tListable; Minus[x, y] (any count other than one) is left unevaluated with Minus::argx.");

    symtab_add_builtin("Internal`SyntacticNegativeQ", builtin_syntacticnegativeq);
    symtab_get_def("Internal`SyntacticNegativeQ")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("Internal`SyntacticNegativeQ",
        "Internal`SyntacticNegativeQ[expr] gives True if expr carries a minus sign "
        "in the leading position of its printed form: a negative number or Rational, "
        "a Complex with a negative imaginary part, or a Times or Plus whose first "
        "argument is one of those (-a - b and -1 + x give True, x - y False).");
}
