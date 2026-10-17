/* Mathilda -- the Struve function H_nu(z).
 *
 *   StruveH[nu, z]   H_nu(z) = (z/2)^(nu+1) (2/(Sqrt[Pi] Gamma[nu+3/2]))
 *                              1F2(1; 3/2, nu+3/2; -z^2/4)
 *
 * Evaluation:
 *   StruveH[nu, 0] = 0            exact, for a real nu > -1
 *   numeric nu and inexact z      via the 1F2 representation (reusing the tested
 *                                 HypergeometricPFQ numeric evaluator)
 *   everything else               stays symbolic (return NULL)
 *
 * The 1F2 identity is also the reduction used by the half-line Mellin transform
 * (reduce_to_hypergeometric in integrate_ramanujan.c), so Integrate[x^(s-1)
 * StruveH[0, a x], {x,0,Infinity}] closes through rec_pfq.
 *
 * Attributes: Listable, NumericFunction, Protected.
 */
#include "struveh.h"

#include "symtab.h"
#include "attr.h"
#include "eval.h"
#include "parse.h"
#include "sym_names.h"
#include "common.h"
#include "numeric.h"
#include "arithmetic.h"         /* is_complex */
#include "flint_num_bridge.h"   /* flint_num_struveh (arb, arbitrary precision) */
#include "sf_machine.h"         /* sf_machine_struve_h (shared double kernel)   */

#include <stdbool.h>

/* Small expr builders (local idiom). */
static Expr* mk_int(int64_t v) { return expr_new_integer(v); }
static Expr* mk_fn1(const char* name, Expr* a) {
    return expr_new_function(expr_new_symbol(name), (Expr*[]){ a }, 1);
}
static Expr* mk_fn2(const char* name, Expr* a, Expr* b) {
    return expr_new_function(expr_new_symbol(name), (Expr*[]){ a, b }, 2);
}
/* evaluate `call` to a fixed point, then free it (evaluate borrows its arg). */
static Expr* sh_eval(Expr* call) {
    Expr* r = evaluate(call);
    expr_free(call);
    return r;
}

/* H_nu(z) via its 1F2 representation, numericised; NULL if the result is not a
 * number (e.g. the pFq evaluator declined). */
static Expr* struveh_numeric(const Expr* nu, const Expr* z) {
    static const char* TMPL =
        "(zz/2)^(nn+1) (2/(Sqrt[Pi] Gamma[nn+3/2])) "
        "HypergeometricPFQ[{1}, {3/2, nn+3/2}, -zz^2/4]";
    Expr* tmpl = parse_expression(TMPL);
    if (!tmpl) return NULL;
    Expr* rules = expr_new_function(expr_new_symbol("List"), (Expr*[]){
        mk_fn2("Rule", expr_new_symbol("nn"), expr_copy((Expr*)nu)),
        mk_fn2("Rule", expr_new_symbol("zz"), expr_copy((Expr*)z)) }, 2);
    Expr* subbed = sh_eval(mk_fn2("ReplaceAll", tmpl, rules));
    Expr* r = sh_eval(mk_fn1("N", subbed));
    if (r && expr_is_numeric_quantity(r)) return r;
    if (r) expr_free(r);
    return NULL;
}

/* Minimum precision (bits) of e's inexact components (machine Real / MPFR,
 * or Complex of those); 0 when e is exact, i.e. N[] has not been applied. */
static long sh_inexact_bits(const Expr* e) {
#ifdef USE_MPFR
    return numeric_min_inexact_bits(e);
#else
    /* numeric_min_inexact_bits is MPFR-only; without MPFR the only inexact
     * leaf is a machine Real (53 bits), possibly inside Complex[...]. */
    Expr *re, *im;
    if (e && e->type == EXPR_REAL) return 53;
    if (e && is_complex((Expr*)e, &re, &im)
            && (re->type == EXPR_REAL || im->type == EXPR_REAL)) return 53;
    return 0;
#endif
}
static bool sh_inexact(const Expr* e) { return sh_inexact_bits(e) > 0; }

Expr* builtin_struveh(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) return NULL;
    Expr* nu = res->data.function.args[0];
    Expr* z  = res->data.function.args[1];

    /* H_nu(0) = 0 for a real nu > -1 (exact). */
    if (z->type == EXPR_INTEGER && z->data.integer == 0) {
        double nuv;
        if (common_machine_real_value(nu, &nuv) && nuv > -1.0) return mk_int(0);
    }

    /* Numeric path: at least one inexact argument, both concrete numbers. */
    if ((sh_inexact(nu) || sh_inexact(z)) &&
        expr_is_numeric_quantity(nu) && expr_is_numeric_quantity(z)) {
        /* Target precision: the min inexact-bit count over the arguments
         * (an exact arg imposes no constraint), floored implicitly at 53. */
        long bnu = sh_inexact_bits(nu);
        long bz  = sh_inexact_bits(z);
        long outb = bnu;
        if (bz > 0 && (outb == 0 || bz < outb)) outb = bz;

        if (outb > 53) {
            /* Arbitrary precision: the rigorous arb kernel (all nu, all C).
             * On a genuine decline (non-finite) leave it symbolic rather than
             * silently downgrade to machine precision. */
            Expr* v = flint_num_struveh(nu, z);
            if (v) return v;
        } else {
            /* Machine precision: the shared double kernel for real nu, z.
             * Bit-identical to the NDArray / Compile surfaces. */
            double nuv, zv;
            if (common_machine_real_value(nu, &nuv) &&
                common_machine_real_value(z, &zv)) {
                double h;
                if (sf_machine_struve_h(nuv, zv, &h)) return expr_new_real(h);
            }
            /* Kernel declined (the hard middle band, or a complex argument):
             * fall back to the symbolic 1F2 template for a real z. */
            if (z->type == EXPR_REAL) {
                double nuv2;
                if (common_machine_real_value(nu, &nuv2)) {
                    Expr* v = struveh_numeric(nu, z);
                    if (v) return v;
                }
            }
        }
    }
    return NULL;
}

void struveh_init(void) {
    symtab_add_builtin("StruveH", builtin_struveh);
    symtab_get_def("StruveH")->attributes |=
        (ATTR_LISTABLE | ATTR_NUMERICFUNCTION | ATTR_PROTECTED);
    /* Docstring lives in info.c (info_init). */
}
