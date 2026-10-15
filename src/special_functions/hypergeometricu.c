/* Mathilda -- HypergeometricU, the Tricomi confluent hypergeometric U(a, b, z).
 *
 *   HypergeometricU[a, b, z]  is the second solution of the confluent
 *   hypergeometric equation z w'' + (b - z) w' - a w = 0, with the branch cut on
 *   z from -infinity to 0 (the z^(1-b) term in its connection to Hypergeometric1F1).
 *
 * Evaluation is layered so each argument shape takes the cheapest accurate route
 * (mirroring Hypergeometric1F1 / WhittakerM):
 *
 *   U(0, b, z) = 1,  U(a, a+1, z) = z^(-a)   -> exact (symbolic or numeric).
 *   a, b, z all machine reals              -> the shared double helper
 *                                             sf_machine_hyperu (bit-for-bit with
 *                                             the NDArray buffer loop below and the
 *                                             Compile[] VM). Declines (z <= 0,
 *                                             cancellation, pole) -> FLINT.
 *   some inexact argument                  -> FLINT's rigorous acb_hypgeom_u via
 *                                             flint_num_hypergeometric_u (MPFR,
 *                                             complex, z <= 0). Kept iff it
 *                                             collapses to a number.
 *   everything else (symbolic)             -> stays symbolic (return NULL).
 *
 * Attributes: Listable, NumericFunction, Protected.
 *
 * Deferred, as the confluent siblings: the general two-branch Series at the
 * origin (z^(1/2 +- ...) / log structure -- U is not analytic at 0) and
 * FunctionExpand. The D[] rule (DLMF 13.3.22) lives in src/calculus/deriv.c.
 */
#include "hypergeometricu.h"

#include <math.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

#include <gmp.h>

#include "arithmetic.h"      /* is_complex */
#include "attr.h"
#include "common.h"          /* builtin_arg_error */
#include "eval.h"            /* eval_and_free */
#include "expr.h"
#include "flint_num_bridge.h"/* flint_num_hypergeometric_u */
#include "ndarray.h"         /* is_ndarray, ndt_get/_set, expr_new_ndarray_like, ... */
#include "numeric.h"         /* numeric_min_inexact_bits */
#include "sf_machine.h"      /* sf_machine_hyperu */
#include "symtab.h"

#ifdef USE_MPFR
#include <mpfr.h>
#endif

/* ------------------------------------------------------------------ */
/* Expression-building shorthands                                      */
/* ------------------------------------------------------------------ */

static Expr* mk_fn2(const char* h, Expr* a, Expr* b) {
    return expr_new_function(expr_new_symbol(h), (Expr*[]){ a, b }, 2);
}

/* ------------------------------------------------------------------ */
/* Argument-shape predicates / extractors (as in whittaker.c)          */
/* ------------------------------------------------------------------ */

static bool hu_is_inexact_leaf(const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_REAL) return true;
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR) return true;
#endif
    return false;
}

static bool hu_is_inexact(const Expr* e) {
    if (hu_is_inexact_leaf(e)) return true;
    Expr *re, *im;
    if (is_complex((Expr*)e, &re, &im))
        return hu_is_inexact_leaf(re) || hu_is_inexact_leaf(im);
    return false;
}

static bool hu_is_number(const Expr* e) {
    if (!e) return false;
    switch (e->type) {
        case EXPR_INTEGER:
        case EXPR_REAL:
        case EXPR_BIGINT:
#ifdef USE_MPFR
        case EXPR_MPFR:
#endif
            return true;
        default: break;
    }
    if (e->type == EXPR_FUNCTION && e->data.function.head &&
        e->data.function.head->type == EXPR_SYMBOL &&
        (strcmp(e->data.function.head->data.symbol.name, "Rational") == 0 ||
         strcmp(e->data.function.head->data.symbol.name, "Complex") == 0)) {
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            if (!hu_is_number(e->data.function.args[i])) return false;
        return true;
    }
    return false;
}

/* Machine-double value of a real/integer leaf only (declines MPFR/BigInt/complex);
 * gates the fast paths, exactly like whittaker.c's whit_leaf_double. */
static bool hu_leaf_double(const Expr* e, double* out) {
    if (!e) return false;
    if (e->type == EXPR_REAL)    { *out = e->data.real; return true; }
    if (e->type == EXPR_INTEGER) { *out = (double)e->data.integer; return true; }
    return false;
}

/* Any real numeric leaf -> double (Integer/Real/BigInt/MPFR). */
static bool hu_num_to_d(const Expr* e, double* out) {
    if (!e) return false;
    switch (e->type) {
        case EXPR_REAL:    *out = e->data.real;              return true;
        case EXPR_INTEGER: *out = (double)e->data.integer;   return true;
        case EXPR_BIGINT:  *out = mpz_get_d(e->data.bigint); return true;
#ifdef USE_MPFR
        case EXPR_MPFR:    *out = mpfr_get_d(e->data.mpfr, MPFR_RNDN); return true;
#endif
        default: return false;
    }
}

/* Demote a numeric FLINT result (real, or Complex of numbers) to machine-precision
 * leaves, so machine-precision inputs give a machine-precision answer rather than a
 * 53-bit MPFR (machine-in / machine-out). Takes ownership of `e`. */
static Expr* hu_demote_machine(Expr* e) {
    double d;
    if (hu_num_to_d(e, &d)) { expr_free(e); return expr_new_real(d); }
    Expr *re, *im;
    if (is_complex(e, &re, &im)) {
        double dr, di;
        if (hu_num_to_d(re, &dr) && hu_num_to_d(im, &di)) {
            expr_free(e);
            return make_complex(expr_new_real(dr), expr_new_real(di));
        }
    }
    return e;
}

/* ------------------------------------------------------------------ */
/* Exact reductions                                                    */
/* ------------------------------------------------------------------ */

/* b == a + 1  (tested by evaluating b - a - 1 to an exact literal 0). */
static bool hu_b_is_a_plus_1(const Expr* a_e, const Expr* b_e) {
    Expr* diff = mk_fn2("Plus", expr_copy((Expr*)b_e),
                        mk_fn2("Plus",
                               mk_fn2("Times", expr_new_integer(-1), expr_copy((Expr*)a_e)),
                               expr_new_integer(-1)));
    Expr* d = eval_and_free(diff);
    bool zero = (d && d->type == EXPR_INTEGER && d->data.integer == 0);
    expr_free(d);
    return zero;
}

/* ------------------------------------------------------------------ */
/* FLINT numeric path                                                  */
/* ------------------------------------------------------------------ */

/* Rigorous MPFR/complex/z<=0 value via acb_hypgeom_u. Returns NULL (stay as the
 * caller left it) when symbolic, a pole, non-finite, or FLINT is absent. */
static Expr* hu_numeric_flint(const Expr* a_e, const Expr* b_e, const Expr* z_e) {
    Expr* out = flint_num_hypergeometric_u(a_e, b_e, z_e);
    if (out && hu_is_number(out)) return out;
    expr_free(out);
    return NULL;
}

/* ------------------------------------------------------------------ */
/* Scalar driver                                                       */
/* ------------------------------------------------------------------ */

static Expr* hypergeometricu_eval(Expr* a_e, Expr* b_e, Expr* z_e) {
    /* Exact reductions (symbolic or numeric). */
    if (a_e->type == EXPR_INTEGER && a_e->data.integer == 0)      /* U(0,b,z) = 1 */
        return expr_new_integer(1);
    if (hu_b_is_a_plus_1(a_e, b_e))                               /* U(a,a+1,z) = z^(-a) */
        return eval_and_free(mk_fn2("Power", expr_copy(z_e),
                                    mk_fn2("Times", expr_new_integer(-1),
                                           expr_copy(a_e))));

    /* Numericise only when some argument is inexact; otherwise stay symbolic. */
    if (!(hu_is_inexact(a_e) || hu_is_inexact(b_e) || hu_is_inexact(z_e)))
        return NULL;

    /* Does the result track machine precision? Yes iff some argument is inexact
     * and every inexact argument is machine (<= 53 bits); an exact Integer/Rational
     * imposes no precision (so HypergeometricU[0.6, 1/3, 1/7.] is machine).
     * numeric_min_inexact_bits lives under USE_MPFR (like every other caller --
     * fourier.c, erfi.c, convolutions.c guard their calls the same way); without
     * MPFR there is no FLINT/arb path to demote below, so machine_prec is moot and
     * the machine double kernel above is the only numeric route. */
#ifdef USE_MPFR
    long gov = 0;
    {
        long bits[3] = { numeric_min_inexact_bits(a_e),
                         numeric_min_inexact_bits(b_e),
                         numeric_min_inexact_bits(z_e) };
        for (int i = 0; i < 3; i++)
            if (bits[i] > 0 && (gov == 0 || bits[i] < gov)) gov = bits[i];
    }
    bool machine_prec = (gov > 0 && gov <= 53);
#else
    bool machine_prec = false;
#endif

    /* Machine real fast path: a, b, z all machine reals. Shares sf_machine_hyperu
     * with the NDArray kernel and the Compile[] VM, so the answer is bit-identical.
     * Declines (z <= 0, pole, cancellation) -> the FLINT path below. */
    {
        double ad, bd, zd;
        if (hu_leaf_double(a_e, &ad) && hu_leaf_double(b_e, &bd) &&
            hu_leaf_double(z_e, &zd)) {
            double v;
            if (sf_machine_hyperu(ad, bd, zd, &v))
                return expr_new_real(v);
        }
    }

    /* Rigorous FLINT path (MPFR precision, complex args, z <= 0, or a machine case
     * the double kernel declined). A machine-precision input must come back
     * machine-precision, so demote the FLINT result in that case. */
    Expr* out = hu_numeric_flint(a_e, b_e, z_e);
    if (out) return machine_prec ? hu_demote_machine(out) : out;
    return NULL;                               /* stay symbolic */
}

/* ------------------------------------------------------------------ */
/* NDArray fast path (machine buffer argument z)                       */
/* ------------------------------------------------------------------ */

/* Thread sf_machine_hyperu -- the SAME helper the scalar machine path and
 * Compile[] use -- over a visible NDArray / packed-List z, holding the scalar
 * parameters a and b constant. Once z IS an NDArray this must NOT return NULL
 * (the scalar driver cannot consume a buffer): it returns a packed result or
 * ndarray_delist_and_reeval(res) so each element takes the accurate scalar route.
 * Mirrors whittaker_m_ndarray. */
static Expr* hypergeometricu_ndarray(const Expr* a_e, const Expr* b_e,
                                     const Expr* z_e, Expr* res) {
    if (!is_ndarray(z_e)) return NULL;             /* scalar driver handles */
    NDType dt = z_e->data.ndarray.dtype;
    double ad, bd;
    if (dt == NDT_INT64 || dt == NDT_BOOL || ndt_is_complex(dt) ||
        !hu_leaf_double(a_e, &ad) || !hu_leaf_double(b_e, &bd))
        return ndarray_delist_and_reeval(res);

    size_t sz = ndarray_size(z_e);
    void* out = malloc(ndt_elem_size(dt) * (sz ? sz : 1));
    if (!out) return ndarray_delist_and_reeval(res);
    bool ok = true;
    for (size_t i = 0; i < sz && ok; i++) {
        double xr, xi, v = 0.0;
        ndt_get(z_e->data.ndarray.data, i, dt, &xr, &xi);
        ok = sf_machine_hyperu(ad, bd, xr, &v);
        if (ok) ndt_set(out, i, dt, v, 0.0);
    }
    if (ok)
        return expr_new_ndarray_like(z_e, z_e->data.ndarray.rank,
                                     z_e->data.ndarray.dims, out, dt);
    free(out);
    return ndarray_delist_and_reeval(res);
}

/* ------------------------------------------------------------------ */
/* Builtin entry point                                                 */
/* ------------------------------------------------------------------ */

Expr* builtin_hypergeometric_u(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    Expr** args = res->data.function.args;
    if (argc != 3)
        return builtin_arg_error("HypergeometricU", argc, 3, 3);
    Expr* nd = hypergeometricu_ndarray(args[0], args[1], args[2], res);
    if (nd) return nd;
    return hypergeometricu_eval(args[0], args[1], args[2]);
}

/* ------------------------------------------------------------------ */
/* Registration                                                        */
/* ------------------------------------------------------------------ */

void hypergeometricu_init(void) {
    symtab_add_builtin("HypergeometricU", builtin_hypergeometric_u);
    symtab_get_def("HypergeometricU")->attributes |=
        (ATTR_LISTABLE | ATTR_NUMERICFUNCTION | ATTR_PROTECTED);
    /* Docstring lives in info.c (info_init). */
}
