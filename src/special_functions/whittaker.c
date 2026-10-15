/* Mathilda -- WhittakerM, the Whittaker M function M_{k,m}(z).
 *
 *   WhittakerM[k, m, z] = e^(-z/2) z^(m+1/2) Hypergeometric1F1[m-k+1/2, 2m+1, z].
 *
 * Evaluation is layered so each argument shape takes the cheapest numeric route
 * (mirroring the LaguerreL / HermiteH confluent-hypergeometric design):
 *
 *   z == 0 (exact), Re(m) > -1/2   ->  exact 0 (M_{k,m} ~ z^(m+1/2) near 0).
 *   k, m, z all machine reals      ->  the shared double helper
 *                                        sf_machine_whittaker_m
 *                                      (bit-for-bit with the NDArray element-wise
 *                                       buffer loop below and the Compile[] VM).
 *                                      Declines (z <= 0, 2m+1 pole, 1F1
 *                                       cancellation, overflow) -> template.
 *   some inexact argument          ->  the confluent-hypergeometric form built as
 *                                        e^(-z/2) z^(m+1/2) Hypergeometric1F1[...],
 *                                      evaluated so it inherits 1F1's machine +
 *                                      MPFR, real + complex numerics. Kept only
 *                                      if it collapses to a number.
 *   everything else (symbolic)     ->  stays symbolic (return NULL).
 *
 * WhittakerM is a single growing PRODUCT (prefactor * 1F1), not a difference, so
 * it adds no cancellation over 1F1 itself and needs no guard digits (unlike
 * HermiteH): bad 1F1 cancellation is declined inside sf_machine_pfq / re-summed
 * in MPFR by the 1F1 evaluator.
 *
 * Attributes: Listable, NumericFunction, Protected.
 *
 * Deferred / best-effort, as the sibling special functions: tight Interval
 * enclosures (the generic NumericFunction + D[] interval threading in eval.c
 * declines here because D[WhittakerM,z] mentions WhittakerM itself -- the
 * self-referential-derivative guard), the general asymptotic Series at infinity
 * (a two-scale Stokes transseries; only the terminating/elementary case is
 * emitted, in calculus/series.c), symbolic FunctionExpand, and general-n
 * SeriesCoefficient. The D[] rule and the origin Series live in
 * src/calculus/deriv.c and src/calculus/series.c respectively.
 */
#include "whittaker.h"

#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include <gmp.h>

#include "arithmetic.h"      /* is_complex */
#include "attr.h"
#include "common.h"          /* builtin_arg_error */
#include "eval.h"            /* eval_and_free */
#include "expr.h"
#include "ndarray.h"         /* is_ndarray, ndt_get/_set, expr_new_ndarray_like, ... */
#include "sf_machine.h"      /* sf_machine_whittaker_m */
#include "symtab.h"

#ifdef USE_MPFR
#include <mpfr.h>
#endif

/* ------------------------------------------------------------------ */
/* Expression-building shorthands                                      */
/* ------------------------------------------------------------------ */

static Expr* mk_fn1(const char* h, Expr* a) {
    return expr_new_function(expr_new_symbol(h), (Expr*[]){ a }, 1);
}
static Expr* mk_fn2(const char* h, Expr* a, Expr* b) {
    return expr_new_function(expr_new_symbol(h), (Expr*[]){ a, b }, 2);
}
static Expr* mk_fn3(const char* h, Expr* a, Expr* b, Expr* c) {
    return expr_new_function(expr_new_symbol(h), (Expr*[]){ a, b, c }, 3);
}
/* 1/2 and -1/2 as exact rationals. */
static Expr* mk_half(void)     { return mk_fn2("Rational", expr_new_integer(1),  expr_new_integer(2)); }
static Expr* mk_neg_half(void) { return mk_fn2("Rational", expr_new_integer(-1), expr_new_integer(2)); }

/* ------------------------------------------------------------------ */
/* Argument-shape predicates / extractors                             */
/* ------------------------------------------------------------------ */

/* True if `e` is an inexact numeric leaf (Real or MPFR). */
static bool whit_is_inexact_leaf(const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_REAL) return true;
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR) return true;
#endif
    return false;
}

/* True if `e` carries an inexact numeric leaf (also inside Complex[..]). */
static bool whit_is_inexact(const Expr* e) {
    if (whit_is_inexact_leaf(e)) return true;
    Expr *re, *im;
    if (is_complex((Expr*)e, &re, &im))
        return whit_is_inexact_leaf(re) || whit_is_inexact_leaf(im);
    return false;
}

/* True iff `e` is a pure numeric quantity (Integer/Real/MPFR/BigInt, or
 * Rational/Complex of such) -- i.e. the result actually collapsed to a number. */
static bool whit_is_number(const Expr* e) {
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
            if (!whit_is_number(e->data.function.args[i])) return false;
        return true;
    }
    return false;
}

/* Machine-double value of a real/integer leaf (declines MPFR/BigInt/complex). */
static bool whit_leaf_double(const Expr* e, double* out) {
    if (!e) return false;
    if (e->type == EXPR_REAL)    { *out = e->data.real; return true; }
    if (e->type == EXPR_INTEGER) { *out = (double)e->data.integer; return true; }
    return false;
}

/* Value of a real numeric leaf as a double (Integer/Real/BigInt/MPFR/Rational). */
static bool whit_real_number_double(const Expr* e, double* out) {
    if (!e) return false;
    switch (e->type) {
        case EXPR_INTEGER: *out = (double)e->data.integer; return true;
        case EXPR_REAL:    *out = e->data.real;            return true;
        case EXPR_BIGINT:  *out = mpz_get_d(e->data.bigint); return true;
#ifdef USE_MPFR
        case EXPR_MPFR:    *out = mpfr_get_d(e->data.mpfr, MPFR_RNDN); return true;
#endif
        default: break;
    }
    if (e->type == EXPR_FUNCTION && e->data.function.head &&
        e->data.function.head->type == EXPR_SYMBOL &&
        strcmp(e->data.function.head->data.symbol.name, "Rational") == 0 &&
        e->data.function.arg_count == 2) {
        double p, q;
        if (whit_real_number_double(e->data.function.args[0], &p) &&
            whit_real_number_double(e->data.function.args[1], &q) && q != 0.0) {
            *out = p / q;
            return true;
        }
    }
    return false;
}

/* Real part of a numeric `e` as a double (for the z == 0 special value).
 * Declines a symbolic or non-numeric argument. */
static bool whit_real_part_double(const Expr* e, double* out) {
    Expr *re, *im;
    if (is_complex((Expr*)e, &re, &im)) return whit_real_number_double(re, out);
    return whit_real_number_double(e, out);
}

/* ------------------------------------------------------------------ */
/* Confluent-hypergeometric template                                   */
/* ------------------------------------------------------------------ */

/* Build e^(-z/2) z^(m+1/2) Hypergeometric1F1[m-k+1/2, 2m+1, z] as an expression
 * and evaluate it, so it inherits 1F1's machine + MPFR, real + complex numerics
 * (and returns a Complex result for z < 0 or complex parameters). Returns NULL
 * if the result did not collapse to a number (symbolic k / m). */
static Expr* whittaker_m_numeric_hyp(const Expr* k_e, const Expr* m_e, const Expr* z_e) {
    /* a = m - k + 1/2 */
    Expr* a = mk_fn3("Plus",
                     expr_copy((Expr*)m_e),
                     mk_fn2("Times", expr_new_integer(-1), expr_copy((Expr*)k_e)),
                     mk_half());
    /* b = 2 m + 1 */
    Expr* b = mk_fn2("Plus",
                     mk_fn2("Times", expr_new_integer(2), expr_copy((Expr*)m_e)),
                     expr_new_integer(1));
    Expr* hyp = mk_fn3("Hypergeometric1F1", a, b, expr_copy((Expr*)z_e));

    /* e^(-z/2) */
    Expr* expfac = mk_fn1("Exp", mk_fn2("Times", mk_neg_half(), expr_copy((Expr*)z_e)));
    /* z^(m+1/2) */
    Expr* zpow = mk_fn2("Power", expr_copy((Expr*)z_e),
                        mk_fn2("Plus", expr_copy((Expr*)m_e), mk_half()));

    Expr* prod = mk_fn3("Times", expfac, zpow, hyp);
    Expr* out = eval_and_free(prod);
    if (out && whit_is_number(out)) return out;
    expr_free(out);
    return NULL;
}

/* ------------------------------------------------------------------ */
/* Scalar driver                                                       */
/* ------------------------------------------------------------------ */

static Expr* whittaker_m_eval(Expr* k_e, Expr* m_e, Expr* z_e) {
    /* Special value: M_{k,m}(0) = 0 for Re(m) > -1/2 (leading power z^(m+1/2)).
     * Only the exact integer 0; a Real 0.0 flows through the template to 0.0. */
    if (z_e && z_e->type == EXPR_INTEGER && z_e->data.integer == 0) {
        double rem;
        if (whit_real_part_double(m_e, &rem) && rem > -0.5)
            return expr_new_integer(0);
        /* symbolic m, or Re(m) <= -1/2 (singular): fall through / stay symbolic */
    }

    /* Numericise only when some argument is inexact; otherwise stay symbolic
     * (so WhittakerM[2,3,x] and all-exact WhittakerM[2,3,4] do not expand). */
    if (!(whit_is_inexact(k_e) || whit_is_inexact(m_e) || whit_is_inexact(z_e)))
        return NULL;

    /* Machine real fast path: k, m, z all machine reals (precision <= 53). Shares
     * sf_machine_whittaker_m with the NDArray kernel and the Compile[] VM, so the
     * answer is bit-identical. Declines (z <= 0, pole, cancellation, overflow) ->
     * the MPFR / complex template below. */
    {
        double kd, md, zd;
        if (whit_leaf_double(k_e, &kd) && whit_leaf_double(m_e, &md) &&
            whit_leaf_double(z_e, &zd)) {
            double v;
            if (sf_machine_whittaker_m(kd, md, zd, &v))
                return expr_new_real(v);
        }
    }

    /* General confluent-hypergeometric template (MPFR precision, complex args,
     * z <= 0, exact-rational k/m with an inexact z). */
    Expr* out = whittaker_m_numeric_hyp(k_e, m_e, z_e);
    if (out) return out;
    return NULL;                               /* stay symbolic */
}

/* ------------------------------------------------------------------ */
/* NDArray fast path (machine buffer argument z)                       */
/* ------------------------------------------------------------------ */

/* Thread sf_machine_whittaker_m -- the SAME helper the scalar machine path and
 * Compile[] use -- over a visible NDArray / packed-List argument `z_e`, holding
 * the scalar parameters k and m constant. This is what makes WhittakerM fast
 * over arrays: a 3-arg call has no arity>=3 REPL element-wise kernel, so without
 * this the buffer would be materialised into boxed Exprs and threaded one scalar
 * at a time.
 *
 * Because the scalar driver (whittaker_m_eval) cannot consume a buffer, once
 * `z_e` IS an NDArray this must NOT return NULL -- it returns a packed result or
 * ndarray_delist_and_reeval(res) (which rebuilds the call with the buffer turned
 * into a nested List and re-evaluates, so Listable threads each element through
 * the accurate scalar / complex route). Declines mid-loop (z <= 0, pole,
 * overflow) take that same delist route so every element is computed correctly.
 */
static Expr* whittaker_m_ndarray(const Expr* k_e, const Expr* m_e,
                                 const Expr* z_e, Expr* res) {
    if (!is_ndarray(z_e)) return NULL;             /* scalar driver handles */
    NDType dt = z_e->data.ndarray.dtype;
    double kd, md;
    /* Real float buffer + machine real k, m only; everything else (int64/bool
     * buffer, complex buffer, non-machine k/m) takes the per-element route. */
    if (dt == NDT_INT64 || dt == NDT_BOOL || ndt_is_complex(dt) ||
        !whit_leaf_double(k_e, &kd) || !whit_leaf_double(m_e, &md))
        return ndarray_delist_and_reeval(res);

    size_t sz = ndarray_size(z_e);
    void* out = malloc(ndt_elem_size(dt) * (sz ? sz : 1));
    if (!out) return ndarray_delist_and_reeval(res);
    bool ok = true;
    for (size_t i = 0; i < sz && ok; i++) {
        double xr, xi, v = 0.0;
        ndt_get(z_e->data.ndarray.data, i, dt, &xr, &xi);
        ok = sf_machine_whittaker_m(kd, md, xr, &v);
        if (ok) ndt_set(out, i, dt, v, 0.0);
    }
    if (ok)
        return expr_new_ndarray_like(z_e, z_e->data.ndarray.rank,
                                     z_e->data.ndarray.dims, out, dt);
    /* A declined element must not reach the scalar driver, which cannot consume
     * the buffer; delist so each element takes the accurate scalar route. */
    free(out);
    return ndarray_delist_and_reeval(res);
}

/* ------------------------------------------------------------------ */
/* Builtin entry point                                                 */
/* ------------------------------------------------------------------ */

Expr* builtin_whittaker_m(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    Expr** args = res->data.function.args;
    if (argc != 3)
        return builtin_arg_error("WhittakerM", argc, 3, 3);
    Expr* nd = whittaker_m_ndarray(args[0], args[1], args[2], res);
    if (nd) return nd;
    return whittaker_m_eval(args[0], args[1], args[2]);
}

/* ------------------------------------------------------------------ */
/* Registration                                                        */
/* ------------------------------------------------------------------ */

void whittaker_init(void) {
    symtab_add_builtin("WhittakerM", builtin_whittaker_m);
    symtab_get_def("WhittakerM")->attributes |=
        (ATTR_LISTABLE | ATTR_NUMERICFUNCTION | ATTR_PROTECTED);
    /* Docstring lives in info.c (info_init). */
}
