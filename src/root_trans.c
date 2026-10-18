/*
 * root_trans.c
 *
 * Transcendental Root objects Root[{f, x0}] and certified real-root isolation
 * of exact univariate expressions.  See root_trans.h for the contract.
 *
 * All arithmetic is arb ball arithmetic through flint_cert_eval_real (the
 * rigorous evaluator behind flint_cert_parts_sign, with the variable bound to
 * a real ball), so every exclusion / monotonicity / existence claim below is a
 * proof, not a numerical estimate.  The only non-arb decision is the
 * recognition of an EXACT rational root, which requires the evaluator to
 * reduce f(q) to the literal integer 0.
 */
#include "root_trans.h"

#include "eval.h"
#include "symtab.h"
#include "attr.h"
#include "sym_names.h"
#include "message.h"

#include <stdlib.h>
#include <string.h>
#include <math.h>

/* ------------------------------------------------------------------ *
 *  Structural helpers (no FLINT needed)                               *
 * ------------------------------------------------------------------ */

static bool rt_is_head(const Expr* e, const char* name) {
    return e && e->type == EXPR_FUNCTION
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == name;
}

static bool rt_is_slot1(const Expr* e) {
    return rt_is_head(e, SYM_Slot) && e->data.function.arg_count == 1
        && e->data.function.args[0]->type == EXPR_INTEGER
        && e->data.function.args[0]->data.integer == 1;
}

static bool rt_is_real_number(const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_INTEGER || e->type == EXPR_REAL || e->type == EXPR_BIGINT)
        return true;
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR) return true;
#endif
    return rt_is_head(e, SYM_Rational);
}

/* Root[{Function[body], x0}]: borrowed body and x0. */
static bool rt_parse(const Expr* e, const Expr** body, const Expr** x0) {
    if (!rt_is_head(e, SYM_Root) || e->data.function.arg_count != 1) return false;
    const Expr* l = e->data.function.args[0];
    if (!rt_is_head(l, SYM_List) || l->data.function.arg_count != 2) return false;
    const Expr* fn = l->data.function.args[0];
    if (!rt_is_head(fn, SYM_Function) || fn->data.function.arg_count != 1) return false;
    if (!rt_is_real_number(l->data.function.args[1])) return false;
    *body = fn->data.function.args[0];
    *x0 = l->data.function.args[1];
    return true;
}

/* Every leaf is a number, a NumericQ constant symbol or #1, and every head a
 * NumericFunction: the body is a function of #1 alone. */
static bool rt_body_ok(const Expr* e, int depth) {
    if (depth > 512) return false;
    if (rt_is_real_number(e)) return true;
    if (rt_is_slot1(e)) return true;
    if (e->type == EXPR_SYMBOL) return expr_is_numeric_quantity(e);
    if (e->type != EXPR_FUNCTION || e->data.function.head->type != EXPR_SYMBOL)
        return false;
    const char* h = e->data.function.head->data.symbol.name;
    if (h == SYM_Complex) return true;
    SymbolDef* def = symtab_get_def(h);
    if (!def || !(def->attributes & ATTR_NUMERICFUNCTION)) return false;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (!rt_body_ok(e->data.function.args[i], depth + 1)) return false;
    return true;
}

bool root_trans_is_object(const Expr* e) {
    const Expr *body, *x0;
    return rt_parse(e, &body, &x0) && rt_body_ok(body, 0);
}

#if defined(USE_FLINT) && defined(USE_MPFR)
/* Structural copy of e with every Slot[1] replaced by `repl` (not descending
 * into a nested Function, whose slots are its own). */
static Expr* rt_subst_slot(const Expr* e, const Expr* repl) {
    if (rt_is_slot1(e)) return expr_copy((Expr*)repl);
    if (e->type != EXPR_FUNCTION || rt_is_head(e, SYM_Function))
        return expr_copy((Expr*)e);
    size_t n = e->data.function.arg_count;
    Expr** a = malloc((n ? n : 1) * sizeof(Expr*));
    for (size_t i = 0; i < n; i++) a[i] = rt_subst_slot(e->data.function.args[i], repl);
    Expr* r = expr_new_function(expr_copy(e->data.function.head), a, n);
    free(a);
    return r;
}
#endif

/* Structural copy of e with the symbol `name` replaced by `repl`. */
static Expr* rt_subst_sym(const Expr* e, const char* name, const Expr* repl) {
    if (e->type == EXPR_SYMBOL && e->data.symbol.name == name) return expr_copy((Expr*)repl);
    if (e->type != EXPR_FUNCTION) return expr_copy((Expr*)e);
    size_t n = e->data.function.arg_count;
    Expr** a = malloc((n ? n : 1) * sizeof(Expr*));
    for (size_t i = 0; i < n; i++) a[i] = rt_subst_sym(e->data.function.args[i], name, repl);
    Expr* r = expr_new_function(rt_subst_sym(e->data.function.head, name, repl), a, n);
    free(a);
    return r;
}

static Expr* rt_call1(const char* head, Expr* a) {
    return eval_and_free(expr_new_function(expr_new_symbol(head), (Expr*[]){ a }, 1));
}

Expr* root_trans_make(const Expr* f, const Expr* x, double approx) {
    /* Present the function as Mathematica does: the numerator of the combined
     * fraction (x/2 - Sin[x] -> x - 2 Sin[x]); it has the same simple roots
     * where f is finite. */
    Expr* g = rt_call1(SYM_Numerator, rt_call1(SYM_Together, expr_copy((Expr*)f)));
    /* Sign convention: the top power of x has a positive coefficient
     * (Cos[x] - x -> x - Cos[x]), so the same root prints the same way
     * whichever side of the relation it came from. */
    Expr* k = eval_and_free(expr_new_function(expr_new_symbol(SYM_Exponent),
        (Expr*[]){ expr_copy(g), expr_copy((Expr*)x) }, 2));
    if (k->type == EXPR_INTEGER && k->data.integer > 0) {
        Expr* c = eval_and_free(expr_new_function(expr_new_symbol(SYM_Coefficient),
            (Expr*[]){ expr_copy(g), expr_copy((Expr*)x), expr_copy(k) }, 3));
        Expr* neg = rt_call1(SYM_Negative, c);
        if (neg->type == EXPR_SYMBOL && neg->data.symbol.name == SYM_True)
            g = rt_call1(SYM_Expand, eval_and_free(expr_new_function(
                    expr_new_symbol(SYM_Times), (Expr*[]){ expr_new_integer(-1), g }, 2)));
        expr_free(neg);
    }
    expr_free(k);
    Expr* slot = expr_new_function(expr_new_symbol(SYM_Slot),
                                   (Expr*[]){ expr_new_integer(1) }, 1);
    Expr* body = rt_subst_sym(g, x->data.symbol.name, slot);
    expr_free(slot);
    expr_free(g);
    Expr* fn = expr_new_function(expr_new_symbol(SYM_Function), (Expr*[]){ body }, 1);
    Expr* lst = expr_new_function(expr_new_symbol(SYM_List),
                                  (Expr*[]){ fn, expr_new_real(approx) }, 2);
    return expr_new_function(expr_new_symbol(SYM_Root), (Expr*[]){ lst }, 1);
}

#if defined(USE_FLINT) && defined(USE_MPFR)

#include "flint_cert.h"
#include <gmp.h>
#include <mpfr.h>
#include <flint/fmpz.h>
#include <flint/fmpq.h>
#include <flint/arf.h>
#include <flint/arb.h>

#define RT_PREC     128     /* isolation working precision (bits)            */
#define RT_MAXD     12      /* highest derivative tried for a multiple root  */
#define RT_RES_BITS 40      /* isolation resolution: 2^-40 relative width    */
#define RT_WIN_BITS 24      /* Root[{f, x0}] certification window, relative  */

/* ------------------------------------------------------------------ *
 *  Number conversions                                                 *
 * ------------------------------------------------------------------ */

static Expr* rt_fmpz_expr(const fmpz_t z) {
    if (fmpz_fits_si(z)) return expr_new_integer((int64_t)fmpz_get_si(z));
    mpz_t m; mpz_init(m);
    fmpz_get_mpz(m, z);
    Expr* r = expr_new_bigint_from_mpz(m);
    mpz_clear(m);
    return r;
}

static Expr* rt_fmpq_expr(const fmpq_t q) {
    Expr* num = rt_fmpz_expr(fmpq_numref(q));
    if (fmpz_is_one(fmpq_denref(q))) return num;
    Expr* den = rt_fmpz_expr(fmpq_denref(q));
    return eval_and_free(expr_new_function(expr_new_symbol(SYM_Times),
        (Expr*[]){ num, expr_new_function(expr_new_symbol(SYM_Power),
                       (Expr*[]){ den, expr_new_integer(-1) }, 2) }, 2));
}

/* A real number Expr (Integer / BigInt / Rational / Real / MPFR) -> arf. */
static bool rt_number_arf(const Expr* e, arf_t out) {
    arb_t b; arb_init(b);
    bool ok = false;
    if (e->type == EXPR_REAL) { arf_set_d(out, e->data.real); ok = isfinite(e->data.real); }
    else if (e->type == EXPR_MPFR) { arf_set_mpfr(out, e->data.mpfr); ok = mpfr_number_p(e->data.mpfr); }
    else if (flint_cert_eval_real(e, NULL, NULL, b, 256)) { arf_set(out, arb_midref(b)); ok = true; }
    arb_clear(b);
    return ok;
}

/* ------------------------------------------------------------------ *
 *  Point signs                                                        *
 * ------------------------------------------------------------------ */

/* Exact test f(q) == 0: the evaluator must reduce it to the literal 0. */
static bool rt_exact_zero(const Expr* f, const char* v, const fmpq_t q) {
    Expr* qe = rt_fmpq_expr(q);
    Expr* val = eval_and_free(rt_subst_sym(f, v, qe));
    expr_free(qe);
    bool z = val->type == EXPR_INTEGER && val->data.integer == 0;
    expr_free(val);
    return z;
}

/* Sign of g at the rational q: -1 / +1 certified by arb, 0 exactly zero
 * (exact evaluation), 2 undecided or outside the supported language. */
static int rt_point_sign(const Expr* g, const char* v, const fmpq_t q) {
    static const slong ladder[] = { 128, 512, 2048 };
    arb_t X, y; arb_init(X); arb_init(y);
    int s = 2;
    bool any = false;
    for (size_t i = 0; i < sizeof ladder / sizeof ladder[0]; i++) {
        arb_set_fmpq(X, q, ladder[i]);
        if (!flint_cert_eval_real(g, v, X, y, ladder[i])) continue;
        any = true;
        if (arb_is_positive(y)) { s = 1; break; }
        if (arb_is_negative(y)) { s = -1; break; }
        if (arb_is_zero(y)) break;
    }
    arb_clear(X); arb_clear(y);
    if (s == 2 && any && rt_exact_zero(g, v, q)) s = 0;
    return s;
}

static int rt_point_sign_arf(const Expr* g, const char* v, const arf_t p) {
    fmpq_t q; fmpq_init(q);
    arf_get_fmpq(q, p);
    int s = rt_point_sign(g, v, q);
    fmpq_clear(q);
    return s;
}

/* ------------------------------------------------------------------ *
 *  Isolation context                                                  *
 * ------------------------------------------------------------------ */

typedef struct {
    const Expr* xsym;
    const char* v;
    Expr*  d[RT_MAXD + 1];       /* d[0] = f, d[j] = j-th derivative (lazy) */
    int    nd;                   /* highest computed index                  */
    fmpq*  ex;  int nex, capex;  /* exact (rational) roots                  */
    arf_struct* ilo; arf_struct* ihi; int niv, capiv;   /* simple-root brackets */
    /* Exclusion zones around certified exact multiple roots: f^(k) != 0 on
     * the zone, so it holds no other root and the search skips it.  Without
     * them the intervals NEXT to a multiple root never certify (ball
     * dependency: f ~ eps^k but its enclosure is ~eps wide). */
    arf_struct* zlo; arf_struct* zhi; int nz, capz;
    int64_t budget;
    int    fail;
} RTCtx;

static void rt_add_zone(RTCtx* c, const arf_t lo, const arf_t hi) {
    if (c->nz == c->capz) {
        int nc = c->capz ? 2 * c->capz : 4;
        c->zlo = realloc(c->zlo, (size_t)nc * sizeof(arf_struct));
        c->zhi = realloc(c->zhi, (size_t)nc * sizeof(arf_struct));
        for (int i = c->capz; i < nc; i++) { arf_init(&c->zlo[i]); arf_init(&c->zhi[i]); }
        c->capz = nc;
    }
    arf_set(&c->zlo[c->nz], lo);
    arf_set(&c->zhi[c->nz], hi);
    c->nz++;
}

static bool rt_in_zone(const RTCtx* c, const arf_t lo, const arf_t hi) {
    for (int i = 0; i < c->nz; i++)
        if (arf_cmp(&c->zlo[i], lo) <= 0 && arf_cmp(hi, &c->zhi[i]) <= 0) return true;
    return false;
}

static const Expr* rt_deriv(RTCtx* c, int j) {
    while (c->nd < j) {
        c->d[c->nd + 1] = eval_and_free(expr_new_function(expr_new_symbol(SYM_D),
            (Expr*[]){ expr_copy(c->d[c->nd]), expr_copy((Expr*)c->xsym) }, 2));
        c->nd++;
    }
    return c->d[j];
}

/* Enclosure of f^(j) over the ball X: the naive interval evaluation
 * intersected with the mean-value form f^(j)(m) + f^(j+1)(X) (X - m).  The
 * naive form alone loses the cancellation in e.g. Cos[x] - 1 near 0 (its
 * width is ~rad(X) while the value is ~x^2), which made the search near a
 * multiple root visit ~1/eps intervals; the centred form is second-order.
 * The derivative factor is itself enclosed the same way, `depth` levels
 * deep (a Taylor-like form), so a root of multiplicity up to ~depth keeps a
 * logarithmic search. */
#define RT_ENCL_DEPTH 6
static bool rt_encl_d(RTCtx* c, int j, const arb_t X, arb_t out, int depth) {
    bool ok = flint_cert_eval_real(rt_deriv(c, j), c->v, X, out, RT_PREC) != 0;
    if (j >= RT_MAXD || depth <= 0 || mag_is_zero(arb_radref(X))) return ok;
    arb_t m, fm, g, t;
    arb_init(m); arb_init(fm); arb_init(g); arb_init(t);
    arb_get_mid_arb(m, X);
    if (flint_cert_eval_real(rt_deriv(c, j), c->v, m, fm, RT_PREC)
        && rt_encl_d(c, j + 1, X, g, depth - 1)) {
        arb_sub(t, X, m, RT_PREC);
        arb_mul(t, t, g, RT_PREC);
        arb_add(t, t, fm, RT_PREC);
        if (!ok) { arb_set(out, t); ok = true; }
        else if (!arb_intersection(out, out, t, RT_PREC)) arb_set(out, t);
    }
    arb_clear(m); arb_clear(fm); arb_clear(g); arb_clear(t);
    return ok;
}

static bool rt_encl(RTCtx* c, int j, const arb_t X, arb_t out) {
    return rt_encl_d(c, j, X, out, RT_ENCL_DEPTH);
}

static void rt_add_exact(RTCtx* c, const fmpq_t q) {
    for (int i = 0; i < c->nex; i++) if (fmpq_equal(&c->ex[i], q)) return;
    if (c->nex == c->capex) {
        int nc = c->capex ? 2 * c->capex : 4;
        c->ex = realloc(c->ex, (size_t)nc * sizeof(fmpq));
        for (int i = c->capex; i < nc; i++) fmpq_init(&c->ex[i]);
        c->capex = nc;
    }
    fmpq_set(&c->ex[c->nex++], q);
}

static void rt_add_bracket(RTCtx* c, const arf_t lo, const arf_t hi) {
    if (c->niv == c->capiv) {
        int nc = c->capiv ? 2 * c->capiv : 4;
        c->ilo = realloc(c->ilo, (size_t)nc * sizeof(arf_struct));
        c->ihi = realloc(c->ihi, (size_t)nc * sizeof(arf_struct));
        for (int i = c->capiv; i < nc; i++) { arf_init(&c->ilo[i]); arf_init(&c->ihi[i]); }
        c->capiv = nc;
    }
    arf_set(&c->ilo[c->niv], lo);
    arf_set(&c->ihi[c->niv], hi);
    c->niv++;
}

static void rt_ctx_clear(RTCtx* c) {
    for (int j = 0; j <= c->nd; j++) expr_free(c->d[j]);
    for (int i = 0; i < c->capex; i++) fmpq_clear(&c->ex[i]);
    free(c->ex);
    for (int i = 0; i < c->capiv; i++) { arf_clear(&c->ilo[i]); arf_clear(&c->ihi[i]); }
    free(c->ilo); free(c->ihi);
    for (int i = 0; i < c->capz; i++) { arf_clear(&c->zlo[i]); arf_clear(&c->zhi[i]); }
    free(c->zlo); free(c->zhi);
}

/* hi - lo below the isolation resolution 2^-RT_RES_BITS * max(1, |lo|, |hi|). */
static bool rt_is_small(const arf_t lo, const arf_t hi) {
    arf_t w, s, t;
    arf_init(w); arf_init(s); arf_init(t);
    arf_sub(w, hi, lo, RT_PREC, ARF_RND_UP);
    arf_abs(s, lo); arf_abs(t, hi);
    if (arf_cmp(t, s) > 0) arf_set(s, t);
    if (arf_cmp_si(s, 1) < 0) arf_one(s);
    arf_mul_2exp_si(s, s, -RT_RES_BITS);
    bool small = arf_cmp(w, s) <= 0;
    arf_clear(w); arf_clear(s); arf_clear(t);
    return small;
}

/* f is strictly monotone on [lo, hi]: at most one root, read off the signs. */
static void rt_monotone(RTCtx* c, const arf_t lo, const arf_t hi) {
    int sl = rt_point_sign_arf(c->d[0], c->v, lo);
    int sh = rt_point_sign_arf(c->d[0], c->v, hi);
    if (sl == 2 || sh == 2) { c->fail = 1; return; }
    fmpq_t q; fmpq_init(q);
    if (sl == 0) { arf_get_fmpq(q, lo); rt_add_exact(c, q); }
    if (sh == 0) { arf_get_fmpq(q, hi); rt_add_exact(c, q); }
    fmpq_clear(q);
    if (sl != 0 && sh != 0 && sl != sh) rt_add_bracket(c, lo, hi);
}

/* A resolution-sized interval on which f' may vanish: the only certifiable
 * case is an exact rational root q of multiplicity k with f^(k) != 0 on the
 * whole interval (then, by Rolle, f has no other root there). */
static void rt_multiple(RTCtx* c, const arf_t lo, const arf_t hi) {
    fmpq_t ql, qh, q; fmpq_init(ql); fmpq_init(qh); fmpq_init(q);
    arf_get_fmpq(ql, lo); arf_get_fmpq(qh, hi);
    fmpq_simplest_between(q, ql, qh);
    if (rt_point_sign(c->d[0], c->v, q) != 0) { c->fail = 1; goto out; }
    for (int j = 1; j <= RT_MAXD; j++) {
        const Expr* dj = rt_deriv(c, j);
        int s = rt_point_sign(dj, c->v, q);
        if (s == 0) continue;
        if (s == 2) break;
        arb_t X, y; arb_init(X); arb_init(y);
        arb_set_interval_arf(X, lo, hi, RT_PREC);
        bool excl = rt_encl(c, j, X, y) && !arb_contains_zero(y);
        if (excl) {
            /* Grow the exclusion zone [q - r, q + r] while f^(j) stays
             * certified nonzero on it: by Rolle, q is f's only root there. */
            arf_t r, zl, zh; arf_init(r); arf_init(zl); arf_init(zh);
            arf_sub(r, hi, lo, RT_PREC, ARF_RND_UP);
            arb_t Z; arb_init(Z);
            for (int g = 0; g < 64; g++) {
                arf_t r2; arf_init(r2);
                arf_mul_2exp_si(r2, r, 1);
                arb_set_fmpq(Z, q, RT_PREC);
                arb_add_error_arf(Z, r2);
                bool more = rt_encl(c, j, Z, y) && !arb_contains_zero(y);
                if (more) arf_set(r, r2);
                arf_clear(r2);
                if (!more) break;
            }
            arb_set_fmpq(Z, q, RT_PREC);
            arb_add_error_arf(Z, r);
            if (rt_encl(c, j, Z, y) && !arb_contains_zero(y)) {
                /* The certified ball contains [q - r, q + r] exactly; the zone
                 * is that interval rounded INWARD (plus [lo, hi], certified
                 * above, which it already contains). */
                arf_set_fmpq(zl, q, RT_PREC, ARF_RND_UP);
                arf_sub(zl, zl, r, RT_PREC, ARF_RND_UP);
                arf_set_fmpq(zh, q, RT_PREC, ARF_RND_DOWN);
                arf_add(zh, zh, r, RT_PREC, ARF_RND_DOWN);
                if (arf_cmp(zl, lo) > 0) arf_set(zl, lo);
                if (arf_cmp(zh, hi) < 0) arf_set(zh, hi);
                rt_add_zone(c, zl, zh);
            }
            arb_clear(Z);
            arf_clear(r); arf_clear(zl); arf_clear(zh);
            arb_clear(X); arb_clear(y);
            rt_add_exact(c, q);
            goto out;
        }
        arb_clear(X); arb_clear(y);
        break;
    }
    c->fail = 1;
out:
    fmpq_clear(ql); fmpq_clear(qh); fmpq_clear(q);
}

static void rt_iso(RTCtx* c, const arf_t lo, const arf_t hi) {
    /* A failure does NOT stop the scan: the rest of the span may still reveal
     * an exclusion zone, after which the driver reruns (see below). */
    if (c->budget < 0) return;
    if (--c->budget < 0) { c->fail = 1; return; }
    if (rt_in_zone(c, lo, hi)) return;
    arb_t X, y; arb_init(X); arb_init(y);
    arb_set_interval_arf(X, lo, hi, RT_PREC);
    bool okF = rt_encl(c, 0, X, y);
    if (okF && !arb_contains_zero(y)) goto done;               /* root-free */
    if (okF) {
        if (rt_encl(c, 1, X, y) && !arb_contains_zero(y)) {
            rt_monotone(c, lo, hi);
            goto done;
        }
    }
    if (rt_is_small(lo, hi)) {
        if (okF) rt_multiple(c, lo, hi);
        else c->fail = 1;                       /* not evaluable: decline */
        goto done;
    }
    {
        arf_t mid; arf_init(mid);
        arf_add(mid, lo, hi, ARF_PREC_EXACT, ARF_RND_DOWN);
        arf_mul_2exp_si(mid, mid, -1);
        rt_iso(c, lo, mid);
        rt_iso(c, mid, hi);
        arf_clear(mid);
    }
done:
    arb_clear(X); arb_clear(y);
}

/* ------------------------------------------------------------------ *
 *  Certification of Root[{f, x0}] (interval Newton)                   *
 * ------------------------------------------------------------------ */

/* Certify a unique simple root of f in the window x0 * (1 +- 2^-RT_WIN_BITS)
 * (absolute for |x0| < 1) and refine it to relative radius 2^-target_bits.
 * On success `out` encloses the root.  Interval Newton: with f' != 0 on X,
 * N(X) = m - f(m)/f'(X) contains every root in X, and N(X) inside X proves
 * that X holds exactly one. */
static bool rt_certify(const Expr* f, const Expr* df, const char* v,
                       const arf_t x0, slong target_bits, arb_t out) {
    slong prec = target_bits + 64;
    arb_t X, m, fm, D, N;
    arb_init(X); arb_init(m); arb_init(fm); arb_init(D); arb_init(N);
    arf_t w; arf_init(w);
    arf_abs(w, x0);
    if (arf_cmp_si(w, 1) < 0) arf_one(w);
    arf_mul_2exp_si(w, w, -RT_WIN_BITS);
    arf_set(arb_midref(X), x0);
    arf_get_mag(arb_radref(X), w);
    bool certified = false, ok = false;
    for (int it = 0; it < 200; it++) {
        arb_get_mid_arb(m, X);
        if (!flint_cert_eval_real(f, v, m, fm, prec)) break;
        if (!flint_cert_eval_real(df, v, X, D, prec) || arb_contains_zero(D)) break;
        arb_div(N, fm, D, prec);
        arb_sub(N, m, N, prec);
        if (arb_contains(X, N)) certified = true;
        mag_t r0; mag_init(r0); mag_set(r0, arb_radref(X));
        if (!arb_intersection(X, X, N, prec)) { mag_clear(r0); break; }   /* no root */
        if (certified && arb_rel_accuracy_bits(X) >= target_bits) {
            mag_clear(r0); ok = true; break;
        }
        if (mag_cmp(arb_radref(X), r0) >= 0) {         /* stalled: more precision */
            prec += 64;
            if (prec > target_bits + 64 + 4096) { mag_clear(r0); break; }
        }
        mag_clear(r0);
    }
    if (ok) arb_set(out, X);
    arb_clear(X); arb_clear(m); arb_clear(fm); arb_clear(D); arb_clear(N);
    arf_clear(w);
    return ok;
}

/* f(t) and f'(t) of a transcendental Root's body, in a private variable. */
static const char* rt_var_name(void) {
    static const char* name = NULL;
    if (!name) {
        Expr* s = expr_new_symbol("RootTrans`t$");
        name = s->data.symbol.name;
        expr_free(s);
    }
    return name;
}

static bool rt_object_fns(const Expr* body, Expr** f, Expr** df) {
    Expr* t = expr_new_symbol(rt_var_name());
    *f = eval_and_free(rt_subst_slot(body, t));
    *df = eval_and_free(expr_new_function(expr_new_symbol(SYM_D),
        (Expr*[]){ expr_copy(*f), t }, 2));
    return true;
}

static bool rt_object_certify(const Expr* e, slong bits, arb_t out) {
    const Expr *body, *x0;
    if (!rt_parse(e, &body, &x0) || !rt_body_ok(body, 0)) return false;
    arf_t a; arf_init(a);
    bool ok = rt_number_arf(x0, a);
    if (ok) {
        Expr *f, *df;
        rt_object_fns(body, &f, &df);
        ok = rt_certify(f, df, rt_var_name(), a, bits, out);
        expr_free(f); expr_free(df);
    }
    arf_clear(a);
    return ok;
}

Expr* root_trans_numericalize(const Expr* e, NumericSpec spec) {
    bool machine = (spec.mode == NUMERIC_MODE_MACHINE);
    slong bits = machine ? 53 : (slong)spec.bits;
    if (bits < 53) bits = 53;
    arb_t r; arb_init(r);
    Expr* out = NULL;
    if (rt_object_certify(e, bits + 8, r)) {
        if (machine) {
            out = expr_new_real(arf_get_d(arb_midref(r), ARF_RND_NEAR));
        } else {
            out = expr_new_mpfr_bits((mpfr_prec_t)bits);
            arf_get_mpfr(out->data.mpfr, arb_midref(r), MPFR_RNDN);
        }
    }
    arb_clear(r);
    return out;
}

Expr* root_trans_canonicalize(const Expr* e, bool* invalid) {
    *invalid = false;
    const Expr *body, *x0;
    if (!rt_parse(e, &body, &x0) || !rt_body_ok(body, 0)) return NULL;
    arb_t r; arb_init(r);
    bool ok = rt_object_certify(e, 96, r);
    double d = ok ? arf_get_d(arb_midref(r), ARF_RND_NEAR) : 0.0;
    arb_clear(r);
    if (!ok) { *invalid = true; return NULL; }
    if (x0->type == EXPR_MPFR) return NULL;          /* keep the given precision */
    if (x0->type == EXPR_REAL) {
        /* Within one ulp counts as canonical: no flip-flop between neighbours. */
        double c = x0->data.real;
        if (c == d || nextafter(c, d) == d) return NULL;
    }
    Expr* lst = expr_new_function(expr_new_symbol(SYM_List),
        (Expr*[]){ expr_new_function(expr_new_symbol(SYM_Function),
                       (Expr*[]){ expr_copy((Expr*)body) }, 1),
                   expr_new_real(d) }, 2);
    return expr_new_function(expr_new_symbol(SYM_Root), (Expr*[]){ lst }, 1);
}

/* ------------------------------------------------------------------ *
 *  Dominance bound                                                    *
 * ------------------------------------------------------------------ */

static bool rt_free_of(const Expr* e, const char* v) {
    if (e->type == EXPR_SYMBOL) return e->data.symbol.name != v;
    if (e->type != EXPR_FUNCTION) return true;
    if (!rt_free_of(e->data.function.head, v)) return false;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (!rt_free_of(e->data.function.args[i], v)) return false;
    return true;
}

static bool rt_is_true(Expr* e) {
    bool t = e->type == EXPR_SYMBOL && e->data.symbol.name == SYM_True;
    expr_free(e);
    return t;
}

/* arg is a polynomial in x with REAL constant coefficients: PolynomialQ, and
 * its enclosure over the real ball [-1, 1] has an exactly-zero imaginary part
 * (a complex coefficient leaves a nonzero imaginary ball). */
static bool rt_real_poly(const Expr* arg, const Expr* x) {
    if (!rt_is_true(eval_and_free(expr_new_function(expr_new_symbol(SYM_PolynomialQ),
            (Expr*[]){ expr_copy((Expr*)arg), expr_copy((Expr*)x) }, 2)))) return false;
    arb_t X, y; arb_init(X); arb_init(y);
    arb_zero(X); mag_one(arb_radref(X));
    int ok = flint_cert_eval_real(arg, x->data.symbol.name, X, y, 64);
    arb_clear(X); arb_clear(y);
    return ok != 0;
}

/* The exponent of an exponential factor (Exp[a] / E^a), borrowed, or NULL. */
static const Expr* rt_exp_arg(const Expr* e) {
    if (rt_is_head(e, SYM_Exp) && e->data.function.arg_count == 1)
        return e->data.function.args[0];
    if (rt_is_head(e, SYM_Power) && e->data.function.arg_count == 2
        && e->data.function.args[0]->type == EXPR_SYMBOL
        && e->data.function.args[0]->data.symbol.name == SYM_E)
        return e->data.function.args[1];
    return NULL;
}

/* Upper bound of |fac| for real x on the half-line side*x >= 0, written to
 * *bnd; false when fac is not bounded there.  Bounded: Sin / Cos / Tanh / Erf
 * / Sech / ArcTan of a real polynomial (everywhere), and Exp[a x + b] with
 * a*side <= 0 (bound Exp[b]).  For a*side < 0, *decay receives a lower bound
 * of |a| (the factor then also absorbs any power of x); else 0. */
static bool rt_bounded_factor(const Expr* fac, const Expr* x, int side, double* bnd,
                              double* decay) {
    *decay = 0.0;
    const Expr* ea = rt_exp_arg(fac);
    if (ea) {
        if (!rt_real_poly(ea, x)) return false;
        Expr* cl = eval_and_free(expr_new_function(expr_new_symbol(SYM_CoefficientList),
            (Expr*[]){ expr_copy((Expr*)ea), expr_copy((Expr*)x) }, 2));
        bool ok = false;
        if (rt_is_head(cl, SYM_List) && cl->data.function.arg_count >= 1
            && cl->data.function.arg_count <= 2) {
            arb_t b, a; arb_init(b); arb_init(a);
            ok = flint_cert_eval_real(cl->data.function.args[0], NULL, NULL, b, RT_PREC) != 0;
            if (ok && cl->data.function.arg_count == 2) {
                ok = flint_cert_eval_real(cl->data.function.args[1], NULL, NULL, a, RT_PREC) != 0;
                if (ok) {
                    if (side < 0) arb_neg(a, a);
                    ok = arb_is_negative(a) || arb_is_zero(a);
                    if (ok && arb_is_negative(a)) {
                        arf_t l; arf_init(l);
                        arb_neg(a, a);
                        arb_get_lbound_arf(l, a, RT_PREC);
                        *decay = arf_get_d(l, ARF_RND_DOWN);
                        arf_clear(l);
                    }
                }
            }
            if (ok) {
                arb_exp(b, b, RT_PREC);
                arf_t u; arf_init(u);
                arb_get_ubound_arf(u, b, RT_PREC);
                *bnd = arf_get_d(u, ARF_RND_UP);
                arf_clear(u);
                ok = isfinite(*bnd) && *bnd < 1e15;
            }
            arb_clear(b); arb_clear(a);
        }
        expr_free(cl);
        return ok;
    }
    if (fac->type != EXPR_FUNCTION || fac->data.function.arg_count != 1
        || fac->data.function.head->type != EXPR_SYMBOL) return false;
    const char* h = fac->data.function.head->data.symbol.name;
    if (h == SYM_Sin || h == SYM_Cos || h == SYM_Tanh || h == SYM_Erf || h == SYM_Sech)
        *bnd = 1.0;
    else if (h == SYM_ArcTan) *bnd = 2.0;              /* Pi/2 < 2 */
    else return false;
    return rt_real_poly(fac->data.function.args[0], x);
}

/* term = c * x^k * B on the half-line side*x >= 0, with c an exact real
 * constant and |B| <= *bnd (B a product of bounded factors; *bounded says
 * whether there was any; *bnd is +Inf for an exponential that grows on this
 * side).  *pos_exp reports that B is a product of exponentials of real
 * polynomials only (so B > 0 everywhere); *sgn is the sign of x^k on the side.
 * An exponential that DECAYS on this side (rate >= d > 0) absorbs the power:
 * |x^k Exp[-d |x|]| <= (k/(d e))^k, so the term is reported with k = 0. */
static bool rt_split_term(const Expr* t, const Expr* x, int side, int64_t* k,
                          arb_t c, double* bnd, bool* bounded, bool* pos_exp,
                          int* sgn) {
    double decay = 0.0;
    const char* v = x->data.symbol.name;
    bool times = rt_is_head(t, SYM_Times);
    size_t n = times ? t->data.function.arg_count : 1;
    *k = 0; *bnd = 1.0; *bounded = false; *pos_exp = true;
    arb_one(c);
    arb_t y; arb_init(y);
    bool ok = true;
    for (size_t i = 0; i < n && ok; i++) {
        const Expr* fac = times ? t->data.function.args[i] : t;
        if (rt_free_of(fac, v)) {
            ok = flint_cert_eval_real(fac, NULL, NULL, y, RT_PREC) != 0;
            if (ok) arb_mul(c, c, y, RT_PREC);
            continue;
        }
        const Expr* base = fac;
        int64_t ex = 1;
        if (rt_is_head(fac, SYM_Power) && fac->data.function.arg_count == 2
            && fac->data.function.args[1]->type == EXPR_INTEGER
            && fac->data.function.args[1]->data.integer > 0) {
            base = fac->data.function.args[0];
            ex = fac->data.function.args[1]->data.integer;
        }
        if (base->type == EXPR_SYMBOL && base->data.symbol.name == v) { *k += ex; continue; }
        double b;
        const Expr* ea = rt_exp_arg(base);
        if (ea && rt_real_poly(ea, x)) {
            *bounded = true;
            double d;
            if (rt_bounded_factor(base, x, side, &b, &d)) {
                *bnd *= pow(b, (double)ex);
                decay += d * (double)ex;
            } else *bnd = INFINITY;                  /* positive but unbounded */
            continue;
        }
        double dd;
        if (!rt_bounded_factor(base, x, side, &b, &dd)) { ok = false; break; }
        *bnd *= pow(b, (double)ex);
        *bounded = true;
        *pos_exp = false;
    }
    arb_clear(y);
    *sgn = (side < 0 && (*k % 2)) ? -1 : 1;
    if (ok && *k > 0 && decay > 0.0 && !isinf(*bnd)) {
        /* max over x >= 0 of x^k e^(-d x) is (k/(d e))^k; pad the double
         * evaluation (e taken as 2.718, slightly low) to stay an upper bound. */
        double kk = (double)*k;
        *bnd *= pow(kk / (decay * 2.718), kk) * (1.0 + 1e-9);
        *k = 0;
    }
    return ok && (isinf(*bnd) || *bnd < 1e15);
}

/* f has the constant nonzero sign *sign for all x with side*x >= *R >= 0.
 * Three certificates, tried in order:
 *   1. sign-definite: every term is c x^(2j) * (exponentials of real
 *      polynomials) with the SAME sign of c, one of them with j = 0;
 *   2. dominance: the top-degree pure power term A x^n beats all lower terms,
 *      whose coefficients are bounded on that side (Cauchy-type bound);
 *   3. monotone tail: f' has constant sign s on side*x >= R1 (recursively),
 *      and f(side R) already has sign s*side at some R >= R1 -- then |f| only
 *      grows from there on (Exp[x] - 3 x on the right).
 * False when none applies (Sin[x], x Sin[x] - 1). */
static bool rt_side_bound(const Expr* f, const Expr* x, int side, int depth,
                          int64_t* R, int* sign) {
    Expr* g = rt_call1(SYM_Expand, expr_copy((Expr*)f));
    bool plus = rt_is_head(g, SYM_Plus);
    size_t nt = plus ? g->data.function.arg_count : 1;
    int64_t* deg = malloc(nt * sizeof(int64_t));
    double* bnd = malloc(nt * sizeof(double));
    bool* bdd = malloc(nt * sizeof(bool));
    bool* pex = malloc(nt * sizeof(bool));
    int* sgf = malloc(nt * sizeof(int));
    arb_struct* coef = malloc(nt * sizeof(arb_struct));
    for (size_t i = 0; i < nt; i++) arb_init(&coef[i]);
    bool ok = true, found = false;
    for (size_t i = 0; i < nt && ok; i++) {
        const Expr* t = plus ? g->data.function.args[i] : g;
        ok = rt_split_term(t, x, side, &deg[i], &coef[i], &bnd[i], &bdd[i], &pex[i], &sgf[i]);
    }
    if (ok) {
        /* 1. sign-definite. */
        int s0 = 0; bool def = true, has_const = false;
        for (size_t i = 0; i < nt && def; i++) {
            /* on the side, x^k has the fixed sign sgf (zero only at 0) */
            int s = arb_is_positive(&coef[i]) ? 1 : (arb_is_negative(&coef[i]) ? -1 : 0);
            s *= sgf[i];
            if (!s || (bdd[i] && !pex[i])) def = false;
            else if (s0 && s != s0) def = false;
            s0 = s;
            if (deg[i] == 0) has_const = true;
        }
        if (def && has_const) { *R = 0; *sign = s0; found = true; }
    }
    if (ok && !found) {
        /* 2. dominance by the top pure power. */
        int64_t n = -1;
        for (size_t i = 0; i < nt; i++) if (!bdd[i] && deg[i] > n) n = deg[i];
        bool dom = n >= 0;
        for (size_t i = 0; i < nt; i++)
            if (bdd[i] && (deg[i] >= n || isinf(bnd[i]))) dom = false;
        arb_t A, M, y, bb; arb_init(A); arb_init(M); arb_init(y); arb_init(bb);
        if (dom) {
            arb_zero(A); arb_zero(M);
            for (size_t i = 0; i < nt; i++) {
                if (!bdd[i] && deg[i] == n) { arb_add(A, A, &coef[i], RT_PREC); continue; }
                arb_abs(y, &coef[i]);
                arb_set_d(bb, bnd[i]);
                arb_mul(y, y, bb, RT_PREC);
                arb_add(M, M, y, RT_PREC);
            }
            int sA = arb_is_positive(A) ? 1 : (arb_is_negative(A) ? -1 : 0);
            arb_abs(A, A);
            if (sA && n == 0) {
                arb_sub(y, A, M, RT_PREC);
                if (arb_is_positive(y)) { *R = 0; *sign = sA; found = true; }
            } else if (sA) {
                /* |x| >= 1 and |x| > M/|A|  =>  |A x^n| > M |x|^(n-1) >= rest. */
                arb_div(y, M, A, RT_PREC);
                arf_t u; arf_init(u);
                arb_get_ubound_arf(u, y, RT_PREC);
                if (arf_cmp_si(u, 1) < 0) arf_one(u);
                if (arf_cmpabs_2exp_si(u, 52) < 0) {
                    *R = (int64_t)floor(arf_get_d(u, ARF_RND_UP)) + 1;
                    *sign = (side < 0 && (n % 2)) ? -sA : sA;
                    found = true;
                }
                arf_clear(u);
            }
        }
        arb_clear(A); arb_clear(M); arb_clear(y); arb_clear(bb);
    }
    for (size_t i = 0; i < nt; i++) arb_clear(&coef[i]);
    free(coef); free(bnd); free(bdd); free(pex); free(sgf); free(deg);
    expr_free(g);
    if (!ok || found) return found;

    /* 3. monotone tail through the derivative. */
    if (depth >= 8) return false;
    Expr* df = eval_and_free(expr_new_function(expr_new_symbol(SYM_D),
        (Expr*[]){ expr_copy((Expr*)f), expr_copy((Expr*)x) }, 2));
    int64_t R1; int s1;
    bool sub = rt_side_bound(df, x, side, depth + 1, &R1, &s1);
    expr_free(df);
    if (!sub) return false;
    int want = s1 * side;
    if (R1 < 1) R1 = 1;
    for (int it = 0; it < 40; it++, R1 *= 2) {
        fmpq_t q; fmpq_init(q);
        fmpq_set_si(q, (slong)(side * R1), 1);
        int s = rt_point_sign(f, x->data.symbol.name, q);
        fmpq_clear(q);
        if (s == want) { *R = R1; *sign = want; return true; }
        if (s == 2) return false;
    }
    return false;
}

Expr* root_trans_dominance_bound(const Expr* f, const Expr* x) {
    if (!x || x->type != EXPR_SYMBOL) return NULL;
    int64_t Rp, Rm; int sp, sm;
    if (!rt_side_bound(f, x, 1, 0, &Rp, &sp)) return NULL;
    if (!rt_side_bound(f, x, -1, 0, &Rm, &sm)) return NULL;
    return expr_new_integer(Rp > Rm ? Rp : Rm);
}

/* ------------------------------------------------------------------ *
 *  Isolation driver                                                   *
 * ------------------------------------------------------------------ */

bool root_trans_isolate(const Expr* f, const Expr* x, const Expr* lo,
                        const Expr* hi, Expr*** out, int* n, int* cap) {
    if (!x || x->type != EXPR_SYMBOL) return false;
    arb_t bl, bh; arb_init(bl); arb_init(bh);
    if (!flint_cert_eval_real(lo, NULL, NULL, bl, RT_PREC)
        || !flint_cert_eval_real(hi, NULL, NULL, bh, RT_PREC)) {
        arb_clear(bl); arb_clear(bh);
        return false;
    }
    arf_t a, b; arf_init(a); arf_init(b);
    arb_get_lbound_arf(a, bl, RT_PREC);
    arb_get_ubound_arf(b, bh, RT_PREC);
    arb_clear(bl); arb_clear(bh);

    RTCtx c;
    memset(&c, 0, sizeof c);
    c.xsym = x;
    c.v = x->data.symbol.name;
    c.d[0] = expr_copy((Expr*)f);
    c.nd = 0;
    c.budget = INT64_C(1) << 17;
    /* Rerun while a failed pass discovered a new exclusion zone: the zone
     * removes exactly the intervals around a multiple root that failed. */
    for (;;) {
        int nz0 = c.nz;
        c.fail = 0;
        c.niv = 0;
        if (arf_cmp(a, b) <= 0) rt_iso(&c, a, b);
        if (!c.fail || c.nz == nz0 || c.budget < 0) break;
    }
    arf_clear(a); arf_clear(b);

    int n0 = *n;
    bool ok = !c.fail;
    /* Exact roots. */
    for (int i = 0; ok && i < c.nex; i++) {
        if (*n == *cap) { *cap = *cap ? 2 * *cap : 8; *out = realloc(*out, (size_t)*cap * sizeof(Expr*)); }
        (*out)[(*n)++] = rt_fmpq_expr(&c.ex[i]);
    }
    /* Simple roots: refine the bracket by bisection, then emit an exact
     * rational if one is a root, else a certified Root[{f, approx}]. */
    arf_t l, h, mid; arf_init(l); arf_init(h); arf_init(mid);
    fmpq_t ql, qh, q; fmpq_init(ql); fmpq_init(qh); fmpq_init(q);
    for (int i = 0; ok && i < c.niv; i++) {
        arf_set(l, &c.ilo[i]); arf_set(h, &c.ihi[i]);
        int sl = rt_point_sign_arf(c.d[0], c.v, l);
        bool exact = false;
        for (int it = 0; it < 400; it++) {
            arf_t w, s; arf_init(w); arf_init(s);
            arf_sub(w, h, l, RT_PREC, ARF_RND_UP);
            arf_abs(s, l);
            if (arf_cmp_si(s, 1) < 0) arf_one(s);
            arf_mul_2exp_si(s, s, -64);
            bool done = arf_cmp(w, s) <= 0;
            arf_clear(w); arf_clear(s);
            if (done) break;
            arf_add(mid, l, h, ARF_PREC_EXACT, ARF_RND_DOWN);
            arf_mul_2exp_si(mid, mid, -1);
            int sm = rt_point_sign_arf(c.d[0], c.v, mid);
            if (sm == 2) { ok = false; break; }
            if (sm == 0) { arf_set(l, mid); arf_set(h, mid); exact = true; break; }
            if (sm == sl) arf_set(l, mid); else arf_set(h, mid);
        }
        if (!ok) break;
        Expr* val = NULL;
        arf_get_fmpq(ql, l); arf_get_fmpq(qh, h);
        fmpq_simplest_between(q, ql, qh);
        if (exact || rt_point_sign(c.d[0], c.v, q) == 0) {
            val = rt_fmpq_expr(q);
        } else {
            arf_add(mid, l, h, ARF_PREC_EXACT, ARF_RND_DOWN);
            arf_mul_2exp_si(mid, mid, -1);
            val = root_trans_make(c.d[0], x, arf_get_d(mid, ARF_RND_NEAR));
            /* The object must name THIS root unambiguously: certify it. */
            arb_t r; arb_init(r);
            bool cert = rt_object_certify(val, 64, r);
            arb_clear(r);
            if (!cert) { expr_free(val); ok = false; break; }
        }
        if (*n == *cap) { *cap = *cap ? 2 * *cap : 8; *out = realloc(*out, (size_t)*cap * sizeof(Expr*)); }
        (*out)[(*n)++] = val;
    }
    arf_clear(l); arf_clear(h); arf_clear(mid);
    fmpq_clear(ql); fmpq_clear(qh); fmpq_clear(q);
    rt_ctx_clear(&c);
    if (!ok) {
        for (int i = n0; i < *n; i++) expr_free((*out)[i]);
        *n = n0;
    }
    return ok;
}

#else /* !(USE_FLINT && USE_MPFR) */

Expr* root_trans_numericalize(const Expr* e, NumericSpec spec) {
    (void)e; (void)spec; return NULL;
}
Expr* root_trans_canonicalize(const Expr* e, bool* invalid) {
    (void)e; *invalid = false; return NULL;
}
Expr* root_trans_dominance_bound(const Expr* f, const Expr* x) {
    (void)f; (void)x; return NULL;
}
bool root_trans_isolate(const Expr* f, const Expr* x, const Expr* lo,
                        const Expr* hi, Expr*** out, int* n, int* cap) {
    (void)f; (void)x; (void)lo; (void)hi; (void)out; (void)n; (void)cap;
    return false;
}

#endif
