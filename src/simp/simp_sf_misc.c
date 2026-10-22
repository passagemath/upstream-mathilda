/*
 * simp_sf_misc.c -- special-function Simplify descriptors: complete
 * elliptic integrals, ProductLog, QPochhammer.
 *
 * Derivations in docs/design/simp_special_functions.md (section 11); every
 * relation's domain probed numerically in all four quadrants and on the
 * cuts.
 *
 * Elliptic (pass): Legendre's relation (parameter convention m = k^2)
 *     E(m) K(1 - m) + E(1 - m) K(m) - K(m) K(1 - m) = Pi/2.
 *   Differentiating the left side in m with the derivative rules
 *   dK/dm = (E - (1-m) K)/(2 m (1-m)), dE/dm = (E - K)/(2m) gives 0, so it is
 *   constant; m -> 0 gives Pi/2.  It held at every probed m, complex and
 *   outside (0, 1) included.  Expressions holding EllipticE[m] next to
 *   EllipticK[1 - m] are reduced modulo it, eliminating E(m) K(1 - m).
 *
 * ProductLog (pass): W = ProductLog[z] satisfies W E^W = z by definition, so
 *     E^(k W) = (z/W)^k  (integer k),   Log[W] = Log[z] - W
 *   (the second held at every probed z, including the cut (-Inf, -1/e)).
 *   Canon: ProductLog[u E^u] = u only for u real and u >= -1 proved (the
 *   principal branch; u = -3/2 gives 0.87, checked).
 *
 * QPochhammer (rank 1 in n): (a; q)_{n+1} = (a; q)_n (1 - a q^n), and the
 *   bulk step (a; q)_{n+k} = (a; q)_n (a q^n; q)_k (the product definition).
 */

#include <stdlib.h>
#include <string.h>
#include "expr.h"
#include "eval.h"
#include "sym_names.h"
#include "simp.h"
#include "simp_internal.h"
#include "simp_sf_engine.h"

/* ---- elliptic ---------------------------------------------------------------- */

typedef struct { const Expr* a[32]; size_t n; } AList;

static void collect_args1(const Expr* e, const char* head, AList* out) {
    if (e->type != EXPR_FUNCTION) return;
    if (sf_is_call(e, head, 1) && out->n < 32) {
        bool dup = false;
        for (size_t i = 0; i < out->n && !dup; i++) dup = expr_eq(out->a[i], e->data.function.args[0]);
        if (!dup) out->a[out->n++] = e->data.function.args[0];
    }
    for (size_t i = 0; i < e->data.function.arg_count; i++) collect_args1(e->data.function.args[i], head, out);
}

static bool same_value(const Expr* a, const Expr* b) {
    Expr* d = eval_and_free(sf_call("Together", (Expr*[]){ sf_minus(sf_cp(a), sf_cp(b)) }, 1));
    bool z = d->type == EXPR_INTEGER && d->data.integer == 0;
    expr_free(d);
    return z;
}

static Expr* pass_legendre_relation(const Expr* e, const Expr* input, const AssumeCtx* ctx) {
    (void)input; (void)ctx;
    AList es = { { NULL }, 0 }, ks = { { NULL }, 0 };
    collect_args1(e, SYM_EllipticE, &es);
    collect_args1(e, SYM_EllipticK, &ks);
    Expr* cur = NULL;
    for (size_t i = 0; i < es.n; i++) {
        const Expr* m = es.a[i];
        /* the complementary parameter, spelled as the expression spells it */
        const Expr* w = NULL;
        Expr* one_m = sf_minus(sf_int(1), sf_cp(m));
        for (size_t j = 0; j < ks.n && !w; j++) if (same_value(ks.a[j], one_m)) w = ks.a[j];
        expr_free(one_m);
        if (!w) continue;
        Expr* Em = sf_call(SYM_EllipticE, (Expr*[]){ sf_cp(m) }, 1);
        Expr* Ew = sf_call(SYM_EllipticE, (Expr*[]){ sf_cp(w) }, 1);
        Expr* Km = sf_call(SYM_EllipticK, (Expr*[]){ sf_cp(m) }, 1);
        Expr* Kw = sf_call(SYM_EllipticK, (Expr*[]){ sf_cp(w) }, 1);
        /* E(m) K(w) + E(w) K(m) - K(m) K(w) - Pi/2; lex E(m) > K(w) > E(w) > K(m) */
        Expr* R = sf_minus(sf_minus(sf_plus(sf_times(sf_cp(Em), sf_cp(Kw)), sf_times(sf_cp(Ew), sf_cp(Km))),
                                    sf_times(sf_cp(Km), sf_cp(Kw))),
                           sf_div(sf_sym(SYM_Pi), sf_int(2)));
        Expr* vars[4] = { Em, Kw, Ew, Km };
        Expr* red = sf_reduce_modulo(cur ? cur : e, R, vars, 4);
        expr_free(R);
        if (red) { if (cur) expr_free(cur); cur = red; }
        expr_free(Em); expr_free(Ew); expr_free(Km); expr_free(Kw);
    }
    return cur;
}

/* ---- ProductLog ---------------------------------------------------------------- */

static bool is_e(const Expr* e) { return e->type == EXPR_SYMBOL && e->data.symbol.name == SYM_E; }

/* E^(k W) -> (z/W)^k for integer k, Log[W] -> Log[z] - W. */
static Expr* plog_walk(const Expr* e, bool* changed) {
    if (e->type != EXPR_FUNCTION) return sf_cp(e);
    size_t n = e->data.function.arg_count;
    Expr** a = malloc(sizeof(Expr*) * (n ? n : 1));
    for (size_t i = 0; i < n; i++) a[i] = plog_walk(e->data.function.args[i], changed);
    Expr* r = eval_and_free(expr_new_function(sf_cp(e->data.function.head), a, n));
    free(a);
    if (sf_is_call(r, SYM_Power, 2) && is_e(r->data.function.args[0])) {
        const Expr* ex = r->data.function.args[1];
        const Expr* W = NULL; int64_t k = 1;
        if (sf_is_call(ex, SYM_ProductLog, 1)) W = ex;
        else if (sf_is_call(ex, SYM_Times, 2) && ex->data.function.args[0]->type == EXPR_INTEGER &&
                 sf_is_call(ex->data.function.args[1], SYM_ProductLog, 1)) {
            W = ex->data.function.args[1]; k = ex->data.function.args[0]->data.integer;
        }
        if (W) {
            Expr* zw = sf_div(sf_cp(W->data.function.args[0]), sf_cp(W));
            Expr* out = sf_pow(zw, sf_int(k));
            expr_free(r);
            *changed = true;
            return out;
        }
    }
    if (sf_is_call(r, SYM_Log, 1) && sf_is_call(r->data.function.args[0], SYM_ProductLog, 1)) {
        const Expr* W = r->data.function.args[0];
        Expr* out = sf_minus(sf_mk1(SYM_Log, sf_cp(W->data.function.args[0])), sf_cp(W));
        expr_free(r);
        *changed = true;
        return out;
    }
    return r;
}

static Expr* pass_productlog(const Expr* e, const Expr* input, const AssumeCtx* ctx) {
    (void)input; (void)ctx;
    bool changed = false;
    Expr* r = plog_walk(e, &changed);
    if (!changed) { expr_free(r); return NULL; }
    return r;
}

/* ProductLog[u E^u] = u for u real, u >= -1 proved */
static Expr* canon_productlog(const Expr* atom, const Expr* root, const AssumeCtx* ctx) {
    (void)root;
    if (!ctx) return NULL;
    const Expr* arg = atom->data.function.args[0];
    if (!sf_is_call(arg, SYM_Times, arg->type == EXPR_FUNCTION ? arg->data.function.arg_count : 0)) return NULL;
    for (size_t i = 0; i < arg->data.function.arg_count; i++) {
        const Expr* f = arg->data.function.args[i];
        if (!(sf_is_call(f, SYM_Power, 2) && is_e(f->data.function.args[0]))) continue;
        const Expr* u = f->data.function.args[1];
        Expr* rest = eval_and_free(sf_div(sf_cp(arg), sf_cp(f)));
        bool match = same_value(rest, u);
        expr_free(rest);
        if (!match || !prov_re(ctx, u)) continue;
        Expr* pred = sf_call(SYM_GreaterEqual, (Expr*[]){ sf_cp(u), sf_int(-1) }, 2);
        int ok = assume_reduce_entails(ctx, pred);
        expr_free(pred);
        if (ok == 1) return sf_cp(u);
    }
    return NULL;
}

/* ---- QPochhammer ---------------------------------------------------------------- */

static bool up_qpoch(Expr* const* a, int d, Expr** M) {
    if (d != 0) return false;
    /* (a; q)_{n+1} = (a; q)_n (1 - a q^n) */
    M[0] = sf_minus(sf_int(1), sf_times(sf_cp(a[0]), sf_pow(sf_cp(a[1]), sf_cp(a[2]))));
    M[1] = sf_int(0);
    return true;
}

static Expr* bulk_qpoch(Expr* const* a, int d, int64_t k) {
    if (d != 0) return NULL;
    Expr* aqn = sf_times(sf_cp(a[0]), sf_pow(sf_cp(a[1]), sf_cp(a[2])));
    if (k > 0) return sf_call(SYM_QPochhammer, (Expr*[]){ aqn, sf_cp(a[1]), sf_int(k) }, 3);
    /* (a; q)_{n+k}/(a; q)_n = 1/(a q^(n+k); q)_{-k} */
    expr_free(aqn);
    Expr* lo = sf_times(sf_cp(a[0]), sf_pow(sf_cp(a[1]), sf_plus(sf_cp(a[2]), sf_int(k))));
    return sf_pow(sf_call(SYM_QPochhammer, (Expr*[]){ lo, sf_cp(a[1]), sf_int(-k) }, 3), sf_int(-1));
}

/* ---- module table ------------------------------------------------------------- */

static const SfDesc MISC_DESCS[] = {
    { &SYM_QPochhammer, 3, 1, 1, { 2 }, up_qpoch, bulk_qpoch, NULL,             NULL, NULL },
    { &SYM_ProductLog,  1, 0, 1, { 0 }, NULL,     NULL,       canon_productlog, NULL, NULL },
    /* registered so the kernel count sees them; Legendre's relation is a pass */
    { &SYM_EllipticK,   1, 0, 1, { 0 }, NULL,     NULL,       NULL,             NULL, NULL },
    { &SYM_EllipticE,   1, 0, 1, { 0 }, NULL,     NULL,       NULL,             NULL, NULL },
};

static const SfPass MISC_PASSES[] = { pass_legendre_relation, pass_productlog };
static const char** MISC_PASS_HEADS[] = { &SYM_EllipticE, &SYM_ProductLog };

const SfModule SF_MODULE_MISC = {
    MISC_DESCS, sizeof MISC_DESCS / sizeof MISC_DESCS[0],
    MISC_PASSES, sizeof MISC_PASSES / sizeof MISC_PASSES[0],
    MISC_PASS_HEADS, sizeof MISC_PASS_HEADS / sizeof MISC_PASS_HEADS[0],
    NULL, NULL,
};
