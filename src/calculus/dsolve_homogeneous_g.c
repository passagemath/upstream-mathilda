/*
 * dsolve_homogeneous_g.c — DSolve`HomogeneousG / DSolve`HomogeneousGImplicit.
 *
 * Maple's "homogeneous class G" (isobaric) first-order ODE: F(x, y, p) == 0
 * (p = y') invariant under the scaling x -> l x, y -> l^k y (so p -> l^(k-1) p)
 * up to an overall factor, F(l x, l^k y, l^(k-1) p) == l^w F(x, y, p).  The
 * substitution y = x^k u, t = Log[x] (so x y' = x^k (k u + u'(t))) maps it to
 *
 *        x^w F~(u, u'(t)) == 0,
 *
 * an AUTONOMOUS ODE in u(t), which the cascade already owns (NthAlgebraic splits
 * a polynomial-in-u' form into quadratures u' == g_i(u)).  Mapping back,
 * u -> y / x^k and t -> Log[x], recovers y.  Works on the UNSPLIT F, so a
 * quadratic-in-y' equation (y'^2 + x y^2 y' + y^3 == 0, k = -1) never passes
 * through the radical branch forms Solve would produce.
 *
 * Weight detection.  Differentiating the invariance at l = 1 gives the Euler
 * identity  x F_x - p F_p + k (y F_y + p F_p) == w F,  linear in the constants
 * (k, w).  It is solved from two exact rational sample points of (x, y, p) --
 * the parameters stay symbolic, so k = -1/n is found for x (a + x y^n) y' + b y
 * -- and then CHECKED identically (Together + zero test); a sample artefact can
 * only decline.  k == 1 is the ordinary homogeneous class (DSolve`Homogeneous).
 *
 * Results.  When every branch of the u(t) solve is explicit the answer is
 * explicit, y = x^k u(Log[x]), returned through dsolve_run (back-substitution
 * verified).  Otherwise every branch is converted to a first integral
 * G(x, y[x]) == C[1] and stashed for the implicit pass, which re-verifies each
 * against the ORIGINAL ODE by the implicit-function rule (dsolve_run_implicit).
 *
 * Runs late in the first-order cascade (it recurses the engine), top-level only,
 * re-entry guarded (the u(t) sub-solve is autonomous and must not re-scale),
 * with a per-top-level decline memo; MATHILDA_NO_HOMOG_G disables it for A/B.
 */
#include "dsolve_common.h"
#include "../core.h"
#include "dsolve.h"
#include "../sym_names.h"
#include "../eval.h"
#include "../sym_intern.h"
#include "../symtab.h"
#include "../attr.h"
#include "../common.h"
#include "../message.h"
#include <stdlib.h>
#include <stdint.h>

#define HG_STASH_TAG   3
#define HG_SUBSOLVE_SEC 5

static int g_hg_active = 0;        /* re-entry guard: no rescaling inside the u(t) solve */

#define HG_MEMO_SLOTS 32
static uint64_t hg_epoch = 0;
static int hg_memo_n = 0;
static uint64_t hg_memo[HG_MEMO_SLOTS];
static void hg_memo_sync(uint64_t tid){ if(tid!=hg_epoch){hg_epoch=tid;hg_memo_n=0;} }
static bool hg_memo_seen(uint64_t h){ for(int i=0;i<hg_memo_n;i++) if(hg_memo[i]==h) return true; return false; }
static void hg_memo_add(uint64_t h){ if(hg_memo_n<HG_MEMO_SLOTS && !hg_memo_seen(h)) hg_memo[hg_memo_n++]=h; }

/* True when e is ALGEBRAIC in the variables x, Y, p: every subexpression that
 * involves one of them is a Plus, Times or Power whose exponent is free of them.
 * An isobaric equation can only be built that way -- Sin[x], E^y, Log[y] have no
 * scaling weight -- and the gate also keeps the weight solve away from sampled
 * transcendental constants (Cot[3/7], ...) whose Together/Simplify can spin
 * (§2.2.53-5238).  e borrowed. */
static bool hg_algebraic(const Expr* e, const char* x, const char* Y, const char* p){
    if (!e || e->type != EXPR_FUNCTION) return true;
    if (ds_free_of(e, x) && ds_free_of(e, Y) && ds_free_of(e, p)) return true;
    if (head_is(e, SYM_Plus) || head_is(e, SYM_Times)) {
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            if (!hg_algebraic(e->data.function.args[i], x, Y, p)) return false;
        return true;
    }
    if (head_is(e, SYM_Power) && e->data.function.arg_count == 2) {
        const Expr* ex = e->data.function.args[1];
        if (!(ds_free_of(ex, x) && ds_free_of(ex, Y) && ds_free_of(ex, p))) return false;
        return hg_algebraic(e->data.function.args[0], x, Y, p);
    }
    return false;
}

/* Rational literal a/b. */
static Expr* hg_q(int64_t a, int64_t b){
    return eval_and_free(expr_new_function(expr_new_symbol(SYM_Rational),
               (Expr*[]){ expr_new_integer(a), expr_new_integer(b) }, 2));
}

/* e at (x, Y, p) = (vx, vY, vp), Together'd; e borrowed, values borrowed. */
static Expr* hg_at(const Expr* e, const char* x, const char* Y, const char* p,
                   const Expr* vx, const Expr* vY, const Expr* vp){
    Expr* r = ds_subst(expr_copy((Expr*)e), expr_new_symbol(x), expr_copy((Expr*)vx));
    r = ds_subst(r, expr_new_symbol(Y), expr_copy((Expr*)vY));
    r = ds_subst(r, expr_new_symbol(p), expr_copy((Expr*)vp));
    return eval_and_free(ds_call1(SYM_Together, r));
}

/* A usable sampled value: free of infinities / Indeterminate. */
static bool hg_finite(const Expr* e){
    return !ds_contains(e, intern_symbol("ComplexInfinity"))
        && !ds_contains(e, intern_symbol("Indeterminate"))
        && !ds_contains(e, intern_symbol("DirectedInfinity"));
}

/* Candidate (k, w) from three sample points pts[i] = (x, Y, p) as rational pairs,
 * or false.  A, B, C borrowed. */
static bool hg_sample_weights(const Expr* A, const Expr* B, const Expr* C,
                              const char* x, const char* Y, const char* p,
                              const int64_t pts[3][6], Expr** kout, Expr** wout){
    Expr *Ai[3], *Bi[3], *Ci[3];
    bool ok = true;
    for (int i = 0; i < 3; i++){
        Expr* vx = hg_q(pts[i][0], pts[i][1]);
        Expr* vY = hg_q(pts[i][2], pts[i][3]);
        Expr* vp = hg_q(pts[i][4], pts[i][5]);
        Ai[i] = hg_at(A, x, Y, p, vx, vY, vp);
        Bi[i] = hg_at(B, x, Y, p, vx, vY, vp);
        Ci[i] = hg_at(C, x, Y, p, vx, vY, vp);
        expr_free(vx); expr_free(vY); expr_free(vp);
        if (!hg_finite(Ai[i]) || !hg_finite(Bi[i]) || !hg_finite(Ci[i])) ok = false;
    }
    Expr* k = NULL; Expr* w = NULL;
    for (int j = 1; j < 3 && ok && !k; j++){
        /* det = C0 Bj - B0 Cj;  k = (A0 Cj - C0 Aj)/det;  w = (A0 Bj - Aj B0)/det */
        Expr* det = eval_and_free(ds_call1(SYM_Together, ds_call2(SYM_Subtract,
                        ds_call2(SYM_Times, expr_copy(Ci[0]), expr_copy(Bi[j])),
                        ds_call2(SYM_Times, expr_copy(Bi[0]), expr_copy(Ci[j])))));
        if (ds_is_zero(det)) { expr_free(det); continue; }
        Expr* kn = ds_call2(SYM_Subtract, ds_call2(SYM_Times, expr_copy(Ai[0]), expr_copy(Ci[j])),
                                          ds_call2(SYM_Times, expr_copy(Ci[0]), expr_copy(Ai[j])));
        Expr* wn = ds_call2(SYM_Subtract, ds_call2(SYM_Times, expr_copy(Ai[0]), expr_copy(Bi[j])),
                                          ds_call2(SYM_Times, expr_copy(Ai[j]), expr_copy(Bi[0])));
        k = ds_simplify(eval_and_free(ds_call1(SYM_Together,
                ds_call2(SYM_Times, kn, ds_call2(SYM_Power, expr_copy(det), expr_new_integer(-1))))));
        w = ds_simplify(eval_and_free(ds_call1(SYM_Together,
                ds_call2(SYM_Times, wn, ds_call2(SYM_Power, det, expr_new_integer(-1))))));
    }
    for (int i = 0; i < 3; i++){ expr_free(Ai[i]); expr_free(Bi[i]); expr_free(Ci[i]); }
    if (k && (!hg_finite(k) || !hg_finite(w)
              || !ds_free_of(k, x) || !ds_free_of(k, Y) || !ds_free_of(k, p)
              || !ds_free_of(w, x) || !ds_free_of(w, Y) || !ds_free_of(w, p))) {
        expr_free(k); expr_free(w); k = w = NULL;
    }
    if (!k) return false;
    *kout = k; *wout = w;
    return true;
}

/* Find the scaling weights (k, w) of R(x, Y, p), checked identically.  Returns
 * true with *k, *w owned.  The first sample set fixes x = Y = 1, so a symbolic
 * exponent (y^n, x^(2n-1)) contributes 1 rather than an uncancellable (3/4)^n,
 * while its derivative factors still carry n; the general set is the fallback. */
static bool hg_weights(const Expr* R, const char* x, const char* Y, const char* p,
                       Expr** kout, Expr** wout){
    Expr* Rx = ds_d(expr_copy((Expr*)R), expr_new_symbol(x));
    Expr* RY = ds_d(expr_copy((Expr*)R), expr_new_symbol(Y));
    Expr* Rp = ds_d(expr_copy((Expr*)R), expr_new_symbol(p));
    /* A = x R_x - p R_p,  B = Y R_Y + p R_p,  C = R:  A + k B - w C == 0 */
    Expr* pRp = eval_and_free(ds_call2(SYM_Times, expr_new_symbol(p), Rp));
    Expr* A = eval_and_free(ds_call2(SYM_Subtract,
                  eval_and_free(ds_call2(SYM_Times, expr_new_symbol(x), Rx)), expr_copy(pRp)));
    Expr* B = eval_and_free(ds_call2(SYM_Plus,
                  eval_and_free(ds_call2(SYM_Times, expr_new_symbol(Y), RY)), pRp));

    static const int64_t unit[3][6] = { {1,1, 1,1, 2,13}, {1,1, 1,1, 11,9}, {1,1, 1,1, 4,17} };
    static const int64_t gen[3][6]  = { {3,7, 5,11, 2,13}, {7,5, 3,4, 11,9}, {5,3, 13,7, 4,17} };
    Expr *k = NULL, *w = NULL;
    bool found = false;
    for (int set = 0; set < 2 && !found; set++){
        if (!hg_sample_weights(A, B, R, x, Y, p, set == 0 ? unit : gen, &k, &w)) continue;
        /* the identity, checked symbolically: A + k B - w C == 0 */
        Expr* Z = eval_and_free(ds_call1(SYM_Together, ds_call2(SYM_Plus, expr_copy(A),
                      ds_call2(SYM_Subtract, ds_call2(SYM_Times, expr_copy(k), expr_copy(B)),
                                             ds_call2(SYM_Times, expr_copy(w), expr_copy((Expr*)R))))));
        found = ds_is_zero(Z);
        expr_free(Z);
        if (!found) { expr_free(k); expr_free(w); k = w = NULL; }
    }
    expr_free(A); expr_free(B);
    if (!found) return false;
    *kout = k; *wout = w;
    return true;
}

/* The core: weights, the autonomous u(t) ODE, and its solve.  On success with an
 * all-explicit answer returns the y-bodies; otherwise (every branch converted)
 * stashes the first-integral LHSs under the residual hash and returns NULL. */
Expr** dsolve_homog_g_try(DSolveProblem* P, size_t* nbranch){
    if (getenv("MATHILDA_NO_HOMOG_G")) return NULL;     /* A/B killswitch */
    if (g_dsolve_depth > 1 || g_hg_active) return NULL;  /* top level, no re-entry */
    if (P->nfun != 1 || P->neq != 1 || P->max_order[0] != 1) return NULL;

    const char* yname = P->fun_names[0];
    const char* xv    = P->ind_names[0];
    uint64_t h = expr_hash(P->eq_residuals[0]);
    hg_memo_sync(eval_toplevel_id());
    if (hg_memo_seen(h)) return NULL;
    hg_memo_add(h);                    /* one attempt per top-level evaluation */

    const char* Yn = intern_symbol("DSolve`hgY");
    const char* Pn = intern_symbol("DSolve`hgP");
    Expr* R = dsolve_algebraic_residual(P, Yn, Pn);
    if (!R) return NULL;
    if (!ds_contains(R, Yn) || !ds_contains(R, Pn) || ds_contains(R, yname)
        || !hg_algebraic(R, xv, Yn, Pn)) { expr_free(R); return NULL; }

    Expr *k = NULL, *w = NULL;
    if (!hg_weights(R, xv, Yn, Pn, &k, &w)) { expr_free(R); return NULL; }
    Expr* km1 = eval_and_free(ds_call2(SYM_Subtract, expr_copy(k), expr_new_integer(1)));
    bool class_a = ds_is_zero(km1);
    expr_free(km1);
    if (class_a) { expr_free(R); expr_free(k); expr_free(w); return NULL; }   /* DSolve`Homogeneous */

    /* F~(U, Q) = x^-w R(x, x^k U, x^(k-1) (k U + Q)), which must be free of x */
    const char* Un = intern_symbol("DSolve`hgU");
    const char* Qn = intern_symbol("DSolve`hgQ");
    Expr* xk   = ds_call2(SYM_Power, expr_new_symbol(xv), expr_copy(k));
    Expr* xk1  = ds_call2(SYM_Power, expr_new_symbol(xv),
                     ds_call2(SYM_Subtract, expr_copy(k), expr_new_integer(1)));
    Expr* S = ds_subst(R, expr_new_symbol(Pn),
                  ds_call2(SYM_Times, xk1, ds_call2(SYM_Plus,
                      ds_call2(SYM_Times, expr_copy(k), expr_new_symbol(Un)), expr_new_symbol(Qn))));
    S = ds_subst(S, expr_new_symbol(Yn), ds_call2(SYM_Times, xk, expr_new_symbol(Un)));
    S = ds_call2(SYM_Times, S, ds_call2(SYM_Power, expr_new_symbol(xv),
                                         ds_call2(SYM_Times, expr_new_integer(-1), w)));
    Expr* Ft = eval_and_free(ds_call1(SYM_Numerator, eval_and_free(ds_call1(SYM_Together,
                   eval_and_free(ds_call1("PowerExpand", eval_and_free(S)))))));
    if (!ds_free_of(Ft, xv) || !ds_contains(Ft, Qn) || !ds_contains(Ft, Un)) {
        expr_free(Ft); expr_free(k); return NULL;
    }

    /* the autonomous ODE  F~(u[t], u'[t]) == 0 */
    const char* un = intern_symbol("DSolve`hgu");
    const char* tn = intern_symbol("DSolve`hgt");
    Ft = ds_subst(Ft, expr_new_symbol(Qn), ds_make_funcapp(un, 1, tn));
    Ft = ds_subst(Ft, expr_new_symbol(Un), ds_make_funcapp(un, 0, tn));
    Expr* ode = expr_new_function(expr_new_symbol(SYM_Equal), (Expr*[]){ Ft, expr_new_integer(0) }, 2);
    Expr* call = expr_new_function(expr_new_symbol(SYM_DSolve),
                     (Expr*[]){ ode, ds_make_funcapp(un, 0, tn), expr_new_symbol(tn) }, 3);
    call = expr_new_function(expr_new_symbol(intern_symbol("TimeConstrained")),
               (Expr*[]){ call, expr_new_integer(HG_SUBSOLVE_SEC),
                          expr_new_symbol(intern_symbol("$Aborted")) }, 3);
    g_hg_active++;
    mth_msg_suppress_push();
    Expr* r = eval_and_free(call);
    mth_msg_suppress_pop();
    g_hg_active--;

    size_t nbr = head_is(r, SYM_List) ? r->data.function.arg_count : 0;
    if (nbr == 0) { expr_free(r); expr_free(k); return NULL; }

    /* back-map: u[t] -> y[x] / x^k, t -> Log[x] */
    Expr* ymap = ds_call2(SYM_Times, ds_make_funcapp(yname, 0, xv),
                     ds_call2(SYM_Power, expr_new_symbol(xv),
                              ds_call2(SYM_Times, expr_new_integer(-1), expr_copy(k))));
    Expr* tmap = ds_call1(SYM_Log, expr_new_symbol(xv));

    size_t nb = 0;
    Expr** ub = dsolve_extract_applied_bodies(r, un, &nb);
    Expr** out = NULL;
    if (ub && nb == nbr) {
        /* all explicit: y = x^k u(Log[x]) */
        out = malloc(nb * sizeof(Expr*));
        for (size_t i = 0; i < nb; i++) {
            Expr* b = ds_subst(ub[i], expr_new_symbol(tn), expr_copy(tmap));
            out[i] = eval_and_free(ds_call2(SYM_Times,
                         ds_call2(SYM_Power, expr_new_symbol(xv), expr_copy(k)), b));
        }
        free(ub);
        *nbranch = nb;
    } else {
        if (ub) { for (size_t i = 0; i < nb; i++) expr_free(ub[i]); free(ub); }
        Expr** rel = NULL; size_t nrel = 0;
        /* An inert Integrate[f(u[t]), u[t]] cannot be back-mapped by substitution:
         * u[t] -> y[x]/x^k would rewrite the integration-variable slot too (a
         * malformed Integrate[.., x^(2/3) y[x]]).  Decline such a family. */
        if (!ds_has_head(r, SYM_Integrate) && !ds_has_inactive_integrate(r)
            && dsolve_result_relations(r, un, tn, &rel, &nrel)) {
            for (size_t i = 0; i < nrel; i++) {
                Expr* G = ds_subst(rel[i], ds_make_funcapp(un, 0, tn), expr_copy(ymap));
                rel[i] = ds_subst(G, expr_new_symbol(tn), expr_copy(tmap));
            }
            dsolve_stash_put(HG_STASH_TAG, h, rel, nrel);
        }
    }
    expr_free(r); expr_free(k); expr_free(ymap); expr_free(tmap);
    return out;
}

/* The implicit pass: consume what the explicit pass stashed for this residual
 * (pinned on its own, run the explicit pass first to populate it). */
static Expr** hg_implicit_impl(DSolveProblem* P, size_t* nbranch, bool may_solve){
    if (P->nfun != 1 || P->neq != 1 || P->max_order[0] != 1) return NULL;
    uint64_t h = expr_hash(P->eq_residuals[0]);
    size_t n = 0;
    Expr** G = dsolve_stash_take(HG_STASH_TAG, h, &n);
    if (!G && may_solve) {
        size_t ne = 0;
        Expr** ex = dsolve_homog_g_try(P, &ne);
        if (ex) { for (size_t i = 0; i < ne; i++) expr_free(ex[i]); free(ex); return NULL; }
        G = dsolve_stash_take(HG_STASH_TAG, h, &n);
    }
    if (!G) return NULL;
    *nbranch = n;
    return G;
}
Expr** dsolve_homog_g_implicit_try(DSolveProblem* P, size_t* nbranch){ return hg_implicit_impl(P, nbranch, false); }
static Expr** hg_implicit_pinned_try(DSolveProblem* P, size_t* nbranch){ return hg_implicit_impl(P, nbranch, true); }

static Expr* builtin_dsolve_homog_g(Expr* res){ return dsolve_method_builtin(res, dsolve_homog_g_try); }
static Expr* builtin_dsolve_homog_g_implicit(Expr* res){ return dsolve_method_builtin_implicit(res, hg_implicit_pinned_try); }

void dsolve_homog_g_init(void){
    tc_register_unwind_counter(&g_hg_active);   /* restored on a timeout unwind */
    symtab_add_builtin("DSolve`HomogeneousG", builtin_dsolve_homog_g);
    symtab_get_def("DSolve`HomogeneousG")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("DSolve`HomogeneousG",
        "DSolve`HomogeneousG[eqn, y, x] solves a first-order ODE F(x, y, y') == 0 that is "
        "isobaric (homogeneous class G: invariant under x -> l x, y -> l^k y up to a "
        "factor) by y = x^k u, t = Log[x], which makes it autonomous in u(t); the "
        "explicit answer y = x^k u(Log[x]).");
    symtab_add_builtin("DSolve`HomogeneousGImplicit", builtin_dsolve_homog_g_implicit);
    symtab_get_def("DSolve`HomogeneousGImplicit")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("DSolve`HomogeneousGImplicit",
        "DSolve`HomogeneousGImplicit[eqn, y, x] is DSolve`HomogeneousG when the autonomous "
        "u(t) equation solves only implicitly: one verified first integral "
        "G[x, y[x]] == C[1] per branch, with u -> y/x^k and t -> Log[x].");
}
