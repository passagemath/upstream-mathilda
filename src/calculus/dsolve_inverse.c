/*
 * dsolve_inverse.c — DSolve`InverseOrientation.
 *
 * A first-order ODE  y'[x] == F(x, y)  that no specialist solves in the y(x)
 * orientation is often trivial in the x(y) orientation: its reciprocal
 *
 *     dx/dy == 1 / F(x, y)
 *
 * — read with x the dependent function of the independent variable y — is a
 * standard class the cascade already handles.  E.g. §2.2.45-4408
 * `y' == 1/(x y + x^3 y^3)` inverts to `x'[y] == x y + x^3 y^3`, a Bernoulli
 * equation in x(y); §2.2.45-4417 inverts to a linear ODE in x(y).  Maple and
 * SymPy both try this orientation; Mathilda's Lie backstop finds a symmetry of
 * the inverse ODE but not the far cheaper direct solve.
 *
 * Method.  Solve for y' to get F, build the reciprocal ODE with x<->y swapped,
 * and RECURSE the cascade on it (TimeConstrained).  An explicit inverse solution
 * x[y] == h(y, C[1]) is turned into a first integral by isolating the constant,
 * Solve[x == h, C[1]] -> C[1] == Phi(x, y); the relation Phi(x, y) == C[1] is the
 * implicit solution in the original orientation, verified by the substrate's
 * implicit-function rule against the ORIGINAL y' (dsolve_run_implicit).  So a
 * mis-inversion can only fail that verify and be dropped, never emit a wrong
 * answer.
 *
 * Runs LATE in the first-order cascade — it recurses the whole engine, so every
 * cheaper specialist (and the Lie backstop, whose cleaner symmetry form is
 * preferred) is tried first.  Re-entry guarded so the inner solve cannot recurse
 * back into the inverse orientation.  No-IVP only (the implicit relation is fitted
 * by the shared substrate when conditions are present, future work).
 */
#include "dsolve_common.h"
#include "../core.h"
#include "dsolve.h"               /* g_dsolve_depth: top-level-only gate */
#include "../sym_names.h"
#include "../sym_intern.h"
#include "../symtab.h"
#include "../attr.h"
#include "../eval.h"
#include "../common.h"            /* head_is */
#include "../message.h"           /* mth_msg_suppress_push/pop: this is a probe */
#include <stdlib.h>

static int g_inv_active = 0;      /* re-entry guard: no inverse inside inverse */

/* decline memo (the evaluator re-invokes a declining builtin ~3x per call) */
#define INV_MEMO_SLOTS 32
static uint64_t inv_epoch = 0;
static int inv_memo_n = 0;
static uint64_t inv_memo[INV_MEMO_SLOTS];
static void inv_memo_sync(uint64_t tid){ if(tid!=inv_epoch){inv_epoch=tid;inv_memo_n=0;} }
static bool inv_memo_seen(uint64_t h){ for(int i=0;i<inv_memo_n;i++) if(inv_memo[i]==h) return true; return false; }
static void inv_memo_add(uint64_t h){ if(inv_memo_n<INV_MEMO_SLOTS && !inv_memo_seen(h)) inv_memo[inv_memo_n++]=h; }

/* TimeConstrained[e, secs, $Aborted], evaluated inside a message-SUPPRESSED region
 * (the inverse solve is a speculative probe: its Solve::nongen / Integrate::nonelem
 * must not leak through Quiet[] or flip Check[], matching the nested-DSolve mute).
 * e consumed. */
static Expr* inv_beval(Expr* e, int secs){
    Expr* g = expr_new_function(expr_new_symbol("TimeConstrained"),
                  (Expr*[]){ e, expr_new_integer(secs), expr_new_symbol(intern_symbol("$Aborted")) }, 3);
    mth_msg_suppress_push();
    Expr* r = eval_and_free(g);
    mth_msg_suppress_pop();
    return r;
}

/* Solve DSolve[eqn, X[t], t] (TimeConstrained) -> the applied-form RHS body h(t),
 * or NULL if the inverse did not solve to a single explicit x[t] = body.  eqn
 * consumed. */
static Expr* inv_run_applied(Expr* eqn, const char* xname, const char* tvar){
    Expr* lhs  = ds_call1(xname, expr_new_symbol(tvar));
    Expr* call = expr_new_function(expr_new_symbol(SYM_DSolve),
                     (Expr*[]){ eqn, lhs, expr_new_symbol(tvar) }, 3);
    Expr* r = inv_beval(call, 5);
    Expr* body = NULL;
    size_t nbody = 0;
    if (head_is(r, SYM_List) && r->data.function.arg_count >= 1){
        Expr* inner = r->data.function.args[0];      /* first branch only */
        if (head_is(inner, SYM_List))
            for (size_t k=0;k<inner->data.function.arg_count;k++){
                Expr* rule = inner->data.function.args[k];
                if (head_is(rule, SYM_Rule) && rule->data.function.arg_count==2){
                    Expr* rl = rule->data.function.args[0];
                    if (rl->type==EXPR_FUNCTION && rl->data.function.head->type==EXPR_SYMBOL
                        && rl->data.function.head->data.symbol.name==xname){
                        if (body) expr_free(body);
                        body = expr_copy(rule->data.function.args[1]);
                        nbody++;
                    }
                }
            }
    }
    expr_free(r);
    if (nbody != 1) { if (body) expr_free(body); return NULL; }   /* want exactly one explicit x[t]=body */
    return body;
}

Expr** dsolve_inverse_try(DSolveProblem* P, size_t* nbranch){
    if (getenv("MATHILDA_NO_INVERSE")) return NULL;   /* A/B killswitch (isolates this method) */
    /* TOP-LEVEL ONLY.  This method recurses the whole engine; letting it fire inside
     * a nested DSolve (a change-of-variable / dAlembert sub-solve) adds its ~several-
     * second recursion to a case another method already answers -- §2.2.44-4323
     * (y'==Sin[x-y]^2) solves in 6 s via nonlinear exactness but, with the inverse
     * firing in its u=x-y sub-solve, timed out at 8 s.  A genuine x(y)-only equation
     * is a top-level problem, so gating here loses nothing. */
    if (g_dsolve_depth > 1) return NULL;
    if (P->nfun != 1 || P->neq != 1) return NULL;
    if (P->max_order[0] != 1) return NULL;
    if (P->ncond > 0) return NULL;            /* no-IVP only */
    if (g_inv_active) return NULL;            /* no inverse inside inverse */

    const char* yname = P->fun_names[0];      /* original dependent function */
    const char* xv    = P->ind_names[0];      /* original independent variable */

    inv_memo_sync(eval_toplevel_id());
    uint64_t h = expr_hash(P->eq_residuals[0]);
    if (inv_memo_seen(h)) return NULL;

    /* F: y'[x] == F(x, y[x]). */
    Expr* F = dsolve_solve_top_derivative(P, 1);
    if (!F) { inv_memo_add(h); return NULL; }

    /* omega(x, Y) with the applied y[x] replaced by a plain symbol Y. */
    const char* Yn = intern_symbol("DSolve`invY");
    Expr* omega = ds_subst(F, ds_make_funcapp(yname, 0, xv), expr_new_symbol(Yn));

    /* Reciprocal ODE for x(y): dx/dy == 1/F, with the two roles swapped --
     *   old dependent value Y  -> the new INDEPENDENT variable (we reuse `yname`),
     *   old independent x       -> the new DEPENDENT function  xv[yname].
     * omega is a function of (xv, Yn); rewrite and reciprocate. */
    Expr* om2 = ds_subst(omega, expr_new_symbol(Yn), expr_new_symbol(yname));       /* Y -> y (new indep) */
    om2 = ds_subst(om2, expr_new_symbol(xv), ds_make_funcapp(xv, 0, yname));          /* x -> x[y] (new dep) */
    Expr* invRHS = eval_and_free(ds_call1(SYM_Together,
                       expr_new_function(expr_new_symbol(SYM_Power),
                           (Expr*[]){ om2, expr_new_integer(-1) }, 2)));              /* 1/omega, one fraction */

    if (ds_contains(invRHS, Yn)) { expr_free(invRHS); inv_memo_add(h); return NULL; }  /* swap left a stray Y */

    /* inverse equation:  x'[y] == invRHS */
    Expr* invEq = expr_new_function(expr_new_symbol(SYM_Equal),
                      (Expr*[]){ ds_make_funcapp(xv, 1, yname), invRHS }, 2);

    g_inv_active = 1;
    Expr* hbody = inv_run_applied(invEq, xv, yname);   /* x[y] = hbody(y, C[1]) */
    g_inv_active = 0;
    if (!hbody) { inv_memo_add(h); return NULL; }

    /* Isolate the generated constant: Solve[x == hbody, C[1]] -> C[1] == Phi(x, y). */
    Expr* Ck = ds_const(1);     /* C[1], the inner DSolve's default generated parameter */
    if (!ds_contains(hbody, intern_symbol("C"))) { expr_free(hbody); expr_free(Ck); inv_memo_add(h); return NULL; }
    Expr* rel = expr_new_function(expr_new_symbol(SYM_Equal),
                    (Expr*[]){ expr_new_symbol(xv), hbody }, 2);         /* x == hbody, consumes hbody */
    Expr* sol = inv_beval(expr_new_function(expr_new_symbol(SYM_Solve),
                    (Expr*[]){ rel, expr_copy(Ck) }, 2), 4);             /* Solve[x == hbody, C[1]] */
    Expr* Phi = NULL;
    if (head_is(sol, SYM_List) && sol->data.function.arg_count >= 1){
        Expr* b0 = sol->data.function.args[0];                          /* first solution branch */
        if (head_is(b0, SYM_List) && b0->data.function.arg_count >= 1){
            Expr* rule = b0->data.function.args[0];
            if (head_is(rule, SYM_Rule) && rule->data.function.arg_count == 2
                && expr_eq(rule->data.function.args[0], Ck))
                Phi = expr_copy(rule->data.function.args[1]);
        }
    }
    expr_free(sol); expr_free(Ck);
    if (!Phi) { inv_memo_add(h); return NULL; }

    /* Map back to the original orientation: the new-indep symbol `yname` is the
     * y-value, so y -> y[x]; `xv` stays the independent variable.  G(x, y[x]) is
     * then a first integral and G == C[1] the implicit solution. */
    Expr* G = ds_subst(Phi, expr_new_symbol(yname), ds_make_funcapp(yname, 0, xv));   /* consumes Phi */
    if (ds_contains(G, "C") || !ds_free_of(G, "DSolve`invY")) { expr_free(G); inv_memo_add(h); return NULL; }
    G = ds_simplify(G);
    if (ds_is_zero(G) || ds_free_of(G, xv)) { expr_free(G); inv_memo_add(h); return NULL; }  /* degenerate */

    inv_memo_add(h);
    Expr** out = malloc(sizeof(Expr*));
    out[0] = G;
    *nbranch = 1;
    return out;
}

static Expr* builtin_dsolve_inverse(Expr* res){
    return dsolve_method_builtin_implicit(res, dsolve_inverse_try);
}

void dsolve_inverse_init(void){
    tc_register_unwind_counter(&g_inv_active);   /* restored on a timeout unwind */
    const char* doc =
        "DSolve`InverseOrientation[eqn, y, x] solves a first-order ODE by its "
        "reciprocal: dx/dy == 1/F for x(y), recursing the cascade (Bernoulli / "
        "linear / separable in x(y)), then returning the implicit first integral "
        "in the original orientation, verified by back-substitution.";
    symtab_add_builtin("DSolve`InverseOrientation", builtin_dsolve_inverse);
    symtab_get_def("DSolve`InverseOrientation")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("DSolve`InverseOrientation", doc);
}
