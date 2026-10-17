/*
 * dsolve_exact.c — DSolve`Exact.
 *
 * Solves  M(x,y) + N(x,y) y' == 0  when it is exact (M_y == N_x), or can be made
 * exact by an integrating factor mu(x) or mu(y):
 *   - (M_y - N_x)/N free of y  ->  mu(x) = Exp[Integrate[(M_y-N_x)/N, x]]
 *   - (N_x - M_y)/M free of x  ->  mu(y) = Exp[Integrate[(N_x-M_y)/M, y]]
 * The potential F with F_x = M, F_y = N is F = Integrate[M, x] + g(y) where
 * g'(y) = N - d/dy Integrate[M, x]; the solution is F(x,y) == C[1], solved for y.
 *
 * Works on the algebraic residual R(x, Y, p) (p = y'), which must be linear in p:
 * N = R_p, M = R|p=0.  Every returned branch is still back-substitution verified
 * by the substrate, so an imperfect potential construction cannot ship a wrong
 * answer — it is dropped.
 */
#include "dsolve_common.h"
#include "../common.h"
#include "integrate.h"          /* g_integrate_quiet */
#include "../sym_names.h"
#include "../eval.h"
#include "../sym_intern.h"
#include "../symtab.h"
#include "../attr.h"
#include <stdlib.h>
#include <time.h>

static Expr* pw_inv(Expr* base) {
    return expr_new_function(expr_new_symbol(SYM_Power), (Expr*[]){ base, expr_new_integer(-1) }, 2);
}

/* True iff `e` is a polynomial in x and Y. */
static bool poly_in_xy(const Expr* e, const char* xvar, const char* Yn) {
    Expr* vars = expr_new_function(expr_new_symbol(SYM_List),
                     (Expr*[]){ expr_new_symbol(xvar), expr_new_symbol(Yn) }, 2);
    Expr* q = eval_and_free(expr_new_function(expr_new_symbol(SYM_PolynomialQ),
                  (Expr*[]){ expr_copy((Expr*)e), vars }, 2));
    bool r = (q->type == EXPR_SYMBOL && q->data.symbol.name == SYM_True);
    expr_free(q);
    return r;
}

/* Propose an integrating factor mu = x^a y^b (constant exponents) and accept it
 * only after a symbolic exactness re-check.  For M + N y' == 0 the exactness
 * condition under mu = x^a y^b is
 *     b (M/y) - a (N/x) + (M_y - N_x) == 0   (identically in x, y),
 * which is linear in (a, b).  We propose (a, b) by sampling that identity at a
 * few generic non-zero points (a heuristic), then require (mu M)_y == (mu N)_x
 * symbolically before returning mu; the substrate's back-substitution verify is
 * the final backstop, so a wrong proposal is dropped, never shipped.  Solves
 * e.g. (x y - 2 x) y' == y - y^2 + 3 x^2 y^3  (mu = x^-2 y^-3).  M, N, diff
 * (= M_y - N_x) are borrowed; returns mu (owned) or NULL. */
static Expr* exact_xayb_factor(const Expr* M, const Expr* N, const Expr* diff,
                               const char* xvar, const char* Yn) {
    /* x^a y^b is a polynomial-ODE tool: gate on M, N being polynomials in x, y so
     * the speculative Solve below never runs on the heavy radical/transcendental
     * systems the AutonomousReduction recursion feeds into Exact (that Solve on
     * Sqrt-coefficient systems is what made a whole test run blow past its
     * budget). */
    if (!poly_in_xy(M, xvar, Yn) || !poly_in_xy(N, xvar, Yn)) return NULL;

    const char* aSym = intern_symbol("DSolve`exA");
    const char* bSym = intern_symbol("DSolve`exB");

    /* cond = a*(-N/x) + b*(M/y) + diff */
    Expr* Pc = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times), (Expr*[]){
                   expr_new_integer(-1), expr_copy((Expr*)N), pw_inv(expr_new_symbol(xvar)) }, 3));
    Expr* Qc = eval_and_free(ds_call2(SYM_Times, expr_copy((Expr*)M), pw_inv(expr_new_symbol(Yn))));
    Expr* cond = eval_and_free(ds_call2(SYM_Plus,
        ds_call2(SYM_Times, expr_new_symbol(aSym), Pc),          /* consumes Pc */
        ds_call2(SYM_Plus,
            ds_call2(SYM_Times, expr_new_symbol(bSym), Qc),      /* consumes Qc */
            expr_copy((Expr*)diff))));

    /* sample at generic non-zero points -> linear system {cond(pt_i)==0} in (a,b) */
    static const int sx[] = { 2, 3, 5 };
    static const int sy[] = { 3, 5, 7 };
    Expr* eqs[3];
    for (int i = 0; i < 3; i++) {
        Expr* ci = ds_subst(
            ds_subst(expr_copy(cond), expr_new_symbol(xvar), expr_new_integer(sx[i])),
            expr_new_symbol(Yn), expr_new_integer(sy[i]));
        eqs[i] = expr_new_function(expr_new_symbol(SYM_Equal),
                     (Expr*[]){ ci, expr_new_integer(0) }, 2);
    }
    expr_free(cond);
    Expr* eqlist  = expr_new_function(expr_new_symbol(SYM_List), eqs, 3);
    Expr* varlist = expr_new_function(expr_new_symbol(SYM_List),
                        (Expr*[]){ expr_new_symbol(aSym), expr_new_symbol(bSym) }, 2);
    Expr* solres = ds_solve(eqlist, varlist);      /* consumes eqlist, varlist */

    /* extract a, b from the first solution {{a -> ..., b -> ...}} */
    Expr* aval = NULL; Expr* bval = NULL;
    if (solres && solres->type == EXPR_FUNCTION && solres->data.function.arg_count >= 1) {
        Expr* first = solres->data.function.args[0];
        if (first && first->type == EXPR_FUNCTION) {
            for (size_t k = 0; k < first->data.function.arg_count; k++) {
                Expr* rule = first->data.function.args[k];
                if (rule && rule->type == EXPR_FUNCTION && rule->data.function.arg_count == 2
                    && rule->data.function.head->type == EXPR_SYMBOL
                    && (rule->data.function.head->data.symbol.name == SYM_Rule
                        || rule->data.function.head->data.symbol.name == SYM_RuleDelayed)) {
                    Expr* lhs = rule->data.function.args[0];
                    Expr* rhs = rule->data.function.args[1];
                    if (lhs->type == EXPR_SYMBOL && lhs->data.symbol.name == aSym && !aval) aval = expr_copy(rhs);
                    else if (lhs->type == EXPR_SYMBOL && lhs->data.symbol.name == bSym && !bval) bval = expr_copy(rhs);
                }
            }
        }
    }
    if (solres) expr_free(solres);
    if (!aval || !bval) { if (aval) expr_free(aval); if (bval) expr_free(bval); return NULL; }

    /* mu = x^a y^b */
    Expr* mu = eval_and_free(ds_call2(SYM_Times,
        eval_and_free(ds_call2(SYM_Power, expr_new_symbol(xvar), aval)),   /* consumes aval */
        eval_and_free(ds_call2(SYM_Power, expr_new_symbol(Yn), bval))));   /* consumes bval */

    /* symbolic exactness re-check: (mu M)_y - (mu N)_x == 0 */
    Expr* mM = eval_and_free(ds_call2(SYM_Times, expr_copy(mu), expr_copy((Expr*)M)));
    Expr* mN = eval_and_free(ds_call2(SYM_Times, expr_copy(mu), expr_copy((Expr*)N)));
    Expr* chk = eval_and_free(ds_call2(SYM_Subtract,
        ds_d(mM, expr_new_symbol(Yn)),          /* consumes mM */
        ds_d(mN, expr_new_symbol(xvar))));      /* consumes mN */
    bool exact = ds_is_zero(chk);
    expr_free(chk);
    if (!exact) { expr_free(mu); return NULL; }
    return mu;
}

/* Clear the residual's common denominator: an INEXACT form  M + N y' == 0  whose
 * coefficients are rational in x, y (a coefficient pole 1/x, 1/y, or the
 * y' == -P/Q spelling) can become an exact  P dx + Q dy  after multiplying
 * through by the common denominator D.  D is then an integrating factor, and the
 * cleared potential is a valid first integral of the original (D == 0 loci are
 * dropped by the substrate's back-substitution verify).  Only reached when the
 * raw form and every mu(x)/mu(y)/x^a y^b candidate already failed, so a
 * genuinely-exact form — which clearing WOULD break — never enters here.  On
 * success writes the owned cleared (Mm, Nn) and returns true; else false.
 * M, N are borrowed. */
/* (Mm, Nn)_exact test: M_y == N_x.  Mm, Nn borrowed. */
static bool exact_pair_is_exact(const Expr* Mm, const Expr* Nn,
                                const char* xvar, const char* Yn) {
    Expr* d = eval_and_free(ds_call2(SYM_Subtract,
                  ds_d(expr_copy((Expr*)Mm), expr_new_symbol(Yn)),
                  ds_d(expr_copy((Expr*)Nn), expr_new_symbol(xvar))));
    bool ex = ds_is_zero(d);
    expr_free(d);
    return ex;
}

/* Multiply (M, N) by candidate integrating factor D (consumed) and accept the
 * cleared pair only if it is exact.  Returns true (writes Mm_out, Nn_out) or
 * false.  A D free of both variables is a no-op factor and rejected. */
static bool exact_apply_factor(Expr* D, const Expr* M, const Expr* N,
                               const char* xvar, const char* Yn,
                               Expr** Mm_out, Expr** Nn_out) {
    bool trivial = ds_free_of(D, xvar) && ds_free_of(D, Yn);
    if (trivial) { expr_free(D); return false; }
    Expr* Mm = eval_and_free(ds_call1("Expand",
                   ds_call2(SYM_Times, expr_copy(D), expr_copy((Expr*)M))));
    Expr* Nn = eval_and_free(ds_call1("Expand",
                   ds_call2(SYM_Times, D, expr_copy((Expr*)N))));   /* consumes D */
    if (exact_pair_is_exact(Mm, Nn, xvar, Yn)) { *Mm_out = Mm; *Nn_out = Nn; return true; }
    expr_free(Mm); expr_free(Nn);
    return false;
}

static bool exact_clear_denominators(const Expr* M, const Expr* N, const char* Pn,
                                     const char* xvar, const char* Yn,
                                     Expr** Mm_out, Expr** Nn_out) {
    /* Candidate A — the SYNTACTIC denominators of M and N.  Handles the y' == P/Q
     * spelling (N == 1, M a ratio) even when P, Q carry a negative exponential
     * like E^(-x), which the Together route (B) mis-factors as a spurious E^x
     * denominator.  Tried first because it needs no Together. */
    Expr* Da = eval_and_free(ds_call2(SYM_Times,
                   ds_call1("Denominator", expr_copy((Expr*)M)),
                   ds_call1("Denominator", expr_copy((Expr*)N))));
    if (exact_apply_factor(Da, M, N, xvar, Yn, Mm_out, Nn_out)) return true;

    /* Candidate B — the common denominator via Together.  Handles summed rational
     * coefficients (1/x, 1/y) that A, seeing a bare Plus, leaves at 1.  Only
     * computed when A fails (Together can be costly on transcendental atoms). */
    Expr* R = eval_and_free(ds_call2(SYM_Plus, expr_copy((Expr*)M),
                  ds_call2(SYM_Times, expr_copy((Expr*)N), expr_new_symbol(Pn))));
    Expr* Db = eval_and_free(ds_call1("Denominator", eval_and_free(ds_call1("Together", R))));
    return exact_apply_factor(Db, M, N, xvar, Yn, Mm_out, Nn_out);
}

/* True when e is algebraic in x and Y: every subexpression involving them is a
 * Plus, Times or Power with an exponent free of them.  Gates mu(u), whose
 * Jacobian zero tests on a trig-rational M, N cost seconds and find nothing
 * (§2.2.53-5238: 7.4 s across the explicit and implicit entries).  e borrowed. */
static bool exact_algebraic(const Expr* e, const char* x, const char* Y) {
    if (!e || e->type != EXPR_FUNCTION) return true;
    if (ds_free_of(e, x) && ds_free_of(e, Y)) return true;
    if (head_is(e, SYM_Plus) || head_is(e, SYM_Times)) {
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            if (!exact_algebraic(e->data.function.args[i], x, Y)) return false;
        return true;
    }
    if (head_is(e, SYM_Power) && e->data.function.arg_count == 2) {
        const Expr* ex = e->data.function.args[1];
        if (!(ds_free_of(ex, x) && ds_free_of(ex, Y))) return false;
        return exact_algebraic(e->data.function.args[0], x, Y);
    }
    return false;
}

/* Integrating factor mu(u) depending on one combination u of x and Y (Kamke /
 * Murphy's classical test).  (mu M)_Y == (mu N)_x gives
 *        mu'/mu == Q,   Q = (M_Y - N_x) / (u_x N - u_Y M),
 * so a factor exists iff Q is a function of u alone -- i.e. the Jacobian
 * Q_x u_Y - Q_Y u_x vanishes identically.  Then phi(u) = Q is read off along the
 * line x = 1 (Y = Y(1, u)), mu = Exp[Integrate[phi, u]] with u -> u(x, Y), and the
 * product pair must be exact (checked).  Candidates: x + Y, x - Y, x Y, x^2 + Y^2,
 * x^2 - Y^2 -- the factors behind (x^2 + y^2 + x) y' == y (mu = 1/(x^2+y^2)) and
 * the 2nd-type Abel rationals of §2.2.52/53.  NOT Y/x: every homogeneous (class A)
 * equation has mu(Y/x) = 1/(x M + Y N), so that candidate would claim the whole
 * class ahead of DSolve`Homogeneous, which owns it (§2.2.52-5120, §2.2.53-5230,
 * §2.2.54-5349 regressed to timeouts).  Returns the factor (owned) or NULL.
 * M, N, diff borrowed. */
static Expr* exact_mu_u_factor(const Expr* M, const Expr* N, const Expr* diff,
                               const char* xvar, const char* Yn) {
    if (getenv("MATHILDA_NO_EXACT_MU_U")) return NULL;   /* A/B killswitch */
    if (!exact_algebraic(M, xvar, Yn) || !exact_algebraic(N, xvar, Yn)) return NULL;
    /* A homogeneous (class A) equation -- M/N of degree 0, i.e. free of x after
     * Y -> v x -- has a mu(u) for several u (x + Y as well as Y/x); it belongs to
     * DSolve`Homogeneous, so leave it there (§2.2.52-5120, §2.2.53-5248). */
    {
        const char* Vn = intern_symbol("DSolve`exV");
        Expr* vx = ds_call2(SYM_Times, expr_new_symbol(Vn), expr_new_symbol(xvar));
        /* sequenced: the two substitutions must not share vx within one call's
         * (unordered) argument evaluation */
        Expr* Mv = ds_subst(expr_copy((Expr*)M), expr_new_symbol(Yn), expr_copy(vx));
        Expr* Nv = ds_subst(expr_copy((Expr*)N), expr_new_symbol(Yn), vx);
        Expr* r = eval_and_free(ds_call1(SYM_Together, ds_call2(SYM_Times, Mv,
                      ds_call2(SYM_Power, Nv, expr_new_integer(-1)))));
        bool homogeneous = ds_free_of(r, xvar);
        expr_free(r);
        if (homogeneous) return NULL;
    }
    const char* Un = intern_symbol("DSolve`exU");
    Expr* X = expr_new_symbol(xvar);
    Expr* Y = expr_new_symbol(Yn);
    Expr* U = expr_new_symbol(Un);
    /* candidate u(x, Y), and Y(1, u) on the line x = 1 */
    Expr* cand[5][2] = {
        { ds_call2(SYM_Plus, expr_copy(X), expr_copy(Y)),
          ds_call2(SYM_Subtract, expr_copy(U), expr_new_integer(1)) },
        { ds_call2(SYM_Subtract, expr_copy(X), expr_copy(Y)),
          ds_call2(SYM_Subtract, expr_new_integer(1), expr_copy(U)) },
        { ds_call2(SYM_Times, expr_copy(X), expr_copy(Y)),
          expr_copy(U) },
        { ds_call2(SYM_Plus, ds_call2(SYM_Power, expr_copy(X), expr_new_integer(2)),
                             ds_call2(SYM_Power, expr_copy(Y), expr_new_integer(2))),
          ds_call1(SYM_Sqrt, ds_call2(SYM_Subtract, expr_copy(U), expr_new_integer(1))) },
        { ds_call2(SYM_Subtract, ds_call2(SYM_Power, expr_copy(X), expr_new_integer(2)),
                                 ds_call2(SYM_Power, expr_copy(Y), expr_new_integer(2))),
          ds_call1(SYM_Sqrt, ds_call2(SYM_Subtract, expr_new_integer(1), expr_copy(U))) },
    };
    expr_free(X); expr_free(Y); expr_free(U);
    Expr* Mu = NULL;
    for (int c = 0; c < 5; c++) {
        Expr* u  = eval_and_free(cand[c][0]);
        Expr* yu = eval_and_free(cand[c][1]);
        if (Mu) { expr_free(u); expr_free(yu); continue; }
        Expr* ux = ds_d(expr_copy(u), expr_new_symbol(xvar));
        Expr* uy = ds_d(expr_copy(u), expr_new_symbol(Yn));
        Expr* den = eval_and_free(ds_call2(SYM_Subtract,
                        ds_call2(SYM_Times, expr_copy(ux), expr_copy((Expr*)N)),
                        ds_call2(SYM_Times, expr_copy(uy), expr_copy((Expr*)M))));
        if (ds_is_zero(den)) { expr_free(den); expr_free(ux); expr_free(uy); expr_free(u); expr_free(yu); continue; }
        Expr* Q = eval_and_free(ds_call1(SYM_Together, ds_call2(SYM_Times, expr_copy((Expr*)diff),
                      ds_call2(SYM_Power, den, expr_new_integer(-1)))));
        /* Jacobian test: Q is a function of u alone */
        Expr* jac = eval_and_free(ds_call1(SYM_Together, ds_call2(SYM_Subtract,
                        ds_call2(SYM_Times, ds_d(expr_copy(Q), expr_new_symbol(xvar)), expr_copy(uy)),
                        ds_call2(SYM_Times, ds_d(expr_copy(Q), expr_new_symbol(Yn)), expr_copy(ux)))));
        bool fn_of_u = ds_is_zero(jac) && !ds_is_zero(Q);
        expr_free(jac); expr_free(ux); expr_free(uy);
        if (fn_of_u) {
            Expr* phi = ds_subst(ds_subst(expr_copy(Q), expr_new_symbol(xvar), expr_new_integer(1)),
                                 expr_new_symbol(Yn), expr_copy(yu));
            phi = eval_and_free(ds_call1(SYM_Together, phi));
            if (ds_free_of(phi, xvar) && ds_free_of(phi, Yn)) {
                g_integrate_quiet++;
                Expr* I = ds_integrate(phi, expr_new_symbol(Un));
                g_integrate_quiet--;
                /* ALGEBRAIC factors only: Integrate[phi, u] must be a log/rational
                 * combination, so mu = Exp[..] is a product of powers.  A
                 * transcendental mu (an ArcTan in the exponent -- the mu(Y/x) of a
                 * homogeneous y y' + a x + b y == 0, §2.2.51-5038) is genuine but
                 * makes the potential integrals below spin, where the methods after
                 * Exact answer in 0.2 s. */
                static const char* const transc[] = { "ArcTan", "ArcTanh", "ArcCot", "ArcCoth",
                    "ArcSin", "ArcCos", "ArcSinh", "ArcCosh", "Erf", "Erfi", "ExpIntegralEi",
                    "RootSum", "Root" };
                bool algebraic = !ds_has_head(I, SYM_Integrate);
                for (size_t t = 0; t < sizeof(transc) / sizeof(transc[0]) && algebraic; t++)
                    if (ds_has_head(I, intern_symbol(transc[t]))) algebraic = false;
                if (algebraic) {
                    Expr* mu = eval_and_free(ds_call1("Exp", ds_subst(I, expr_new_symbol(Un), expr_copy(u))));
                    Expr* Mm = eval_and_free(ds_call2(SYM_Times, expr_copy(mu), expr_copy((Expr*)M)));
                    Expr* Nn = eval_and_free(ds_call2(SYM_Times, expr_copy(mu), expr_copy((Expr*)N)));
                    if (exact_pair_is_exact(Mm, Nn, xvar, Yn)) Mu = mu;
                    else expr_free(mu);
                    expr_free(Mm); expr_free(Nn);
                } else expr_free(I);
            } else expr_free(phi);
        }
        expr_free(Q); expr_free(u); expr_free(yu);
    }
    return Mu;
}

/* Integrate[e, v] -- unbounded as before, or (`deadline` > 0) inside what is left
 * of that wall-clock deadline, a timeout reading as non-elementary (a decline).
 * The bound is used only for a mu(u) factor, whose potentials can be far harder
 * than the equation (§2.2.54-5320 spun past the corpus wall where the methods
 * after Exact answer in 0.5 s).  e consumed. */
#define EXACT_MU_U_BUDGET_SEC 2
static Expr* exact_integrate(Expr* e, const char* var, time_t deadline) {
    if (deadline <= 0) return ds_integrate(e, expr_new_symbol(var));
    long left = (long)(deadline - time(NULL));
    if (left < 1) { expr_free(e); return NULL; }
    Expr* call = expr_new_function(expr_new_symbol(SYM_Integrate),
                     (Expr*[]){ e, expr_new_symbol(var) }, 2);
    Expr* r = eval_and_free(expr_new_function(expr_new_symbol(SYM_TimeConstrained),
                  (Expr*[]){ call, expr_new_integer(left),
                             expr_new_symbol(intern_symbol("$Aborted")) }, 3));
    if (r && r->type == EXPR_SYMBOL && r->data.symbol.name == intern_symbol("$Aborted")) {
        expr_free(r); return NULL;
    }
    return r;
}
/* True when an exact_integrate result is unusable (timed out or non-elementary). */
static bool exact_int_failed(const Expr* r) { return !r || ds_has_head(r, SYM_Integrate); }

/* Build the exact potential F(x,Y) with F_x = Mu*M, F_y = Mu*N (Mu an
 * integrating factor: 1, mu(x), mu(y), or x^a y^b) for M + N y' == 0.  Returns
 * F (owned), or NULL when the residual is not linear in p, cannot be made
 * exact, or a potential integral is non-elementary.  Shared by the explicit
 * Exact entry (which then Solves F == C[1] for Y) and the implicit entry (which
 * returns F(x, y[x]) == C[1] for a transcendental, non-invertible F). */
/* mu_u_only selects the LATE pass (DSolve`ExactMuU): skip the classical factors
 * (1, mu(x), mu(y), x^a y^b -- the plain pass's question) and try only mu(u).  It is
 * a separate cascade slot because an exact-by-mu(u) relation is often only an
 * implicit answer where a later specialist closes explicitly: p' == (x + p)^2
 * (u = x + p) gives ArcTan[x + p] - x == C[1] where the explicit p = -x - Tan[..]
 * follows, and a second-order reduction (y'' == (x + y')^2, §2.2.2-154) needs the
 * explicit one. */
static Expr* exact_potential(DSolveProblem* P, const char* xvar, const char* Yn,
                             const char* Pn, bool mu_u_only) {
    Expr* R = dsolve_algebraic_residual(P, Yn, Pn);
    if (!R) return NULL;

    /* R = M + N p, N = R_p free of p, M = R|p=0 */
    Expr* N = ds_d(expr_copy(R), expr_new_symbol(Pn));
    if (!ds_free_of(N, Pn)) { expr_free(N); expr_free(R); return NULL; }
    Expr* M = ds_subst(expr_copy(R), expr_new_symbol(Pn), expr_new_integer(0));
    Expr* recon = eval_and_free(ds_call2(SYM_Plus, expr_copy(M),
                        ds_call2(SYM_Times, expr_copy(N), expr_new_symbol(Pn))));
    Expr* lin = eval_and_free(ds_call2(SYM_Subtract, expr_copy(R), recon));
    bool linear = ds_is_zero(lin);
    expr_free(lin); expr_free(R);
    if (!linear) { expr_free(M); expr_free(N); return NULL; }

    /* diff = M_y - N_x */
    Expr* My = ds_d(expr_copy(M), expr_new_symbol(Yn));
    Expr* Nx = ds_d(expr_copy(N), expr_new_symbol(xvar));
    Expr* diff = eval_and_free(ds_call2(SYM_Subtract, My, Nx));

    Expr* Mu = NULL;
    time_t pot_deadline = 0;           /* bounded potential integrals (mu(u) only) */
    if (mu_u_only) {
        if (!ds_is_zero(diff)) {
            Mu = exact_mu_u_factor(M, N, diff, xvar, Yn);
            if (Mu) pot_deadline = time(NULL) + EXACT_MU_U_BUDGET_SEC;
        }
    } else if (ds_is_zero(diff)) {
        Mu = expr_new_integer(1);
    } else {
        /* mu(x): (M_y - N_x)/N free of Y */
        Expr* r1 = eval_and_free(ds_call2(SYM_Times, expr_copy(diff), pw_inv(expr_copy(N))));
        if (ds_free_of(r1, Yn)) {
            Expr* r1int = ds_integrate(expr_copy(r1), expr_new_symbol(xvar));
            if (!ds_has_head(r1int, SYM_Integrate)) Mu = eval_and_free(ds_call1("Exp", r1int));
            else expr_free(r1int);
        }
        expr_free(r1);
        /* mu(y): (N_x - M_y)/M free of x.  Tried whenever mu(x) yielded no factor
         * -- its free-of test can mis-decide on a trig-rational derivative, or its
         * integral come back non-elementary, even when mu(y) succeeds cleanly
         * (e.g. E^x + (E^x Cot y + 2 y Csc y) y' == 0, mu = Sin y). */
        if (!Mu) {
            Expr* r2 = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times), (Expr*[]){
                expr_new_integer(-1), expr_copy(diff), pw_inv(expr_copy(M)) }, 3));
            if (ds_free_of(r2, xvar)) {
                Expr* r2int = ds_integrate(expr_copy(r2), expr_new_symbol(Yn));
                if (!ds_has_head(r2int, SYM_Integrate)) Mu = eval_and_free(ds_call1("Exp", r2int));
                else expr_free(r2int);
            }
            expr_free(r2);
        }
    }
    /* fallback: mu = x^a y^b (constant exponents) when neither mu(x) nor mu(y) fits */
    if (!Mu && !mu_u_only) Mu = exact_xayb_factor(M, N, diff, xvar, Yn);
    expr_free(diff);

    Expr* Mm, * Nn;
    if (Mu) {
        Mm = eval_and_free(ds_call2(SYM_Times, expr_copy(Mu), M));   /* consumes M */
        Nn = eval_and_free(ds_call2(SYM_Times, expr_copy(Mu), N));   /* consumes N */
        expr_free(Mu);
    } else {
        /* last resort: clear a rational form's common denominator (integrating
         * factor mu = D), recognising y' == -P/Q and 1/x, 1/y coefficient poles.
         * Restricted to the condition-free general solve: with an initial
         * condition the cleared potential is a Root-form explicit whose C[1] does
         * not inverse-fit, so an IVP that reaches here is left to whichever method
         * already fits it (the explicit Root would otherwise shadow it, unfitted). */
        bool ok = (P->ncond == 0) &&
                  exact_clear_denominators(M, N, Pn, xvar, Yn, &Mm, &Nn);
        expr_free(M); expr_free(N);
        if (!ok) return NULL;
    }

    /* Build the potential from whichever coefficient integrates cleanly.
     *  Path 1:  F = ∫Mm dx + g(Y),  g'(Y) = Nn - ∂_Y ∫Mm dx   (free of x, exact)
     *  Path 2:  F = ∫Nn dY + h(x),  h'(x) = Mm - ∂_x ∫Nn dY   (free of Y, exact)
     * The two are equal potentials, but one antiderivative can come back far
     * messier than the other — e.g. the E^(x y) exact family, where ∫Mm dx lands
     * in a Tan half-angle form that leaves g'(Y) only *cancelling* to a function
     * of Y (not syntactically free of x), so integrating it would hang.  Take the
     * path whose remainder derivative is already free of the other variable, with
     * NO Simplify on the hot path. */
    Expr* Fpot = NULL;
    Expr* Fx = exact_integrate(expr_copy(Mm), xvar, pot_deadline);
    if (!exact_int_failed(Fx)) {
        Expr* gp = eval_and_free(ds_call2(SYM_Subtract, expr_copy(Nn),
                       ds_d(expr_copy(Fx), expr_new_symbol(Yn))));
        if (ds_free_of(gp, xvar)) {
            Expr* g = exact_integrate(expr_copy(gp), Yn, pot_deadline);
            if (!exact_int_failed(g))
                Fpot = eval_and_free(ds_call2(SYM_Plus, expr_copy(Fx), g));
            else if (g) expr_free(g);
        }
        expr_free(gp);
    }
    if (Fx) expr_free(Fx);
    if (!Fpot) {                                    /* Path 2: integrate Nn dY */
        Expr* Fy = exact_integrate(expr_copy(Nn), Yn, pot_deadline);
        if (!exact_int_failed(Fy)) {
            Expr* hp = eval_and_free(ds_call2(SYM_Subtract, expr_copy(Mm),
                           ds_d(expr_copy(Fy), expr_new_symbol(xvar))));
            if (ds_free_of(hp, Yn)) {
                Expr* h = exact_integrate(expr_copy(hp), xvar, pot_deadline);
                if (!exact_int_failed(h))
                    Fpot = eval_and_free(ds_call2(SYM_Plus, expr_copy(Fy), h));
                else if (h) expr_free(h);
            }
            expr_free(hp);
        }
        if (Fy) expr_free(Fy);
    }
    expr_free(Mm); expr_free(Nn);
    return Fpot;                                     /* NULL if both paths messy */
}

static Expr** exact_explicit(DSolveProblem* P, size_t* nbranch, bool mu_u_only) {
    if (P->nfun != 1 || P->neq != 1) return NULL;
    if (P->max_order[0] != 1) return NULL;
    const char* xvar = P->ind_names[0];
    const char* Yn = intern_symbol("DSolve`Y");
    const char* Pn = intern_symbol("DSolve`p");

    Expr* Fpot = exact_potential(P, xvar, Yn, Pn, mu_u_only);
    if (!Fpot) return NULL;

    /* Gate the explicit inversion to a potential that is a RATIONAL function of Y.
     * Solve on a rational-in-Y equation reduces to a terminating polynomial solve,
     * but a potential carrying FRACTIONAL powers of Y (Y^(4/3), Y^(3/2) — the
     * class-G / radical exact ODEs, e.g. 2.2.3-204:
     *   6 x^(3/2) y^(4/3) - 10 x^(6/5) y^(3/2) == C) sends Solve into a
     * non-terminating radical-inversion search.  Those are exactly the potentials
     * the implicit entry returns verbatim as F(x, y[x]) == C[1] (as Maple and
     * Mathematica do), so decline here and let the cascade's
     * dsolve_exact_implicit_try own them. */
    if (!ds_is_rational_in(Fpot, Yn)) { expr_free(Fpot); return NULL; }

    /* explicit solution F(x, Y) == C[1], solved for Y */
    Expr* eq = expr_new_function(expr_new_symbol(SYM_Equal), (Expr*[]){ Fpot, ds_const(1) }, 2);
    Expr* solres = ds_solve(eq, expr_new_symbol(Yn));
    size_t nb = 0;
    Expr** bodies = dsolve_extract_solutions(solres, Yn, &nb);
    if (solres) expr_free(solres);
    if (!bodies) return NULL;                 /* transcendental F -> implicit entry */
    *nbranch = nb;
    return bodies;
}

/* Implicit first-integral entry: when Solve cannot invert F(x, Y) == C[1] for Y
 * (a transcendental potential, e.g. x + E^(x y) + y^2 for problem 140), return
 * the relation F(x, y[x]) == C[1] via dsolve_run_implicit, which verifies it by
 * the implicit-function rule y' == -F_x/F_y == -M/N (the ODE itself). */
static Expr** exact_implicit(DSolveProblem* P, size_t* nbranch, bool mu_u_only) {
    if (P->nfun != 1 || P->neq != 1) return NULL;
    if (P->max_order[0] != 1) return NULL;
    const char* xvar = P->ind_names[0];
    const char* yname = P->fun_names[0];
    const char* Yn = intern_symbol("DSolve`Y");
    const char* Pn = intern_symbol("DSolve`p");

    Expr* Fpot = exact_potential(P, xvar, Yn, Pn, mu_u_only);
    if (!Fpot) return NULL;
    Expr* G = ds_subst(Fpot, expr_new_symbol(Yn), ds_make_funcapp(yname, 0, xvar));
    Expr** out = malloc(sizeof(Expr*));
    out[0] = G;
    *nbranch = 1;
    return out;
}

Expr** dsolve_exact_try(DSolveProblem* P, size_t* nbranch) {
    return exact_explicit(P, nbranch, false);
}
Expr** dsolve_exact_implicit_try(DSolveProblem* P, size_t* nbranch) {
    return exact_implicit(P, nbranch, false);
}
/* The late mu(u) pass (see exact_potential); MATHILDA_NO_EXACT_MU_U disables it. */
Expr** dsolve_exact_mu_u_try(DSolveProblem* P, size_t* nbranch) {
    return exact_explicit(P, nbranch, true);
}
Expr** dsolve_exact_mu_u_implicit_try(DSolveProblem* P, size_t* nbranch) {
    return exact_implicit(P, nbranch, true);
}

static Expr* builtin_dsolve_exact_mu_u(Expr* res) {
    Expr* r = dsolve_method_builtin(res, dsolve_exact_mu_u_try);
    if (!r) r = dsolve_method_builtin_implicit(res, dsolve_exact_mu_u_implicit_try);
    return r;
}

static Expr* builtin_dsolve_exact(Expr* res) {
    Expr* r = dsolve_method_builtin(res, dsolve_exact_try);
    if (!r) r = dsolve_method_builtin_implicit(res, dsolve_exact_implicit_try);
    return r;
}

void dsolve_exact_init(void) {
    symtab_add_builtin("DSolve`Exact", builtin_dsolve_exact);
    symtab_get_def("DSolve`Exact")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("DSolve`Exact",
        "DSolve`Exact[eqn, y, x] solves M(x,y) + N(x,y) y' == 0 when exact "
        "(M_y == N_x), or made exact by an integrating factor mu(x), mu(y) or x^a y^b.");
    symtab_add_builtin("DSolve`ExactMuU", builtin_dsolve_exact_mu_u);
    symtab_get_def("DSolve`ExactMuU")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("DSolve`ExactMuU",
        "DSolve`ExactMuU[eqn, y, x] solves M + N y' == 0 made exact by an integrating "
        "factor mu(u) of one combination u in {x+y, x-y, x y, x^2+y^2, x^2-y^2} with an "
        "algebraic integral; a late cascade pass after DSolve`Exact.");
}
