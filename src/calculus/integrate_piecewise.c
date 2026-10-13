/* integrate_piecewise.c
 *
 * Definite integration of piecewise integrands.  See integrate_piecewise.h for
 * the contract.  The integrand is a combination of the branch-defining heads
 * Abs / Sign / HeavisideTheta / UnitStep / Min / Max / Piecewise whose
 * arguments depend on the integration variable.  We find the breakpoints where
 * a branch switches, split [a, b] at the interior ones, reduce each sub-interval
 * to a single smooth branch (chosen by a *proved* sign of each governing
 * argument at an interior test point), integrate that with the ordinary
 * cascade, and sum.  Anything undecidable yields NULL (left unevaluated).
 */

#include "integrate_piecewise.h"
#include "expr.h"
#include "eval.h"
#include "symtab.h"
#include "attr.h"
#include "sym_names.h"

#include <string.h>
#include <stdbool.h>
#include <stdlib.h>

/* ---- small construction / evaluation helpers ---------------------------- */

static Expr* cp(const Expr* e) { return expr_copy((Expr*)e); }

static Expr* mk_fn(const char* head, Expr** args, size_t n) {
    return expr_new_function(expr_new_symbol(head), args, n);
}
static Expr* mk_fn1(const char* head, Expr* a) {
    Expr* args[1] = { a }; return mk_fn(head, args, 1);
}
static Expr* mk_fn2(const char* head, Expr* a, Expr* b) {
    Expr* args[2] = { a, b }; return mk_fn(head, args, 2);
}
static Expr* neg_of(Expr* g) { return mk_fn2("Times", expr_new_integer(-1), g); }

/* Evaluate `call` (taking ownership) and return the result. */
static Expr* eval_take(Expr* call) {
    Expr* r = evaluate(call);   /* evaluate does not free its input */
    expr_free(call);
    return r;
}

static bool head_name_is(const Expr* e, const char* name) {
    return e && e->type == EXPR_FUNCTION &&
           e->data.function.head->type == EXPR_SYMBOL &&
           strcmp(e->data.function.head->data.symbol.name, name) == 0;
}
static bool sym_is(const Expr* e, const char* name) {
    return e && e->type == EXPR_SYMBOL && strcmp(e->data.symbol.name, name) == 0;
}

static bool contains_symbol(const Expr* e, const Expr* x) {
    if (!e) return false;
    if (e->type == EXPR_SYMBOL)
        return x->type == EXPR_SYMBOL && e->data.symbol.name == x->data.symbol.name;
    if (e->type != EXPR_FUNCTION) return false;
    if (contains_symbol(e->data.function.head, x)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (contains_symbol(e->data.function.args[i], x)) return true;
    return false;
}
static bool contains_head(const Expr* e, const char* name) {
    if (!e) return false;
    if (head_name_is(e, name)) return true;
    if (e->type != EXPR_FUNCTION) return false;
    if (contains_head(e->data.function.head, name)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (contains_head(e->data.function.args[i], name)) return true;
    return false;
}

/* ---- infinities (bounds may be +/-Infinity) ----------------------------- */

static bool is_pos_infinity(const Expr* b) {
    if (sym_is(b, "Infinity")) return true;
    return head_name_is(b, "DirectedInfinity") && b->data.function.arg_count == 1 &&
           b->data.function.args[0]->type == EXPR_INTEGER &&
           b->data.function.args[0]->data.integer == 1;
}
static bool is_neg_infinity(const Expr* a) {
    if (head_name_is(a, "DirectedInfinity") && a->data.function.arg_count == 1 &&
        a->data.function.args[0]->type == EXPR_INTEGER &&
        a->data.function.args[0]->data.integer == -1)
        return true;
    if (head_name_is(a, "Times") && a->data.function.arg_count == 2) {
        Expr* p = a->data.function.args[0];
        Expr* q = a->data.function.args[1];
        if (p->type == EXPR_INTEGER && p->data.integer == -1 && is_pos_infinity(q))
            return true;
    }
    return false;
}

/* ---- assumption-aware predicate proving --------------------------------- */

/* True iff Refine[pred, as] proves `pred` to the symbol True.  Consumes pred. */
static bool refine_true(Expr* pred, const Expr* as) {
    Expr* call = as ? mk_fn2("Refine", pred, cp(as)) : mk_fn1("Refine", pred);
    Expr* r = eval_take(call);
    bool t = sym_is(r, "True");
    if (r) expr_free(r);
    return t;
}

/* ReplaceAll[e, x -> val], evaluated. */
static Expr* subst(const Expr* e, const Expr* x, const Expr* val) {
    Expr* rule = mk_fn2("Rule", cp(x), cp(val));
    return eval_take(mk_fn2("ReplaceAll", cp(e), rule));
}

/* Is `lo` < `c` provable (lo, c may be +/-Infinity)? */
static bool lt_bound(const Expr* lo, const Expr* c, const Expr* as) {
    if (is_neg_infinity(lo)) return true;
    if (is_pos_infinity(lo)) return false;
    if (is_pos_infinity(c))  return true;
    if (is_neg_infinity(c))  return false;
    return refine_true(mk_fn2("Less", cp(lo), cp(c)), as);
}

/* ---- breakpoint set ----------------------------------------------------- */

typedef struct { Expr** pts; int n, cap; } BPSet;

static void bp_init(BPSet* s) { s->pts = NULL; s->n = 0; s->cap = 0; }
static void bp_free(BPSet* s) {
    for (int i = 0; i < s->n; i++) expr_free(s->pts[i]);
    free(s->pts);
    s->pts = NULL; s->n = 0; s->cap = 0;
}
/* Add `p` (takes ownership); drop an exact structural duplicate. */
static void bp_add(BPSet* s, Expr* p) {
    for (int i = 0; i < s->n; i++)
        if (expr_eq(s->pts[i], p)) { expr_free(p); return; }
    if (s->n == s->cap) {
        s->cap = s->cap ? s->cap * 2 : 4;
        s->pts = realloc(s->pts, (size_t)s->cap * sizeof(Expr*));
    }
    s->pts[s->n++] = p;
}

/* Real roots of g == 0 in x, appended to s.  Complex roots and roots that
 * still mention x are skipped. */
static void solve_roots(const Expr* g, const Expr* x, BPSet* s) {
    Expr* eq  = mk_fn2("Equal", cp(g), expr_new_integer(0));
    /* Quiet: Solve is an internal probe here; its Solve::ifun note on a
     * transcendental breakpoint argument must not surface from Integrate. */
    Expr* sol = eval_take(mk_fn1("Quiet", mk_fn2("Solve", eq, cp(x))));
    if (head_name_is(sol, "List")) {
        for (size_t i = 0; i < sol->data.function.arg_count; i++) {
            Expr* si = sol->data.function.args[i];
            if (!head_name_is(si, "List")) continue;
            for (size_t j = 0; j < si->data.function.arg_count; j++) {
                Expr* rule = si->data.function.args[j];
                if (!head_name_is(rule, "Rule") ||
                    rule->data.function.arg_count != 2) continue;
                if (!expr_eq(rule->data.function.args[0], (Expr*)x)) continue;
                Expr* root = rule->data.function.args[1];
                if (contains_head(root, "Complex") || contains_symbol(root, x))
                    continue;
                bp_add(s, cp(root));
            }
        }
    }
    if (sol) expr_free(sol);
}

/* Breakpoints contributed by a condition `cond` (relational boundaries). */
static void collect_condition_boundaries(const Expr* cond, const Expr* x, BPSet* s) {
    if (!cond || cond->type != EXPR_FUNCTION) return;
    const char* h = cond->data.function.head->type == EXPR_SYMBOL
                      ? cond->data.function.head->data.symbol.name : "";
    if ((strcmp(h, "Less") == 0 || strcmp(h, "LessEqual") == 0 ||
         strcmp(h, "Greater") == 0 || strcmp(h, "GreaterEqual") == 0 ||
         strcmp(h, "Equal") == 0 || strcmp(h, "Unequal") == 0) &&
        cond->data.function.arg_count == 2) {
        Expr* diff = mk_fn2("Plus", cp(cond->data.function.args[0]),
                            neg_of(cp(cond->data.function.args[1])));
        if (contains_symbol(diff, x)) solve_roots(diff, x, s);
        expr_free(diff);
        return;
    }
    /* And / Or / Inequality / Not: recurse into the pieces. */
    for (size_t i = 0; i < cond->data.function.arg_count; i++)
        collect_condition_boundaries(cond->data.function.args[i], x, s);
}

/* Walk the integrand, collecting every branch-switch breakpoint in x. */
static void collect_breakpoints(const Expr* e, const Expr* x, BPSet* s) {
    if (!e || e->type != EXPR_FUNCTION) return;
    const char* h = e->data.function.head->type == EXPR_SYMBOL
                      ? e->data.function.head->data.symbol.name : "";
    size_t n = e->data.function.arg_count;
    if ((strcmp(h, "Abs") == 0 || strcmp(h, "Sign") == 0 ||
         strcmp(h, "HeavisideTheta") == 0 || strcmp(h, "UnitStep") == 0) &&
        n == 1 && contains_symbol(e->data.function.args[0], x)) {
        solve_roots(e->data.function.args[0], x, s);
    } else if ((strcmp(h, "Min") == 0 || strcmp(h, "Max") == 0) && n >= 2) {
        for (size_t i = 0; i < n; i++)
            for (size_t j = i + 1; j < n; j++) {
                Expr* d = mk_fn2("Plus", cp(e->data.function.args[i]),
                                 neg_of(cp(e->data.function.args[j])));
                if (contains_symbol(d, x)) solve_roots(d, x, s);
                expr_free(d);
            }
    } else if (strcmp(h, "Piecewise") == 0 && n >= 1 &&
               head_name_is(e->data.function.args[0], "List")) {
        Expr* clauses = e->data.function.args[0];
        for (size_t i = 0; i < clauses->data.function.arg_count; i++) {
            Expr* cl = clauses->data.function.args[i];
            if (head_name_is(cl, "List") && cl->data.function.arg_count == 2)
                collect_condition_boundaries(cl->data.function.args[1], x, s);
        }
    }
    /* Recurse into head + args to catch nested / product forms. */
    collect_breakpoints(e->data.function.head, x, s);
    for (size_t i = 0; i < n; i++)
        collect_breakpoints(e->data.function.args[i], x, s);
}

/* Does `e` contain any branch-defining head with an x-dependent argument? */
static bool has_target(const Expr* e, const Expr* x) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    const char* h = e->data.function.head->type == EXPR_SYMBOL
                      ? e->data.function.head->data.symbol.name : "";
    size_t n = e->data.function.arg_count;
    if ((strcmp(h, "Abs") == 0 || strcmp(h, "Sign") == 0 ||
         strcmp(h, "HeavisideTheta") == 0 || strcmp(h, "UnitStep") == 0 ||
         strcmp(h, "Min") == 0 || strcmp(h, "Max") == 0 ||
         strcmp(h, "Piecewise") == 0) && contains_symbol(e, x))
        return true;
    if (has_target(e->data.function.head, x)) return true;
    for (size_t i = 0; i < n; i++)
        if (has_target(e->data.function.args[i], x)) return true;
    return false;
}

/* ---- branch resolution at an interior test point ------------------------ */

/* Compare ei vs ej at x = m under `as`: -1 (<), 0 (==), 1 (>), 2 undecidable. */
static int cmp_at(const Expr* ei, const Expr* ej, const Expr* x,
                  const Expr* m, const Expr* as) {
    Expr* im = subst(ei, x, m);
    Expr* jm = subst(ej, x, m);
    int r;
    if (refine_true(mk_fn2("Less", cp(im), cp(jm)), as)) r = -1;
    else if (refine_true(mk_fn2("Greater", cp(im), cp(jm)), as)) r = 1;
    else if (refine_true(mk_fn2("Equal", cp(im), cp(jm)), as)) r = 0;
    else r = 2;
    expr_free(im); expr_free(jm);
    return r;
}

static Expr* resolve_at(const Expr* e, const Expr* x, const Expr* m,
                        const Expr* as, bool* ok);

/* Abs/Sign/HeavisideTheta/UnitStep of `g` (already resolved; consumed). */
static Expr* resolve_unary(const char* h, Expr* g, const Expr* x,
                           const Expr* m, const Expr* as, bool* ok) {
    Expr* gm = subst(g, x, m);
    bool pos = refine_true(mk_fn2("Greater", cp(gm), expr_new_integer(0)), as);
    bool neg = pos ? false
                   : refine_true(mk_fn2("Less", cp(gm), expr_new_integer(0)), as);
    expr_free(gm);
    if (!pos && !neg) { expr_free(g); *ok = false; return NULL; }
    if (strcmp(h, "Abs") == 0) return pos ? g : neg_of(g);
    expr_free(g);
    if (strcmp(h, "Sign") == 0) return expr_new_integer(pos ? 1 : -1);
    /* HeavisideTheta / UnitStep */
    return expr_new_integer(pos ? 1 : 0);
}

/* Piecewise[clauses, default?] active branch at x = m. */
static Expr* resolve_piecewise(const Expr* e, const Expr* x, const Expr* m,
                               const Expr* as, bool* ok) {
    Expr* clauses = e->data.function.args[0];
    for (size_t i = 0; i < clauses->data.function.arg_count; i++) {
        Expr* cl = clauses->data.function.args[i];
        if (!head_name_is(cl, "List") || cl->data.function.arg_count != 2) {
            *ok = false; return NULL;
        }
        Expr* cond_m = subst(cl->data.function.args[1], x, m);
        bool t = refine_true(cp(cond_m), as);
        bool f = t ? false : refine_true(mk_fn1("Not", cp(cond_m)), as);
        expr_free(cond_m);
        if (t) return resolve_at(cl->data.function.args[0], x, m, as, ok);
        if (!f) { *ok = false; return NULL; }   /* condition undecidable here */
    }
    if (e->data.function.arg_count >= 2)
        return resolve_at(e->data.function.args[1], x, m, as, ok);
    return expr_new_integer(0);   /* Mathematica's Piecewise default */
}

/* Rewrite `e` to the single smooth branch active on the sub-interval whose
 * interior test point is `m`.  *ok is cleared when a branch cannot be proved. */
static Expr* resolve_at(const Expr* e, const Expr* x, const Expr* m,
                        const Expr* as, bool* ok) {
    if (!*ok) return NULL;
    if (!e || e->type != EXPR_FUNCTION) return cp(e);
    const char* h = e->data.function.head->type == EXPR_SYMBOL
                      ? e->data.function.head->data.symbol.name : "";
    size_t n = e->data.function.arg_count;

    if ((strcmp(h, "Abs") == 0 || strcmp(h, "Sign") == 0 ||
         strcmp(h, "HeavisideTheta") == 0 || strcmp(h, "UnitStep") == 0) &&
        n == 1 && contains_symbol(e->data.function.args[0], x)) {
        Expr* g = resolve_at(e->data.function.args[0], x, m, as, ok);
        if (!*ok) { if (g) expr_free(g); return NULL; }
        return resolve_unary(h, g, x, m, as, ok);
    }
    if ((strcmp(h, "Min") == 0 || strcmp(h, "Max") == 0) && n >= 1 &&
        contains_symbol(e, x)) {
        bool want_max = (strcmp(h, "Max") == 0);
        size_t best = 0;
        for (size_t j = 1; j < n; j++) {
            int c = cmp_at(e->data.function.args[j], e->data.function.args[best],
                           x, m, as);
            if (c == 2) { *ok = false; return NULL; }
            if (( want_max && c > 0) || (!want_max && c < 0)) best = j;
        }
        return resolve_at(e->data.function.args[best], x, m, as, ok);
    }
    if (strcmp(h, "Piecewise") == 0 && n >= 1 &&
        head_name_is(e->data.function.args[0], "List")) {
        return resolve_piecewise(e, x, m, as, ok);
    }

    /* Generic node: resolve head + args and rebuild. */
    Expr* rhead = resolve_at(e->data.function.head, x, m, as, ok);
    Expr** rargs = malloc(n * sizeof(Expr*));
    size_t built = 0;
    for (size_t i = 0; i < n && *ok; i++)
        rargs[built++] = resolve_at(e->data.function.args[i], x, m, as, ok);
    if (!*ok) {
        if (rhead) expr_free(rhead);
        for (size_t i = 0; i < built; i++) if (rargs[i]) expr_free(rargs[i]);
        free(rargs);
        return NULL;
    }
    Expr* out = expr_new_function(rhead, rargs, n);
    free(rargs);
    return out;
}

/* ---- the integrator ----------------------------------------------------- */

/* Integrate[piece, {x, p, q}, Assumptions -> as] -> value, or NULL if it does
 * not close to an x-free, Integrate-free result. */
static Expr* integrate_piece(const Expr* piece, const Expr* x,
                             const Expr* p, const Expr* q, const Expr* as) {
    Expr* spec = mk_fn("List", (Expr*[]){ cp(x), cp(p), cp(q) }, 3);
    Expr* call;
    if (as) {
        Expr* opt = mk_fn2("Rule", expr_new_symbol("Assumptions"), cp(as));
        call = mk_fn("Integrate", (Expr*[]){ cp(piece), spec, opt }, 3);
    } else {
        call = mk_fn("Integrate", (Expr*[]){ cp(piece), spec }, 2);
    }
    Expr* v = eval_take(call);
    if (!v) return NULL;
    if (contains_head(v, "Integrate") || contains_symbol(v, x)) {
        expr_free(v); return NULL;
    }
    return v;
}

Expr* integrate_piecewise_try(Expr* f, Expr* x, Expr* a, Expr* b,
                              Expr* assumptions) {
    if (!f || !x || !a || !b || x->type != EXPR_SYMBOL) return NULL;
    if (!has_target(f, x)) return NULL;

    /* 1. Collect breakpoints, keep the strictly-interior ones. */
    BPSet all; bp_init(&all);
    collect_breakpoints(f, x, &all);
    BPSet inside; bp_init(&inside);
    for (int i = 0; i < all.n; i++) {
        if (lt_bound(a, all.pts[i], assumptions) &&
            lt_bound(all.pts[i], b, assumptions))
            bp_add(&inside, cp(all.pts[i]));
    }
    bp_free(&all);

    /* 2. Sort interior breakpoints (insertion sort, proved order). */
    for (int i = 1; i < inside.n; i++) {
        Expr* key = inside.pts[i];
        int j = i - 1;
        while (j >= 0) {
            if (refine_true(mk_fn2("Less", cp(key), cp(inside.pts[j])),
                            assumptions)) {
                inside.pts[j + 1] = inside.pts[j]; j--;
            } else if (refine_true(mk_fn2("Less", cp(inside.pts[j]), cp(key)),
                                   assumptions) ||
                       refine_true(mk_fn2("Equal", cp(inside.pts[j]), cp(key)),
                                   assumptions)) {
                break;
            } else {
                bp_free(&inside); return NULL;   /* order undecidable */
            }
        }
        inside.pts[j + 1] = key;
    }

    /* 3. Endpoints [a, c1, ..., cn, b]; sum the resolved pieces. */
    int npts = inside.n + 2;
    Expr** ep = malloc((size_t)npts * sizeof(Expr*));
    ep[0] = cp(a);
    for (int i = 0; i < inside.n; i++) ep[i + 1] = cp(inside.pts[i]);
    ep[npts - 1] = cp(b);
    bp_free(&inside);

    Expr** terms = malloc((size_t)(npts - 1) * sizeof(Expr*));
    int nterms = 0;
    bool ok = true;
    for (int i = 0; i < npts - 1 && ok; i++) {
        Expr* p = ep[i];
        Expr* q = ep[i + 1];
        /* interior test point m of [p, q] (p, q may be +/-Infinity) */
        Expr* m;
        bool pinf = is_pos_infinity(p) || is_neg_infinity(p);
        bool qinf = is_pos_infinity(q) || is_neg_infinity(q);
        if (!pinf && !qinf)
            m = eval_take(mk_fn2("Divide", mk_fn2("Plus", cp(p), cp(q)),
                                 expr_new_integer(2)));           /* (p+q)/2 */
        else if (!pinf && qinf)
            m = eval_take(mk_fn2("Plus", cp(p), expr_new_integer(1)));   /* p+1 */
        else if (pinf && !qinf)
            m = eval_take(mk_fn2("Plus", cp(q), expr_new_integer(-1)));  /* q-1 */
        else
            m = expr_new_integer(0);   /* whole line: test at 0 */
        Expr* piece = resolve_at(f, x, m, assumptions, &ok);
        expr_free(m);
        if (!ok || !piece) { if (piece) expr_free(piece); ok = false; break; }
        Expr* val = integrate_piece(piece, x, p, q, assumptions);
        expr_free(piece);
        if (!val) { ok = false; break; }
        terms[nterms++] = val;
    }
    for (int i = 0; i < npts; i++) expr_free(ep[i]);
    free(ep);

    if (!ok) {
        for (int i = 0; i < nterms; i++) expr_free(terms[i]);
        free(terms);
        return NULL;
    }

    Expr* sum;
    if (nterms == 1) { sum = terms[0]; free(terms); }
    else {
        sum = expr_new_function(expr_new_symbol(SYM_Plus), terms, (size_t)nterms);
        free(terms);
    }
    /* Clean up under the assumptions (a>0 etc. tidy the closed form). */
    Expr* simp;
    if (assumptions)
        simp = eval_take(mk_fn2("Simplify", sum, cp(assumptions)));
    else
        simp = eval_take(mk_fn1("Simplify", sum));
    return simp;
}

/* ---- builtin ------------------------------------------------------------ */

Expr* builtin_integrate_piecewise(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 2) return NULL;
    Expr* f    = res->data.function.args[0];
    Expr* spec = res->data.function.args[1];
    if (!head_name_is(spec, "List") || spec->data.function.arg_count != 3)
        return NULL;
    Expr* x = spec->data.function.args[0];
    Expr* a = spec->data.function.args[1];
    Expr* b = spec->data.function.args[2];
    if (x->type != EXPR_SYMBOL) return NULL;

    Expr* assumptions = NULL;
    for (size_t t = 2; t < argc; t++) {
        Expr* opt = res->data.function.args[t];
        if (opt->type == EXPR_FUNCTION && opt->data.function.arg_count == 2 &&
            opt->data.function.head->type == EXPR_SYMBOL &&
            (opt->data.function.head->data.symbol.name == SYM_Rule ||
             opt->data.function.head->data.symbol.name == SYM_RuleDelayed)) {
            Expr* lhs = opt->data.function.args[0];
            if (lhs->type == EXPR_SYMBOL &&
                strcmp(lhs->data.symbol.name, "Assumptions") == 0) {
                assumptions = opt->data.function.args[1];
                continue;
            }
        }
        return NULL;
    }
    return integrate_piecewise_try(f, x, a, b, assumptions);
}

void integrate_piecewise_init(void) {
    symtab_add_builtin("Integrate`Piecewise", builtin_integrate_piecewise);
    symtab_get_def("Integrate`Piecewise")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("Integrate`Piecewise",
        "Integrate`Piecewise[f, {x, a, b}] evaluates a definite integral whose "
        "integrand is built from Abs, Sign, HeavisideTheta, UnitStep, Min, Max, "
        "or Piecewise.  It splits [a, b] at the interior breakpoints where a "
        "branch switches (decided assumption-aware with Refine), reduces each "
        "sub-interval to its single active branch, integrates that, and sums.  "
        "Returns unevaluated when a breakpoint position or branch sign cannot be "
        "proved, or a piece does not close.");
}
