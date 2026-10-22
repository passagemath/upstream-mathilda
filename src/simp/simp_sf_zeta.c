/*
 * simp_sf_zeta.c -- special-function Simplify descriptors: PolyGamma,
 * HarmonicNumber, the zeta family (Zeta, HurwitzZeta, LerchPhi) and PolyLog.
 *
 * Derivations in docs/design/simp_special_functions.md (sections PolyGamma,
 * Zeta, PolyLog); every identity is checked numerically by
 * tests/test_simp_sf_zeta.c at real and complex points.
 *
 * Recurrences (rank 1, affine; up-matrix row given):
 *
 *   PolyGamma[n, z+1] = PolyGamma[n, z] + (-1)^n n! z^(-n-1)  [[1, (-1)^n n! z^(-n-1)]]
 *       n-th derivative of psi(z+1) = psi(z) + 1/z (log-derivative of the
 *       Gamma recurrence).  Integer n >= 0 only: at a non-integer order the
 *       fractional-derivative PolyGamma does not obey it (checked).
 *   HarmonicNumber[n+1]    = HarmonicNumber[n] + 1/(n+1)         [[1, 1/(n+1)]]
 *   HarmonicNumber[n+1, r] = HarmonicNumber[n, r] + (n+1)^-r     [[1, (n+1)^-r]]
 *       from H(n) = psi(n+1) + EulerGamma, resp. the defining partial sums
 *       continued analytically.
 *   HurwitzZeta[s, a+1] = HurwitzZeta[s, a] - a^-s               [[1, -a^-s]]
 *       Sum_{k>=0} (k+a)^-s, first term split off; analytic continuation.
 *   Zeta[s, a+1]        = Zeta[s, a] - (a^2)^(-s/2)              [[1, -(a^2)^(-s/2)]]
 *       Zeta[s, a] sums ((k+a)^2)^(-s/2) (the convention for Re a < 0;
 *       = (k+a)^-s when Re a > 0), so the split-off term is (a^2)^(-s/2).
 *   LerchPhi[z, s, a+1] = (LerchPhi[z, s, a] - a^-s)/z           [[1/z, -a^-s/z]]
 *       Phi = Sum z^k (k+a)^-s = a^-s + z Phi(z, s, a+1).
 *
 * Reflection / multiplication (canon hooks, applied when the partner atom is
 * present and the atom is the less canonical side, as for Gamma):
 *
 *   PolyGamma[n, 1 - w] = (-1)^n (PolyGamma[n, w] + Pi D^n_w Cot[Pi w])
 *       n-th derivative of psi(1-w) - psi(w) = Pi Cot[Pi w] (log-derivative
 *       of Euler's reflection).
 *   PolyGamma[0, m z] = Log[m] + (1/m) Sum_j PolyGamma[0, z + j/m]
 *   PolyGamma[n, m z] = m^(-n-1) Sum_j PolyGamma[n, z + j/m]        (n >= 1)
 *       log-derivative (then n derivatives) of Gauss multiplication.
 *   Zeta[1-w] = 2^(1-w) Pi^-w Cos[Pi w/2] Gamma[w] Zeta[w]  (Riemann, at s = 1-w)
 *   HarmonicNumber[n] = PolyGamma[0, n+1] + EulerGamma  (when PolyGamma[0, _]
 *       is present to meet it).
 *
 *   PolyLog (order n, same in the pair):
 *     reflection  Li2[1-u] = Pi^2/6 - Log[u] Log[1-u] - Li2[u]
 *     Landen      Li2[u/(u-1)] = -Log[1-u]^2/2 - Li2[u]       (u < 1 proved)
 *     duplication Li_n[u^2] = 2^(n-1) (Li_n[u] + Li_n[-u])
 *     inversion   Li2[1/u] = -Pi^2/6 - Log[-u]^2/2 - Li2[u]   (u < 0 proved)
 *   Reflection and duplication hold for every complex u, cuts included
 *   (checked in all four quadrants and on (1, Inf) and (-Inf, 0)).  Landen's
 *   two sides disagree on (1, Inf) and inversion's on (0, Inf) -- a real
 *   half-line, not a removable set -- so those need u < 1 and u < 0 proved.
 */

#include <stdlib.h>
#include <string.h>
#include "expr.h"
#include "eval.h"
#include "sym_names.h"
#include "simp.h"
#include "simp_internal.h"
#include "simp_sf_engine.h"

#define PG_MULT_MAX 16

/* Literal machine integer >= 0? */
static bool nonneg_int(const Expr* e, int64_t* v) {
    if (e->type != EXPR_INTEGER || e->data.integer < 0) return false;
    *v = e->data.integer;
    return true;
}

static Expr* factorial_i(int64_t n) { return sf_mk1(SYM_Factorial, sf_int(n)); }

/* ---- recurrences ---------------------------------------------------------- */

static bool up_polygamma(Expr* const* a, int d, Expr** M) {
    int64_t n;
    if (d != 0 || !nonneg_int(a[0], &n)) return false;
    M[0] = sf_int(1);
    M[1] = sf_times(sf_times(sf_int((n % 2) ? -1 : 1), factorial_i(n)),
                    sf_pow(sf_cp(a[1]), sf_int(-n - 1)));
    return true;
}

static bool up_harmonic(Expr* const* a, int d, Expr** M) {
    if (d != 0) return false;
    M[0] = sf_int(1);
    M[1] = sf_pow(sf_plus(sf_cp(a[0]), sf_int(1)), sf_int(-1));
    return true;
}

static bool up_harmonic2(Expr* const* a, int d, Expr** M) {
    if (d != 0) return false;
    M[0] = sf_int(1);
    M[1] = sf_pow(sf_plus(sf_cp(a[0]), sf_int(1)), sf_times(sf_int(-1), sf_cp(a[1])));
    return true;
}

static bool up_hurwitz(Expr* const* a, int d, Expr** M) {
    if (d != 0) return false;
    M[0] = sf_int(1);
    M[1] = sf_times(sf_int(-1), sf_pow(sf_cp(a[1]), sf_times(sf_int(-1), sf_cp(a[0]))));
    return true;
}

static bool up_zeta2(Expr* const* a, int d, Expr** M) {
    if (d != 0) return false;
    M[0] = sf_int(1);
    M[1] = sf_times(sf_int(-1), sf_pow(sf_pow(sf_cp(a[1]), sf_int(2)),
                                       sf_times(sf_rat(-1, 2), sf_cp(a[0]))));
    return true;
}

static bool up_lerch(Expr* const* a, int d, Expr** M) {
    if (d != 0) return false;
    Expr* zinv = sf_pow(sf_cp(a[0]), sf_int(-1));
    M[0] = sf_cp(zinv);
    M[1] = sf_times(sf_int(-1), sf_times(zinv, sf_pow(sf_cp(a[2]), sf_times(sf_int(-1), sf_cp(a[1])))));
    return true;
}

/* ---- shared class helpers -------------------------------------------------- */

static bool class_eq(const Expr* a, const Expr* b) {
    Expr *sa, *sb; int64_t na, da, nb, db;
    if (!sf_split(a, &sa, &na, &da)) return false;
    if (!sf_split(b, &sb, &nb, &db)) { if (sa) expr_free(sa); return false; }
    bool eq = da == db && ((na - nb) % da) == 0 &&
              ((!sa && !sb) || (sa && sb && expr_eq(sa, sb)));
    if (sa) expr_free(sa);
    if (sb) expr_free(sb);
    return eq;
}

/* Is there a call head[..., arg_i, ...] in e (other than `skip`) whose
 * argument `ai` lies in the class of mu and whose argument `oi` (if >= 0)
 * equals `other`? */
static bool has_atom_in_class(const Expr* e, const char* head, size_t arity, size_t ai,
                              const Expr* mu, int oi, const Expr* other, const Expr* skip) {
    if (e->type != EXPR_FUNCTION) return false;
    if (sf_is_call(e, head, arity) && !(skip && expr_eq(e, skip)) &&
        class_eq(e->data.function.args[ai], mu) &&
        (oi < 0 || expr_eq(e->data.function.args[oi], other)))
        return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (has_atom_in_class(e->data.function.args[i], head, arity, ai, mu, oi, other, skip)) return true;
    return false;
}

/* Is there a call head[.., arg ai with symbolic part s, ..] (other arg oi equal)? */
static bool has_symbolic_part(const Expr* e, const char* head, size_t arity, size_t ai,
                              const Expr* s, int oi, const Expr* other) {
    if (e->type != EXPR_FUNCTION) return false;
    if (sf_is_call(e, head, arity) && (oi < 0 || expr_eq(e->data.function.args[oi], other))) {
        Expr* sp = NULL; int64_t num, den;
        bool hit = sf_split(e->data.function.args[ai], &sp, &num, &den) && sp && expr_eq(sp, s);
        if (sp) expr_free(sp);
        if (hit) return true;
    }
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (has_symbolic_part(e->data.function.args[i], head, arity, ai, s, oi, other)) return true;
    return false;
}

static int64_t floor_div(int64_t a, int64_t b) {
    int64_t q = a / b;
    return (a % b != 0 && ((a < 0) != (b < 0))) ? q - 1 : q;
}

/* Is mu the less canonical side of the reflection pair (mu, 1 - mu)? */
static bool less_canonical(const Expr* mu) {
    Expr* s = NULL; int64_t num, den;
    if (!sf_split(mu, &s, &num, &den)) return false;
    int64_t k = floor_div(num, den), r = num - k * den;
    bool lc = s ? sf_reads_negative(s) : (2 * r > den);
    if (s) expr_free(s);
    return lc;
}

static Expr* expand1(Expr* e) { return sf_mk1("Expand", e); }

/* ---- PolyGamma ---------------------------------------------------------- */

/* D[Cot[Pi u], {u, n}] /. u -> w (owned), via a private dummy symbol. */
static Expr* cot_derivative(int64_t n, const Expr* w) {
    Expr* u = sf_sym("SfEngine$u");
    Expr* f = sf_call(SYM_Cot, (Expr*[]){ sf_times(sf_sym(SYM_Pi), sf_cp(u)) }, 1);
    Expr* dn = eval_and_free(sf_call("D", (Expr*[]){ f, sf_call(SYM_List, (Expr*[]){ sf_cp(u), sf_int(n) }, 2) }, 2));
    Expr* rule = sf_call(SYM_Rule, (Expr*[]){ u, sf_cp(w) }, 2);
    return eval_and_free(sf_call(SYM_ReplaceAll, (Expr*[]){ dn, rule }, 2));
}

static Expr* polygamma_reflect(const Expr* atom, const Expr* root) {
    int64_t n;
    const Expr* ord = atom->data.function.args[0];
    const Expr* mu = atom->data.function.args[1];
    if (!nonneg_int(ord, &n) || !less_canonical(mu)) return NULL;
    Expr* w = expand1(sf_minus(sf_int(1), sf_cp(mu)));
    if (!has_atom_in_class(root, SYM_PolyGamma, 2, 1, w, 0, ord, atom)) { expr_free(w); return NULL; }
    Expr* pg = sf_call(SYM_PolyGamma, (Expr*[]){ sf_int(n), sf_cp(w) }, 2);
    Expr* body = sf_plus(pg, sf_times(sf_sym(SYM_Pi), cot_derivative(n, w)));
    expr_free(w);
    return (n % 2) ? sf_times(sf_int(-1), body) : body;
}

static int64_t gcd64(int64_t a, int64_t b) {
    if (a < 0) a = -a;
    if (b < 0) b = -b;
    while (b) { int64_t t = a % b; a = b; b = t; }
    return a;
}

static int64_t integer_content(const Expr* s) {
    const Expr* const* terms = &s;
    size_t nt = 1;
    if (s->type == EXPR_FUNCTION && s->data.function.head->type == EXPR_SYMBOL &&
        s->data.function.head->data.symbol.name == SYM_Plus) {
        terms = (const Expr* const*)s->data.function.args;
        nt = s->data.function.arg_count;
    }
    int64_t g = 0;
    for (size_t i = 0; i < nt; i++) {
        const Expr* t = terms[i];
        int64_t c = 1;
        if (t->type == EXPR_FUNCTION && t->data.function.head->type == EXPR_SYMBOL &&
            t->data.function.head->data.symbol.name == SYM_Times && t->data.function.arg_count >= 2) {
            const Expr* c0 = t->data.function.args[0];
            if (c0->type == EXPR_INTEGER) c = c0->data.integer;
            else if (sf_is_call(c0, SYM_Rational, 2)) return 1;
        }
        g = gcd64(g, c);
    }
    return g ? g : 1;
}

static Expr* polygamma_multiply(const Expr* atom, const Expr* root) {
    int64_t n;
    const Expr* ord = atom->data.function.args[0];
    if (!nonneg_int(ord, &n)) return NULL;
    Expr* mu = expand1(sf_cp(atom->data.function.args[1]));
    Expr* s = NULL; int64_t num, den;
    bool ok = sf_split(mu, &s, &num, &den);
    expr_free(mu);
    if (!ok || !s) { if (s) expr_free(s); return NULL; }
    int64_t m = integer_content(s);
    if (m < 2 || m > PG_MULT_MAX) { expr_free(s); return NULL; }
    Expr* s0 = expand1(sf_times(sf_rat(1, m), s));
    if (!has_symbolic_part(root, SYM_PolyGamma, 2, 1, s0, 0, ord)) { expr_free(s0); return NULL; }
    Expr* z = sf_plus(s0, sf_rat(num, den * m));
    Expr* sum = sf_int(0);
    for (int64_t j = 0; j < m; j++)
        sum = sf_call(SYM_Plus, (Expr*[]){ sum, sf_call(SYM_PolyGamma,
                     (Expr*[]){ sf_int(n), sf_plus(sf_cp(z), sf_rat(j, m)) }, 2) }, 2);
    expr_free(z);
    if (n == 0)
        return sf_plus(sf_mk1(SYM_Log, sf_int(m)), sf_times(sf_rat(1, m), sum));
    return sf_times(sf_pow(sf_int(m), sf_int(-n - 1)), sum);
}

static Expr* canon_polygamma(const Expr* atom, const Expr* root, const AssumeCtx* ctx) {
    (void)ctx;
    Expr* r = polygamma_multiply(atom, root);
    if (r) return r;
    return polygamma_reflect(atom, root);
}

/* PolyGamma[z] (one argument) is PolyGamma[0, z]; the evaluator already
 * rewrites it, so no descriptor is needed for arity 1. */

static Expr* canon_harmonic(const Expr* atom, const Expr* root, const AssumeCtx* ctx) {
    (void)ctx;
    if (!sf_contains_head(root, SYM_PolyGamma)) return NULL;
    return sf_plus(sf_call(SYM_PolyGamma, (Expr*[]){ sf_int(0),
                   sf_plus(sf_cp(atom->data.function.args[0]), sf_int(1)) }, 2),
                   sf_sym(SYM_EulerGamma));
}

/* ---- Riemann zeta functional equation ------------------------------------ */

/* Cos[Pi w/2] with an integer offset k in w = k + t folded by quarter
 * turns: k mod 4 = 0, 1, 2, 3 -> Cos[Pi t/2], -Sin[Pi t/2], -Cos[Pi t/2],
 * Sin[Pi t/2] (the evaluator leaves Cos[Pi/2 + Pi s/2] alone). */
static Expr* cos_half_pi(const Expr* w) {
    Expr* t = NULL; int64_t num, den;
    if (sf_split(w, &t, &num, &den) && t && den == 1) {
        int64_t k = ((num % 4) + 4) % 4;
        Expr* arg = sf_times(sf_times(sf_sym(SYM_Pi), t), sf_rat(1, 2));
        Expr* f = sf_mk1((k % 2) ? SYM_Sin : SYM_Cos, arg);
        return (k == 1 || k == 2) ? sf_times(sf_int(-1), f) : f;
    }
    if (t) expr_free(t);
    return sf_mk1(SYM_Cos, sf_times(sf_times(sf_sym(SYM_Pi), sf_cp(w)), sf_rat(1, 2)));
}

static Expr* canon_zeta1(const Expr* atom, const Expr* root, const AssumeCtx* ctx) {
    (void)ctx;
    const Expr* mu = atom->data.function.args[0];
    if (!less_canonical(mu)) return NULL;
    Expr* w = expand1(sf_minus(sf_int(1), sf_cp(mu)));     /* partner argument */
    bool partner = false;
    /* Zeta[w] present (exactly; the functional equation is not a lattice) */
    Expr* zw = eval_and_free(sf_call(SYM_Zeta, (Expr*[]){ sf_cp(w) }, 1));
    partner = sf_is_call(zw, SYM_Zeta, 1) && (sf_contains_head(root, SYM_Zeta)) &&
              has_atom_in_class(root, SYM_Zeta, 1, 0, w, -1, NULL, atom);
    if (!partner) { expr_free(zw); expr_free(w); return NULL; }
    /* Zeta[1 - w] = 2^(1-w) Pi^-w Cos[Pi w/2] Gamma[w] Zeta[w]: the Riemann
     * form at s = 1 - w, with Sin[Pi (1 - w)/2] = Cos[Pi w/2], written in the
     * partner's own variable */
    Expr* r = sf_times(sf_times(sf_pow(sf_int(2), sf_minus(sf_int(1), sf_cp(w))),
                                sf_pow(sf_sym(SYM_Pi), sf_times(sf_int(-1), sf_cp(w)))),
                       sf_times(cos_half_pi(w),
                                sf_times(sf_call(SYM_Gamma, (Expr*[]){ sf_cp(w) }, 1), zw)));
    expr_free(w);
    return r;
}

/* ---- PolyLog --------------------------------------------------------------- */

/* Collect the PolyLog[n, w] atoms of e with order n (borrowed, excluding
 * any equal to `skip`), deduplicated. */
typedef struct { const Expr* a[64]; size_t n; } AtomList;

static void collect_polylogs(const Expr* e, const Expr* n, const Expr* skip, AtomList* out) {
    if (e->type != EXPR_FUNCTION) return;
    if (sf_is_call(e, SYM_PolyLog, 2) && expr_eq(e->data.function.args[0], n) &&
        !(skip && expr_eq(e, skip))) {
        bool dup = false;
        for (size_t i = 0; i < out->n && !dup; i++) dup = expr_eq(out->a[i], e);
        if (!dup && out->n < sizeof out->a / sizeof out->a[0]) out->a[out->n++] = e;
    }
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        collect_polylogs(e->data.function.args[i], n, skip, out);
}

/* Is a - b identically zero as a rational function (Together)? Borrowed. */
static bool same_value(const Expr* a, const Expr* b) {
    Expr* d = eval_and_free(sf_call("Together", (Expr*[]){ sf_minus(sf_cp(a), sf_cp(b)) }, 1));
    bool z = d->type == EXPR_INTEGER && d->data.integer == 0;
    expr_free(d);
    return z;
}

/* Strict total order on arguments picking the side of a pair to rewrite:
 * the larger leaf count, ties by canonical expression order. */
static bool rewrite_this_side(const Expr* mine, const Expr* partner) {
    size_t a = score_with_func(mine, NULL), b = score_with_func(partner, NULL);
    if (a != b) return a > b;
    return expr_compare(mine, partner) > 0;
}

static bool is_two(const Expr* e) { return e->type == EXPR_INTEGER && e->data.integer == 2; }

static Expr* pi2_6(void) { return sf_times(sf_rat(1, 6), sf_pow(sf_sym(SYM_Pi), sf_int(2))); }

static Expr* li(const Expr* n, Expr* u) { return sf_call(SYM_PolyLog, (Expr*[]){ sf_cp(n), u }, 2); }

/* w < c under the assumptions (Reduce/CAD entailment: 2 x < 0 and
 * x - 1 < 0 from x < 0). */
static bool proved_less(const AssumeCtx* ctx, const Expr* w, int64_t c) {
    Expr* pred = sf_call(SYM_Less, (Expr*[]){ sf_cp(w), sf_int(c) }, 2);
    int r = assume_reduce_entails(ctx, pred);
    expr_free(pred);
    return r == 1;
}

/* w < 0: the cheap prover first, then the entailment. */
static bool proved_negative(const AssumeCtx* ctx, const Expr* w) {
    return prov_neg(ctx, w) || proved_less(ctx, w, 0);
}

/* PolyLog[n, v]: rewrite onto a present partner PolyLog[n, w], matched by
 * VALUE (Together), and written with the partner's own argument w. */
static Expr* canon_polylog(const Expr* atom, const Expr* root, const AssumeCtx* ctx) {
    const Expr* n = atom->data.function.args[0];
    const Expr* v = atom->data.function.args[1];
    /* one spelling per argument value: PolyLog[2, -(1/10 + z)] and the
     * engine's PolyLog[2, -1/10 - z] must be one kernel */
    if (v->type == EXPR_FUNCTION) {
        Expr* xv = sf_mk1("Expand", sf_cp(v));
        if (!expr_eq(xv, v)) return li(n, xv);
        expr_free(xv);
    }
    AtomList pl = { { NULL }, 0 };
    collect_polylogs(root, n, atom, &pl);
    if (!pl.n) return NULL;
    for (size_t i = 0; i < pl.n; i++) {
        const Expr* w = pl.a[i]->data.function.args[1];
        /* duplication: v = w^2  ->  2^(n-1) (Li_n[w] + Li_n[-w]) */
        Expr* w2 = sf_pow(sf_cp(w), sf_int(2));
        bool dup = same_value(v, w2);
        expr_free(w2);
        if (dup)
            return sf_times(sf_pow(sf_int(2), sf_minus(sf_cp(n), sf_int(1))),
                            sf_plus(li(n, sf_cp(w)), li(n, sf_mk1("Expand", sf_times(sf_int(-1), sf_cp(w))))));
        if (!is_two(n) || !rewrite_this_side(v, w)) continue;
        /* reflection: v = 1 - w */
        Expr* t = sf_minus(sf_int(1), sf_cp(w));
        bool refl = same_value(v, t);
        expr_free(t);
        if (refl)
            return sf_minus(sf_minus(pi2_6(), sf_times(sf_mk1(SYM_Log, sf_cp(w)),
                                                       sf_mk1(SYM_Log, sf_minus(sf_int(1), sf_cp(w))))),
                            li(n, sf_cp(w)));
        /* Landen: v = w/(w - 1), only for w < 1 proved: the two sides
         * disagree on the real half-line (1, Inf) (checked: 5.575 i apart at
         * w = 17/10), and a Simplify identity must hold for every value */
        t = sf_div(sf_cp(w), sf_minus(sf_cp(w), sf_int(1)));
        bool landen = ctx && same_value(v, t) && proved_less(ctx, w, 1);
        expr_free(t);
        if (landen)
            return sf_minus(sf_times(sf_rat(-1, 2),
                                     sf_pow(sf_mk1(SYM_Log, sf_minus(sf_int(1), sf_cp(w))), sf_int(2))),
                            li(n, sf_cp(w)));
        /* inversion: v = 1/w, w < 0 proved */
        t = sf_pow(sf_cp(w), sf_int(-1));
        bool inv = ctx && same_value(v, t) && proved_negative(ctx, w);
        expr_free(t);
        if (inv)
            return sf_minus(sf_minus(sf_times(sf_int(-1), pi2_6()),
                                     sf_times(sf_rat(1, 2), sf_pow(sf_mk1(SYM_Log,
                                              sf_times(sf_int(-1), sf_cp(w))), sf_int(2)))),
                            li(n, sf_cp(w)));
    }
    return NULL;
}

/* ---- module table -------------------------------------------------------- */

static const SfDesc ZETA_DESCS[] = {
    { &SYM_PolyGamma,      2, 1, 1, { 1 }, up_polygamma, NULL, canon_polygamma, NULL, NULL },
    { &SYM_HarmonicNumber, 1, 1, 1, { 0 }, up_harmonic,  NULL, canon_harmonic,  NULL, NULL },
    { &SYM_HarmonicNumber, 2, 1, 1, { 0 }, up_harmonic2, NULL, NULL,            NULL, NULL },
    { &SYM_HurwitzZeta,    2, 1, 1, { 1 }, up_hurwitz,   NULL, NULL,            NULL, NULL },
    { &SYM_Zeta,           2, 1, 1, { 1 }, up_zeta2,     NULL, NULL,            NULL, NULL },
    { &SYM_Zeta,           1, 0, 1, { 0 }, NULL,         NULL, canon_zeta1,     NULL, NULL },
    { &SYM_LerchPhi,       3, 1, 1, { 2 }, up_lerch,     NULL, NULL,            NULL, NULL },
    { &SYM_PolyLog,        2, 0, 1, { 0 }, NULL,         NULL, canon_polylog,   NULL, NULL },
};

const SfModule SF_MODULE_ZETA = {
    ZETA_DESCS, sizeof ZETA_DESCS / sizeof ZETA_DESCS[0], NULL, 0, NULL, 0
};
