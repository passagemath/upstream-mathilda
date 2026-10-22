/*
 * simp_sf_orthopoly.c -- special-function Simplify descriptors: Legendre P
 * and Q (and associated P), Chebyshev T and U, Hermite H, Laguerre L (and
 * generalized L^a).
 *
 * Derivations in docs/design/simp_special_functions.md (section 8); every
 * recurrence below was checked numerically at NON-INTEGER and complex
 * degree (the functions are continued in the degree, and the recurrences
 * hold there), and each case again by tests/test_simp_sf_orthopoly.c.
 *
 * Three-term recurrences, rank 2, basis v(n) = (C[n], C[n+1], 1); the up
 * row gives C[n+2] = A C[n] + B C[n+1]:
 *
 *   LegendreP, LegendreQ   (n+1) C[n+1] = (2n+1) x C[n] - n C[n-1]
 *       A = -(n+1)/(n+2),        B = (2n+3) x/(n+2)
 *   LegendreP[n, m, x]     (n-m+1) P[n+1] = (2n+1) x P[n] - (n+m) P[n-1]
 *       A = -(n+m+1)/(n-m+2),    B = (2n+3) x/(n-m+2)
 *   ChebyshevT, ChebyshevU C[n+1] = 2x C[n] - C[n-1]
 *       (cos((n+1)t) + cos((n-1)t) = 2 cos t cos(n t), and likewise for
 *        sin((n+1)t)/sin t)              A = -1, B = 2x
 *   HermiteH               H[n+1] = 2x H[n] - 2n H[n-1]
 *       (from the generating function e^(2xt - t^2))   A = -2(n+1), B = 2x
 *   LaguerreL[n, x]        (n+1) L[n+1] = (2n+1-x) L[n] - n L[n-1]
 *       A = -(n+1)/(n+2),        B = (2n+3-x)/(n+2)
 *   LaguerreL[n, a, x]     (n+1) L[n+1] = (2n+1+a-x) L[n] - (n+a) L[n-1]
 *       A = -(n+1+a)/(n+2),      B = (2n+3+a-x)/(n+2)
 *     and in a:  x L^(a+1)_n = (n+a+1) L^a_n - (n+1) L^a_{n+1}, which with the
 *     n-recurrence gives the a-direction matrix
 *       [[(n+a+1)/x, -(n+1)/x], [(n+a+1)/x, (x-n-1)/x]].
 *
 * Pair: LegendreP/LegendreQ Casoratian P[n+1] Q[n] - P[n] Q[n+1] = 1/(n+1)
 *   (Abel: (1-x^2) W{P_n, Q_n} is constant; the recurrence turns W into the
 *    Casoratian; checked numerically inside and off (-1, 1)).
 *
 * Cross-family canon:
 *   ChebyshevT[n] = (ChebyshevU[n] - ChebyshevU[n-2])/2   (every n; when a U
 *                                            of the same x and class is present)
 *   LaguerreL[n, x] = LaguerreL[n, 0, x]    (when a generalized L of x is
 *                                            present; refolded afterwards)
 *
 * Reflections (canon):
 *   ChebyshevT[-n] = ChebyshevT[n]          (cos is even; every n)
 *   ChebyshevU[-n] = -ChebyshevU[n - 2]     (sin((1-n)t) = -sin((n-1)t))
 *   LegendreP[-n-1] = LegendreP[n]          (the ODE depends on n(n+1) only,
 *                                            and P is the solution regular
 *                                            at x = 1; every n), with a
 *                                            class-preserving canonical sign
 *   parity C[n, -x] = (-1)^n C[n, x] for P, T, U, H at INTEGER n proved.
 */

#include <stdlib.h>
#include <string.h>
#include "expr.h"
#include "eval.h"
#include "sym_names.h"
#include "simp.h"
#include "simp_internal.h"
#include "simp_sf_engine.h"

/* M rows 0..1: [[0, 1, 0], [A, B, 0]] for the order n = a[0]. */
static void three_term(Expr** M, Expr* A, Expr* B) {
    M[0] = sf_int(0); M[1] = sf_int(1); M[2] = sf_int(0);
    M[3] = A;         M[4] = B;         M[5] = sf_int(0);
}

static bool up_legendre(Expr* const* a, int d, Expr** M) {
    if (d != 0) return false;
    const Expr* n = a[0]; const Expr* x = a[1];
    /* A and B in separate statements: B consumes den, and C leaves the
     * order of a call's argument evaluation unspecified */
    Expr* den = sf_plus(sf_cp(n), sf_int(2));
    Expr* A = sf_times(sf_int(-1), sf_div(sf_plus(sf_cp(n), sf_int(1)), sf_cp(den)));
    Expr* B = sf_div(sf_times(sf_plus(sf_times(sf_int(2), sf_cp(n)), sf_int(3)), sf_cp(x)), den);
    three_term(M, A, B);
    return true;
}

static bool up_legendre3(Expr* const* a, int d, Expr** M) {
    if (d != 0) return false;
    const Expr* n = a[0]; const Expr* m = a[1]; const Expr* x = a[2];
    Expr* den = sf_plus(sf_minus(sf_cp(n), sf_cp(m)), sf_int(2));
    Expr* A = sf_times(sf_int(-1), sf_div(sf_plus(sf_plus(sf_cp(n), sf_cp(m)), sf_int(1)), sf_cp(den)));
    Expr* B = sf_div(sf_times(sf_plus(sf_times(sf_int(2), sf_cp(n)), sf_int(3)), sf_cp(x)), den);
    three_term(M, A, B);
    return true;
}

static bool up_chebyshev(Expr* const* a, int d, Expr** M) {
    if (d != 0) return false;
    three_term(M, sf_int(-1), sf_times(sf_int(2), sf_cp(a[1])));
    return true;
}

static bool up_hermite(Expr* const* a, int d, Expr** M) {
    if (d != 0) return false;
    three_term(M, sf_times(sf_int(-2), sf_plus(sf_cp(a[0]), sf_int(1))), sf_times(sf_int(2), sf_cp(a[1])));
    return true;
}

static bool up_laguerre(Expr* const* a, int d, Expr** M) {
    if (d != 0) return false;
    const Expr* n = a[0]; const Expr* x = a[1];
    Expr* den = sf_plus(sf_cp(n), sf_int(2));
    Expr* A = sf_times(sf_int(-1), sf_div(sf_plus(sf_cp(n), sf_int(1)), sf_cp(den)));
    Expr* B = sf_div(sf_minus(sf_plus(sf_times(sf_int(2), sf_cp(n)), sf_int(3)), sf_cp(x)), den);
    three_term(M, A, B);
    return true;
}

static bool up_laguerre3(Expr* const* a, int d, Expr** M) {
    const Expr* n = a[0]; const Expr* al = a[1]; const Expr* x = a[2];
    if (d == 0) {
        Expr* den = sf_plus(sf_cp(n), sf_int(2));
        Expr* A = sf_times(sf_int(-1), sf_div(sf_plus(sf_plus(sf_cp(n), sf_int(1)), sf_cp(al)), sf_cp(den)));
        Expr* B = sf_div(sf_minus(sf_plus(sf_plus(sf_times(sf_int(2), sf_cp(n)), sf_int(3)), sf_cp(al)), sf_cp(x)), den);
        three_term(M, A, B);
        return true;
    }
    if (d == 1) {
        Expr* na1 = sf_plus(sf_plus(sf_cp(n), sf_cp(al)), sf_int(1));
        M[0] = sf_div(sf_cp(na1), sf_cp(x));
        M[1] = sf_times(sf_int(-1), sf_div(sf_plus(sf_cp(n), sf_int(1)), sf_cp(x)));
        M[2] = sf_int(0);
        M[3] = sf_div(na1, sf_cp(x));
        M[4] = sf_div(sf_minus(sf_minus(sf_cp(x), sf_cp(n)), sf_int(1)), sf_cp(x));
        M[5] = sf_int(0);
        return true;
    }
    return false;
}

/* P[n+1] Q[n] - P[n] Q[n+1] - 1/(n+1) */
static Expr* rel_PQ(Expr* const* a, const Expr* c0, const Expr* c1, const Expr* d0, const Expr* d1) {
    Expr* w = sf_pow(sf_plus(sf_cp(a[0]), sf_int(1)), sf_int(-1));
    return sf_minus(sf_minus(sf_times(sf_cp(c1), sf_cp(d0)), sf_times(sf_cp(c0), sf_cp(d1))), w);
}

/* ---- reflections ------------------------------------------------------------ */


static Expr* neg_expand(const Expr* e) { return sf_mk1("Expand", sf_times(sf_int(-1), sf_cp(e))); }

/* C[n, -x] = (-1)^n C[n, x] at integer n proved (P, T, U, H). */
static Expr* integer_parity(const Expr* atom, const AssumeCtx* ctx) {
    const Expr* n = atom->data.function.args[0];
    const Expr* x = atom->data.function.args[1];
    if (!ctx || !sf_reads_negative(x) || !sf_proved_integer(ctx, n)) return NULL;
    const char* h = atom->data.function.head->data.symbol.name;
    /* n = c + s: (-1)^n written (-1)^c (-1)^s, so it meets the (-1)^s the
     * assumption rules leave on the input ((-1)^(1 + n) does not fold) */
    Expr* sp = NULL; int64_t num, den;
    Expr* sign;
    if (sf_split(n, &sp, &num, &den) && den == 1) {
        sign = sf_int((num % 2 == 0) ? 1 : -1);
        if (sp) sign = sf_times(sign, sf_pow(sf_int(-1), sf_cp(sp)));
    } else {
        sign = sf_pow(sf_int(-1), sf_cp(n));
    }
    if (sp) expr_free(sp);
    return sf_times(sign, sf_call(h, (Expr*[]){ sf_cp(n), neg_expand(x) }, 2));
}

/* Do two orders lie in one class (same symbolic part, rational part mod 1)? */
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

/* Is there a head[nu, ..., x] (argument `xi` == x) with nu in n's class? */
static bool has_kin(const Expr* e, const char* head, size_t arity, size_t xi,
                    const Expr* n, const Expr* x) {
    if (e->type != EXPR_FUNCTION) return false;
    if (sf_is_call(e, head, arity) && expr_eq(e->data.function.args[xi], x) &&
        class_eq(e->data.function.args[0], n))
        return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (has_kin(e->data.function.args[i], head, arity, xi, n, x)) return true;
    return false;
}

/* T[n] = (U[n] - U[n-2])/2 for every n (cos(n t) = (sin((n+1)t) -
 * sin((n-1)t))/(2 sin t)), applied when a U of the same argument and order
 * class is present -- D[ChebyshevT[n, x], x] = n ChebyshevU[n - 1, x] puts
 * both in the Chebyshev ODE. */
static Expr* canon_chebyshevT(const Expr* atom, const Expr* root, const AssumeCtx* ctx) {
    const Expr* n = atom->data.function.args[0];
    const Expr* x = atom->data.function.args[1];
    if (sf_reflection_flips(n, 0))          /* T[-n] = T[n], class-preserving */
        return sf_call(SYM_ChebyshevT, (Expr*[]){ neg_expand(n), sf_cp(x) }, 2);
    if (has_kin(root, SYM_ChebyshevU, 2, 1, n, x)) {
        Expr* u0 = sf_call(SYM_ChebyshevU, (Expr*[]){ sf_cp(n), sf_cp(x) }, 2);
        Expr* u2 = sf_call(SYM_ChebyshevU, (Expr*[]){ sf_mk1("Expand", sf_minus(sf_cp(n), sf_int(2))), sf_cp(x) }, 2);
        return sf_times(sf_rat(1, 2), sf_minus(u0, u2));
    }
    return integer_parity(atom, ctx);
}

/* LaguerreL[n, x] = LaguerreL[n, 0, x], applied when a generalized LaguerreL
 * of the same x is present (D[LaguerreL[n, x], x] = -LaguerreL[n-1, 1, x]
 * puts both in the Laguerre ODE), so the two meet on the (n, a) lattice. */
static bool has_laguerre3_of(const Expr* e, const Expr* x) {
    if (e->type != EXPR_FUNCTION) return false;
    if (sf_is_call(e, SYM_LaguerreL, 3) && expr_eq(e->data.function.args[2], x)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (has_laguerre3_of(e->data.function.args[i], x)) return true;
    return false;
}

static Expr* canon_laguerre2(const Expr* atom, const Expr* root, const AssumeCtx* ctx) {
    (void)ctx;
    const Expr* x = atom->data.function.args[1];
    if (!has_laguerre3_of(root, x)) return NULL;
    return sf_call(SYM_LaguerreL, (Expr*[]){ sf_cp(atom->data.function.args[0]), sf_int(0), sf_cp(x) }, 3);
}

/* Refold LaguerreL[n, 0, x] -> LaguerreL[n, x] when the input used the
 * two-argument form. */
static Expr* refold_laguerre_walk(const Expr* e, bool* changed) {
    if (e->type != EXPR_FUNCTION) return sf_cp(e);
    size_t n = e->data.function.arg_count;
    Expr** a = malloc(sizeof(Expr*) * (n ? n : 1));
    for (size_t i = 0; i < n; i++) a[i] = refold_laguerre_walk(e->data.function.args[i], changed);
    Expr* r = expr_new_function(sf_cp(e->data.function.head), a, n);
    free(a);
    if (sf_is_call(r, SYM_LaguerreL, 3) && r->data.function.args[1]->type == EXPR_INTEGER &&
        r->data.function.args[1]->data.integer == 0) {
        Expr* f = sf_call(SYM_LaguerreL, (Expr*[]){ sf_cp(r->data.function.args[0]), sf_cp(r->data.function.args[2]) }, 2);
        expr_free(r);
        *changed = true;
        return f;
    }
    return r;
}

static bool has_laguerre2(const Expr* e) {
    if (e->type != EXPR_FUNCTION) return false;
    if (sf_is_call(e, SYM_LaguerreL, 2)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (has_laguerre2(e->data.function.args[i])) return true;
    return false;
}

static Expr* pass_refold_laguerre(const Expr* e, const Expr* input, const AssumeCtx* ctx) {
    (void)ctx;
    if (!has_laguerre2(input)) return NULL;
    bool changed = false;
    Expr* r = refold_laguerre_walk(e, &changed);
    if (!changed) { expr_free(r); return NULL; }
    return r;
}

static Expr* canon_chebyshevU(const Expr* atom, const Expr* root, const AssumeCtx* ctx) {
    (void)root;
    const Expr* n = atom->data.function.args[0];
    /* U[-n] = -U[n - 2]: -n - 2 lies in the class of -n, so the flip uses
     * the class-preserving rule (k = 2: U[-1] = -U[-1] = 0, U[-3/2] =
     * -U[-1/2]) */
    if (sf_reflection_flips(n, 2)) {
        Expr* flipped = sf_mk1("Expand", sf_minus(neg_expand(n), sf_int(2)));
        return sf_times(sf_int(-1), sf_call(SYM_ChebyshevU, (Expr*[]){ flipped, sf_cp(atom->data.function.args[1]) }, 2));
    }
    return integer_parity(atom, ctx);
}

/* LegendreP[mu] = LegendreP[-mu - 1].  -mu - 1 lies in the class of -mu
 * (mod 1), so the flip is class-preserving (sf_reflection_flips with k = 1:
 * P[-1] = P[0], P[-3/2] = P[1/2]). */

static Expr* canon_legendreP(const Expr* atom, const Expr* root, const AssumeCtx* ctx) {
    (void)root;
    const Expr* n = atom->data.function.args[0];
    if (sf_reflection_flips(n, 1)) {
        Expr* m = sf_mk1("Expand", sf_minus(neg_expand(n), sf_int(1)));
        return sf_call(SYM_LegendreP, (Expr*[]){ m, sf_cp(atom->data.function.args[1]) }, 2);
    }
    return integer_parity(atom, ctx);
}

static Expr* canon_parity_only(const Expr* atom, const Expr* root, const AssumeCtx* ctx) {
    (void)root; return integer_parity(atom, ctx);
}

static const SfDesc ORTHO_DESCS[] = {
    { &SYM_LegendreP,  2, 1, 2, { 0 },    up_legendre,  NULL, canon_legendreP,   &SYM_LegendreQ, rel_PQ },
    { &SYM_LegendreQ,  2, 1, 2, { 0 },    up_legendre,  NULL, NULL,              NULL, NULL },
    { &SYM_LegendreP,  3, 1, 2, { 0 },    up_legendre3, NULL, NULL,              NULL, NULL },
    { &SYM_ChebyshevT, 2, 1, 2, { 0 },    up_chebyshev, NULL, canon_chebyshevT,  NULL, NULL },
    { &SYM_ChebyshevU, 2, 1, 2, { 0 },    up_chebyshev, NULL, canon_chebyshevU,  NULL, NULL },
    { &SYM_HermiteH,   2, 1, 2, { 0 },    up_hermite,   NULL, canon_parity_only, NULL, NULL },
    { &SYM_LaguerreL,  2, 1, 2, { 0 },    up_laguerre,  NULL, canon_laguerre2,   NULL, NULL },
    { &SYM_LaguerreL,  3, 2, 2, { 0, 1 }, up_laguerre3, NULL, NULL,              NULL, NULL },
};

static const SfPass ORTHO_PASSES[] = { pass_refold_laguerre };
static const char** ORTHO_PASS_HEADS[] = { &SYM_LaguerreL };

const SfModule SF_MODULE_ORTHOPOLY = {
    ORTHO_DESCS, sizeof ORTHO_DESCS / sizeof ORTHO_DESCS[0],
    ORTHO_PASSES, sizeof ORTHO_PASSES / sizeof ORTHO_PASSES[0],
    ORTHO_PASS_HEADS, sizeof ORTHO_PASS_HEADS / sizeof ORTHO_PASS_HEADS[0],
};
