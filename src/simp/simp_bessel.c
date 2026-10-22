/*
 * simp_bessel.c -- Simplify pass: Bessel functions of contiguous order.
 *
 * Bessel functions of one kind and argument whose orders differ by integers
 * are not independent.  The three-term recurrences (DLMF 10.6.1, 10.29.1)
 *
 *     J, Y:  C[m-1] + C[m+1] =  (2m/z) C[m]
 *     I:     I[m-1] - I[m+1] =  (2m/z) I[m]
 *     K:     K[m-1] - K[m+1] = -(2m/z) K[m]
 *
 * all have the shape  C[m+1] = p(m) C[m] + q C[m-1]  with q = -1 (J, Y) or
 * q = +1 (I, K), and p(m) = 2m/z (J, Y, K) or -2m/z (I).  Iterating it
 * expresses every C[b+k] as  A_k C[b] + B_k C[b+1]  with A_k, B_k rational in
 * z (and in the order, when it is symbolic) -- the Lommel polynomials.  So a
 * polynomial in {C[b+k]} rewrites to a polynomial in the basis pair
 * {C[b], C[b+1]} over Q(m)(z).
 *
 * That rewrite is a canonical form, not just a heuristic: for an order
 * m not in 1/2 + Z, C[m] and C[m+1] (equivalently C[m] and C[m]') are
 * algebraically independent over C(z) (Siegel), so after the rewrite the
 * Together numerator vanishes iff the expression does.  This is what makes
 *
 *     Simplify[D[Integrate[x^2 BesselJ[1/3, x]^3, x], x]]
 *
 * -- a polynomial in BesselJ at orders -2/3, 1/3, 4/3, 7/3 -- collapse back to
 * x^2 BesselJ[1/3, x]^3.  (Half-integer orders are elementary; the rewrite is
 * still an identity there, merely not a complete zero test.)
 *
 * Grouping needs no evaluation: an order is split as  s + c  with c its
 * rational numeric addend, and two kernels belong to one family iff they share
 * head, argument z, symbolic part s, and c mod 1 -- i.e. their orders differ by
 * an integer.  Symbolic orders (n - 1, n, n + 1) group the same way.
 *
 * Contract (mirrors simp_radical_rational): input borrowed; returns a fresh
 * Expr only when the rewritten form scores STRICTLY below the input (then
 * dispatched further), else NULL and simp_dispatch continues unchanged.
 */

#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include "expr.h"
#include "eval.h"
#include "sym_names.h"
#include "simp.h"
#include "simp_internal.h"

/* One Bessel kernel occurrence class: head, argument, symbolic order part,
 * and the rational order constant c = num/den (den shared by the family). */
typedef struct {
    const char*  head;       /* SYM_BesselJ / Y / I / K (interned)           */
    const Expr*  z;          /* borrowed argument                            */
    const Expr*  s;          /* borrowed symbolic order part, or NULL        */
    int64_t      den;        /* denominator of c (>= 1)                      */
    int64_t      cls;        /* num mod den: the c mod 1 class               */
    /* members: offsets (num - base_num)/den relative to base_num, plus the
     * original kernel atom each offset was read from */
    int64_t      base_num;
    int64_t*     off;
    const Expr** atom;
    size_t       n, cap;
} BesselFamily;

typedef struct {
    BesselFamily* f;
    size_t n, cap;
} FamilySet;

static bool is_bessel_head(const char* h) {
    return h == SYM_BesselJ || h == SYM_BesselY || h == SYM_BesselI || h == SYM_BesselK;
}

/* Read a machine rational: Integer or Rational[Integer, Integer]. */
static bool read_rational(const Expr* e, int64_t* num, int64_t* den) {
    if (e->type == EXPR_INTEGER) { *num = e->data.integer; *den = 1; return true; }
    if (e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL &&
        e->data.function.head->data.symbol.name == SYM_Rational &&
        e->data.function.arg_count == 2 &&
        e->data.function.args[0]->type == EXPR_INTEGER &&
        e->data.function.args[1]->type == EXPR_INTEGER) {
        *num = e->data.function.args[0]->data.integer;
        *den = e->data.function.args[1]->data.integer;
        return *den > 0;
    }
    return false;
}

/* Split an order mu as s + num/den.  `s` is NULL for a numeric order; for
 * Plus[c, rest...] with a single numeric addend the remaining terms must form
 * exactly one argument (the evaluator puts the number first), otherwise the
 * whole order is symbolic with c = 0.  Returns false on a bignum constant. */
static bool split_order(const Expr* mu, const Expr** s, int64_t* num, int64_t* den) {
    if (read_rational(mu, num, den)) { *s = NULL; return true; }
    if (mu->type == EXPR_INTEGER || mu->type == EXPR_BIGINT) return false;
    if (mu->type == EXPR_FUNCTION && mu->data.function.head->type == EXPR_SYMBOL &&
        mu->data.function.head->data.symbol.name == SYM_Plus &&
        mu->data.function.arg_count == 2 &&
        read_rational(mu->data.function.args[0], num, den)) {
        *s = mu->data.function.args[1];
        return true;
    }
    *s = mu; *num = 0; *den = 1;
    return true;
}

static int64_t mod_pos(int64_t a, int64_t m) { int64_t r = a % m; return r < 0 ? r + m : r; }

static void family_add(BesselFamily* f, int64_t num, const Expr* atom) {
    int64_t k = (num - f->base_num) / f->den;
    for (size_t i = 0; i < f->n; i++) if (f->off[i] == k) return;
    if (f->n == f->cap) {
        f->cap = f->cap ? f->cap * 2 : 4;
        f->off  = realloc(f->off,  sizeof(int64_t) * f->cap);
        f->atom = realloc((void*)f->atom, sizeof(Expr*) * f->cap);
    }
    f->off[f->n] = k; f->atom[f->n] = atom; f->n++;
}

static void collect_families(const Expr* e, FamilySet* fs) {
    if (e->type != EXPR_FUNCTION) return;
    const Expr* h = e->data.function.head;
    if (h->type == EXPR_SYMBOL && is_bessel_head(h->data.symbol.name) &&
        e->data.function.arg_count == 2) {
        const Expr* s; int64_t num, den;
        if (split_order(e->data.function.args[0], &s, &num, &den)) {
            const Expr* z = e->data.function.args[1];
            int64_t cls = mod_pos(num, den);
            for (size_t i = 0; i < fs->n; i++) {
                BesselFamily* f = &fs->f[i];
                if (f->head == h->data.symbol.name && f->den == den && f->cls == cls &&
                    expr_eq((Expr*)f->z, (Expr*)z) &&
                    ((!f->s && !s) || (f->s && s && expr_eq((Expr*)f->s, (Expr*)s)))) {
                    family_add(f, num, e);
                    goto args;
                }
            }
            if (fs->n == fs->cap) {
                fs->cap = fs->cap ? fs->cap * 2 : 4;
                fs->f = realloc(fs->f, sizeof(BesselFamily) * fs->cap);
            }
            BesselFamily* f = &fs->f[fs->n++];
            memset(f, 0, sizeof *f);
            f->head = h->data.symbol.name; f->z = z; f->s = s;
            f->den = den; f->cls = cls; f->base_num = num;
            family_add(f, num, e);
        }
    }
args:
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        collect_families(e->data.function.args[i], fs);
}

static void families_free(FamilySet* fs) {
    for (size_t i = 0; i < fs->n; i++) { free(fs->f[i].off); free((void*)fs->f[i].atom); }
    free(fs->f);
}

/* A family is worth rewriting when it has >= 2 orders that are not exactly
 * one adjacent pair {m, m+1} (that pair is already a basis). */
static bool family_triggers(const BesselFamily* f) {
    if (f->n < 2) return false;
    if (f->n == 2 && (f->off[0] - f->off[1] == 1 || f->off[1] - f->off[0] == 1)) return false;
    return true;
}

/* ---- expression builders (all evaluate) ---------------------------------- */

static Expr* mk2(const char* head, Expr* a, Expr* b) {
    return eval_and_free(expr_new_function(expr_new_symbol(head), (Expr*[]){ a, b }, 2));
}

/* The order at offset k: s + (base_num + k den)/den. */
static Expr* order_at(const BesselFamily* f, int64_t k) {
    Expr* c = eval_and_free(expr_new_function(expr_new_symbol(SYM_Rational),
        (Expr*[]){ expr_new_integer(f->base_num + k * f->den), expr_new_integer(f->den) }, 2));
    return f->s ? mk2(SYM_Plus, expr_copy((Expr*)f->s), c) : c;
}

/* p(mu) for the family's kind, mu the order at offset k: +-2 mu / z. */
static Expr* rec_p(const BesselFamily* f, int64_t k) {
    int64_t sign = (f->head == SYM_BesselI) ? -2 : 2;
    Expr* zinv = mk2(SYM_Power, expr_copy((Expr*)f->z), expr_new_integer(-1));
    return mk2(SYM_Times, mk2(SYM_Times, expr_new_integer(sign), order_at(f, k)), zinv);
}
static int64_t rec_q(const BesselFamily* f) {
    return (f->head == SYM_BesselJ || f->head == SYM_BesselY) ? -1 : 1;
}

static Expr* together(Expr* e) {
    Expr* t = call_unary_copy("Together", e);
    if (!t) return e;
    expr_free(e);
    return t;
}

/* lin = q*a + p*b  (p owned, a/b borrowed), Together'd. */
static Expr* lin_comb(int64_t q, const Expr* a, Expr* p, const Expr* b) {
    Expr* t1 = mk2(SYM_Times, expr_new_integer(q), expr_copy((Expr*)a));
    Expr* t2 = mk2(SYM_Times, p, expr_copy((Expr*)b));
    return together(mk2(SYM_Plus, t1, t2));
}

/* ---- rewrite one family onto basis {b, b+1} ------------------------------ */

typedef struct { const Expr** from; Expr** to; size_t n; } SubstMap;

static Expr* subst_apply(const Expr* e, const SubstMap* m) {
    for (size_t i = 0; i < m->n; i++)
        if (expr_eq((Expr*)e, (Expr*)m->from[i])) return expr_copy(m->to[i]);
    if (e->type != EXPR_FUNCTION) return expr_copy((Expr*)e);
    size_t n = e->data.function.arg_count;
    Expr** args = malloc(sizeof(Expr*) * (n ? n : 1));
    for (size_t i = 0; i < n; i++) args[i] = subst_apply(e->data.function.args[i], m);
    Expr* r = expr_new_function(expr_copy(e->data.function.head), args, n);
    free(args);
    return r;
}

/* Returns the evaluated rewrite of `e` with every member of `f` expressed in
 * C[b], C[b+1]; NULL on allocation failure. */
static Expr* rewrite_family(const Expr* e, const BesselFamily* f, int64_t b) {
    int64_t lo = b, hi = b + 1;
    for (size_t i = 0; i < f->n; i++) {
        if (f->off[i] < lo) lo = f->off[i];
        if (f->off[i] > hi) hi = f->off[i];
    }
    size_t span = (size_t)(hi - lo + 1);
    Expr** A = calloc(span, sizeof(Expr*));
    Expr** B = calloc(span, sizeof(Expr*));
    if (!A || !B) { free(A); free(B); return NULL; }
    int64_t q = rec_q(f);
#define IX(k) ((size_t)((k) - lo))
    A[IX(b)] = expr_new_integer(1);     B[IX(b)] = expr_new_integer(0);
    A[IX(b + 1)] = expr_new_integer(0); B[IX(b + 1)] = expr_new_integer(1);
    /* upward: C[k+1] = p(k) C[k] + q C[k-1] */
    for (int64_t k = b + 1; k < hi; k++) {
        A[IX(k + 1)] = lin_comb(q, A[IX(k - 1)], rec_p(f, k), A[IX(k)]);
        B[IX(k + 1)] = lin_comb(q, B[IX(k - 1)], rec_p(f, k), B[IX(k)]);
    }
    /* downward: C[k-1] = q C[k+1] - q p(k) C[k]   (q = +-1, so 1/q = q) */
    for (int64_t k = b; k > lo; k--) {
        Expr* p1 = mk2(SYM_Times, expr_new_integer(-q), rec_p(f, k));
        Expr* p2 = mk2(SYM_Times, expr_new_integer(-q), rec_p(f, k));
        A[IX(k - 1)] = lin_comb(q, A[IX(k + 1)], p1, A[IX(k)]);
        B[IX(k - 1)] = lin_comb(q, B[IX(k + 1)], p2, B[IX(k)]);
    }
    Expr* cb  = expr_new_function(expr_new_symbol(f->head),
                    (Expr*[]){ order_at(f, b), expr_copy((Expr*)f->z) }, 2);
    Expr* cb1 = expr_new_function(expr_new_symbol(f->head),
                    (Expr*[]){ order_at(f, b + 1), expr_copy((Expr*)f->z) }, 2);
    SubstMap m = { f->atom, malloc(sizeof(Expr*) * f->n), f->n };
    for (size_t i = 0; i < f->n; i++) {
        int64_t k = f->off[i];
        m.to[i] = mk2(SYM_Plus, mk2(SYM_Times, expr_copy(A[IX(k)]), expr_copy(cb)),
                                mk2(SYM_Times, expr_copy(B[IX(k)]), expr_copy(cb1)));
    }
#undef IX
    Expr* out = eval_and_free(subst_apply(e, &m));
    for (size_t i = 0; i < f->n; i++) expr_free(m.to[i]);
    free(m.to);
    for (size_t i = 0; i < span; i++) { if (A[i]) expr_free(A[i]); if (B[i]) expr_free(B[i]); }
    free(A); free(B);
    expr_free(cb); expr_free(cb1);
    return out;
}

/* Lowest-scoring of the plain rewrite, its Together form and its Expand form
 * (owned in, owned out).  Together cancels over a common denominator; Expand
 * is needed when the coefficients carry symbolic powers -- x^n (n/x) Y[n] and
 * n x^(n-1) Y[n] only cancel once the product is distributed, which Together
 * does not do for a non-integer exponent. */
static Expr* best_form(Expr* r, const Expr* cf, size_t* score) {
    static const char* const forms[] = { "Together", "Expand" };
    size_t best = score_with_func(r, cf);
    for (size_t i = 0; i < sizeof forms / sizeof forms[0]; i++) {
        Expr* t = call_unary_copy(forms[i], r);
        if (!t) continue;
        size_t st = score_with_func(t, cf);
        if (st < best) { expr_free(r); r = t; best = st; }
        else expr_free(t);
    }
    *score = best;
    return r;
}

Expr* simp_bessel_contiguous(const Expr* input, const AssumeCtx* ctx,
                             const Expr* complexity_func) {
    FamilySet fs = { NULL, 0, 0 };
    collect_families(input, &fs);
    bool any = false;
    for (size_t i = 0; i < fs.n && !any; i++) any = family_triggers(&fs.f[i]);
    if (!any) { families_free(&fs); return NULL; }

    size_t in_score = score_with_func(input, complexity_func);
    Expr* cur = expr_copy((Expr*)input);
    size_t cur_score = in_score;
    bool changed = false;

    /* Families are independent: choose each basis greedily.  Candidate bases
     * {b, b+1} touch a present order (b or b+1 is one) inside the family span;
     * the choice only affects how a non-collapsing remainder is presented. */
    for (size_t fi = 0; fi < fs.n; fi++) {
        const BesselFamily* f = &fs.f[fi];
        if (!family_triggers(f)) continue;
        int64_t lo = f->off[0], hi = f->off[0];
        for (size_t i = 1; i < f->n; i++) {
            if (f->off[i] < lo) lo = f->off[i];
            if (f->off[i] > hi) hi = f->off[i];
        }
        Expr* fam_best = NULL;
        size_t fam_score = 0;
        for (int64_t b = lo; b < hi; b++) {
            bool touches = false;
            for (size_t i = 0; i < f->n && !touches; i++)
                touches = (f->off[i] == b || f->off[i] == b + 1);
            if (!touches) continue;
            /* f->atom borrows kernels of the original input; they are still
             * present verbatim in `cur`, since an earlier family's rewrite
             * only replaces kernels with a different family key. */
            Expr* r = rewrite_family(cur, f, b);
            if (!r) continue;
            size_t s;
            r = best_form(r, complexity_func, &s);
            if (!fam_best || s < fam_score) {
                if (fam_best) expr_free(fam_best);
                fam_best = r; fam_score = s;
            } else {
                expr_free(r);
            }
        }
        if (fam_best) {
            expr_free(cur);
            cur = fam_best; cur_score = fam_score;
            changed = true;
        }
    }
    families_free(&fs);

    if (!changed || cur_score >= in_score) { expr_free(cur); return NULL; }
    if (simp_debug_enabled()) simp_debug_log("BesselContiguous", input, cur, 0.0);

    /* Every family is now on an adjacent basis pair, so this recursion does
     * not re-enter the rewrite. */
    Expr* fin = simp_dispatch(cur, ctx, complexity_func);
    if (fin && score_with_func(fin, complexity_func) <= cur_score) {
        expr_free(cur);
        return fin;
    }
    if (fin) expr_free(fin);
    return cur;
}
