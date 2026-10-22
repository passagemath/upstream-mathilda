/*
 * simp_bessel.c -- special-function Simplify descriptors: Bessel J, Y, I, K.
 *
 * Bessel functions of one kind and argument whose orders differ by integers
 * are not independent.  The three-term recurrences (DLMF 10.6.1, 10.29.1)
 *
 *     J, Y:  C[m-1] + C[m+1] =  (2m/z) C[m]
 *     I:     I[m-1] - I[m+1] =  (2m/z) I[m]
 *     K:     K[m-1] - K[m+1] = -(2m/z) K[m]
 *
 * all have the shape  C[m+1] = p(m) C[m] + q C[m-1]  with q = -1 (J, Y) or
 * q = +1 (I, K), and p(m) = 2m/z (J, Y, K) or -2m/z (I).  In the engine's
 * basis v(m) = (C[m], C[m+1], 1) the up step is
 *
 *     v(m+1) = [[0, 1, 0], [q, p(m+1), 0], [0, 0, 1]] v(m),
 *
 * so every C[b+k] becomes A_k C[b] + B_k C[b+1] with A_k, B_k rational in z
 * (and in the order, when it is symbolic) -- the Lommel polynomials.
 *
 * That rewrite is a canonical form, not just a heuristic: for an order m not
 * in 1/2 + Z, C[m] and C[m+1] (equivalently C[m] and C[m]') are algebraically
 * independent over C(z) (Siegel), so after the rewrite the Together numerator
 * vanishes iff the expression does.  This is what makes
 *
 *     Simplify[D[Integrate[x^2 BesselJ[1/3, x]^3, x], x]]
 *
 * collapse back to x^2 BesselJ[1/3, x]^3.  (Half-integer orders are
 * elementary; the rewrite is still an identity there, merely not a complete
 * zero test.)
 *
 * Pairs.  J and Y (I and K) obey the same recurrence, so their Casoratian is
 * a recurrence invariant up to the factor the recurrence carries; evaluating
 * it from the Wronskian and the recurrence gives (DLMF 10.5.5, 10.28.2)
 *
 *     J[m+1] Y[m] - J[m] Y[m+1] = 2/(Pi z)
 *     I[m+1] K[m] + I[m] K[m+1] = 1/z
 *
 * Once both families sit on one basis {m, m+1}, the expression is reduced
 * modulo that relation (eliminating the monomial C[m+1] D[m]), which makes
 * the four generators canonical.
 *
 * Reflection.  K is even in its order for every order (from the integral
 * K_nu(z) = Int_0^Inf e^{-z cosh t} cosh(nu t) dt), so BesselK[mu] and
 * BesselK[-mu] are one value; the canonical sign keeps families intact (see
 * K_flips).  For an integer order n, J[-n] = (-1)^n J[n], Y[-n] = (-1)^n
 * Y[n] and I[-n] = I[n] (from the generating function); those are applied
 * only when the order is PROVABLY an integer under the assumptions.
 */

#include "expr.h"
#include "eval.h"
#include "sym_names.h"
#include "simp.h"
#include "simp_internal.h"
#include "simp_sf_engine.h"

/* p(mu) = +-2 mu / z for the kind; q = -1 (J, Y) / +1 (I, K). */
static Expr* rec_p(const char* head, const Expr* mu, const Expr* z) {
    int64_t sign = (head == SYM_BesselI) ? -2 : 2;
    return sf_div(sf_times(sf_int(sign), sf_cp(mu)), sf_cp(z));
}

static bool bessel_up_kind(const char* head, Expr* const* a, int d, Expr** M) {
    if (d != 0) return false;
    int64_t q = (head == SYM_BesselJ || head == SYM_BesselY) ? -1 : 1;
    Expr* m1 = sf_plus(sf_cp(a[0]), sf_int(1));
    M[0] = sf_int(0); M[1] = sf_int(1);        M[2] = sf_int(0);
    M[3] = sf_int(q); M[4] = rec_p(head, m1, a[1]); M[5] = sf_int(0);
    expr_free(m1);
    return true;
}

static bool up_J(Expr* const* a, int d, Expr** M) { return bessel_up_kind(SYM_BesselJ, a, d, M); }
static bool up_Y(Expr* const* a, int d, Expr** M) { return bessel_up_kind(SYM_BesselY, a, d, M); }
static bool up_I(Expr* const* a, int d, Expr** M) { return bessel_up_kind(SYM_BesselI, a, d, M); }
static bool up_K(Expr* const* a, int d, Expr** M) { return bessel_up_kind(SYM_BesselK, a, d, M); }

/* J[b+1] Y[b] - J[b] Y[b+1] - 2/(Pi z) */
static Expr* rel_JY(Expr* const* a, const Expr* c0, const Expr* c1,
                    const Expr* d0, const Expr* d1) {
    Expr* w = sf_div(sf_int(2), sf_times(sf_sym(SYM_Pi), sf_cp(a[1])));
    return sf_minus(sf_minus(sf_times(sf_cp(c1), sf_cp(d0)), sf_times(sf_cp(c0), sf_cp(d1))), w);
}

/* I[b+1] K[b] + I[b] K[b+1] - 1/z */
static Expr* rel_IK(Expr* const* a, const Expr* c0, const Expr* c1,
                    const Expr* d0, const Expr* d1) {
    Expr* w = sf_pow(sf_cp(a[1]), sf_int(-1));
    return sf_minus(sf_plus(sf_times(sf_cp(c1), sf_cp(d0)), sf_times(sf_cp(c0), sf_cp(d1))), w);
}

/* Does the order read negative (its symbolic part, or the number itself)? */
static bool order_reads_negative(const Expr* mu) {
    Expr* s; int64_t num, den;
    if (!sf_split(mu, &s, &num, &den)) return false;
    bool neg = s ? sf_reads_negative(s) : num < 0;
    if (s) expr_free(s);
    return neg;
}

/* -mu, distributed (the evaluator leaves -(1 - n) as Times[-1, Plus[...]],
 * which would not match the order -1 + n of the family it joins). */
static Expr* neg_order(const Expr* mu) {
    return sf_mk1("Expand", sf_times(sf_int(-1), sf_cp(mu)));
}

/* Do two orders lie in one family class (same symbolic part, same rational
 * part mod 1)? */
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

/* Is there a BesselI[nu, z] in `e` with nu in mu's class? */
static bool has_I_in_class(const Expr* e, const Expr* z, const Expr* mu) {
    if (e->type != EXPR_FUNCTION) return false;
    if (sf_is_call(e, SYM_BesselI, 2) && expr_eq(e->data.function.args[1], z) &&
        class_eq(e->data.function.args[0], mu))
        return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (has_I_in_class(e->data.function.args[i], z, mu)) return true;
    return false;
}

/* Canonical sign for K's even order.  K must stay in the class of a BesselI
 * of the same argument when one is present (the I/K pair relation needs both
 * families on one basis), so: keep mu when an I shares its class, flip when
 * an I shares -mu's class.  Otherwise a symbolic order flips when its
 * symbolic part reads negative, and a numeric order mu = num/den flips
 * exactly when -mu has the smaller class: frac(mu) > 1/2, or frac(mu) in
 * {0, 1/2} with mu < 0.  Then K[-2/3] (class 1/3) stays beside K[1/3],
 * K[4/3], and K[2/3] (class 2/3) joins them as K[-2/3]. */
static bool K_flips(const Expr* mu, const Expr* z, const Expr* root) {
    if (has_I_in_class(root, z, mu)) return false;
    Expr* nm = neg_order(mu);
    bool partner_neg = has_I_in_class(root, z, nm);
    expr_free(nm);
    if (partner_neg) return true;
    Expr* s; int64_t num, den;
    if (!sf_split(mu, &s, &num, &den)) return false;
    if (s) { bool neg = sf_reads_negative(s); expr_free(s); return neg; }
    int64_t r = num % den; if (r < 0) r += den;     /* frac = r/den */
    if (2 * r > den) return true;
    if (2 * r == den || r == 0) return num < 0;
    return false;
}

static Expr* canon_K(const Expr* atom, const Expr* root, const AssumeCtx* ctx) {
    (void)ctx;
    const Expr* mu = atom->data.function.args[0];
    const Expr* z = atom->data.function.args[1];
    if (!K_flips(mu, z, root)) return NULL;
    return sf_call(SYM_BesselK, (Expr*[]){ neg_order(mu), sf_cp(z) }, 2);
}

/* Integer-order reflection for J, Y (sign (-1)^n) and I (no sign). */
static Expr* canon_int_reflect(const Expr* atom, const AssumeCtx* ctx, bool signed_) {
    const Expr* mu = atom->data.function.args[0];
    if (!ctx || !order_reads_negative(mu) || !sf_proved_integer(ctx, mu)) return NULL;
    const char* h = atom->data.function.head->data.symbol.name;
    Expr* flipped = sf_call(h, (Expr*[]){ neg_order(mu), sf_cp(atom->data.function.args[1]) }, 2);
    if (!signed_) return flipped;
    /* mu = c + s (c integer): (-1)^-mu = (-1)^c (-1)^-s, written with the
     * positive-reading symbolic part -s so it combines with a (-1)^n the
     * input already carries ((-1)^(n-1) would not fold against (-1)^n). */
    Expr* sp = NULL; int64_t num, den;
    if (!sf_split(mu, &sp, &num, &den) || den != 1) { expr_free(flipped); if (sp) expr_free(sp); return NULL; }
    Expr* sign = sf_int((num % 2 == 0) ? 1 : -1);
    if (sp) sign = sf_times(sign, sf_pow(sf_int(-1), neg_order(sp)));
    if (sp) expr_free(sp);
    return sf_call(SYM_Times, (Expr*[]){ sign, flipped }, 2);
}

static Expr* canon_J(const Expr* atom, const Expr* root, const AssumeCtx* ctx) {
    (void)root; return canon_int_reflect(atom, ctx, true);
}
static Expr* canon_I(const Expr* atom, const Expr* root, const AssumeCtx* ctx) {
    (void)root; return canon_int_reflect(atom, ctx, false);
}

static const SfDesc BESSEL_DESCS[] = {
    { &SYM_BesselJ, 2, 1, 2, { 0 }, up_J, NULL, canon_J, &SYM_BesselY, rel_JY },
    { &SYM_BesselY, 2, 1, 2, { 0 }, up_Y, NULL, canon_J, NULL, NULL },
    { &SYM_BesselI, 2, 1, 2, { 0 }, up_I, NULL, canon_I, &SYM_BesselK, rel_IK },
    { &SYM_BesselK, 2, 1, 2, { 0 }, up_K, NULL, canon_K, NULL, NULL },
};

const SfModule SF_MODULE_BESSEL = {
    BESSEL_DESCS, sizeof BESSEL_DESCS / sizeof BESSEL_DESCS[0], NULL, 0, NULL, 0
};
