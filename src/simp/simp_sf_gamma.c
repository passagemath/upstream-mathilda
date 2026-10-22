/*
 * simp_sf_gamma.c -- special-function Simplify descriptors: the Gamma family.
 *
 * Heads: Gamma[z], Gamma[a, z] (upper incomplete), LogGamma, BarnesG, and the
 * Gamma-ratio heads Pochhammer, Beta, Binomial, Factorial, which are
 * canonicalised onto Gamma.  Every identity below is derived in
 * docs/design/simp_special_functions.md (section Gamma) and checked
 * numerically by tests/test_simp_sf_gamma.c.
 *
 * Recurrences (rank 1, basis v(z) = (C[z], 1), up matrix rows given):
 *
 *   Gamma[z+1]    = z Gamma[z]                       [[z, 0]]
 *       Gamma(z+1) = Int_0^Inf t^z e^-t dt, integrate by parts once.
 *       Bulk: Gamma[z+k]/Gamma[z] = Pochhammer[z, k]  (k > 0),
 *             1/Pochhammer[z+k, -k]                   (k < 0).
 *   Gamma[a+1, z] = a Gamma[a, z] + z^a E^-z         [[a, z^a E^-z]]
 *       Gamma(a+1, z) = Int_z^Inf t^a e^-t dt, by parts:
 *       [-t^a e^-t]_z^Inf + a Int_z^Inf t^(a-1) e^-t dt.
 *   LogGamma[z+1] = LogGamma[z] + Log[z]             [[1, Log[z]]]
 *       Both sides are analytic on C \ (-Inf, 0] and agree for z > 0
 *       (log of the Gamma recurrence), hence agree on the whole slit plane:
 *       the principal-branch identity holds wherever LogGamma is defined.
 *   BarnesG[z+1]  = Gamma[z] BarnesG[z]              [[Gamma[z], 0]]
 *       The defining functional equation of the Barnes G-function.
 *
 * Canonical generators (applied when the expression holds >= 2 Gamma-type
 * atoms, so a lone Beta[a, b] is never expanded for nothing):
 *
 *   Factorial[n]     = Gamma[n + 1]
 *   Pochhammer[a, n] = Gamma[a + n] / Gamma[a]
 *   Beta[a, b]       = Gamma[a] Gamma[b] / Gamma[a + b]     (Euler)
 *   Binomial[n, k]   = Gamma[n + 1] / (Gamma[k + 1] Gamma[n - k + 1])
 *
 *   (as meromorphic identities; at the removable cases where Binomial is
 *   defined by a limit -- negative integer n -- the rewrite, like every
 *   Simplify rewrite, holds generically).
 *
 * Reflection (Euler): Gamma[w] Gamma[1 - w] = Pi / Sin[Pi w].  Gamma[mu] is
 * reflected onto Gamma[1 - mu] only when an atom of the class of 1 - mu is
 * present, and mu is the less canonical of the two (its symbolic part reads
 * negative; numerically, frac(mu) > 1/2).  For mu = k + r - t (k integer,
 * 0 <= r < 1, t the positive-reading symbolic part) the sine phase is folded:
 * Sin[Pi mu] = -(-1)^k Sin[Pi (t - r)], so the family of t stays intact.
 *
 * Multiplication (Gauss): Gamma[m z] = (2 Pi)^((1-m)/2) m^(m z - 1/2)
 * Prod_{j=0}^{m-1} Gamma[z + j/m], applied to Gamma[m s + c] when a Gamma of
 * symbolic part s is present (m = 2 is Legendre's duplication formula).
 *
 * Refold passes put Gamma ratios back into the heads the input was written
 * with (Factorial, Pochhammer, Beta, Binomial) when that scores lower.
 */

#include <stdlib.h>
#include <string.h>
#include "expr.h"
#include "eval.h"
#include "sym_names.h"
#include "simp.h"
#include "simp_internal.h"
#include "simp_sf_engine.h"

/* Gauss multiplication is a heuristic rewrite (one Gamma becomes m), not a
 * decision procedure: it is attempted for multipliers up to this bound. */
#define GAMMA_MULT_MAX 16

/* ---- recurrences ---------------------------------------------------------- */

static bool up_gamma(Expr* const* a, int d, Expr** M) {
    if (d != 0) return false;
    M[0] = sf_cp(a[0]); M[1] = sf_int(0);
    return true;
}

static Expr* bulk_gamma(Expr* const* a, int d, int64_t k) {
    if (d != 0) return NULL;
    if (k > 0) return sf_call("Pochhammer", (Expr*[]){ sf_cp(a[0]), sf_int(k) }, 2);
    Expr* lo = sf_plus(sf_cp(a[0]), sf_int(k));
    return sf_pow(sf_call("Pochhammer", (Expr*[]){ lo, sf_int(-k) }, 2), sf_int(-1));
}

static bool up_gamma_inc(Expr* const* a, int d, Expr** M) {
    if (d != 0) return false;
    M[0] = sf_cp(a[0]);
    M[1] = sf_times(sf_pow(sf_cp(a[1]), sf_cp(a[0])),
                    sf_pow(sf_sym(SYM_E), sf_times(sf_int(-1), sf_cp(a[1]))));
    return true;
}

static bool up_loggamma(Expr* const* a, int d, Expr** M) {
    if (d != 0) return false;
    M[0] = sf_int(1); M[1] = sf_mk1(SYM_Log, sf_cp(a[0]));
    return true;
}

static bool up_barnesg(Expr* const* a, int d, Expr** M) {
    if (d != 0) return false;
    M[0] = sf_mk1(SYM_Gamma, sf_cp(a[0])); M[1] = sf_int(0);
    return true;
}

/* ---- Gamma-type atoms and their implied Gamma arguments -------------------- */

/* The Gamma arguments an atom stands for (owned into out[], count returned;
 * 0 for a non-Gamma-type atom).  out must hold 3. */
static size_t gamma_args_of(const Expr* e, Expr** out) {
    if (e->type != EXPR_FUNCTION || e->data.function.head->type != EXPR_SYMBOL) return 0;
    const char* h = e->data.function.head->data.symbol.name;
    size_t n = e->data.function.arg_count;
    Expr* const* a = e->data.function.args;
    if (h == SYM_Gamma && n == 1) { out[0] = sf_cp(a[0]); return 1; }
    if (h == SYM_Factorial && n == 1) { out[0] = sf_plus(sf_cp(a[0]), sf_int(1)); return 1; }
    if (h == SYM_Pochhammer && n == 2) {
        out[0] = sf_cp(a[0]); out[1] = sf_plus(sf_cp(a[0]), sf_cp(a[1])); return 2;
    }
    if (h == SYM_Beta && n == 2) {
        out[0] = sf_cp(a[0]); out[1] = sf_cp(a[1]);
        out[2] = sf_plus(sf_cp(a[0]), sf_cp(a[1])); return 3;
    }
    if (h == SYM_Binomial && n == 2) {
        out[0] = sf_plus(sf_cp(a[0]), sf_int(1));
        out[1] = sf_plus(sf_cp(a[1]), sf_int(1));
        out[2] = sf_mk1("Expand", sf_plus(sf_minus(sf_cp(a[0]), sf_cp(a[1])), sf_int(1)));
        return 3;
    }
    return 0;
}

/* Number of implied Gamma arguments in e: a lone Beta[a, 1 - a] counts 3,
 * since its own arguments can interact (reflection, a + b = 1). */
static size_t count_gamma_type(const Expr* e) {
    if (e->type != EXPR_FUNCTION) return 0;
    Expr* tmp[3];
    size_t k = gamma_args_of(e, tmp);
    for (size_t i = 0; i < k; i++) expr_free(tmp[i]);
    size_t c = k;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        c += count_gamma_type(e->data.function.args[i]);
    return c;
}

/* Memo of count_gamma_type(root) across one canon walk (root is stable). */
static const Expr* g_count_root = NULL;
static size_t      g_count_val  = 0;
static size_t root_gamma_count(const Expr* root) {
    if (root != g_count_root) { g_count_root = root; g_count_val = count_gamma_type(root); }
    return g_count_val;
}

/* ---- class arithmetic ---------------------------------------------------- */

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

/* Is some implied Gamma argument in `root` in the class of mu?  `skip` is
 * the atom being canonicalised (its own arguments do not count). */
static bool root_has_class(const Expr* e, const Expr* mu, const Expr* skip) {
    if (e->type != EXPR_FUNCTION) return false;
    if (!(skip && expr_eq(e, skip))) {
        Expr* g[3];
        size_t k = gamma_args_of(e, g);
        bool hit = false;
        for (size_t i = 0; i < k; i++) {
            if (!hit && class_eq(g[i], mu)) hit = true;
            expr_free(g[i]);
        }
        if (hit) return true;
    }
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (root_has_class(e->data.function.args[i], mu, skip)) return true;
    return false;
}

/* Is there an implied Gamma argument whose symbolic part is exactly s? */
static bool symbolic_part_is(const Expr* g, const Expr* s) {
    Expr* sp; int64_t num, den;
    if (!sf_split(g, &sp, &num, &den)) return false;
    bool eq = sp && expr_eq(sp, s);
    if (sp) expr_free(sp);
    return eq;
}

static bool root_has_symbolic_part(const Expr* e, const Expr* s) {
    if (e->type != EXPR_FUNCTION) return false;
    Expr* g[3];
    size_t k = gamma_args_of(e, g);
    bool hit = false;
    for (size_t i = 0; i < k; i++) {
        if (!hit && symbolic_part_is(g[i], s)) hit = true;
        expr_free(g[i]);
    }
    if (hit) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (root_has_symbolic_part(e->data.function.args[i], s)) return true;
    return false;
}

static int64_t floor_div(int64_t a, int64_t b) {
    int64_t q = a / b;
    return (a % b != 0 && ((a < 0) != (b < 0))) ? q - 1 : q;
}

/* ---- reflection and multiplication ---------------------------------------- */

/* Gamma[mu] -> Pi / (Sin[Pi mu] Gamma[1 - mu]) with the sine phase folded,
 * when the class of 1 - mu is present and mu is the less canonical side. */
static Expr* gamma_reflect(const Expr* atom, const Expr* root) {
    const Expr* mu = atom->data.function.args[0];
    Expr* s; int64_t num, den;
    if (!sf_split(mu, &s, &num, &den)) return NULL;
    int64_t k = floor_div(num, den), r = num - k * den;    /* c = k + r/den */
    bool less_canonical = s ? sf_reads_negative(s) : (2 * r > den);
    if (!less_canonical) { if (s) expr_free(s); return NULL; }
    Expr* one_minus = sf_mk1("Expand", sf_minus(sf_int(1), sf_cp(mu)));
    bool partner = root_has_class(root, one_minus, atom);
    if (!partner) { expr_free(one_minus); if (s) expr_free(s); return NULL; }
    /* Sin[Pi mu] = Sin[Pi (1 - mu)]: write the sine of the partner's own
     * (positive-reading) argument, so it meets the input's trig verbatim.
     * In the half-integer class 1 - mu = k + 1/2 + t the sine is
     * (-1)^k Cos[Pi t]. */
    Expr* sine;
    Expr* ps; int64_t pn, pd;
    if (s && sf_split(one_minus, &ps, &pn, &pd) && ps && pd == 2) {
        int64_t pk = floor_div(pn, 2);        /* 1 - mu = pk + 1/2 + ps */
        sine = sf_times(sf_int((pk % 2 == 0) ? 1 : -1),
                        sf_mk1(SYM_Cos, sf_times(sf_sym(SYM_Pi), ps)));
    } else {
        if (s && ps) expr_free(ps);
        sine = sf_mk1(SYM_Sin, sf_times(sf_sym(SYM_Pi), sf_cp(one_minus)));
    }
    if (s) expr_free(s);
    Expr* den_e = sf_times(sine, sf_call(SYM_Gamma, (Expr*[]){ one_minus }, 1));
    return sf_div(sf_sym(SYM_Pi), den_e);
}

static int64_t gcd64(int64_t a, int64_t b) {
    if (a < 0) a = -a;
    if (b < 0) b = -b;
    while (b) { int64_t t = a % b; a = b; b = t; }
    return a;
}

/* Numeric coefficient of one term (1 when it has none; 0 if not integer). */
static int64_t term_coeff(const Expr* t) {
    if (t->type == EXPR_FUNCTION && t->data.function.head->type == EXPR_SYMBOL &&
        t->data.function.head->data.symbol.name == SYM_Times && t->data.function.arg_count >= 2) {
        const Expr* c = t->data.function.args[0];
        if (c->type == EXPR_INTEGER) return c->data.integer;
        if (sf_is_call(c, SYM_Rational, 2)) return 0;
    }
    return 1;
}

/* gcd of the integer coefficients of s's terms (1 for a non-integer one). */
static int64_t integer_content(const Expr* s) {
    if (s->type == EXPR_FUNCTION && s->data.function.head->type == EXPR_SYMBOL &&
        s->data.function.head->data.symbol.name == SYM_Plus) {
        int64_t g = 0;
        for (size_t i = 0; i < s->data.function.arg_count; i++) {
            int64_t c = term_coeff(s->data.function.args[i]);
            if (c == 0) return 1;
            g = gcd64(g, c);
        }
        return g ? g : 1;
    }
    int64_t c = term_coeff(s);
    return c > 0 ? c : (c < 0 ? -c : 1);
}

/* Gamma[m s + c] -> Gauss product over Gamma[s + c/m + j/m] when a Gamma of
 * symbolic part s is present. */
static Expr* gamma_multiply(const Expr* atom, const Expr* root) {
    Expr* mu = sf_mk1("Expand", sf_cp(atom->data.function.args[0]));   /* 2 (x+1) -> 2 + 2 x */
    Expr* s = NULL; int64_t num, den;
    bool ok = sf_split(mu, &s, &num, &den);
    expr_free(mu);
    if (!ok || !s) { if (s) expr_free(s); return NULL; }
    /* m = integer content of the symbolic part (2 x -> 2; 2 a + 2 b -> 2) */
    int64_t m = integer_content(s);
    if (m < 2 || m > GAMMA_MULT_MAX) { expr_free(s); return NULL; }
    Expr* s0 = sf_mk1("Expand", sf_times(sf_rat(1, m), sf_cp(s)));
    expr_free(s);
    s = s0;                          /* freed below as `s` */
    if (!root_has_symbolic_part(root, s0)) { expr_free(s); return NULL; }
    /* z = mu / m = s0 + num/(den m).  Prefactor (2 Pi)^((1-m)/2) m^(m z - 1/2)
     * with m = 2^e m' (m' odd) written as
     *     2^((1-m)/2 + e (m z - 1/2)) m'^(m z - 1/2) Pi^((1-m)/2),
     * one power per base -- the evaluator does not merge 2^u 4^v itself, and
     * the input's own prefactor normalises to the same per-base form. */
    Expr* z = sf_plus(sf_cp(s0), sf_rat(num, den * m));
    int64_t e2 = 0, modd = m;
    while (modd % 2 == 0) { modd /= 2; e2++; }
    Expr* mz_half = sf_mk1("Expand", sf_minus(sf_times(sf_int(m), sf_cp(z)), sf_rat(1, 2)));
    Expr* prod = sf_times(sf_pow(sf_int(2), sf_mk1("Expand", sf_plus(sf_rat(1 - m, 2),
                                                    sf_times(sf_int(e2), sf_cp(mz_half))))),
                          sf_pow(sf_sym(SYM_Pi), sf_rat(1 - m, 2)));
    if (modd > 1) prod = sf_times(prod, sf_pow(sf_int(modd), sf_cp(mz_half)));
    expr_free(mz_half);
    for (int64_t j = 0; j < m; j++) {
        Expr* arg = sf_plus(sf_cp(z), sf_rat(j, m));
        prod = sf_call(SYM_Times, (Expr*[]){ prod, sf_call(SYM_Gamma, (Expr*[]){ arg }, 1) }, 2);
    }
    expr_free(z);
    expr_free(s);
    return prod;
}

static Expr* canon_gamma(const Expr* atom, const Expr* root, const AssumeCtx* ctx) {
    (void)ctx;
    if (root_gamma_count(root) < 2) return NULL;
    Expr* r = gamma_multiply(atom, root);
    if (r) return r;
    return gamma_reflect(atom, root);
}

/* ---- canonical generators for the ratio heads ------------------------------ */

static Expr* gamma_of(Expr* a) { return sf_call(SYM_Gamma, (Expr*[]){ a }, 1); }
static Expr* inv(Expr* a) { return sf_call(SYM_Power, (Expr*[]){ a, sf_int(-1) }, 2); }

/* Does e hold a Gamma-type head other than Factorial? */
static bool has_non_factorial_gamma_head(const Expr* e) {
    return sf_contains_head(e, SYM_Gamma) || sf_contains_head(e, SYM_Pochhammer) ||
           sf_contains_head(e, SYM_Beta) || sf_contains_head(e, SYM_Binomial);
}

/* Does the argument have a non-integer rational part (Factorial[n + 1/2])? */
static bool fractional_class(const Expr* mu) {
    Expr* s = NULL; int64_t num, den;
    bool frac = sf_split(mu, &s, &num, &den) && den > 1;
    if (s) expr_free(s);
    return frac;
}

/* Factorial-only expressions in the integer class belong to simp_factorial
 * (the FactorialRules seed), which groups factorials by integer offset with
 * or without assumptions and knows Factorial2.  The engine converts
 * Factorial when it must meet another Gamma-type head, or when the argument
 * lies in a fractional class (Factorial[n + 1/2]), which simp_factorial
 * does not group. */
static Expr* canon_factorial(const Expr* atom, const Expr* root, const AssumeCtx* ctx) {
    (void)ctx;
    if (root_gamma_count(root) < 2) return NULL;
    if (!has_non_factorial_gamma_head(root) && !fractional_class(atom->data.function.args[0]))
        return NULL;
    return gamma_of(sf_plus(sf_cp(atom->data.function.args[0]), sf_int(1)));
}

static Expr* canon_pochhammer(const Expr* atom, const Expr* root, const AssumeCtx* ctx) {
    (void)ctx;
    if (root_gamma_count(root) < 2) return NULL;
    const Expr* a = atom->data.function.args[0];
    const Expr* n = atom->data.function.args[1];
    return sf_call(SYM_Times, (Expr*[]){ gamma_of(sf_plus(sf_cp(a), sf_cp(n))),
                                         inv(gamma_of(sf_cp(a))) }, 2);
}

static Expr* canon_beta(const Expr* atom, const Expr* root, const AssumeCtx* ctx) {
    (void)ctx;
    if (root_gamma_count(root) < 2) return NULL;
    const Expr* a = atom->data.function.args[0];
    const Expr* b = atom->data.function.args[1];
    return sf_call(SYM_Times, (Expr*[]){ gamma_of(sf_cp(a)), gamma_of(sf_cp(b)),
                                         inv(gamma_of(sf_plus(sf_cp(a), sf_cp(b)))) }, 3);
}

static Expr* canon_binomial(const Expr* atom, const Expr* root, const AssumeCtx* ctx) {
    (void)ctx;
    if (root_gamma_count(root) < 2) return NULL;
    const Expr* n = atom->data.function.args[0];
    const Expr* k = atom->data.function.args[1];
    Expr* nk1 = sf_mk1("Expand", sf_plus(sf_minus(sf_cp(n), sf_cp(k)), sf_int(1)));
    return sf_call(SYM_Times, (Expr*[]){ gamma_of(sf_plus(sf_cp(n), sf_int(1))),
                                         inv(gamma_of(sf_plus(sf_cp(k), sf_int(1)))),
                                         inv(gamma_of(nk1)) }, 3);
}

/* ---- refold passes ------------------------------------------------------- */

/* Factors of a Times (or the single factor), borrowed. */
static size_t factors_of(const Expr* e, const Expr*** out) {
    if (e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL &&
        e->data.function.head->data.symbol.name == SYM_Times) {
        *out = (const Expr**)e->data.function.args;
        return e->data.function.arg_count;
    }
    *out = NULL;
    return 0;
}

static bool is_gamma1(const Expr* e) { return sf_is_call(e, SYM_Gamma, 1); }
static bool is_inv_gamma1(const Expr* e) {
    return sf_is_call(e, SYM_Power, 2) && is_gamma1(e->data.function.args[0]) &&
           e->data.function.args[1]->type == EXPR_INTEGER &&
           e->data.function.args[1]->data.integer == -1;
}

static bool expands_to_zero(Expr* d) {
    Expr* x = sf_mk1("Expand", d);
    bool z = x->type == EXPR_INTEGER && x->data.integer == 0;
    expr_free(x);
    return z;
}

/* In one Times, fold Gamma ratios into `head` (Pochhammer / Beta / Binomial).
 * Returns an owned replacement Times or NULL. */
static Expr* refold_times(const Expr* t, const char* head) {
    const Expr** f;
    size_t n = factors_of(t, &f);
    if (n < 2) return NULL;
    bool* used = calloc(n, sizeof(bool));
    Expr* folded = NULL;
    for (size_t i = 0; i < n && !folded; i++) {
        if (!is_gamma1(f[i])) continue;
        const Expr* u = f[i]->data.function.args[0];
        if (head == SYM_Pochhammer) {
            for (size_t j = 0; j < n && !folded; j++) {
                if (!is_inv_gamma1(f[j])) continue;
                const Expr* v = f[j]->data.function.args[0]->data.function.args[0];
                Expr* d = sf_mk1("Expand", sf_minus(sf_cp(u), sf_cp(v)));
                if (d->type == EXPR_INTEGER) { expr_free(d); continue; } /* evaluator's job */
                folded = sf_call(SYM_Pochhammer, (Expr*[]){ sf_cp(v), d }, 2);
                used[i] = used[j] = true;
            }
        } else if (head == SYM_Beta) {
            for (size_t j = i + 1; j < n && !folded; j++) {
                if (!is_gamma1(f[j])) continue;
                const Expr* w = f[j]->data.function.args[0];
                for (size_t k = 0; k < n && !folded; k++) {
                    if (!is_inv_gamma1(f[k])) continue;
                    const Expr* c = f[k]->data.function.args[0]->data.function.args[0];
                    if (expands_to_zero(sf_minus(sf_plus(sf_cp(u), sf_cp(w)), sf_cp(c)))) {
                        folded = sf_call(SYM_Beta, (Expr*[]){ sf_cp(u), sf_cp(w) }, 2);
                        used[i] = used[j] = used[k] = true;
                    }
                }
            }
        } else if (head == SYM_Binomial) {
            /* Gamma[N] / (Gamma[K] Gamma[M]) with N - 1 = (K - 1) + (M - 1) */
            for (size_t j = 0; j < n && !folded; j++) {
                if (!is_inv_gamma1(f[j])) continue;
                const Expr* K = f[j]->data.function.args[0]->data.function.args[0];
                for (size_t k = j + 1; k < n && !folded; k++) {
                    if (!is_inv_gamma1(f[k])) continue;
                    const Expr* Mm = f[k]->data.function.args[0]->data.function.args[0];
                    if (expands_to_zero(sf_minus(sf_plus(sf_cp(u), sf_int(1)),
                                                 sf_plus(sf_cp(K), sf_cp(Mm))))) {
                        folded = sf_call(SYM_Binomial,
                                         (Expr*[]){ sf_minus(sf_cp(u), sf_int(1)),
                                                    sf_minus(sf_cp(K), sf_int(1)) }, 2);
                        used[i] = used[j] = used[k] = true;
                    }
                }
            }
        }
    }
    Expr* out = NULL;
    if (folded) {
        Expr** a = malloc(sizeof(Expr*) * (n + 1));
        size_t m = 0;
        a[m++] = folded;
        for (size_t i = 0; i < n; i++) if (!used[i]) a[m++] = sf_cp(f[i]);
        out = sf_call(SYM_Times, a, m);
        free(a);
    }
    free(used);
    return out;
}

static Expr* refold_walk(const Expr* e, const char* head, bool* changed) {
    if (e->type != EXPR_FUNCTION) return sf_cp(e);
    size_t n = e->data.function.arg_count;
    Expr** a = malloc(sizeof(Expr*) * (n ? n : 1));
    for (size_t i = 0; i < n; i++) a[i] = refold_walk(e->data.function.args[i], head, changed);
    Expr* r = eval_and_free(expr_new_function(sf_cp(e->data.function.head), a, n));
    free(a);
    /* repeat: a Times may hold several foldable ratios */
    for (int guard = 0; guard < 8; guard++) {
        Expr* t = refold_times(r, head);
        if (!t) break;
        expr_free(r);
        r = eval_and_free(t);
        *changed = true;
    }
    return r;
}

/* Gamma[u] -> Factorial[u - 1] where that is no larger (u = c + s, c >= 1). */
static Expr* refold_factorial_walk(const Expr* e, bool* changed) {
    if (e->type != EXPR_FUNCTION) return sf_cp(e);
    size_t n = e->data.function.arg_count;
    Expr** a = malloc(sizeof(Expr*) * (n ? n : 1));
    for (size_t i = 0; i < n; i++) a[i] = refold_factorial_walk(e->data.function.args[i], changed);
    Expr* r = eval_and_free(expr_new_function(sf_cp(e->data.function.head), a, n));
    free(a);
    if (is_gamma1(r)) {
        Expr* s = NULL; int64_t num, den;
        if (sf_split(r->data.function.args[0], &s, &num, &den) && den == 1 && num >= 1 && s) {
            Expr* f = sf_mk1(SYM_Factorial, sf_minus(sf_cp(r->data.function.args[0]), sf_int(1)));
            expr_free(r);
            r = f;
            *changed = true;
        }
        if (s) expr_free(s);
    }
    return r;
}

static Expr* refold_pass(const Expr* e, const Expr* input, const char* head) {
    if (!sf_contains_head(input, head) || !sf_contains_head(e, SYM_Gamma)) return NULL;
    bool changed = false;
    Expr* r = head == SYM_Factorial ? refold_factorial_walk(e, &changed)
                                    : refold_walk(e, head, &changed);
    if (!changed) { expr_free(r); return NULL; }
    return r;
}

static Expr* pass_refold_pochhammer(const Expr* e, const Expr* in, const AssumeCtx* c) {
    (void)c; return refold_pass(e, in, SYM_Pochhammer);
}
static Expr* pass_refold_beta(const Expr* e, const Expr* in, const AssumeCtx* c) {
    (void)c; return refold_pass(e, in, SYM_Beta);
}
static Expr* pass_refold_binomial(const Expr* e, const Expr* in, const AssumeCtx* c) {
    (void)c; return refold_pass(e, in, SYM_Binomial);
}
static Expr* pass_refold_factorial(const Expr* e, const Expr* in, const AssumeCtx* c) {
    (void)c; return refold_pass(e, in, SYM_Factorial);
}

/* ---- module table ------------------------------------------------------------ */

static const SfDesc GAMMA_DESCS[] = {
    { &SYM_Gamma,      1, 1, 1, { 0 }, up_gamma,     bulk_gamma, canon_gamma,      NULL, NULL },
    { &SYM_Gamma,      2, 1, 1, { 0 }, up_gamma_inc, NULL,       NULL,             NULL, NULL },
    { &SYM_LogGamma,   1, 1, 1, { 0 }, up_loggamma,  NULL,       NULL,             NULL, NULL },
    { &SYM_BarnesG,    1, 1, 1, { 0 }, up_barnesg,   NULL,       NULL,             NULL, NULL },
    { &SYM_Factorial,  1, 0, 1, { 0 }, NULL,         NULL,       canon_factorial,  NULL, NULL },
    { &SYM_Pochhammer, 2, 0, 1, { 0 }, NULL,         NULL,       canon_pochhammer, NULL, NULL },
    { &SYM_Beta,       2, 0, 1, { 0 }, NULL,         NULL,       canon_beta,       NULL, NULL },
    { &SYM_Binomial,   2, 0, 1, { 0 }, NULL,         NULL,       canon_binomial,   NULL, NULL },
};

static const SfPass GAMMA_PASSES[] = {
    pass_refold_pochhammer, pass_refold_beta, pass_refold_binomial, pass_refold_factorial,
};

static const char** GAMMA_PASS_HEADS[] = { &SYM_Gamma };

const SfModule SF_MODULE_GAMMA = {
    GAMMA_DESCS, sizeof GAMMA_DESCS / sizeof GAMMA_DESCS[0],
    GAMMA_PASSES, sizeof GAMMA_PASSES / sizeof GAMMA_PASSES[0],
    GAMMA_PASS_HEADS, sizeof GAMMA_PASS_HEADS / sizeof GAMMA_PASS_HEADS[0],
    NULL, NULL,
};
