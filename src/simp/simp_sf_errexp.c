/*
 * simp_sf_errexp.c -- special-function Simplify descriptors: error functions
 * (Erf, Erfc, Erfi, the generalized Erf[z0, z1], InverseErf, InverseErfc),
 * the exponential-integral family (ExpIntegralE, ExpIntegralEi,
 * LogIntegral, SinIntegral, CosIntegral, SinhIntegral, CoshIntegral) and
 * the Fresnel integrals.
 *
 * Derivations in docs/design/simp_special_functions.md (section 6); each
 * identity checked numerically by tests/test_simp_sf_errexp.c, and its
 * domain probed in all four quadrants and on the cuts (domain gates below
 * are the outcome of those probes, not assumptions).
 *
 * Canonical generators (unconditional identities):
 *   Erfc[z]      = 1 - Erf[z]                 (definition)
 *   Erf[z0, z1]  = Erf[z1] - Erf[z0]          (definition)
 *   Erfi[z]      = -I Erf[I z]                (substitute t -> i t in the
 *                                              integral), when an Erf/Erfc is
 *                                              present to meet it
 *   InverseErfc[x] = InverseErf[1 - x]        (Erfc = 1 - Erf), when an
 *                                              InverseErf is present
 *   Erf[InverseErf[x]] = x,  Erfc[InverseErfc[x]] = x   (inverse functions)
 *   Erf[InverseErfc[x]] = 1 - x               (Erf = 1 - Erfc)
 *   InverseErf[Erf[x]] = x                    (x real proved: Erf is not
 *                                              injective on C)
 *   LogIntegral[x] = ExpIntegralEi[Log[x]]    (li = Ei o log; holds on the
 *                                              whole plane, cuts included --
 *                                              checked), when an
 *                                              ExpIntegralEi is present
 *   FresnelS[I z] = -I FresnelS[z],  FresnelC[I z] = I FresnelC[z]
 *                                             (substitute t -> i t; the
 *                                              integrands are even in t^2)
 *   SinIntegral[I z] = I SinhIntegral[z],  SinhIntegral[I z] = I SinIntegral[z]
 *                                             (substitute t -> i t), when the
 *                                              partner head is present
 *
 * Gated (fail on an open set or a real half-line otherwise -- checked):
 *   CosIntegral[-x]  = CosIntegral[x] + I Pi     (x > 0 proved)
 *   CoshIntegral[-x] = CoshIntegral[x] + I Pi    (x > 0 proved)
 *   ExpIntegralE[n, z] = z^(n-1) Gamma[1 - n, z] (z > 0 proved), when an
 *                                                 incomplete Gamma of z is
 *                                                 present
 *
 * Recurrence (rank 1, affine, in the order n):
 *   ExpIntegralE[n+1, z] = (E^-z - z ExpIntegralE[n, z]) / n
 *       E_{n+1}(z) = Int_1^Inf e^(-z t) t^(-n-1) dt; by parts with
 *       dv = t^(-n-1) dt: [-e^(-z t) t^-n / n]_1^Inf - (z/n) Int e^(-z t) t^-n dt
 *       = e^-z/n - (z/n) E_n(z).  Holds for every n != 0 and every z by
 *       continuation (checked at real, negative and complex z).
 *
 * Refold pass: Erf[z] -> 1 - Erfc[z] for each z the input wrote with Erfc,
 * kept when it scores no worse (the vocabulary of the input).
 */

#include <stdlib.h>
#include <string.h>
#include "expr.h"
#include "eval.h"
#include "sym_names.h"
#include "simp.h"
#include "simp_internal.h"
#include "simp_sf_engine.h"

static Expr* I_unit(void) { return sf_sym(SYM_I); }

/* Is e = I w with w free of I?  The evaluator spells I x as
 * Times[Complex[0, 1], x], and the engine's normal form expands
 * I (1/3 + x) to I/3 + I x, so test the value: w = Expand[-I e] must be
 * free of Complex, and e itself must hold one.  Then *w is owned. */
static bool is_i_times(const Expr* e, Expr** w) {
    if (!sf_contains_head(e, SYM_Complex)) return false;
    Expr* q = sf_mk1("Expand", sf_times(sf_times(sf_int(-1), I_unit()), sf_cp(e)));
    if (sf_contains_head(q, SYM_Complex)) { expr_free(q); return false; }
    *w = q;
    return true;
}

static bool reads_negative_arg(const Expr* x) {
    Expr* s = NULL; int64_t num, den;
    if (!sf_split(x, &s, &num, &den)) return false;
    bool neg = s ? sf_reads_negative(s) : num < 0;
    if (s) expr_free(s);
    return neg;
}

/* ---- error functions ------------------------------------------------------ */

static Expr* canon_erfc(const Expr* atom, const Expr* root, const AssumeCtx* ctx) {
    (void)root; (void)ctx;
    const Expr* z = atom->data.function.args[0];
    if (sf_is_call(z, SYM_InverseErfc, 1)) return sf_cp(z->data.function.args[0]);
    return sf_minus(sf_int(1), sf_mk1(SYM_Erf, sf_cp(z)));
}

static Expr* canon_erf(const Expr* atom, const Expr* root, const AssumeCtx* ctx) {
    (void)root; (void)ctx;
    const Expr* z = atom->data.function.args[0];
    if (sf_is_call(z, SYM_InverseErf, 1)) return sf_cp(z->data.function.args[0]);
    /* Erf[InverseErfc[x]] = 1 - Erfc[InverseErfc[x]] = 1 - x */
    if (sf_is_call(z, SYM_InverseErfc, 1)) return sf_minus(sf_int(1), sf_cp(z->data.function.args[0]));
    return NULL;
}

static Expr* canon_erf2(const Expr* atom, const Expr* root, const AssumeCtx* ctx) {
    (void)root; (void)ctx;
    return sf_minus(sf_mk1(SYM_Erf, sf_cp(atom->data.function.args[1])),
                    sf_mk1(SYM_Erf, sf_cp(atom->data.function.args[0])));
}

static Expr* canon_erfi(const Expr* atom, const Expr* root, const AssumeCtx* ctx) {
    (void)ctx;
    if (!sf_contains_head(root, SYM_Erf) && !sf_contains_head(root, SYM_Erfc)) return NULL;
    return sf_times(sf_times(sf_int(-1), I_unit()),
                    sf_mk1(SYM_Erf, sf_times(I_unit(), sf_cp(atom->data.function.args[0]))));
}

static Expr* canon_inverseerfc(const Expr* atom, const Expr* root, const AssumeCtx* ctx) {
    (void)ctx;
    if (!sf_contains_head(root, SYM_InverseErf)) return NULL;
    return sf_mk1(SYM_InverseErf, sf_minus(sf_int(1), sf_cp(atom->data.function.args[0])));
}

static Expr* canon_inverseerf(const Expr* atom, const Expr* root, const AssumeCtx* ctx) {
    (void)root;
    const Expr* x = atom->data.function.args[0];
    if (sf_is_call(x, SYM_Erf, 1) && ctx && prov_re(ctx, x->data.function.args[0]))
        return sf_cp(x->data.function.args[0]);
    return NULL;
}

/* Erf[z] -> 1 - Erfc[z] for every z the input wrote with Erfc. */
typedef struct { const Expr* a[64]; size_t n; } ArgList;

static void collect_args(const Expr* e, const char* head, ArgList* out) {
    if (e->type != EXPR_FUNCTION) return;
    if (sf_is_call(e, head, 1) && out->n < sizeof out->a / sizeof out->a[0]) {
        bool dup = false;
        for (size_t i = 0; i < out->n && !dup; i++) dup = expr_eq(out->a[i], e->data.function.args[0]);
        if (!dup) out->a[out->n++] = e->data.function.args[0];
    }
    for (size_t i = 0; i < e->data.function.arg_count; i++) collect_args(e->data.function.args[i], head, out);
}

static Expr* subst_erf(const Expr* e, const ArgList* zs, bool* changed) {
    if (e->type != EXPR_FUNCTION) return sf_cp(e);
    if (sf_is_call(e, SYM_Erf, 1))
        for (size_t i = 0; i < zs->n; i++)
            if (expr_eq(e->data.function.args[0], zs->a[i])) {
                *changed = true;
                return sf_minus(sf_int(1), sf_call(SYM_Erfc, (Expr*[]){ sf_cp(zs->a[i]) }, 1));
            }
    size_t n = e->data.function.arg_count;
    Expr** a = malloc(sizeof(Expr*) * (n ? n : 1));
    for (size_t i = 0; i < n; i++) a[i] = subst_erf(e->data.function.args[i], zs, changed);
    Expr* r = expr_new_function(sf_cp(e->data.function.head), a, n);
    free(a);
    return r;
}

static Expr* pass_refold_erfc(const Expr* e, const Expr* input, const AssumeCtx* ctx) {
    (void)ctx;
    if (!sf_contains_head(input, SYM_Erfc) || !sf_contains_head(e, SYM_Erf)) return NULL;
    ArgList zs = { { NULL }, 0 };
    collect_args(input, SYM_Erfc, &zs);
    bool changed = false;
    Expr* r = eval_and_free(subst_erf(e, &zs, &changed));
    if (!changed) { expr_free(r); return NULL; }
    return r;
}

/* ---- exponential integrals ------------------------------------------------ */

static bool up_expintegrale(Expr* const* a, int d, Expr** M) {
    if (d != 0) return false;
    if (a[0]->type == EXPR_INTEGER && a[0]->data.integer == 0) return false;   /* n = 0: E_1 = (e^-z - 0 E_0)/0 */
    Expr* ninv = sf_pow(sf_cp(a[0]), sf_int(-1));
    M[0] = sf_times(sf_times(sf_int(-1), sf_cp(a[1])), sf_cp(ninv));
    M[1] = sf_times(sf_pow(sf_sym(SYM_E), sf_times(sf_int(-1), sf_cp(a[1]))), ninv);
    return true;
}

/* Does e hold an incomplete Gamma[_, z]? */
static bool has_gamma_inc_of(const Expr* e, const Expr* z) {
    if (e->type != EXPR_FUNCTION) return false;
    if (sf_is_call(e, SYM_Gamma, 2) && expr_eq(e->data.function.args[1], z)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (has_gamma_inc_of(e->data.function.args[i], z)) return true;
    return false;
}

static Expr* canon_expintegrale(const Expr* atom, const Expr* root, const AssumeCtx* ctx) {
    const Expr* n = atom->data.function.args[0];
    const Expr* z = atom->data.function.args[1];
    if (!ctx || !has_gamma_inc_of(root, z) || !prov_pos(ctx, z)) return NULL;
    return sf_times(sf_pow(sf_cp(z), sf_minus(sf_cp(n), sf_int(1))),
                    sf_call(SYM_Gamma, (Expr*[]){ sf_minus(sf_int(1), sf_cp(n)), sf_cp(z) }, 2));
}

/* Does e contain exactly `sub`? */
static bool contains_expr(const Expr* e, const Expr* sub) {
    if (expr_eq(e, sub)) return true;
    if (e->type != EXPR_FUNCTION) return false;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (contains_expr(e->data.function.args[i], sub)) return true;
    return false;
}

/* LogIntegral[x] -> ExpIntegralEi[Log[x]] only to meet that exact atom. */
static Expr* canon_logintegral(const Expr* atom, const Expr* root, const AssumeCtx* ctx) {
    (void)ctx;
    Expr* ei = sf_mk1(SYM_ExpIntegralEi, sf_mk1(SYM_Log, sf_cp(atom->data.function.args[0])));
    if (!contains_expr(root, ei)) { expr_free(ei); return NULL; }
    return ei;
}

/* Ci[-x] = Ci[x] + I Pi and Chi[-x] = Chi[x] + I Pi, for x > 0 proved. */
static Expr* reflect_plus_ipi(const Expr* atom, const AssumeCtx* ctx) {
    const Expr* a = atom->data.function.args[0];
    if (!ctx || !reads_negative_arg(a)) return NULL;
    Expr* x = sf_mk1("Expand", sf_times(sf_int(-1), sf_cp(a)));
    if (!prov_pos(ctx, x)) { expr_free(x); return NULL; }
    const char* h = atom->data.function.head->data.symbol.name;
    return sf_plus(sf_call(h, (Expr*[]){ x }, 1), sf_times(I_unit(), sf_sym(SYM_Pi)));
}

static Expr* canon_cosintegral(const Expr* atom, const Expr* root, const AssumeCtx* ctx) {
    (void)root; return reflect_plus_ipi(atom, ctx);
}

/* Si[I z] = I Shi[z] / Shi[I z] = I Si[z], when the partner head is present. */
static Expr* canon_sin_sinh(const Expr* atom, const Expr* root, const char* partner) {
    if (!sf_contains_head(root, partner)) return NULL;
    Expr* w = NULL;
    if (!is_i_times(atom->data.function.args[0], &w)) return NULL;
    return sf_times(I_unit(), sf_call(partner, (Expr*[]){ w }, 1));
}

static Expr* canon_sinintegral(const Expr* atom, const Expr* root, const AssumeCtx* ctx) {
    (void)ctx; return canon_sin_sinh(atom, root, SYM_SinhIntegral);
}
static Expr* canon_sinhintegral(const Expr* atom, const Expr* root, const AssumeCtx* ctx) {
    (void)ctx; return canon_sin_sinh(atom, root, SYM_SinIntegral);
}

/* ---- Fresnel ---------------------------------------------------------------- */

static Expr* canon_fresnel(const Expr* atom, int sign) {
    Expr* w = NULL;
    if (!is_i_times(atom->data.function.args[0], &w)) return NULL;
    const char* h = atom->data.function.head->data.symbol.name;
    return sf_times(sf_times(sf_int(sign), I_unit()), sf_call(h, (Expr*[]){ w }, 1));
}

static Expr* canon_fresnels(const Expr* atom, const Expr* root, const AssumeCtx* ctx) {
    (void)root; (void)ctx; return canon_fresnel(atom, -1);
}
static Expr* canon_fresnelc(const Expr* atom, const Expr* root, const AssumeCtx* ctx) {
    (void)root; (void)ctx; return canon_fresnel(atom, 1);
}

/* ---- module table -------------------------------------------------------- */

static const SfDesc ERREXP_DESCS[] = {
    { &SYM_Erfc,          1, 0, 1, { 0 }, NULL,            NULL, canon_erfc,         NULL, NULL },
    { &SYM_Erf,           1, 0, 1, { 0 }, NULL,            NULL, canon_erf,          NULL, NULL },
    { &SYM_Erf,           2, 0, 1, { 0 }, NULL,            NULL, canon_erf2,         NULL, NULL },
    { &SYM_Erfi,          1, 0, 1, { 0 }, NULL,            NULL, canon_erfi,         NULL, NULL },
    { &SYM_InverseErf,    1, 0, 1, { 0 }, NULL,            NULL, canon_inverseerf,   NULL, NULL },
    { &SYM_InverseErfc,   1, 0, 1, { 0 }, NULL,            NULL, canon_inverseerfc,  NULL, NULL },
    { &SYM_ExpIntegralE,  2, 1, 1, { 0 }, up_expintegrale, NULL, canon_expintegrale, NULL, NULL },
    { &SYM_LogIntegral,   1, 0, 1, { 0 }, NULL,            NULL, canon_logintegral,  NULL, NULL },
    /* no relations of its own; registered so the kernel count sees it */
    { &SYM_ExpIntegralEi, 1, 0, 1, { 0 }, NULL,            NULL, NULL,               NULL, NULL },
    { &SYM_CosIntegral,   1, 0, 1, { 0 }, NULL,            NULL, canon_cosintegral,  NULL, NULL },
    { &SYM_CoshIntegral,  1, 0, 1, { 0 }, NULL,            NULL, canon_cosintegral,  NULL, NULL },
    { &SYM_SinIntegral,   1, 0, 1, { 0 }, NULL,            NULL, canon_sinintegral,  NULL, NULL },
    { &SYM_SinhIntegral,  1, 0, 1, { 0 }, NULL,            NULL, canon_sinhintegral, NULL, NULL },
    { &SYM_FresnelS,      1, 0, 1, { 0 }, NULL,            NULL, canon_fresnels,     NULL, NULL },
    { &SYM_FresnelC,      1, 0, 1, { 0 }, NULL,            NULL, canon_fresnelc,     NULL, NULL },
};

static const SfPass ERREXP_PASSES[] = { pass_refold_erfc };
static const char** ERREXP_PASS_HEADS[] = { &SYM_Erf };

const SfModule SF_MODULE_ERREXP = {
    ERREXP_DESCS, sizeof ERREXP_DESCS / sizeof ERREXP_DESCS[0],
    ERREXP_PASSES, sizeof ERREXP_PASSES / sizeof ERREXP_PASSES[0],
    ERREXP_PASS_HEADS, sizeof ERREXP_PASS_HEADS / sizeof ERREXP_PASS_HEADS[0],
    NULL, NULL,
};
