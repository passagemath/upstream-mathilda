/*
 * simp_sf_airy.c -- special-function Simplify descriptors: Airy functions,
 * StruveH, Whittaker M and W, Fibonacci and LucasL (numbers and
 * polynomials).
 *
 * Derivations in docs/design/simp_special_functions.md (section 9); each
 * relation checked numerically at non-integer and complex parameters, and
 * case by case by tests/test_simp_sf_airy.c.
 *
 * Airy (pass): the Wronskian  Ai Bi' - Ai' Bi = 1/Pi.  Ai and Bi solve
 *   w'' = z w, which has no w' term, so their Wronskian is constant (Abel);
 *   at z = 0, Ai(0) = 3^(-2/3)/Gamma(2/3), Ai'(0) = -3^(-1/3)/Gamma(1/3),
 *   Bi = Sqrt[3] Ai, Bi' = -Sqrt[3] Ai' give 2 Sqrt[3] 3^-1 / (Gamma(1/3)
 *   Gamma(2/3)) = 1/Pi.  Every expression holding AiryAiPrime[z] and
 *   AiryBi[z] (same z) is reduced modulo the relation, eliminating the
 *   monomial Ai' Bi -- a canonical form for the four kernels.  (Ai'' = z Ai
 *   is the derivative rule the evaluator already applies.)
 *
 * Recurrences (rank 2, basis (C[k], C[k+1], 1); up row gives C[k+2]):
 *   StruveH   H[v-1] + H[v+1] = (2v/z) H[v] + (z/2)^v/(Sqrt[Pi] Gamma[v + 3/2])
 *       (the series; DLMF 11.4.23) -- inhomogeneous:
 *       H[v+2] = -H[v] + 2(v+1)/z H[v+1] + (z/2)^(v+1)/(Sqrt[Pi] Gamma[v + 5/2])
 *   WhittakerW (in k)  W[k+1] + (2k - z) W[k] - (m - k + 1/2)(m + k - 1/2) W[k-1] = 0
 *   WhittakerM (in k)  (2m + 2k + 1) M[k+1] + 2 (z - 2k) M[k] - (2m - 2k + 1) M[k-1] = 0
 *       (contiguous relations of the confluent hypergeometric functions
 *        M = z^(m+1/2) e^(-z/2) 1F1(m - k + 1/2; 1 + 2m; z), U likewise;
 *        both checked numerically at complex k, m, z)
 *   Fibonacci[n] / LucasL[n]        C[n+2] = C[n+1] + C[n]
 *   Fibonacci[n, x] / LucasL[n, x]  C[n+2] = x C[n+1] + C[n]
 *       (Binet: F(v) = (phi^v - cos(Pi v) phi^-v)/Sqrt[5], and
 *        cos(Pi (v+2)) = cos(Pi v), so the recurrence holds at every v)
 *
 * Pair: Fibonacci/LucasL  F[n+1] L[n] - F[n] L[n+1] = 2 Cos[Pi n]
 *   (Binet for both; holds at non-integer n -- checked).
 */

#include <stdlib.h>
#include <string.h>
#include "expr.h"
#include "eval.h"
#include "sym_names.h"
#include "simp.h"
#include "simp_internal.h"
#include "simp_sf_engine.h"

static void three_term_affine(Expr** M, Expr* A, Expr* B, Expr* c) {
    M[0] = sf_int(0); M[1] = sf_int(1); M[2] = sf_int(0);
    M[3] = A;         M[4] = B;         M[5] = c;
}

static bool up_struve(Expr* const* a, int d, Expr** M) {
    if (d != 0) return false;
    const Expr* v = a[0]; const Expr* z = a[1];
    Expr* v1 = sf_plus(sf_cp(v), sf_int(1));
    /* (z/2)^(v+1) spelled 2^-(v+1) z^(v+1), as D writes it */
    Expr* inh = sf_div(sf_times(sf_pow(sf_int(2), sf_times(sf_int(-1), sf_cp(v1))), sf_pow(sf_cp(z), sf_cp(v1))),
                       sf_times(sf_pow(sf_sym(SYM_Pi), sf_rat(1, 2)),
                                sf_call(SYM_Gamma, (Expr*[]){ sf_plus(sf_cp(v), sf_rat(5, 2)) }, 1)));
    three_term_affine(M, sf_int(-1), sf_div(sf_times(sf_int(2), v1), sf_cp(z)), inh);
    return true;
}

static bool up_whittakerW(Expr* const* a, int d, Expr** M) {
    if (d != 0) return false;
    const Expr* k = a[0]; const Expr* m = a[1]; const Expr* z = a[2];
    /* W[k+2] = (m - k - 1/2)(m + k + 1/2) W[k] + (z - 2(k+1)) W[k+1] */
    Expr* A = sf_times(sf_minus(sf_minus(sf_cp(m), sf_cp(k)), sf_rat(1, 2)),
                       sf_plus(sf_plus(sf_cp(m), sf_cp(k)), sf_rat(1, 2)));
    Expr* B = sf_minus(sf_cp(z), sf_times(sf_int(2), sf_plus(sf_cp(k), sf_int(1))));
    three_term_affine(M, A, B, sf_int(0));
    return true;
}

static bool up_whittakerM(Expr* const* a, int d, Expr** M) {
    if (d != 0) return false;
    const Expr* k = a[0]; const Expr* m = a[1]; const Expr* z = a[2];
    /* (2m + 2k + 3) M[k+2] = (2m - 2k - 1) M[k] - 2 (z - 2k - 2) M[k+1] */
    Expr* den = sf_plus(sf_times(sf_int(2), sf_plus(sf_cp(m), sf_cp(k))), sf_int(3));
    Expr* A = sf_div(sf_minus(sf_times(sf_int(2), sf_minus(sf_cp(m), sf_cp(k))), sf_int(1)), sf_cp(den));
    Expr* B = sf_div(sf_times(sf_int(-2), sf_minus(sf_minus(sf_cp(z), sf_times(sf_int(2), sf_cp(k))), sf_int(2))), den);
    three_term_affine(M, A, B, sf_int(0));
    return true;
}

static bool up_fib(Expr* const* a, int d, Expr** M) {
    (void)a;
    if (d != 0) return false;
    three_term_affine(M, sf_int(1), sf_int(1), sf_int(0));
    return true;
}

static bool up_fib_poly(Expr* const* a, int d, Expr** M) {
    if (d != 0) return false;
    three_term_affine(M, sf_int(1), sf_cp(a[1]), sf_int(0));
    return true;
}

/* F[n+1] L[n] - F[n] L[n+1] - 2 Cos[Pi n] */
static Expr* rel_FL(Expr* const* a, const Expr* c0, const Expr* c1, const Expr* d0, const Expr* d1) {
    Expr* w = sf_times(sf_int(2), sf_mk1(SYM_Cos, sf_times(sf_sym(SYM_Pi), sf_cp(a[0]))));
    return sf_minus(sf_minus(sf_times(sf_cp(c1), sf_cp(d0)), sf_times(sf_cp(c0), sf_cp(d1))), w);
}

/* ---- Airy Wronskian pass -------------------------------------------------- */

typedef struct { const Expr* z[32]; size_t n; } ZList;

static void collect_airy_args(const Expr* e, ZList* out) {
    if (e->type != EXPR_FUNCTION) return;
    if (sf_is_call(e, SYM_AiryAiPrime, 1) && out->n < 32) {
        bool dup = false;
        for (size_t i = 0; i < out->n && !dup; i++) dup = expr_eq(out->z[i], e->data.function.args[0]);
        if (!dup) out->z[out->n++] = e->data.function.args[0];
    }
    for (size_t i = 0; i < e->data.function.arg_count; i++) collect_airy_args(e->data.function.args[i], out);
}

static bool contains_expr(const Expr* e, const Expr* sub) {
    if (expr_eq(e, sub)) return true;
    if (e->type != EXPR_FUNCTION) return false;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (contains_expr(e->data.function.args[i], sub)) return true;
    return false;
}

static Expr* pass_airy_wronskian(const Expr* e, const Expr* input, const AssumeCtx* ctx) {
    (void)input; (void)ctx;
    ZList zs = { { NULL }, 0 };
    collect_airy_args(e, &zs);
    Expr* cur = NULL;
    for (size_t i = 0; i < zs.n; i++) {
        const Expr* z = zs.z[i];
        Expr* ai  = sf_call(SYM_AiryAi,      (Expr*[]){ sf_cp(z) }, 1);
        Expr* aip = sf_call(SYM_AiryAiPrime, (Expr*[]){ sf_cp(z) }, 1);
        Expr* bi  = sf_call(SYM_AiryBi,      (Expr*[]){ sf_cp(z) }, 1);
        Expr* bip = sf_call(SYM_AiryBiPrime, (Expr*[]){ sf_cp(z) }, 1);
        const Expr* base = cur ? cur : e;
        if (contains_expr(base, bi)) {
            /* Ai Bi' - Ai' Bi - 1/Pi; vars lex Ai' > Bi > Ai > Bi' eliminate Ai' Bi */
            Expr* R = sf_minus(sf_minus(sf_times(sf_cp(ai), sf_cp(bip)), sf_times(sf_cp(aip), sf_cp(bi))),
                               sf_pow(sf_sym(SYM_Pi), sf_int(-1)));
            Expr* vars[4] = { aip, bi, ai, bip };
            Expr* red = sf_reduce_modulo(base, R, vars, 4);
            expr_free(R);
            if (red) { if (cur) expr_free(cur); cur = red; }
        }
        expr_free(ai); expr_free(aip); expr_free(bi); expr_free(bip);
    }
    return cur;
}

/* ---- module table ------------------------------------------------------------ */

static const SfDesc AIRY_DESCS[] = {
    { &SYM_StruveH,     2, 1, 2, { 0 }, up_struve,     NULL, NULL, NULL,       NULL },
    { &SYM_WhittakerW,  3, 1, 2, { 0 }, up_whittakerW, NULL, NULL, NULL,       NULL },
    { &SYM_WhittakerM,  3, 1, 2, { 0 }, up_whittakerM, NULL, NULL, NULL,       NULL },
    { &SYM_Fibonacci,   1, 1, 2, { 0 }, up_fib,        NULL, NULL, &SYM_LucasL, rel_FL },
    { &SYM_LucasL,      1, 1, 2, { 0 }, up_fib,        NULL, NULL, NULL,       NULL },
    { &SYM_Fibonacci,   2, 1, 2, { 0 }, up_fib_poly,   NULL, NULL, NULL,       NULL },
    { &SYM_LucasL,      2, 1, 2, { 0 }, up_fib_poly,   NULL, NULL, NULL,       NULL },
    /* registered so the kernel count sees them; the Wronskian is a pass */
    { &SYM_AiryAi,      1, 0, 1, { 0 }, NULL,          NULL, NULL, NULL,       NULL },
    { &SYM_AiryAiPrime, 1, 0, 1, { 0 }, NULL,          NULL, NULL, NULL,       NULL },
    { &SYM_AiryBi,      1, 0, 1, { 0 }, NULL,          NULL, NULL, NULL,       NULL },
    { &SYM_AiryBiPrime, 1, 0, 1, { 0 }, NULL,          NULL, NULL, NULL,       NULL },
};

static const SfPass AIRY_PASSES[] = { pass_airy_wronskian };
static const char** AIRY_PASS_HEADS[] = { &SYM_AiryAiPrime };

const SfModule SF_MODULE_AIRY = {
    AIRY_DESCS, sizeof AIRY_DESCS / sizeof AIRY_DESCS[0],
    AIRY_PASSES, sizeof AIRY_PASSES / sizeof AIRY_PASSES[0],
    AIRY_PASS_HEADS, sizeof AIRY_PASS_HEADS / sizeof AIRY_PASS_HEADS[0],
    NULL, NULL,
};
