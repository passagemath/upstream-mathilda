/*
 * simp_sf_engine.h -- table-driven Simplify engine for special functions.
 *
 * A special function C[a1, ..., an] whose "shift" arguments differ by
 * integers belongs to a FAMILY: a lattice of values tied together by exact
 * recurrences.  The engine rewrites every member of a family onto one basis
 * -- {C[b]} for a rank-1 family (Gamma: C[b+1] = b C[b]), {C[b], C[b+e0]} for
 * a rank-2 family (Bessel, Legendre, 2F1: three-term recurrences) -- after
 * which Together/Expand decides whether the expression collapses.
 *
 * Every recurrence is stated as an AFFINE SHIFT MATRIX.  With the basis
 * vector at lattice point p
 *
 *     v(p) = (C[p], 1)                 rank 1
 *     v(p) = (C[p], C[p+e0], 1)        rank 2
 *
 * the descriptor's `up` builder supplies M(p, d) with v(p + e_d) = M(p, d) v(p).
 * The trailing 1 carries inhomogeneous terms (PolyGamma, HurwitzZeta,
 * ExpIntegralE, Struve).  Down steps use the inverse of M(p - e_d, d), which
 * the engine computes, so a descriptor only ever states the raising relation.
 *
 * Optional hooks:
 *   canon      -- rewrite one atom onto its canonical generator (Erfc -> 1-Erf,
 *                 Beta -> Gamma ratio, BesselK[-nu] -> BesselK[nu], Gamma
 *                 reflection when the reflected family is present, ...).
 *   bulk       -- rank-1 homogeneous closed form of C[p + k e_d] / C[p]
 *                 (Gamma: Pochhammer[p, k]) so large offsets cost O(1).
 *   partner +
 *   pair_relation -- a rank-2 family D sharing C's recurrence: after both
 *                 are on the same basis b, the expression is reduced modulo
 *                 their Casoratian relation (C[b+e0] D[b] - C[b] D[b+e0] =
 *                 w(b) up to sign), which makes the four generators
 *                 canonical (Bessel J/Y and I/K Wronskians).
 *
 * Extra passes (SfPass) are whole-expression candidates (Airy Wronskian,
 * Legendre duplication, PolyLog reflection, ...), each kept only when it
 * strictly lowers the score.
 *
 * Contract of simp_special_functions (mirrors simp_radical_rational): input
 * borrowed; returns a fresh Expr only when the final form scores STRICTLY
 * below the input (after a further simp_dispatch), else NULL.
 *
 * Every identity a descriptor encodes is derived in
 * docs/design/simp_special_functions.md and checked numerically by the
 * tests/test_simp_sf_*.c suites.
 */
#ifndef SIMP_SF_ENGINE_H
#define SIMP_SF_ENGINE_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include "expr.h"
#include "simp.h"

#define SF_MAX_SHIFT 3
#define SF_MAX_RANK  2

typedef struct SfDesc {
    const char** head;                 /* &SYM_X (interned at init)          */
    uint8_t      arity;                /* exact argument count               */
    uint8_t      nshift;               /* lattice coordinates (0 = none)     */
    uint8_t      rank;                 /* 1 or 2                             */
    uint8_t      shift_arg[SF_MAX_SHIFT]; /* arg index per coordinate; [0] is
                                           the basis direction e0          */
    /* Fill the (rank+1)x(rank+1) row-major matrix M (owned, unevaluated Exprs
     * allowed; the engine evaluates) for the up step from the point whose
     * arguments are `a` in direction d.  The last row is (0,...,0,1) and is
     * filled by the engine: the builder writes rows 0..rank-1 only.  Return
     * false to decline (the family is then left alone). */
    bool  (*up)(Expr* const* a, int d, Expr** M);
    /* Optional: C[p + k e_d] / C[p] for a rank-1 homogeneous family, k != 0.
     * Owned, unevaluated allowed.  NULL to fall back to stepping. */
    Expr* (*bulk)(Expr* const* a, int d, int64_t k);
    /* Optional canonicaliser: an owned replacement for `atom` (a call of this
     * head), or NULL to keep it.  `root` is the whole expression, so a hook
     * can make its choice depend on which other atoms are present. */
    Expr* (*canon)(const Expr* atom, const Expr* root, const AssumeCtx* ctx);
    /* Optional rank-2 partner family D (same recurrence) and the quadratic
     * relation tying the pair together at basis point b (arguments `a` of
     * C[b]): an owned R with R == 0, polynomial in c0 = C[b], c1 = C[b+e0],
     * d0 = D[b], d1 = D[b+e0] whose lex-leading monomial (order c1 > d0 > c0
     * > d1) is c1*d0 -- e.g. the Bessel Wronskian
     * J[b+1] Y[b] - J[b] Y[b+1] - 2/(Pi z). */
    const char** partner;
    Expr* (*pair_relation)(Expr* const* a, const Expr* c0, const Expr* c1,
                           const Expr* d0, const Expr* d1);
} SfDesc;

/* A whole-expression candidate transformation on the current form `e`
 * (owned result or NULL).  `input` is the expression Simplify handed the
 * engine, so a pass can e.g. fold Gamma ratios back into Pochhammer only
 * when the input was written with Pochhammer. */
typedef Expr* (*SfPass)(const Expr* e, const Expr* input, const AssumeCtx* ctx);

/* Per-module tables (one per simp_sf_*.c file). */
typedef struct SfModule {
    const SfDesc*  descs;
    size_t         ndesc;
    const SfPass*  passes;
    size_t         npass;
    /* heads that make the module's passes worth running */
    const char***  pass_heads;
    size_t         npass_heads;
} SfModule;

/* ---- shared helpers for descriptor files ------------------------------- */

/* Build and evaluate head[a, b] (both owned). */
Expr* sf_mk2(const char* head, Expr* a, Expr* b);
/* Build and evaluate head[a] (owned). */
Expr* sf_mk1(const char* head, Expr* a);
/* Evaluated Times[a, b] / Plus[a, b] / Power[a, b] (owned). */
Expr* sf_times(Expr* a, Expr* b);
Expr* sf_plus(Expr* a, Expr* b);
Expr* sf_pow(Expr* a, Expr* b);
/* a - b and a / b (owned). */
Expr* sf_minus(Expr* a, Expr* b);
Expr* sf_div(Expr* a, Expr* b);
Expr* sf_int(int64_t n);
/* Rational n/d, evaluated. */
Expr* sf_rat(int64_t n, int64_t d);
Expr* sf_sym(const char* name);
/* Borrowed-argument copy. */
Expr* sf_cp(const Expr* e);
/* Call head[args...] with `n` owned arguments, unevaluated. */
Expr* sf_call(const char* head, Expr** args, size_t n);

/* True when `e` is a call of interned `head` with `n` arguments. */
bool sf_is_call(const Expr* e, const char* head, size_t n);

/* Split mu = s + num/den (s owned or NULL for a numeric mu).  Plus with a
 * leading rational addend contributes it to num/den; every other shape is
 * wholly symbolic.  Returns false on a bignum constant. */
bool sf_split(const Expr* mu, Expr** s, int64_t* num, int64_t* den);

/* True when the symbolic expression s "reads negative": a negative number,
 * or a Times whose numeric coefficient is negative, or a Plus whose first
 * (canonically ordered) term reads negative.  Used to pick the canonical
 * sign of a reflected argument. */
bool sf_reads_negative(const Expr* s);

/* Is mu provably an integer under ctx?  mu = s + c with c an integer and s
 * absent or provably integer (prov_int alone does not see n + 1). */
bool sf_proved_integer(const AssumeCtx* ctx, const Expr* mu);

/* Class-preserving canonical sign for an even/odd-in-order reflection
 * mu <-> -mu - k (k an integer, k = 0 for BesselK and ChebyshevT): does
 * the flipped order have the smaller class?  A symbolic order flips when
 * its symbolic part reads negative; a numeric mu = num/den flips when
 * frac(mu) > 1/2, or frac(mu) in {0, 1/2} and mu < -k/2. */
bool sf_reflection_flips(const Expr* mu, int64_t k);

/* Reduce e modulo the polynomial relation R == 0 in the generator atoms
 * vars[0..nv-1] (lex order, vars[0] highest): numerator remainder over the
 * common denominator.  Owned result or NULL (e borrowed, R/vars borrowed). */
Expr* sf_reduce_modulo(const Expr* e, const Expr* R, Expr* const* vars, size_t nv);

/* Does `e` contain a call of `head` (cheap tree walk)? */
bool sf_contains_head(const Expr* e, const char* head);

/* ---- entry point -------------------------------------------------------- */

Expr* simp_special_functions(const Expr* input, const AssumeCtx* ctx,
                             const Expr* complexity_func);

/* Module tables, defined by the descriptor files. */
extern const SfModule SF_MODULE_BESSEL;
extern const SfModule SF_MODULE_GAMMA;
extern const SfModule SF_MODULE_ZETA;
extern const SfModule SF_MODULE_ERREXP;
extern const SfModule SF_MODULE_ORTHOPOLY;
extern const SfModule SF_MODULE_AIRY;

#endif
