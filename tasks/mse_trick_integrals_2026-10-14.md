---
title: MSE "simple techniques" definite-integral campaign
status: COMPLETE — A+B (v0.404), C narrow stage (v0.405)
started: 2026-10-14
baseline: 17/23 CORRECT, 6 UNEVALUATED, 0 WRONG (v0.402, commit 30964d9)
now: corpus 15/15 CORRECT, 0 WRONG (all three missing general methods landed)
corpus: tests/MseTrickIntegralsCorpus.m; gate mse_trick_corpus_tests (WRONG==0, CORRECT>=15)
plan: ~/.claude/plans/we-should-look-at-enumerated-lark.md
git: A+B committed+tagged v0.404 (e548b0ea); C (v0.405) commit pending
---

## Outcome

- **B. Series term-by-term (v0.403)** — DONE. `src/calculus/integrate_series.{c,h}`.
  Closes #7 (−Catalan), #16 (π³/16), +5 more corpus cases, + bonus cyclotomic
  (1+x^3, 1+x^2+x^4). Correct-by-construction + NIntegrate backstop. 0 regressions.
- **A. General King reflection (v0.404)** — DONE. `src/calculus/integrate_king.{c,h}`.
  Closes #1 (π/4) and the whole 1/(1+Tan^s) family. Numeric-only, midpoint-constant
  + numeric screen + NIntegrate backstop.
- **C. Feynman on a present parameter (#13) — DONE (v0.405, narrow stage)**.
  `src/calculus/integrate_feynparam.{c,h}`, `Integrate\`FeynmanParameter`. Closes
  the Frullani-log family `∫₀¹ (Σ cᵢ x^{eᵢ})/Log[x] = Σ cᵢ Log[eᵢ+1]` (needs
  `Σcᵢ=0`, each `eᵢ+1>0`): `(x^a-1)/Log[x] → Log[a+1]`,
  `(x^a-x^b)/Log[x] → Log[(a+1)/(b+1)]`, `(2x^4-x-1)/Log[x] → Log[25/2]`.
  **Side-steps the two blockers below** by computing the inner x-power integral
  directly as `1/(e+1)` (no `Beta`/`Gamma` form → no reduction needed) and doing no
  long re-entrant `Integrate` (stays out of the crash territory). Correct-by-
  construction + NIntegrate backstop; declines `1/Log[x]` (divergent) and non-power
  integrands. Not "fully general" Feynman — the two items below remain true and
  would be needed for the unrestricted method:
  1. **Gamma-simplification gap** ([[project_gamma_ratio_simplification_gap]]):
     `Simplify`/`FunctionExpand` do not reduce `Beta[1+a,1]` / `Gamma[1+a]/Gamma[2+a]`
     to `1/(1+a)`, so a general integrate-back over the parameter cannot close.
  2. **Crash-prone territory** ([[project_parametric_integrate_stage_crash_gate]]):
     heavy parametric runs trip a PRE-EXISTING latent crash (see below).

## PRE-EXISTING latent crash found (NOT introduced here; worked around)

`Integrate[(x^a-1)/Log[x], {x,0,1}]` under Automatic SIGSEGVs after ~10 s. valgrind:
use-after-free / 1-slot args overflow in `optionvalue_inject_context`
(src/options_builtin.c:322-323) during a deep recursive rebuild, written past by
`evaluate` (eval.c:2796) — fired when the cumulative cascade runs long enough that
the async integrate-TimeConstraint abort lands mid-malloc (cf. memory
`project_timeconstrained_async_malloc_lock`). It is NOT reached on a clean v0.402
(bisected: both new stages off → clean). King's heavy re-entrant Simplify/N on the
parametric integrand was pushing the cascade over that threshold; **fixed** by the
King numeric-domain precondition (declines parametric integrands before any heavy
work). The underlying options/async-TC bug remains a separate latent issue worth a
dedicated fix.

# MSE tough-integrals campaign

Goal: fully general, algorithmic definite-integration methods closing the 6
UNEVALUATED cases from the technique battery. Invariant: ZERO wrong answers —
each stage correct-by-construction + NIntegrate backstop, else decline (NULL).

## Methods (landing order: B → A → C)

- [ ] **B. Series expand → term-by-term → resum** (`integrate_series.c`)
      closes #7 (−Catalan), #16 (π³/16). Alternating sibling of ratlog01.
  - [ ] Verify Sum[] dependency (DONE: Catalan, π³/32, π²/12, ¾ζ3, π²/48 all close)
  - [ ] integrate_series.{c,h}
  - [ ] wire into integrate.c (enum, method_from_string, include, cascade @1378, init, N-L exclusion)
  - [ ] docstring + attributes (Protected) + register Integrate`SeriesExpansion
  - [ ] corpus rows + gate green
  - [ ] docs/spec/builtins/calculus.md + changelog 2026-10-12.md
  - [ ] version bump + tag
- [ ] **A. General King reflection** (`integrate_king.c`) closes #1 (π/4)
- [ ] **C. Feynman on present parameter** (extend integrate_diffunderint.c) closes #11 (Log[a+1])

## Corpus (closed forms, NIntegrate-verified by hand)
Method B new: ∫₀¹ Log x/(1+x²)→−Catalan; Log²x/(1+x²)→π³/16; Log x/(1+x)→−π²/12;
  Log²x/(1+x)→3/2 ζ3; Log³x/(1+x)→−7π⁴/120; x Log x/(1+x²)→−π²/48;
  x² Log x/(1+x²)→Catalan−1 (poly-div path).
Method B regression (ratlog01 owns, must not perturb): Log x/(1−x)→−π²/6;
  Log²x/(1−x)→2ζ3; Log³x/(1−x)→−π⁴/15.
Method A: ∫₀^{π/2} 1/(1+Tan^√2)→π/4; ∫₀^π x Sin x/(1+Cos²x)→π²/4.
Method C: ∫₀¹ (xᵃ−1)/Log x→Log[a+1] (a>−1, inst a→2).
Negative controls: ∫₀¹ Log[1−x]/x (arg≠x→fall through); ∫₀¹ 1/(1+x²)→π/4 (k=0→decline).

## Follow-ons (NOT this campaign)
- #3 ∫₀¹ ln(1+x)/(1+x²)→π ln2/8 : substitution x=tanθ then King (composite).
- #5 ∫₀^{π/2} x cot x→(π/2)ln2 : definite integration-by-parts → #4.
