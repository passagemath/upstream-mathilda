---
title: Assumption-Dependent Definite-Integral Stress Campaign
date_started: 2026-10-13
status: in-progress
---

# Assumed definite-integral stress campaign (50 cases)

Battery: `assumedIntegralStressTest` — 50 definite integrals whose value depends on
`Assumptions` (convergence `a>0`, sign `a<0`, domain `n ∈ Integers`). Each case is
`Integrate[...] - Expected`.

## Harness (gold-standard: one process per case)

A single long process accumulating many `TimeConstrained` aborts **segfaults**
(EXIT=139) — a state-dependent TimeConstrained async fragility, NOT any single
integral (every case run alone exits 0). So the authoritative harness runs **one
Mathilda process per case** (`scratchpad/alldefs.m` + a bash loop), capturing exit
codes. Two metrics per case:

- **ASSUM** (authoritative): `Simplify[raw - expected, assumptions]` — "is the value
  correct under the stated assumptions?"
- **LITERAL** (battery as written): `Simplify[raw - expected]` — penalises
  branch-cut-sensitive FORM differences the assumptions would resolve (e.g.
  `Sqrt[Pi]/Sqrt[a]` vs `Sqrt[Pi/a]`). Secondary; even Mathematica would "fail" some
  of these literally.

## Results (2026-10-13, pre-fix baseline, v0.376)

| Metric | PASS | WRONG | UNEVAL | ABORTED |
|---|---|---|---|---|
| **ASSUM** (authoritative) | **35** | 5 | 6 | 4 |
| LITERAL | 25 | 14 | 7 | 4 |

**40/50 compute the correct value** (35 ASSUM-PASS + 5 "WRONG" that are correct but
Simplify-blocked). **No true wrong answers.** Only **10 cases have a real integration
gap.**

### Correct value, blocked only by a Simplify gap (ASSUM=WRONG) — 5
| # | Integrand | Raw (correct) | Blocking Simplify gap |
|---|---|---|---|
| 8 | `Exp[-a x]Cos[b x]` | `ConditionalExpression[1/(√(a²+b²)√(1+b²/a²)), …]` | ConditionalExpression discharge + nested radical |
| 9 | `Exp[-a x]Sin[b x]` | `ConditionalExpression[(a b √(1+b²/a²))/(a²+b²)^(3/2), …]` | same |
| 18 | `Sin[x]^n` [0,π/2] | `ConditionalExpression[Beta[(1+n)/2,1/2]/2, …]` | ConditionalExpression + Beta→Gamma |
| 26 | `x^a/(1+x²)` | `-π Csc[π(1+a)] Sin[πa/2]` | trig periodicity `Csc[π(1+a)]→−Csc[πa]` |
| 37 | `Exp[-a x]BesselJ[0,b x]` | `1/(a √(1+b²/a²))` | nested radical |

Four independent Simplify gaps, each confirmed un-closed by FullSimplify-under-assumptions:
1. **ConditionalExpression discharge** when assumptions entail the condition (8,9,18).
2. **Nested radical** `Sqrt[1+p/q] → Sqrt[(q+p)/q]` pull-out under positivity (8,9,37).
3. **Trig periodicity** `Csc/Sec/Sin/Cos[π(k+a)]`, integer k (26).
4. **`Beta[x,y] → Gamma`** expansion (FunctionExpand returns unevaluated) (18).

### Correct value, un-reduced FORM (ASSUM=PASS, LITERAL fails) — 10
10 (`√π/√a`), 12–17 (trig orthogonality: `2Sin[πn]/n` etc.), 31 (forced via Simplify),
48,49 (Fresnel `√π/(2√2√a)`). All fixed by a **driver post-simplify under Assumptions**
(reduces `Sin[πn]→0` for integer n, etc.) — Mathematica returns the reduced form.

### Genuine integration gaps (value NOT obtained) — 10
| # | Integrand | Status | Theme |
|---|---|---|---|
| 43 | `HeavisideTheta[a-x]` [0,∞) | UNEVAL | piecewise |
| 44 | `Sign[x-a]` [0,2a] | UNEVAL | piecewise |
| 45 | `Abs[x-a]` [0,2a] | UNEVAL | piecewise |
| 46 | `Min[x,a]` [0,2a] | UNEVAL | piecewise |
| 47 | `Max[x,a]` [0,2a] | UNEVAL | piecewise |
| 27 | `1/√(a²-x²)` [0,a] | UNEVAL | improper endpoint (ArcSin) |
| 20 | `Sin[a x]/x` [0,∞), a<0 | ABORTED | sign reflection (odd in a) |
| 32 | `Log[x]/(x²+a²)` [0,∞) | ABORTED | log-rational half-line |
| 38 | `Log[a²-2ab Cos x+b²]` [0,π], a>b | ABORTED | Poisson log |
| 39 | same, b>a | ABORTED | Poisson log |

## Plan (supersedes the pre-run new-handler-heavy plan; trig-orthogonality handler DROPPED)

**Theme 1 — Result cleanliness under Assumptions (Mathematica parity).**
- 1a. Driver post-simplify: `integrate_definite` returns `Simplify[r, assumptions]` when
  assumptions given (gated → zero blast radius on no-assumption calls). Fixes 12–17 +
  cleans 10,31,48,49.
- 1b. ConditionalExpression discharge under entailed assumptions (8,9,18).
- 1c. Nested-radical combination under positivity (8,9,37).
- 1d. Beta→Gamma expansion (18).
- 1e. Trig periodicity shift (26).

**Theme 2 — Genuine integration gaps.**
- 2a. Piecewise handler `integrate_piecewise.c` (43–47). CONFIRMED needed.
- 2b. Improper endpoint / ArcSin (27).
- 2c. Sign reflection `Sin[ax]/x`, a<0 (20).
- 2d. Log-rational half-line `Log[x]/(x²+a²)` (32).
- 2e. Poisson log on [0,π] (38,39).

**Theme 3 — Robustness.** TimeConstrained-accumulation segfault: characterise;
likely pre-existing TC async issue (memory: timeconstrained_async_malloc_lock).

## Soundness
Zero true wrong answers at baseline. Every fix must preserve that; re-run the full
battery (both metrics, per-case isolated) after each change — no regression allowed.

## Verification artifacts
- `scratchpad/alldefs.m` (reusable 50-case driver) → promote to
  `tests/AssumedDefiniteIntegralsCorpus.m` + runner (model `tests/test_crc_corpus.c`).
