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

## Progress (live)

| Version | Fix | Cases | Status |
|---|---|---|---|
| v0.377 | Piecewise handler (`integrate_piecewise.c`) | 43–47 | DONE, shipped |
| v0.378 | `Sin[r x]^k/x^m` negative frequency | 20 | DONE, shipped |
| v0.379 | Refine radical reduction (products/quotients) + assumption-aware Integrate result | 37 | DONE |
| v0.380 | Simplify trig periodicity shift `trig[Pi·n+r] → (-1)^n trig[r]` | 26 | DONE, shipped |
| v0.381 | Real `ArcSin` antiderivative for `1/Sqrt[c−x²]` symbolic c (prerequisite for 27) | 27 (indef) | DONE, shipped |

ASSUM pass-count: baseline 35/50 → **43/50** (no regression, 0 wrong answers).
Remaining 7 (8, 9, 18 WRONG; 27 UNEVAL; 32, 38, 39 ABORTED) — each needs COUPLED
delicate changes; sharpened diagnoses below.
- **27**: indefinite now gives `ArcSin[x/Sqrt[a²]]` (v0.381). The *definite* still
  declines because Newton–Leibniz runs its **pole placement** (`nl_collect_poles`:
  can't decide x=±a vs (0,a) without a>0 → undecidable → numeric scan fails on symbolic
  a → decline) AND its **boundary limit** (`nl_limit`: `Limit[ArcSin[x/√a²], x→a]`
  without a>0 → `ArcSin[a/√a²]`, unresolved) **without Assumptions**. Fix = thread
  `Assumptions` into `integrate_newton_leibniz_try_pv` → both. `Limit[…, a>0]` already
  gives Pi/2 and the driver-Refine already cleans `ArcSin[a/√a²]→Pi/2`, so the threading
  alone closes it — but it is a core-FTC change with broad regression surface.
Remaining 7 — each needs COUPLED delicate changes (no single-change closure):
- **8, 9**: the value is already clean (`a/(a²+b²)`); blocked by a spurious `b>0` the
  **OscillatoryPower / general Ramanujan Mellin strip** attaches (cos/sin even/odd in b →
  valid all real b). Pinned `DiffUnderInt` returns it unconditioned. Fix = widen that
  Mellin convergence strip (delicate: 51 ramanujan cases) OR make the osc/Mellin path
  decline the `e^{-ax} trig(bx)` family so Automatic falls through to DiffUnderInt.
- **18**: needs BOTH (a) ConditionalExpression discharge — `Refine[CE[v, (1+n)/2>0], n>-1]`
  does NOT currently collapse to `v` even though entailed; AND (b) `Beta↔Gamma`
  normalization (`FullSimplify` keeps the shorter `Beta` side; `FunctionExpand[Beta]` is
  unimplemented) so `Beta[(1+n)/2,1/2]/2` equals the expected Gamma form.
- **27** complex-log antiderivative not `ArcSin`; **32** log-rational complex-pole Mellin;
  **38, 39** Poisson recognizer + radicand factoring + `Log[Abs]` composite-Simplify gap.

## Remaining 9 — precise root causes & fix recipes

**Keystone (unblocks 37, and 38/39 after their recognizer): assumption-aware radical
reduction in `src/refine.c` `deep_positivity_walk`.** The block at `refine.c:479`
already does `Sqrt[p^2] -> ±p` via the CAD-backed `deep_sign`; extend it to a
**product/quotient of squared factors**: `Sqrt[Times[.., base^(2k), ..]] ->
|base|^k Sqrt[remainder]` (and `Together` a `Plus` radicand first to expose a
denominator, so `Sqrt[1+b^2/a^2] = Sqrt[(a^2+b^2)/a^2] -> Sqrt[a^2+b^2]/a`). Each
pulled `|base|` is `deep_sign`-resolved, so it stays sound. **RISK: shared Refine
surface (Simplify / integration gates / zero-tests / Solve) — do with a full
refine+simplify+zero_test+integrate suite run, revert on any regression.** Gate the
new path to `Times`/`Plus` radicands so the trusted single-square path is untouched.

| # | Integrand | Root cause | Recipe |
|---|---|---|---|
| 37 | `e^{-ax}J0(bx)` | `Sqrt[1+b^2/a^2]` not reduced | **FIXED v0.379** (keystone radical + driver Refine) |
| 8,9 | `e^{-ax}cos/sin(bx)` | **After v0.379 the value is already clean**: raw is `ConditionalExpression[a/(a^2+b^2), a>0 && b>0]`. The spurious `b>0` is the only blocker (cos/sin even/odd in b → valid for all real b). Pinned `Method->"DiffUnderInt"` returns `a/(a^2+b^2)` with NO condition; under Automatic the Ramanujan Mellin strip (`ramanujan.c` ~2419, `ConditionalExpression[res, cs]`) over-constrains to `b>0`. | widen that Mellin convergence strip for the `e^{-ax} trig(bx)` family to `Re[a]>0` (b real), OR have the stage decline so Automatic falls through to DiffUnderInt. Delicate (strip derivation). |
| 18 | `sin(x)^n` [0,π/2] | `Beta`↔`Gamma`: FullSimplify keeps `Beta` on one side, `Gamma` on the other, so the difference doesn't cancel | register the reverse `Gamma[a]Gamma[b]/Gamma[a+b] -> Beta[a,b]` transform (or have `integrate_beta`/`trigpower` return the Gamma form) so both sides normalise alike |
| 26 | `x^a/(1+x^2)` | `Csc[π(1+a)]` periodicity not applied by Simplify | add a trig-periodicity shift `Csc/Sec/Sin/Cos[π(k+a)]`, integer k, to the trig-Pi simp pass (`src/simp/simp_trig_pi.c`) |
| 27 | `1/√(a²−x²)` [0,a] | indefinite integral returns a **complex-log** antiderivative, not `ArcSin[x/a]`; N–L can't resolve its boundary limits at the singular endpoint x=a (numeric case works; `Limit[ArcSin[x/a],x→a]` works) | either emit the real `ArcSin` antiderivative for `(a²−x²)^(−1/2)`, or add a Beta recognizer for `(a²−x²)^s` on [0,a] (`x=a t` → `a^{2s+1}/2 · B(1/2,s+1)`) |
| 32 | `log(x)/(x²+a²)` | `ratlogpow` only covers R with **negative-real-axis** poles; `1/(x²+a²)` has complex poles ±ia | extend the Mellin `∫₀^∞ x^{s-1}/(x²+a²) dx = (π/2) a^{s−2}/sin(πs/2)`, then d/ds at s=1 for the log weight |
| 38,39 | Poisson `log(a²−2ab cos x+b²)` [0,π] | no recognizer | new recognizer `∫₀^π log(A+B cos x) dx = π log((A+√(A²−B²))/2)`, A,B free of x, A≥|B|. BUT the result's `√((a²+b²)²−(2ab)²)` is an *unfactored* Plus (`(a²+b²)²−4a²b²`), which the v0.379 radical fix does NOT reduce (it handles Times/Together'd-quotient radicands, not polynomial-perfect-square factoring); and `π log(a²) − 2π log(a)` leaves a `log(Abs[a])` residual Simplify won't drop in that composite. So this needs (a) factor-the-radicand-before-extract, and (b) the `log(Abs[a])→log(a)` composite-Simplify gap closed. |

## Soundness
Zero true wrong answers at baseline AND after both shipped fixes. Every fix must
preserve that; re-run the full battery (both metrics, per-case isolated) after each
change — no regression allowed. The remaining 9 are all **correct-value or additive**:
the 5 Simplify-blocked cases already return the right value (they fail only the
form-match), and the 4 integration cases currently abort/uneval (an additive
recognizer cannot regress a passing case). So no outstanding wrong answer exists; the
risk to manage is a shared-Refine/Simplify change silently breaking an unrelated path.

## Verification artifacts
- `scratchpad/alldefs.m` (reusable 50-case driver) → promote to
  `tests/AssumedDefiniteIntegralsCorpus.m` + runner (model `tests/test_crc_corpus.c`).
