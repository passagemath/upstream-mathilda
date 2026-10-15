# DSolve — next 500 of Nasser's corpus (§2.2.41–§2.2.45)

Plan: `~/.claude/plans/let-s-continue-our-implementation-sleepy-rain.md`
Start version: 0.411 (after M67 / §2.2.40 / problem 4000).

Per block: fetch → convert → indvar-audit → register gate+scoreboard → measure
baseline → diagnose dominant gap → general root-cause fix → re-measure + lower
baseline → regression A/B → track (STATUS/README/reports/changelog/PLAN) →
version bump + git tag.

## Block 1 — §2.2.41 (problems 4001–4100), milestone M69, target v0.412
- [x] Fetch `Ch2.S1.SS41.htm`
- [x] Convert → `DE_examples_2241.m` (label 2.2.41) — 100 scalar, 6 IVP
- [x] `make check-corpus-indvar` green + eyeball spot-check
- [x] Add `dsolve_corpus_2_2_41_tests` to CMakeLists (baseline 25) + STATUS/README blocks
- [x] Measure honest baseline → 71/100 (29 non-PASS); dominant gap = regular-singular 2nd-order linear
- [x] Land general fix: confluent-Whittaker recogniser in normal-form pre-pass → 75/100
- [x] No masked-wrong-answer (fully 1F1-verified; degenerate 2μ∈ℤ correctly declines)
- [x] Regression A/B on 7 prior sections: 0 regress, §2.2.20 −4 bonus; specialform unit-test exprs byte-identical HEAD
- [x] version bump 0.411→0.412; changelog + PLAN M69 + report written
- [x] gates: §2.2.41 ctest 25≤25 PASS; check-c99 rc=0; check-messages OK; dsolve_tests aborts on PRE-EXISTING exact-ODE test (byte-identical HEAD, not mine)
- [x] commit + tag v0.412

## Block 2 — §2.2.42 (4101–4200), M70, v0.413 ✅
- [x] Fetch + convert → DE_examples_2242.m (88 scalar/20 IVP, 12 systems); indvar green
- [x] Baseline 98/100 (2 non-PASS: 4175 coupled system, 4183 integer-c Jacobi)
- [x] General fix: coupled leading-matrix LinearSolve fallback in dsolve_linsys_extract_Ab → 99/100
- [x] A/B systems-heavy (§2223/2229/2239/2240): 0 regress, 0 FAIL; valgrind no new leaks (baseline-identical)
- [x] gate baseline 1; STATUS/README/PLAN/changelog; version 0.412→0.413
- [x] commit + tag v0.413

## Block 3 — §2.2.43 (4201–4300), M71, v0.414 ✅
- [x] Fetch + convert → DE_examples_2243.m (100 scalar/15 IVP); indvar green
- [x] Baseline 92/100 (8 non-PASS; dominant = 5-case first-order _linear cluster)
- [x] Root cause: ∫Sec·Tan² (≡Sin²/Cos³) spins 12s COLD in Integrate Linearity stage (TrigReduce can't clear trig denom)
- [x] General fix: Linearity declines a pure Sin^a Cos^b monomial w/ trig denom → Weierstrass closes (0.03s). integrate.c → 93/100
- [x] Caught + fixed own over-broad-guard regression (§2.2.35-3461); narrowed to pure monomials
- [x] integrate+risch ctest 100% (33/33); A/B §2230/35/37 no regression; 3497 pre-existing flaky
- [x] gate baseline 7; STATUS/README/PLAN/changelog; version 0.413→0.414; check-c99/messages green
- [x] commit + tag v0.414

## Block 4 — §2.2.44 (4301–4400), M72, v0.415 ✅
- [x] Fetch + convert → DE_examples_2244.m (100 scalar/2 IVP; 29 sympy=False); indvar green
- [x] Baseline 91/100 (9 non-PASS; one = 4369 non-interruptible HANG)
- [x] Root cause: 4369 (trivial linear ODE) needs ∫Sin²/(1+Cos) — a NON-monomial rational trig that thrashes Linearity (M71 monomial guard missed it)
- [x] General fix: generalize guard to is_rational_trig_in_x (Weierstrass domain); 4369 hang → 0s PASS → 92/100
- [x] integrate+risch ctest 100%; A/B §2230/35/36/37 identical non-PASS; §2.2.43 still 7; 3461 still PASS
- [x] gate baseline 8; STATUS/README/PLAN/changelog; version 0.414→0.415; c99/messages green
- [x] commit + tag v0.415

## Block 5 — §2.2.45 (4401–4500), M73, v0.416 ✅
- [x] Fetch + convert → DE_examples_2245.m (100 scalar/0 IVP; 17 sympy=False); indvar green
- [x] Baseline 93/100 (7 non-PASS; heterogeneous hard/gated gaps)
- [x] General fix: NEW DSolve`InverseOrientation method (reciprocal dx/dy==1/F, recurse, implicit first integral) → 4408, 94/100
- [x] Caught+fixed own latency regression (4323 via nested recursion) → gated top-level-only; killswitch A/B 0 regression
- [x] Bonus: fixed pre-existing solverad.c container-array leak (valgrind back to baseline)
- [x] gate baseline 6; STATUS/README/PLAN/changelog/calculus.md; version 0.415→0.416; c99/messages green
- [x] commit + tag v0.416

## Review — the next 500 (§2.2.41–§2.2.45, M69–M73, v0.412–v0.416)
All five 100-problem blocks (problems 4001–4500) landed as gated corpus sections,
each with a general root-cause fix, 0 FAIL, 0 crash, A/B-confirmed no regression:
- §2.2.41 (M69, v0.412): 71→75. Confluent-Whittaker recogniser in the normal-form pre-pass (1F1 closed forms for y'-carrying regular-singular 2nd-order linear).
- §2.2.42 (M70, v0.413): 98→99. Coupled leading-form first-order linear SYSTEM (leading-matrix LinearSolve fallback).
- §2.2.43 (M71, v0.414): 92→93. Trig-rational ∫Sin^a/Cos^b monomial no longer spins in the Integrate Linearity stage (→ Weierstrass).
- §2.2.44 (M72, v0.415): 91→92. Generalised that guard to any rational-trig-of-x (removes a non-interruptible hang; 4369).
- §2.2.45 (M73, v0.416): 93→94. New DSolve`InverseOrientation (x(y) reciprocal) + a Solve-subsystem leak fix.
Each general fix validated: integrate+risch ctest 100% (33/33) for the Integrate
changes; killswitch / HEAD A/B for the DSolve changes; valgrind clean; all prior
gates held (the pre-existing M68 drift in older §2.2.x and the pre-existing
dsolve_tests exact-ODE abort were confirmed HEAD-identical, not introduced here).
5 git tags pushed (v0.412–v0.416).
