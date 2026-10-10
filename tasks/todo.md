# DSolve corpus wave M67 — §2.2.40 (Problems 3901–4000)

Plan: `/Users/user/.claude/plans/let-s-continue-our-implementation-cozy-minsky.md`

## Phase A — land & audit the corpus
- [x] Fetch upstream §2.1.40 `Ch2.S1.SS40.htm` (browser-UA curl) → scratchpad (HTTP 200)
- [x] Convert → `DSolve_test_status/DE_examples_2240.m` (label 2.2.40)
- [x] `make check-corpus-indvar` green
- [x] 100 records: 74 scalar (58 IVP) + 26 systems — Laplace-transform-style chapter

## Phase B — measure baseline
- [x] Build `dsolve_corpus_tests`; baseline run → 100/100, 0 FAIL/UNEVAL/crash
- [x] Staff-engineer spot-check (NOT vacuous): forced FIRST-ORDER family (3955–3961,
      3970–3975) are **masked wrong answers** — inert `∫μq` + malformed `Integrate[.,0]`
      from the IC-fit; wrong past the step/impulse (3955 y(2) garbage; 3972 drops impulse)

## Phase C — diagnose & fix (correctness: masked FAIL class)
- [x] Root cause: `DSolve`LinearFirstOrder` (slot 4) owns step/impulse IVPs with the
      indefinite integrating-factor form; 2nd-order already correct (PiecewiseForcing / VoP)
- [x] Fix: route first-order — numeric-breakpoint step → PiecewiseForcing; impulse /
      symbolic-breakpoint → new `dsolve_linear_factor_solve_causal` (definite convolution);
      shared `linfac_mu` refactor. `dsolve_linear1.c` + `dsolve_common.{c,h}`
- [x] §2.2.40 → 100/100, 0 FAIL; all 13 rerouted cases now genuinely correct (verified)
- [x] Regression (7 sections): 0 FAIL; §2.2.1/6/16 +2 each, others unchanged
- [x] A/B vs HEAD: 753/2973 flutters + dsolve_tests `t_m17_normalform_bessel` abort all
      PRE-EXISTING (identical on HEAD binary) — not this wave

## Phase D — lock in & document
- [x] Gate `dsolve_corpus_2_2_40_tests` baseline 0; `ctest` passes (32 s)
- [x] STATUS.md §2.2.40 block + README Contents row; DSOLVE_PLAN.md M67 block
- [x] Changelog `docs/spec/changelog/2026-10-05.md`; spec `calculus.md` LinearFirstOrder row
- [x] Bump `src/version.h` → 0.344
- [x] `make check-c99` ✓ / `make check-messages` ✓
- [x] Regenerate `reports/2.2.40.{tsv,md}`
- [ ] Commit + tag v0.344 + push --follow-tags

## Review
**Outcome.** §2.2.40 lands **100/100, 0 FAIL, 0 crash, 0 UNEVAL**. The count was 100/100
out of the box, so the wave is a *correctness* fix the harness could not see: 13 forced
first-order IVPs (UnitStep/Piecewise/DiracDelta) were returning functions wrong past the
step/impulse, masked because the harness samples the ODE residual (FTC-zero on the inert
integral) and checks the IC only at t=0.

**Fix (one root cause, 3 files).** `dsolve_linear1_try` stops building the broken inert
integrating-factor form for a distributional/step IVP and routes to a correct method:
numeric-breakpoint step → the proven `DSolve`PiecewiseForcing`; DiracDelta impulse or
symbolic breakpoint → the causal definite convolution `y=μ⁻¹(C[1]+∫_0^x μq ds)` (new
`dsolve_linear_factor_solve_causal`, definite DiracDelta sifts to HeavisideTheta). Gated to
IVPs; smooth / no-IC path behaviour-identical (shared `linfac_mu`).

**Verification.** 13 rerouted cases now clean Piecewise/HeavisideTheta forms, residual-
verified at points past the breakpoint. Side-effect: §2.2.1/6/16 +2 each. 0 regression,
A/B-confirmed. Residue: none.

**Pre-existing, out of scope (flagged):** `dsolve_tests` aborts at `t_m17_normalform_bessel`
(SpecialFunctionForm declines the spherical-Bessel eqn `y''+(2/x)y'+y==0`) — identical on
HEAD; §2.2.8-753 (21 s) / §2.2.30-2973 (>30 s) latency drift, identical on HEAD.
