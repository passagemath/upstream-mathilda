# COMPLETE campaign: Extreme definite integrals — battery #2 (2026-10-13)

Findings: `tasks/extreme_definite_integrals_2026-10-13.md`
Scope (user): FULL — close all 24 open scored cases.  **DONE: 48/48 scored correct, 0 wrong.**
Baseline (v0.387, verified oracles): 26 CORRECT / 0 WRONG / 19 UNEVAL / 4 TLE / 1 EXCLUDED(37).
Final (v0.400): 48 CORRECT / 0 WRONG / 0 UNEVAL / 0 TLE / 2 EXCLUDED (19, 37).
Recognizers added v0.388–v0.400 (13 bumps): Frullani-cos+CosIntegral-series, hyperbolic-halfline
table, Weber-Schafheitlin Bessel, keyhole-log symbolic-param hang fix (+battery#1 48/49),
rational-halfline PV, branch-power Mellin PV, whole-line log-quad (+even-fold), Raabe Log-Gamma,
ratlog01 Hurwitz-zeta, half-Beta (1-x^2)^e, Gaussian-moment difference, Mellin PV 1-x^d,
whole-line Log[Abs]. Battery #1 bonus: 42 -> 44.
Corpus `tests/ExtremeDefiniteIntegralsCorpus.m`; gate `extreme_definite_corpus_tests`;
baseline driver (scratchpad) `run_extreme.sh`.
4 corrected battery oracles: 6 (×2), 24 (wrong formula), 36 (sign), 37 (excluded).
2 repointed instantiations: 12 (a→1,b→3), 19 (a→1/3).

## Phase 0 — corpus + committed gate  [DONE, contributor-facing, no bump]
- [x] corpus (50 records, 4 corrected oracles, 2 repointed instantiations)
- [x] committed fork-per-case gate (WRONG baseline 0; ConditionalExpression-strip; case-37 exclude)
- [x] CMake registration; build + run verified EXIT=0 (26 correct / 0 wrong)

## Phase 1+ — recognizers (correctness-first, additive clusters before PV-core)
Clusters (see findings): Frullani-cos+Gaussian-diff (2,6); hyperbolic kernels (24,25,29,30);
Bessel Weber-Schafheitlin (26,27,28); log half-line/Raabe (41,46,47,50); whole-line log (38,44);
half-line rational PV on-axis pole (31,32,33,34); Mellin+log/√ PV (12,19,35,36); Fresnel form (21).

---

# PRIOR campaign (COMPLETE): Singular / FTC-boundary definite integrals (battery #1)

Plan: `~/.claude/plans/pasted-content-id-a30c-this-50-expressi-dazzling-sunset.md`
Findings: `tasks/singular_definite_integrals_2026-10-13.md`
Scope (user): FULL — close all 50.
Baseline (v0.381, corrected oracles): 20 CORRECT / 4 WRONG / 23 UNEVAL / 3 TIMEOUT.
Final: 42/50 correct, 0 wrong (v0.382–v0.387, all pushed).
Corpus `tests/SingularDefiniteIntegralsCorpus.m`; gate `singular_definite_corpus_tests`;
driver `scratchpad/run_corpus.sh`.

## Phase 0 — corpus + harness  [DONE, contributor-facing, no version bump]
- [x] corpus file (50 records, corrected oracles 15,41)
- [x] committed fork-per-case gate (WRONG baseline 5, ratchets down)
- [x] CMake registration + build/run verified (EXIT=0)
- [x] fast one-process-per-case driver

## Phase 1 — soundness (eliminate the wrong answers)  [WRONG 5 → 0, v0.382, verifying/release]
- [x] 1a. √(perfect square)→Abs normalization + periodic breakpoints (18) — piecewise
- [x] 1b. Tan/Cot interior-pole detection (trig heads) + periodic enum → PV=0 / non-PV idiv (24,25)
- [x] 1c. reflect nonpositive sub-interval x→-x (6,14); dispatcher rejects x-leaking result (12 sound UNEVAL)
- [ ] 1d. 48/49 `Log^n x/(a^2+x^2)` half-line hang — DEFERRED (pre-existing residue loop on
      symbolic param; Ramanujan solves it instantly → needs ordering/early-decline fix)
- [x] N-L defensive nets: radical-branch-hazard crosscheck; complex-value→idiv; Quiet pole-solve
- [x] committed gate WRONG baseline 5 → 0 (assert-empty); corpus 24 correct / 24 uneval / 2 timeout
- [x] regression: 16 integration/limit/refine suites green; check-messages, check-c99 pass
- [x] docs: changelog v0.382 + calculus.md (piecewise, N-L pole detection, continuity guard); version bump

## Phase 2 — structural keystones
- [x] 2a. assumption-aware Newton–Leibniz (unlocks 4,40 + more) [v0.383]
      `assumptions` threaded into `_try_pv`; `nl_classify_symbolic` (Refine-proved interior/
      exterior for symbolic root/bound); boundary Limit+Refine under assumptions; internal
      Refine probes Quiet-wrapped. Corpus 24→31 correct, 0 wrong; 17 suites green.
- [~] 2b. PrincipalValue symbolic bounds — DONE for concrete-pole cases 21,22,27,28,29
      (+24,25 from P1). `nl_principal_value` accepts symbolic bounds (local odd-order probe +
      assumption-aware boundary eval); numeric realness gate rejects a spurious complex PV;
      diffunderint skipped under PV; driver Refine guarded against INTRODUCING complexity
      (was stripping a realness-protecting Re[...]). Corpus 31→36 correct, 0 wrong.
      STILL OPEN: 23,26 (symbolic pole x=a — needs symbolic-pole PV placement); 30 (whole-line
      PV, residue indented contour). The symmetry-idiv-leak under PV (25) remains cosmetic.

## Phase 3 — additive recognizers  [corpus 39/50 correct, 0 wrong, v0.385]
- [x] 3a. Log[Sin]/Log[Cos] (9,10) — `integrate_logtrig_try` in integrate_beta.c [v0.385]
- [x] 3c. reciprocal u=1/x on [a,Inf) (41,42,43) — free from 2a (assumption-aware boundary limits)
- [x] 12. symbolic |x|^-b — N-L solves the base of a Power denominator when Solve[x^b==0]
      is undecidable; transient 0^(1-b) Quiet-gated in the probes [v0.385]
- [x] 3e-partial. (a^2-x^2)^s via 2a: 4,5,11(no — 11 still uneval),40 closed
- [x] 3b. Poisson log Log[1-2a Cos x+a^2] (36,37) — integrate_poissonlog_try [v0.386]
- [x] 3d. King-reflection x Sin x/(1+Cos^2) (35) — integrate_kingreflect_try [v0.387]
- [x] 11. 1/Sqrt[Abs[x]] on [-a,a] — free from the Power-base pole fix [v0.385]
- [ ] 3f. Tan^2/Sin^2 over a period, suppress Solve::ifun (38,39) — OPEN
- [ ] 16. 1/CubeRoot[x] odd — needs CubeRoot[-x]=-CubeRoot[x] + indefinite ∫ (CubeRoot machinery) — OPEN
- [ ] hard PV: 23,26 (symbolic pole x=a), 30 (whole-line indented contour) — OPEN
- [ ] 48,49 residue hang (pre-existing; Ramanujan solves instantly) — OPEN

## STATE: corpus 42/50 correct, 0 wrong (v0.382–v0.387 all committed+pushed to main).
Open 8: 16, 23, 26, 30, 38, 39 (uneval) + 48, 49 (residue hang).

## Phase 4 — form cleanup
- [ ] 4a. Beta[k,1]→1/k, Beta[1,k]→1/k, PolyGamma weight (1,2,3)
- [ ] 4b. Sqrt[(1-a^2)^2] reduction under assumptions (33,34)

## Release
- Local work + verify + version bump per substantive change; commit/tag/push AWAIT USER GO.
- Each change: full corpus (no new WRONG), touched unit suites, check-messages, check-c99,
  valgrind new C; docs/spec/builtins/calculus.md + changelog/2026-10-12.md.
