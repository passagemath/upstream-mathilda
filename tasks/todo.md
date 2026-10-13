# Active campaign: Singular / FTC-boundary definite integrals (2026-10-13)

Plan: `~/.claude/plans/pasted-content-id-a30c-this-50-expressi-dazzling-sunset.md`
Findings: `tasks/singular_definite_integrals_2026-10-13.md`
Scope (user): FULL — close all 50.
Baseline (v0.381, corrected oracles): 20 CORRECT / 4 WRONG / 23 UNEVAL / 3 TIMEOUT.
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
- [ ] 3b. Poisson log Log[1-2a Cos x+a^2] (36,37)
- [ ] 3d. King-reflection x Sin x/(1+Cos^2) (35)
- [ ] 3f. Tan^2/Sin^2 over a period, suppress Solve::ifun (38,39)
- [ ] 16. 1/CubeRoot[x] odd — needs CubeRoot[-x]=-CubeRoot[x] + indefinite ∫ (CubeRoot machinery)
- [ ] 11. 1/Sqrt[Abs[x]] on [-a,a] — piecewise reflect + Sqrt endpoint (recheck)
- [ ] hard PV: 23,26 (symbolic pole x=a), 30 (whole-line indented contour)
- [ ] 48,49 residue hang (pre-existing)

## Phase 4 — form cleanup
- [ ] 4a. Beta[k,1]→1/k, Beta[1,k]→1/k, PolyGamma weight (1,2,3)
- [ ] 4b. Sqrt[(1-a^2)^2] reduction under assumptions (33,34)

## Release
- Local work + verify + version bump per substantive change; commit/tag/push AWAIT USER GO.
- Each change: full corpus (no new WRONG), touched unit suites, check-messages, check-c99,
  valgrind new C; docs/spec/builtins/calculus.md + changelog/2026-10-12.md.
