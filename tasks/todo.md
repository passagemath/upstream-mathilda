# ACTIVE: auto-pack complex arrays (RandomComplex + complex Lists) (2026-10-14)

Surfaced by the HermiteH efficiency review (below): complex machine arithmetic is
25x faster than numpy and 1000x faster than Wolfram **when the data is a packed
complex64 buffer**, but `RandomComplex`/complex `List`s do NOT auto-pack, so the
common workflow sits at 0.74 M/s (slower than numpy) instead of 115 M/s — a 156x
self-gap. System-wide (Sin/Exp equally slow on a complex list), not HermiteH-specific.

Root cause (pack.c:80-90, deliberate+documented): `leaf_class` declines Complex
because a packed complex element with im=0 would materialise as `Complex[re,0]`,
which the evaluator folds to Real → "would not round-trip."

Why it is now safe to lift the exclusion (verified empirically 2026-10-14):
  - Mathilda folds `Complex[x,0.]`->Real AND `3.===Complex[3.,0.]` is True, so a
    list can be all-Complex-headed ONLY if every element has im!=0 (any zero-im
    Complex already folded to Real, making the list mixed). => packing an
    all-Complex list is inherently lossless with NO materialisation fold.
  - Wolfram also refuses a mixed {Complex,Real} list (PackedArrayQ False) => the
    right rule is "homogeneous-Complex only", matching the existing no-mix policy
    for Int/Real.

## Plan  [IMPLEMENTED — v0.402]
### Part 1 — pack.c substrate (makes ALL complex Lists auto-pack)  [DONE]
- [x] `PackClass`: added `PK_COMPLEX`.
- [x] `leaf_class`: Complex[Real,Real] -> PK_COMPLEX (exact Gaussian declines).
- [x] `pack_fold_class`: PK_COMPLEX unifies with PK_COMPLEX; coerce widens Int/Real->complex.
- [x] `pack_sniff` packed-row branch: complex dtype -> PK_COMPLEX.
- [x] `pack_build`: cls==PK_COMPLEX -> NDT_COMPLEX64; dropped the complex-want refusal.
- [x] `pack_flatten` + `pack_complex_parts` helper: complex dtype writes (re,im) per leaf.
- [x] NO change to `ndarray_buffer_element_to_expr` (fold unnecessary — verified).
### Part 2 — random.c (RandomComplex builds a packed buffer directly)  [DONE]
- [x] 1D + nD machine RandomComplex fills a complex64 buffer via ndbuild_open;
      degenerate zero-imaginary range keeps the Real-valued List path; re-then-im
      draw order preserved (seeded stream unchanged, A/B-confirmed).
### Verify  [DONE except fastpath-sweep in-flight]
- [x] check-array-exactness PASS (0 MIXED); check-nd-surfaces PASS; check-compile-coverage
      PASS (no new gaps). check-packed-aware: only the pre-existing StandardForm red
      (print.c; my diff adds no AWARE entry/dispatch — provably unchanged).
- [ ] check-fastpath-sweep (running; expect only pre-existing VertexQ/WeightedGraphQ/Xor).
- [x] check-c99 PASS; check-messages PASS (BASELINE empty).
- [x] focused blast-radius A/B (complex/ndarray/packed/random/linalg/fourier, 22 suites):
      IDENTICAL FAIL set my-build vs HEAD -> NO regression. Pre-existing fails only.
- [x] valgrind: definitely-lost = 13,440 B / 420 blk (identical fixed startup baseline);
      32-byte blocks are macOS dyld/objc, not pack/random/ndarray -> NO new leak.
- [x] functional: RandomComplex 1D+nD packs; complex Lists (Table) auto-pack; degenerate
      {0,1} stays Real; seeded values stable; HermiteH complex values unchanged.
- [x] bench: HermiteH[20, RandomComplex[1e6]] 1.29s->0.017s (74x); Sin 2.68s->0.015s (175x).
- [x] docs: changelog v0.402 + packed_arrays.md §6/intro updated; version.h -> 0.402.
- [ ] commit + tag v0.402 (AWAIT USER GO). NB: HermiteH (v0.401) still uncommitted too.

## HermiteH efficiency review verdict (2026-10-14) — ALL GREEN except the above
Correctness vs mpmath oracle (Wolfram + scipy cross-checked): machine real/complex
rel.err <=9.5e-16 (= Wolfram, ULP-level); MPFR meets requested precision every case
(50-dig 4.7e-52, cancellation 40-dig 1.3e-41); n=1000 overflow correctly falls to
MPFR (scipy nan). Exact polys incl. deg-101 coeff list identical to Wolfram; Series
0/x0/Inf + SeriesCoefficient + derivative algebraically identical to Wolfram.
Throughput (elem/s): real 135M (Mathilda) vs 25M scipy vs 0.3M Wolfram; complex
*packed* 115M vs 4.6M numpy vs 0.11M Wolfram; MPFR+non-integer: Mathilda fastest of
the three. => HermiteH itself needs no change; the one actionable item is the
complex auto-packing above.

---

# DONE (awaiting commit): HermiteH[n, x] — physicists' Hermite polynomials (2026-10-14)

Plan: `~/.claude/plans/pasted-content-id-629c-in-src-special-f-zesty-neumann.md`
Template: ChebyshevU (integer coeffs, FLINT `fmpz_poly_hermite_h`); reuse `chebyshev_common.{c,h}`.

## Phase 1 — core evaluator + numerics + fast paths  [DONE]
- [x] `src/special_functions/hermite_h.{c,h}` (driver, poly, ncpx recurrence, H_n(0), 1F1 form w/ Re(x^2) cancellation guard)
- [x] `sf_machine.{h,c}`: `sf_machine_hermite_h` + `_c`
- [x] `sym_names.{h,c}`: SYM_HermiteH (extern + =NULL + intern)
- [x] `core.c`: include + `hermite_h_init()`
- [x] `info.c`: docstring
- [x] `ndkernels.c`: `ndk_HermiteH_c` + `REG_B(HermiteH, ...)`
- [x] build clean; REPL smoke all match spec (poly, machine 940., complex -731.082+940.003I, MPFR 50dig, argrx); Compile scalar+array True; ND real+complex on buffer
## Phase 2 — derivative rule  [DONE]
- [x] `deriv.c`: `SYM_HermiteH` 2-arg block → `2 n HermiteH[-1+n,z]`; chain rule; free Series at 0/x0/1/2 all match spec
## Phase 3 — bespoke Series-at-Infinity + symbolic-m SeriesCoefficient  [DONE]
- [x] `series.c`: `try_series_hermiteh_at_infinity` ((2x)^n 2F0) + dispatch; numeric-verified ~1e-16
- [x] `series.c`: symbolic-m SeriesCoefficient Piecewise (x0=0); verified ~1e-49 vs concrete
## Phase 4 — tests + docs + changelog + version  [DONE]
- [x] `tests/test_hermite_h.c` (13 groups) + CMake (COMMON_SRC + target) — all pass, 0 FAIL
- [x] `docs/spec/builtins/special-functions.md` (## HermiteH + intro); changelog 2026-10-12.md; version.h -> 0.401
## Verify  [in progress]
- [x] clean build; check-c99 PASS; check-messages PASS
- [x] hermite_h_tests PASS (0 FAIL); valgrind: identical 13440B/420blk baseline as cheb/laguerre → NO new leak
- [x] series/deriv regression: 7 PRE-EXISTING series_tests FAILs (printer `1/(360 x^3)`, ParallelMixed Integrate head, Sqrt-assumptions) — none HermiteH; all other 7 suites 0 FAIL
- [x] 5 packed/compile audits: array-exactness PASS (0 MIXED), nd-surfaces PASS, compile-coverage PASS ("no new gaps", HermiteH compiles) — HermiteH clean on all three. check-packed-aware + fastpath-sweep fail ONLY on pre-existing unrelated heads (StandardForm/print.c; VertexQ/WeightedGraphQ/Xor) — my diff touches none of those files, HermiteH not flagged.
- [ ] commit + tag v0.401 (AWAIT USER GO)

## RESULT: HermiteH complete & verified. Builtin + numerics (machine/MPFR, real/complex,
## order & arg) + NDArray/packed/Compile + derivative + Series(0/x0/Infinity) + SeriesCoefficient.
## All spec examples reproduced. 3/5 surface audits green (HermiteH clean on all 5); 2 reds pre-existing.
## No new leak. Not committed — awaiting user go for commit + `git tag v0.401`.

---

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
