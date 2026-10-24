# SparseArray through the functional and procedural builtins (2026-10-24)

Plan: ~/.claude/plans/let-s-continue-our-implementation-clever-frog.md. Spec: docs/spec/builtins/data-structures.md ("Functional and procedural heads").

- [x] Reference ground truth (wolframscript) for ~120 functional/procedural call shapes
- [x] W1 level specs: native Map (leaf-level O(nnz), shape-preserving resparse), MapIndexed/Apply/Scan with a level spec (src/sparse/sparse_func.c)
- [x] W2 iterators: Table/Do/Sum/Product over {i, s} (iter.c, one fix point)
- [x] W3 ordering/set/fold heads answered on Normal (Sort/SortBy kept sparse; Fold keeps a sparse initial value); atom ::normal for Sort/SortBy/Union/Intersection/Complement/Split/DeleteDuplicates (builtin_normal_error)
- [x] W4 AppendTo/PrependTo on sparse; s[[All,k]] = sparse; shape-breaking assignment leaves a List
- [x] W5 Thread over sparse (incl. Thread[s]); Flatten DEEP; Outer sparse only for Times/one-arg; Times with dense arrays sparse; lower-rank sparse broadcast (m + v)
- [x] W6 Clip/Chop element-wise; resparse background = commonest element when the leading one vanished
- [x] listable_sparse: densify cap + Thread::tdlen settle (Range[0] bigv no longer allocates 10^9)
- [x] Dense fixes found: Map/Apply/Scan atomic Rational/Complex; Map All; packed w[[All]] = 0
- [x] tools/sparse_head_sweep.py (sparse vs dense vs reference, 958 heads x 10 shapes); fuzzer func/proc families
- [x] tests/test_sparse_func.c
- [ ] Follow-up (pre-existing, dense): Order/Sort compare Lists from the LAST element (Plus/Times canonical order) -- Sort[{{0,x},{1,0}}] gives {{1,0},{0,x}}; reference {{0,x},{1,0}}. Core expr_compare change, large blast radius.
- [ ] Follow-up (pre-existing, dense): ~700 dense-side differences in the head sweep (DENSE class), e.g. Floor[0, 0] -> 0 (ref Indeterminate), AnyTrue with a non-Boolean test stays unevaluated (ref builds Or), Thread with unequal lengths gives no Thread::tdlen. Re-run: python3 tools/sparse_head_sweep.py
- [ ] Follow-up: Map[f, 1/2, {0, Infinity}] returns 1/2 (ref f[1/2])
- [ ] Follow-up: First[x, s] / Last[x, s] with a sparse default densify the default (head only)
- [ ] Follow-up: Table[i, {i, v}] is a List; the reference re-sparsifies a one-dimensional Table (documented)

# SparseArray stress test, fixes, and beating scipy.sparse (2026-10-24)

Plan: ~/.claude/plans/let-s-thoroughly-stress-test-indexed-hopcroft.md.

- [x] Phase 0: head sweep (1247 heads x 22 call shapes) -> gate modes; `tools/sparse_stress.py` differential fuzzer
- [x] Phase 1: transparency gate (`src/sparse/sparse_gate.c`, `SymbolDef.sparse_mode`, eval pre/post hooks, 2^27 cap + SparseArray::dense, MATHILDA_SPARSE_DIAG=gate)
- [x] Fuzz-found fixes: packed 0^0, `{0.,1.}^0` -> 1, `2.5^0` -> 1., -0.0 printing, symmetric Eigenvalues tie order, linalg FullForm messages, double message on declining builtins (evaluator settled step)
- [ ] Full suite A/B vs HEAD worktree (core evaluator change) before push
- [ ] Phase 2: native kernels (Dot SpMV/SpMM/SpGEMM, Transpose, structural ops, reductions, IdentityMatrix/DiagonalMatrix/KroneckerProduct sparse constructors, sparse LinearSolve direct+Krylov, Eigenvalues[s, k] Lanczos/Arnoldi)
- [ ] Phase 3: benchmark 120 broadened, 3-way, every row <= 1.0x
- [ ] Phase 4: docs, valgrind on the fuzz corpus, `make check-sparse`
- [ ] Follow-up (pre-existing, not sparse): `x = Total[t]` re-emits a message after Set (eval clock churn invalidates the stamp)
- [ ] Follow-up (pre-existing): `Eigensystem` does not exist
- [ ] Follow-up (pre-existing): non-symmetric machine Eigenvalues tie order: {{1.,2.},{0.,-1.}} -> {1.,-1.}, expected {-1.,1.}
- [ ] Follow-up (pre-existing): `ReplaceAll[s, 1]` (a non-rule) returns s; it should stay unevaluated

# SparseArray (2026-10-23)

Plan: ~/.claude/plans/pasted-content-id-bdcb-let-s-plan-smooth-shore.md. Spec: docs/spec/builtins/data-structures.md.

- [x] Semantics probes (canonical form, raw-object matching, construction, merge rule, Band, Part, properties, messages)
- [x] EXPR_SPARSE node type (src/sparse/): payload, shared pattern, canonical kinds; expr/eval/match/print/sort/precision plumbing
- [x] Constructor (rules, {pos}->{vals}, patterns, Band incl. blocks, dense/packed, existing arrays, standard form) + messages
- [x] COO -> CSR: composite-key LSD radix (12-bit digits), first-wins dedup, background drop; packed-rule direct keys
- [x] Normal (packed when one machine type), Length/Dimensions/ArrayDepth/ArrayQ/VectorQ/MatrixQ, Part read + assign, properties, ArrayRules
- [x] Listable threading (union rule, Times intersection), typed Plus/Times/Power incl. fused merges, ND-kernel reuse, generic tier, N
- [x] Tests: tests/test_sparsearray.c (constructors, messages, raw object, structure, Part, arithmetic, functions, fast-path cross-checks); valgrind clean (startup baseline only)
- [x] Audits: check-c99, check-messages, check-packed-aware, check-compile-coverage, check-nd-surfaces, check-array-exactness
- [x] Docs: data-structures.md, changelog v0.572, docstrings (SparseArray, ArrayRules, ArrayQ, ArrayDepth, Band, Normal)
- [x] benchmarks/120-sparse-arrays vs scipy.sparse: 13 cases, all at parity or ahead (v0.573)
- [ ] Follow-up session: linear algebra (Dot/SpMV/SpGEMM, Transpose, Total, LinearSolve, Det, Inverse, eigen), Map/Flatten, graph heads returning SparseArray
- [ ] Follow-up (pre-existing, not sparse): Sqrt[-4.] gives 1.22e-16 + 2. I (scalar and packed); the exact answer is 0. + 2. I
- [ ] Follow-up (pre-existing): standalone `r = pos -> vals` unpacks both sides (Rule is not packed-aware by the no-nesting invariant)

## Review
- Element-wise results equal the operation on Normal[...] element by element (no whole-operand Real coercion) -- documented in data-structures.md.
- Performance at or ahead of scipy.sparse on every measured op (build 40 vs 41 ms; A+B 5.8 vs 10; A*B 3.0 vs 7.2; Sin 2.8 vs 16).
- Pre-existing failures seen in neighbouring suites (packed_list ToNDArray mixed, numeric N[E^(I Pi/4)-(-1)^(1/4),32], series Integrate[Series]) reproduce identically on a HEAD build.

# Special-function Simplify engine (2026-10-22)

Plan: ~/.claude/plans/similar-to-our-simp-bessel-merry-pudding.md. Derivations: docs/design/simp_special_functions.md.

- [x] W0 engine (`simp_sf_engine.{c,h}`) + Bessel port; Casoratians, K evenness, integer reflection
- [x] W1 Gamma / incomplete Gamma / LogGamma / BarnesG / Pochhammer / Beta / Binomial / Factorial; reflection, multiplication
- [x] W2 PolyGamma / HarmonicNumber / HurwitzZeta / Zeta[s,a] / LerchPhi / Riemann FE / PolyLog
- [x] W3 Erf family, InverseErf, ExpIntegralE, LogIntegral/Ei, Si/Ci/Shi/Chi, Fresnel
- [x] W4 orthogonal polynomials, Airy, StruveH, Whittaker, Fibonacci/LucasL
- [x] W5 hypergeometric (PFQ canonical order, contiguous lattices, Kummer/Euler)
- [x] W6 elliptic K/E (Legendre relation), ProductLog, QPochhammer
- [ ] Follow-up: N[LaguerreL[nu, -1, x]] stays unevaluated (numeric generalized Laguerre at a = -1)
- [ ] Follow-up: Simplify bottom-up spends minutes in TrigReduce on `E^-w (poly in w^(2a+j))` (no special function involved); see test_simp_sf_gamma.c GammaIncomplete note

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

---

# DSolve — the next 500 (§2.2.46–§2.2.50, Problems 4501–5000, M74–M78)

Plan: `~/.claude/plans/let-s-continue-our-implementation-twinkling-reef.md`. Executed as ONE
batch (fixes are shared across sections) in two commits: core CAS fixes (v0.434) + DSolve
batch with the five corpora and gates (v0.435).

- [x] Fetch + convert all five (§2.2.46 has 55 systems); indvar green; spot-audit →
      converter fix `\operatorname{arccot}` → ArcCot (4924, §2.1.2-756)
- [x] Baseline 89/83/76/85/90 = 423/500, 0 FAIL
- [x] General fixes (see DSOLVE_PLAN M74–M78): Together/Cancel binder + Log semantics,
      Factor SIGFPE, ∫Secⁿ, quiet Risch probes, inert integrating factor, Riccati
      particular, ExactODE inert, Kummer U, separable split + branch check, UC screen,
      system trig 2nd pass
- [x] Full per-record A/B vs HEAD (51 sections) → found + fixed 8 regressions incl. 2 FAIL,
      1 wrong-as-UNEVAL, 2 crashes; final: 0 lost, +41 old; new 94/89/85/88/92 = 448/500
- [x] Unit suite: only pre-existing HEAD failures; ratcanon test updated; c99/messages green
- [x] Gates 6/11/16/12/8; reports; STATUS/README/PLAN/changelog/calculus.md/algebra.md
- [x] valgrind on the new paths
- [x] commit + tag v0.434, v0.435

## Review — next 500 (§2.2.46–50)
423 → 448 / 500, 0 FAIL, 0 crash; older sections +41 with zero per-record losses. The
biggest wins were not DSolve methods but core correctness: Together/Cancel were silently
value-changing (binder substitution, Log expansion), Factor could crash the kernel, and
∫Secⁿ hung. Ten older gates are red on HEAD already (drift since M73) — flagged, not
rebaselined.

# DSolve — the next 500 (§2.2.51–§2.2.55, Problems 5001–5500, M79–M83)

Plan: `~/.claude/plans/let-s-continue-our-dsolve-snoopy-biscuit.md`. One batch, as M74–M78.

- [x] Fetch + convert all five (all 100 scalar, 0 IVP, 0 systems; Kamke-ordered first-order:
      Abel 2nd type / homogeneous / dAlembert / Clairaut / quadrature-in-y'); indvar green
- [x] Baseline (HEAD v0.445 binary snapshot) → reports
- [x] Triage + general fixes
- [x] Full per-record A/B vs HEAD on all 56 sections; 0 FAIL / 0 crash / 0 lost
- [ ] Unit suites, c99/messages, valgrind
- [x] Gates, reports, STATUS/README/PLAN/changelog/builtins docs
- [ ] Version bump(s), commit, tag, push

## Review — next 500 (§2.2.51–55)

- **Result:** 431 → **462 / 500**, 0 FAIL. Gates 18 / 4 / 5 / 9 / 7 (gap 17 / 4 / 4 / 8 / 5 + a named load margin). Three of HEAD's 431 were
  wrong answers (5005/5007/5008).
- **A/B** vs HEAD v0.449, all 51 older sections, same load: **+22** (passed by neither HEAD
  run). The 4 losses are HEAD wrong answers (2.1.2-777, 2.2.48-4777/4779/4792; residuals
  checked), plus 2 flutterers that pass on rerun. §2.2.48 gate 16 → 20 (+2 for the slow new gains 4774/4727), §2.2.35 2 → 1,
  §2.2.36 7 → 5.
- **The late-slot lesson:** the first full A/B showed 14 real losses. Bernoulli-on-`Expand`
  and Exact-by-μ(u) sat *inside* their early methods and claimed equations the
  specialists closed. Moving both to late slots, plus the linear integrating-factor
  constant drop, cleared all 14.
- **Open:** intermittent SIGILL under load. A TimeConstrained async jump lands inside a
  non-GMP malloc (`symmemo_lookup`, Expr nodes). It is pre-existing and occurs at the
  same rate on HEAD. The fix belongs to a TC design pass (all-cooperative jumps, or
  guarded allocation of Expr nodes).
- **Residue 38:** Abel 2nd type, Riccati/Chini with symbolic exponents, elliptic
  separables, symmetry classes, `y = G(x,y')`, lincoeff `F(L1/L2)`, the `u = E^(2y)`
  substitution.
