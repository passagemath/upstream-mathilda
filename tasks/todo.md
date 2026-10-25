# FeynCalc H1: named characters, test harness, first baseline (2026-10-24)

Plan: ~/.claude/plans/let-s-continue-our-feyncalc-humble-karp.md. Gap doc: packages/feyncalc/FEYNCALC_GAPS.md (local).

- [x] P1 named characters in source: WL classification table, lexer (letters/operators/constants/space), tests, docs ; v0.587 (full-suite A/B vs HEAD clean; match_stress_corpus 2000-deep case overflows only the unoptimised main test build)
- [x] P2 harness: tools/feyncalc/run_fc_tests.m, tools/check_feyncalc.py, make check-feyncalc
- [x] T2 Check/Quiet nesting, Quiet/Check message specs, named emission of legacy helpers ; v0.588 (A/B Release vs Release: one regression found and fixed -- Inverse::sing via an anonymous helper)
- [x] T3 assignment LHS pattern-free parts evaluated (h[y_ - 4]) ; v0.589 (Contract 67 -> 92 / 119)
- [x] T7 Cycles/PermutationGroup/GroupElements/Permute over a group ; v0.590
- [x] P3 first baseline (11362 pass, ~85% of runnable) + three triage agents (NCA+Tables, Dirac+Feynman, Shared+LoopIntegrals) -> FEYNCALC_GAPS.md §7
- [x] Batch 1 ; v0.591: -(x+y) sign distribution + flat-chain parse provenance (+ Root sign, negative coefficient under radicals, FactorList sign, series/FunctionExpand cancelling pairs), HoldPattern LHS, ::write, N-values, indexed iterators, valued pattern heads, OptionsPattern guard/binding, 2ep lexer. Full A/B Release vs Release clean except pre-existing Fit WorkingPrecision sequence bug
- [ ] Follow-up (pre-existing): Fit[..., WorkingPrecision -> 24] returns the exact fit after a WorkingPrecision -> Infinity call in the same session
- [ ] Follow-up: positive non-reducing coefficient under radicals (Sqrt[12 a] -> 2 Sqrt[3] Sqrt[a] in the reference)
- [x] WL reference run of the harness (--per-dir --binary wolframscript): 13412/13412 runnable pass, 235 s -- the harness is faithful; per-dir is the comparison mode (FeynCalc caches across tests)
- [x] Batch 2 ; v0.592: Sequence splice only at pattern-variable sites, nested option lists (OptionValue/OptionsPattern; fixed the LeviCivita loop v0.591 exposed), Coefficient negative orders / denominator terms / expand only in the form; TimeConstrained unwinds only at a safe point (grace for the evaluator poll, then a frame-pointer check that no frame is inside libc) -- the Lorentz directory died with SIGILL after a timeout, losing 880 tests (now 1738/1812)
- [ ] O1 canonical order (src/sort.c expr_compare rewrite, worktree o1wt): numbers by value (complex: re, |im|, sign), Times > Power > Plus compared as polynomial terms from the last factor, atoms before compounds, strings before symbols, symbol collation, head/length/args. Random corpus 400/400, held-out 567/569 (2 misses are evaluation differences). Full-suite A/B and expectation triage against the reference pending
- [ ] Follow-up: DummyIndexFreeQ on a 50000-term sum: 25x the reference -- re-evaluation churn (Block bumps the eval clock, the big Plus is re-walked per call) and matcher speed
- [ ] Follow-up: Infinity / ComplexInfinity / -Infinity are symbols / Times[-1, Infinity] here, DirectedInfinity[...] in the reference
- [ ] Follow-up (pre-existing, wrong answer): PossibleZeroQ[D[F[x+1],x] - Cosh[x+1] Cosh[2x+2] + 10^-2] is True (reference False; F the Weierstrass antiderivative of Cosh[x] Cosh[2x]) -- every sample point lands at large |x|, where the operand scale (~10^40) hides the constant
- [x] v0.592 per-dir FeynCalc run: 12248/13412 (91.3%; Lorentz 877 -> 1738 after the TimeConstrained fix); baseline rewritten from it (per-dir mode)
- [ ] Follow-up: 32 ids pass per-file but fail per-dir (state left by earlier tests; reference passes them per-dir) -- mostly DiracTrick*FreeIndices-ID8, GluonVertex-ID1, FCFeynmanPrepare-ID51/53 (list: scratchpad fc_lost_v592.txt)
- [x] v0.593 PossibleZeroQ exact-point ladder; v0.594 Times numeric-base power merge; v0.595 assumption-rule pattern capture (Refine[1/Sqrt[r^2], r > 0] was -2)
- [x] v0.596 PossibleZeroQ climbs every rung on exact huge-scale points; v0.597 canonical order (O1) + Solve/SolveAlways ordering + Meijer-G convolution side + Simplify power-distribution polish
- [ ] O2 power normal forms (prototype in scratchpad o1fix: power.c numeric-factor extraction for any non-integer exponent, (1/d)^e -> d^-e; times.c rational-coefficient absorption 12 2^n -> 3 2^(2+n)). All match the reference, but they break recognisers that read the old forms: Simplify half-angle on Sqrt[(1 - Cos[x])/2] (now Sqrt[2 - 2 Cos[x]]/2), infinite Product families on (1/3)^(2^k), Gamma duplication Simplify, a Mellin Bessel integral, Goursat/Residue integrals. Needs those recognisers updated with it
- [ ] Follow-up: the symbolic-n Struve ODE identity takes ~22 s (1.2 s before the canonical order) -- the special-function engine families churn Factor/Together over an expansion
- [ ] Follow-up: numeric_sign (simp_assume.c) reads only Integer/BigInt/Real, so x > 1/2 does not imply x > 0 in fact_implies_strict_positive
- [ ] Follow-up: Sqrt[(2 I) x] -> Sqrt[2] Sqrt[I x] in the reference (complex coefficient)
- [x] v0.597 per-dir FeynCalc run: 12644/13412 (94.3%); 405 gained, 9 lost (DiracOrder-ID54 is dummy-index relabelling of an equal answer; FCFeynmanParametrize-ID28..30 and FCToTeXPreviewTermOrder-ID3 need the Times-level power split below; TID-ID42/44/65 and SetMandelstam-ID7 not yet traced)
- [x] v0.598 Multinomial + Binomial negative/half-integer closed forms (PaVeUVPart 137/139); v0.599 FormatValues (SharedObjectsUpValues 92/97); v0.600 Return from a definition (13 FCLoopBasis aborts); v0.6001 linear box syntax (48 ToExpression::sntx)
- [x] v0.6002 Infinity ordering, v0.6003 PolynomialGCD/Cancel over compound variables (x[i], Pair[..]), v0.6004 NValues + N[x, {p, a}] + N not Listable
- [x] v0.6004 per-dir FeynCalc run: 12951/13412 (96.6%); 317 gained, 10 lost vs v0.597: 6 FCLoopIntegralToGraph (per-dir state; per-file v0.6004 fails 8 where v0.597 fails 12), 4 timeouts under A/B load (DiracTrace-ID29, TdecGrc-ID3, TID-ID13/14 -- being timed)
- [x] v0.6005 O3 power-of-product split + Mellin fullsimp2 shifts only partnered Gammas
- [ ] v0.6006 exponent generators (Variables/Factor/FactorList over b^t, expsplit.c), Simplify exponent/argument normal form + cross-term coprime basis, O4 Times coprime basis for numeric powers -- combined full A/B running
- [ ] v0.6009 number atoms (rules/Cases/Count/Position/DeleteCases never enter Rational/Complex; FreeQ/MemberQ see their head) -- A/B running; v0.6010 DownValue order ignores conditions (FromGFAD 28), v0.6011 OptionsPattern in any position (SUNSimplify 10), v0.6012 SyntacticNegativeQ sums/complexes (SUNSimplify 2) -- written, A/B next
- [ ] Matcher: in a Flat+Orderless head, rest_ s[a_] binds the specific pattern to the EARLIEST factors in the reference (Replace[s[1] s[2] s[3], rest_ s[a_] :> {rest, a}] is {s[2] s[3], 1}; here {s[1] s[2], 3}) -- SUNSimplify ID83-86
- [x] DownValue order: a later repeated-name pattern p[a_, a_] goes ahead of p[a1_, b1_] defined earlier, as in the reference (v0.6022, generality order)
- [ ] Simplify[(18 - 6 n^2 + n^4)/(96 n^2)] keeps the input in the reference (here the expanded 1/96 (-6 + 18/n^2 + n^2)) -- SUNSimplify ID23
- [ ] Risch real-trig reconstruction (rt_realify): 1/(3 + 5 Cos[x]) (real zeros) comes back as complex-argument logs; the test's FreeQ[g, Complex] was vacuous until v0.6009
- [x] v0.6013 D of a compound head + Factor num/den via FLINT; v0.6014 ReleaseHold in heads; v0.6015 ToExpression non-string; v0.6016 Times O4 loop fix + merged-exponent sign (last term, complex) + numeric opposite exponents; v0.6017 Together/Cancel indexed variables (29 s -> 0.02 s); v0.6018 Factor Gaussian content (DiracTrace 53/61); v0.6019 contagion only numeric summands; v0.6020 pattern-headed call single operand in Flat -- all pushed 2026-10-26; FeynCalc per-dir measure pending
- [x] v0.6020 per-dir FeynCalc run: 13132/13412 (97.9%); 165 gained, 12 lost -- 9 FCFeynmanPrepare (v0.6009 FreeQ-on-Rational let Factor2 send 2 - D/2 inside Gamma to Together, which rewrote kernel arguments: fixed v0.6021, Together kernels opaque as in the reference), DiracTrickInsideDiracTrace-ID17/29 (v0.6010 exposed score-based DownValue order: rule order by pattern generality, v0.6022), FCFeynmanFindDivergences-ID24 timeout (not yet traced)
- [x] Perf: memoised definitions were quadratic (f[n_] := f[n] = n over 20000 values ~38 s); per-symbol hash index of the pattern-free prefix -> 0.024 s ; v0.6023 (targeted rule suites + reference diff of Unset/Clear/Block/DownValues=/conditional literals; full-suite A/B still to run)
- [ ] v0.6022 per-dir FeynCalc run in progress (scratchpad fc6022): Dirac 1881 (+7), Shared 4621 (+5), SUN 162 (+2), Pauli 225 (+1) so far; then --update-baseline (local)
- [ ] FCFeynmanFindDivergences-ID24: 15 s alone on v0.6008 and v0.6022 (passes); 120 s TIMEOUT only in the per-dir run -- per-dir state or load, not a kernel regression
- [ ] DownValue redefinition: the reference also replaces equivalent-but-not-alpha-equal patterns (x_ then _ ; x_h then h[___]; x:(_Integer|_Real) then x:(_Real|_Integer)) -- one rule there, two here
- [ ] Apart: splits over a linear factor with radical or I constants only when they are symbolic -- Apart[1/((x - Sqrt[2])(x - Sqrt[3])), x], Apart[1/((x - (-1)^(1/3))(x + (-1)^(1/3))), x], Apart[1/((-i Sqrt[3] + x)(i Sqrt[3] + x)), x] are recombined here, split in the reference; Apart[g[x]/((x - 1)(x + 1)), x] (non-polynomial numerator) declines here; Apart[(z^2 - 1)/z, y] (variable absent) splits here, the reference returns the factored Together -- ToFCPartialFractionForm 6
- [ ] Solve[1 + z^2 + z^4 == 0, z] gives nested Sqrt[1/2 (-1 - I Sqrt[3])] here, (-1)^(1/3)-style roots in the reference; Roots[...] unevaluated here
- [ ] (2026-10-26 finding) the reference's multivariate factor sign is NOT a monomial order: 2-variable cases are lex (x > y), but b - a c / a - b c / a b - c and b c - a d all keep the term without the alphabetically largest variable positive, and no lex/grlex/grevlex permutation or numeric-leading-coefficient heuristic beats plain lex on 559 random factors (493/559) -- an artifact of its factoring algorithm; deprioritised (3 DiracTrace tests). Unit content does differ simply: Factor[4 I a - 4 I b] is (4 I) (a - b), Factor[-I a + I b] is -I (a - b)
- [ ] Factor over Z[i] when the input has complex coefficients: Factor[I x^2 + I] is (1 + I x) (I + x) in the reference (I (1 + x^2) here); unit choice for general complex content ((-4 + 2 I) a + (6 - 2 I) b is (2 I) ((1 + 2 I) a - (1 + 3 I) b) there); Factor[(I x + I)^2] gives -(-1 - x)^2 (pre-existing); Factor of a fraction with Gaussian content (num factored internally, not via builtin)
- [ ] Simplify's RadicalCanon re-dispatch (simp_search.c ~1308, `<=` gate) recurses without bound whenever radical_canon and Times disagree on a canonical form (M11338 stack overflow, fixed at the Times side in v0.6016) -- consider a cycle guard
- [ ] Plus/Times: none open from contagion; Flat ReplaceList singleton wraps (fl[a]) still differ (documented)
- [ ] Simplify tie: (-3 + d)(-2 + d) and 6 - 5 d + d^2 both score 9; update_best keeps the first (strict <); the reference returns the factored one (LeafCount 7 vs 8?) -- unverified tiebreak, wide blast radius (DiracTrace ID71)
- [ ] Factor sign/unit convention: Factor[8 b c - 8 a d] is -8 (-b c + a d) here, 8 (b c - a d) in the reference (FactorList already agrees); I-content not pulled: Factor[4 I a - 4 I b] is (4 I) (a - b) there -- DiracTrace ID45/53/61
- [ ] Factor with Gaussian-integer coefficients is very slow multivariate (18 vars: 8 s vs 0 s) -- DiracTrace ID46/29 timeouts
- [ ] (f -> Times sub-product part fixed v0.6020; remaining: the reference enumerates subject SUBSETS in order and assigns pattern elements within each, Mathilda is pattern-element-major, so q[x_, z_] orientation is retried per f there, not here) Matcher: an UpValue on a Times pattern q[z_, x_] f_[y__] (q Orderless) binds f_ to a sub-product of the Flat Times and tries orientations in a different order (q[a, b] e[a] g[b]: e[a] r[b, g] here, g[b] r[a, e] there); condition evaluations 2/16/48/320 for n = 1/2/3/5 vs 1/4/9/25, plus a second full pass -- DiracTrace ID70 and ID29 (PairContract)
- [ ] Together of a sum with kernel factors does not combine: Together[f[a] (-p - q) + f[a] (p + q)] is 0 in the reference (unchanged here; with a symbol it works); the symbol path in turn expands a single product (Together[x (y + 1)] is x + x y here, x (1 + y) there) -- see feedback_together_no_expand: wide blast radius -- DiracTrace ID61#2
- [ ] Simplify should prefer the factored form: Simplify[(2 - d)(3 - d)] is (-3 + d)(-2 + d) in the reference (6 - 5 d + d^2 here) -- DiracTrace ID71
- [ ] AbsoluteTime[] has 1 s resolution; Unique[{"a","b","c"}] gives all three the same number
- [ ] Factor of a large rational expression is slow (FactorList2-ID10: Factor[q] 51 s while Together 1 s + numerator/denominator factoring 0.4 s): the classical exact_poly_div path in builtin_factor's rational branch; FactorList2 ID5/ID10 time out
- [ ] Coefficient[x^(1+e) + 3 x^(2+e), x] is 0 in the reference (x^e here)
- [ ] FactorList order with a radical generator: {Sqrt[x+y], 1} before {(x+y)^e, 1} in the reference
- [x] O3 Times-level power split (v0.6005), measured against the reference: in a product, Power[Times[c, rest], e] with c a positive number splits into c^e rest^e (y (2x)^b is 2^b x^b y; Pi^a (2 Pi)^b is 2^b Pi^(a+b); (-2x)^b x^a is 2^b (-x)^b x^a), while a lone (2x)^b stays. Needed by FCFeynmanParametrize (Pi^(D/2) (2 Pi)^-D) and FCFeynmanFindDivergences
- [x] Variables/Factor/FactorList split b^(n + k e) into b^n (b^e)^k (FactorList[(x+y)^(-3+3e)] is {{1,1},{(x+y)^e,3},{x+y,-3}} in the reference) -- FactorList2 ID7/8/11, FCFeynmanFindDivergences
- [x] Simplify[a^(3(-1+e)) - a^(-3+3e)] is 0 in the reference (exponent normal form); Simplify[la^2 (la^2)^(e-1), la > 0] leaves la^(2 + 2(-1+e)) (assumption power rule result not re-simplified) -- FactorList2 aborts, FCFeynmanProjectiveQ
- [ ] Together does not pull out integer content: Together[-2 a x6 + 2 b x8] is -2 (a x6 - b x8) in the reference
- [ ] Indeterminate < 1 is False here (IEEE-unordered, deliberate) but unevaluated in the reference
- [ ] Parser: Hold[_Symbol - 1] keeps Times[-1, 1] (the reference folds the literal to -1)
- [ ] P4 acceptance loop: fix largest root-cause groups in the kernel, ratchet

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
