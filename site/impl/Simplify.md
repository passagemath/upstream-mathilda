---
references:
  - "Joel S. Cohen, *Computer Algebra and Symbolic Computation: Mathematical Methods* (A K Peters, 2003)."
source:
  - src/simp/simp_builtins.c
  - src/simp/simp_search.c
  - src/simp/simp_assume.c
  - src/simp/simp_assume_rewrite.c
---
**Algorithm.** Simplify is a complexity-weighted, memoized candidate-set search
over the existing battery of algebraic transforms. `builtin_simplify`
(`simp_builtins.c`) parses options — a positional or `Assumptions ->`
assumption (combined into the `$Assumptions` default with `And`), a
`ComplexityFunction`, and `TransformationFunctions` — builds an `AssumeCtx` of
normalised facts, then drives the search. Inexact inputs are rationalised on
entry and numericalised on exit (`internal_rationalize_then_numericalize`).
`Equal`/`Less`/`And`/... heads are threaded manually (with a relational
rebalancing candidate). The head carries only `Protected`. A `List`
argument is simplified as one expression, except inexact input, which
`simplify_thread_list` threads element by element so that numeric
contagion does not cross list boundaries.

The core driver is `simp_bottomup` (`simp_bottomup.c`), which descends into
every `Plus`/`Times`/`Power` child (memoizing each result in a `SimpMemo` hash
table) and dispatches each node through `simp_dispatch` → `simp_search`.
`simp_classify` routes pure-polynomial/rational shapes to dedicated pipelines
(`simp_pipeline_polynomial`/`_rational`/`_logexp`) and there are top-level fast
paths: a `SHAPE_RATIONAL` shortcut that runs `Together`/`Cancel`/`Factor` once
at the top, an algebraic-tower `Together[expr, Extension -> Automatic]` collapse,
and a `simp_trig_rational` substitution that maps trig/opaque subtrees to ground
symbols and works in the quotient ring.

`simp_search` is the heuristic engine. It seeds a `CandSet` with the input plus
the output of a long list of correctness-preserving rewriters (assumption rules,
log/exp identities, `SimpLogRules`, trig roundtrip, ExpToTrig, Pythagorean
square-completion / reduction / canonicalisation, trig-at-rational-Pi, tan-
addition, half-angle, radical product combine, sqrt/cube-root/algebraic
denesting, roots-of-unity, factorial decomposition, per-variable Collect). It
then runs `SIMP_ROUNDS` (= 2) rounds in which every seed is fed through the
`SIMP_TRANSFORMS` table — `Together, Cancel, Expand, ExpandNumerator,
ExpandDenominator, Factor, FactorSquareFree, FactorTerms, Apart, TrigExpand,
TrigFactor, TrigReduce, TrigToExp` — plus chained Pythagorean/radical/trig-Pi
passes. Each transform call is gated by `transform_can_fire` (cheap precondition
check) and the running best is updated by `update_best` against the complexity
score; candidates strictly worse than their parent are dropped from the next
round (with a loosened `2*parent + 8` bound for TrigExpand and the seed-phase
blow-up guard). A short-circuit (`simp_best_is_zero`) exits as soon as a literal
`0` is reached. Final polish passes apply `simp_lift_common_factor`,
`transform_pythag_reduce`, and `canon_negate_pairs`.

After the search, `builtin_simplify` runs a fixed tail of exact folds, each kept
only when it does not raise the score: a constant-algebraic fold
(`flint_qqbar_is_constant_algebraic` / `flint_qqbar_reduce_coeffs`, which treats
trig at any rational multiple of `Pi` as an algebraic number through the root of
unity `exp(i pi r)`, so `Tan[Pi/7] Tan[2 Pi/7] Tan[3 Pi/7]` becomes `Sqrt[7]`);
`simp_const_cpow_fold` (principal value of a constant power with a
non-positive-real base, so `I^I - E^(-Pi/2)` goes to `0`); `simp_invtrig_combo_is_zero` /
`simp_invhyp_combo_is_zero` (a Z-combination of inverse trig/hyperbolic
functions of constants is proved a multiple of `2 Pi` by checking `exp(i e) == 1`
in qqbar, then a numeric screen picks the branch, as in Machin's formula); and, under
assumptions, `simp_invhyp_addition_is_zero` and `simp_region_identity_is_zero`,
which prove an inverse-function identity on a box region from
derivative constancy. Every partial derivative must be exactly 0, every
branch-sensitive subterm is shown by `Reduce` to stay off its cut, and the
value at a `FindInstance` point must be exactly 0. Last comes a confluence
post-pass. When the result still holds a rounding or `Mod` head, it runs
`apply_assumption_rules` to a fixpoint, so `Simplify[Floor[n], n ∈ Integers]`
agrees with `Refine`. A `transform_trigexp_vanish` fast path ahead of the
bottom-up descent proves rational functions of trig/exp kernels identically
zero (`E^(k I x) -> t^k`, exact numerator test), so they never reach the
per-node search.

**Assumption engine.** `AssumeCtx` facts feed two layers. The first is the fast
structural provers (`prov_pos`, `prov_int`, `prov_re`, …) and
`assume_structural_rewrite` (`simp_assume_rewrite.c`). That file holds the
integer-assumption rules: the Boolean-ring parity reduction `simp_reduce_mod2`
and `canon_neg1_pow` route every root of unity `(-1)^Y`, `I^Y` and
`E^(I Pi Y)` to a reduced sign. `trig_pi_shift_symbolic` gives
`Sin[t + k Pi] -> (-1)^k Sin[t]` for a provably-integer `k`. It also holds the
floor-division and `Mod` reductions, `GCD`/`LCM`/`LegendreP[n,-1]` folds,
Gamma/Binomial/Pochhammer ratios via `simp_factorial`, and inverse-of-forward
collapses on a known branch: `invtrig_branch_reduce` takes a two-sided numeric
interval and verifies it lies in a single branch. `invtrig_multiangle_reduce`
emits `ArcCos[1-2x^2] -> 2 ArcSin[x]` only on a `Reduce` proof. The same file
also holds the Schwarz reflection `Conjugate[H[z]] -> H[Conjugate[z]]` for
branch-cut-free `H`. The second layer
is a `Reduce`/CAD entailment oracle (`assume_reduce_entails`, as used by
`Refine`). A fact `P` holds if `Reduce[A && !P]` is `False`. Behind it,
`assume_sign_nonneg_deep` / `assume_sign_nonpos_deep` decide the sign of an
arbitrary subterm, splitting `Tan`/`Sec`/`Csc` into `Sin`/`Cos` signs first. That
is what resolves `Abs[g] -> ±g` and `Sqrt[g^2] -> ±g` (`simp_logexp_abs.c`,
with `sqrt_prepare_radicand` turning `1 + Cos[2 u]` into `2 Cos[u]^2` inside a
radical). Assumption-gated power rules (`transform_powbase_to_exp_assuming`,
`simp_radical_combine_assuming`, `(a^p)^q -> a^(p q)` only when `p` is
provably real) complete the set. Every gated rule declines off its
assumption, and no assumption path uses `PossibleZeroQ`.

**TimeConstraint.** A scalar `TimeConstraint -> t` (or `{tLoc, tTot}`) arms a
whole-call deadline, `g_simp_call_deadline` (`simp_arm_call_deadline`). It is
checked in `simp_dispatch`, the specialised pipelines, `simp_bottomup`, the seed
phase and the round loop. A nested `Simplify` can only tighten the deadline. The
check is synchronous: it cannot interrupt a single running kernel call. The poly
layer's degree-monotonicity guards (`pseudo_rem_standard`) cover that case.

**Complexity scorer.** Candidates are ranked by `score_with_func`: with no user
`ComplexityFunction` this is `simp_default_complexity` (the SimplifyCount metric
— LeafCount with integers contributing their decimal-digit count) plus a
`nested_radical_penalty` (+3 per truly-nested radical). A user function is
evaluated as `f[candidate]`; an Integer result is used directly, a BigInt scores
`SIMP_SCORE_INF`, anything else falls back to the default. Lower wins; ties
favour the form that reached the score plateau first (force-take semantics let
assumption/log/factorial rewrites win even on a tie).

**Data structures.** `CandSet` (dynamic `Expr*` array, dedup via `expr_eq`, capped
at `SIMP_CAND_CAP` = 12); `SimpMemo` (256-bucket chained hash of input→best for
the bottom-up driver); `AssumeCtx` (flat fact array from `assume_ctx_from_expr`);
a per-Simplify `FactorMemo` shared by Factor/Trig* so duplicate subexpression
work hits the cache. Each transform is invoked as `f[candidate]` through the real
evaluator (`traced_call_unary`), so Simplify composes the same builtins exposed
in the REPL.

**Limits.** Simplify is not a complete decision procedure. The `Reduce`
entailment is capped by variable count, and a decline leaves the subterm
alone, so wins under assumptions depend on what the provers and CAD can
settle. Because the search is scored, a correct rewrite that does not lower
`SimplifyCount` can be dropped. On `Pi/2 < x < Pi`, `ArcSin[Sin[x]]` stays
as it is, while `Refine` returns `Pi - x`. Simplify is symbolic and has no
NDArray or `Compile[]` path.
