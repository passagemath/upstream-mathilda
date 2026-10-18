---
references:
  - "D. Gruntz, *On Computing Limits in a Symbolic Manipulation System*, PhD thesis, ETH Zürich, 1996."
source:
  - src/calculus/limit.c
  - src/calculus/gruntz.c
  - src/calculus/limit_osc.c
---
**Algorithm.** `builtin_limit` first rationalises any inexact coefficients,
because the symbolic machinery only handles rational coefficients, and
numericalises the value afterwards. It then mutes the transient
`Power::infy`/`Infinity::indet` warnings its probes would raise. Next,
`builtin_limit_impl` normalises the three calling forms:

- `Limit[f, x -> a]`
- the nested `Limit[f, {x1 -> a1, ...}]`, innermost limit last
- the joint `Limit[f, {x, y} -> {a, b}]`, which uses a polar substitution and
  directional probes

The same step parses three options:

- `Direction`, which takes `"FromAbove"`/`"FromBelow"`/`±1`, `Reals`,
  `Complexes`, or a complex direction such as `I`.
- `Method`.
- `Assumptions`, built into an `AssumeCtx` by `limit_effective_assumptions`.
  It falls back to the ambient `$Assumptions` or `Assuming[...]`.

The core is `compute_limit`, run under a `LimitCtx {x, point, direction, depth,
method, assume}`. It begins with pre-passes:

- One-sided step functions at integers (`Floor`, `Ceiling`, `Round`,
  `IntegerPart`, `FractionalPart`), bounded opaque factors times vanishing
  factors, and `Mod[a, x]` at infinity.
- A refusal on undefined heads `f[...x...]`.
- A reciprocal-trig rewrite.
- At `±Infinity`, a hyperbolic→exponential rewrite plus `Expand`.

Then the layer cascade runs in this fixed order. The first layer that returns a
value wins.

1. `layer_param_power`: assumption-gated `x^n` and `a^x`.
2. `layer_special_fn_asymptotic`. This layer is the Gamma firewall. It handles
   Gamma ratios and Stirling's formula, Hurwitz `Zeta` with two-term tails,
   `BesselJ`/`BesselY` at small argument, `Erfc`, and incomplete `Gamma`. It runs
   before substitution can produce husks such as `BesselJ[a, 0] 0^(-a)`.
3. `layer1_fast_paths`: numeric and continuous substitution through `Together`.
4. `layer_abs_rewrite`.
5. `layer_compose_at_infinity`: a head's own value at infinity, as in
   `Erf[Infinity] = 1`.
6. `layer_exp_decay_re`: the sign of `Re[rate]` under the assumptions.
7. `layer3_rational`: leading-degree comparison for `P/Q`.
8. The Log layers: log of a finite inner limit, the Gruntz-lite dominant-summand
   `Log[sum]`, and Log+linear merge.
9. Term-wise `Plus`.
10. The `f^g` → `Exp[g Log f]` reduction.
11. The bounded layers: bounded envelope (squeeze), assumption-pinned dominant
    products, and `Max`/`Min` against a bounded oscillation.
12. `layer_oscillatory`, from `limit_osc.c`. It rewrites a sum at `±Infinity`
    into a normal form `c0 + Σ cj E^(I θj)` with pairwise-distinct phases. A sum
    with no dominant phase is proved to have no limit, and the layer returns
    `Indeterminate`.
13. `layer2_series`, the workhorse. It calls `Series` at orders 4 → 32 and reads
    the leading term.
14. `layer_gruntz`, from `gruntz.c`. This is a port of Gruntz's mrv algorithm
    with the structure of SymPy's `gruntz.py`: `mrv` sets, `compare`, `rewrite`
    to `w → 0+`, `mrv_leadterm` and a `sign` oracle. It reduces every point and
    direction to `x → +∞`. Its budgets are a depth of 80, series orders 2–16 and
    6000 `mrv_leadterm` calls. If it hits any cap, or an oracle is undecided, the
    limit is left unevaluated.
15. `layer5_lhospital`, with at most 12 iterations and a leaf-growth cap.
16. `layer6_bounded`, which returns an `Interval`.
17. Late fallbacks, in this order: a convergent-summand `Plus` split,
    atom substitution for a `b^e(x)` subterm, a parametric-monomial substitution
    `t = x^a`, small-argument linearisation `H[g] ~ g`, constant-factor
    linearity, and a one-sided-disagreement probe that returns `Indeterminate`.

**Method dispatch.** The `TRY(group, layer)` macro tags every layer with one of
eight groups: `"Substitution"`, `"RationalFunction"`, `"Asymptotic"`,
`"Bounded"`, `"Oscillatory"`, `"Series"`, `"Gruntz"`, `"LHospital"`. The default
`Automatic` runs all of them. A named `Method` restricts **only the outermost
call** (`depth == 1`), so recursive sub-limits always see the full cascade. A
named method that does not apply leaves the `Limit` unevaluated. Each method is
also exposed as a head of its own. For example, ``Limit`Gruntz[f, x -> a]`` is
the same as `Limit[f, x -> a, Method -> "Gruntz"]`.

**Assumptions** (v0.366). With no informative assumption the context is `NULL`,
and the legacy path runs byte-for-byte. Under assumptions, the following layers
switch on:

- the signs of symbolic exponents, rates and coefficients (`literal_sign_ctx`,
  `exponent_sign_ctx`)
- dominant products: `x^p Log[x] → 0` and vanishing×bounded
- `x^a → 0+` for `a < 0`
- `Re[w] > 0` decay
- symbolic dominant-log summands
- an `Abs`-of-complex pre-pass through `ComplexExpand` plus `Refine`

**Data structures.** Subexpressions are `Expr*` trees. Most of the analysis
calls the evaluator symbolically (`Series`, `D`, `Together`, `Expand`,
`ComplexExpand`) rather than through direct C calls. Gruntz keeps a `SubsSet`
of mrv elements mapped to fresh dummies, plus an 8192-bucket memo.

**Complexity / limits.** Recursion is capped at `LIMIT_MAX_DEPTH = 24`. `Series`
covers most `0/0` and `∞/∞` forms. Gruntz covers nested exp-log cancellation.
L'Hospital is kept for shapes `Series` cannot expand. Discontinuous and
undefined heads are refused, never evaluated from one side. A bounded
oscillation such as `Sin[x]` at infinity gives `Indeterminate`.
