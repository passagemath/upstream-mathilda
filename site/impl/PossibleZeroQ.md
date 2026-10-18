---
references:
  - "J. T. Schwartz, \"Fast probabilistic algorithms for verification of polynomial identities\", JACM 27 (1980)."
  - "R. Zippel, \"Probabilistic algorithms for sparse polynomials\", EUROSAM 1979."
source: src/zero_test.c
---
**Algorithm.** `builtin_possible_zero_q` (`src/zero_test.c`) threads a `List`
first argument by hand and broadcasts the remaining arguments unchanged, so a
list of assumptions is never threaded against the expressions. The head is
not `Listable`. Next it forms the effective assumption: an
`Assumptions -> X` option overrides `$Assumptions`, and a positional
assumption is `And`-ed with it. With no usable facts it calls
`zero_test_decide`; otherwise it builds an `AssumeCtx` and calls
`zero_test_decide_assuming`. Both entries first apply
`zt_normalize_exp_kernels`, which combines exponentials so that
Gaussian × `Erf` style products cannot hang the later stages. They then run
a staged pipeline that exits on the first definite verdict:

- *Stage 0, structural* (`decide_structural`): literal zeros, `Complex[0,0]`
  and nonzero constants are settled at once. `provably_nonzero` then
  certifies shapes such as `c E^f`, which can never vanish, so `-E^(-a x)`
  is `False` even though sampling sees only underflow.
- *Guards.* A free symbol mixed with an algebraic constant, or a symbolic
  exponent, skips straight to sampling, because `Together ∘ Cancel` over
  `Q(α)` can fail to terminate.
- *Stage 1, rational* (`decide_rational`): `Together ∘ Cancel` and `Expand`
  decide identities in `Q(x_1, …, x_n)`. `False` is trusted only on the pure
  rational-function path.
- *Stage 1b, exact trig/exp* (`trigexp_rational_is_zero`): a rational
  function of one kernel `t = E^(i x)` is decided by exact evaluation on a
  Nullstellensatz grid. A constant phase `k x + c Pi` is angle-expanded
  first (`expand_pi_phase`), so `Tan[Pi/2 - x] - Cot[x]` comes out `True`
  rather than a spurious `False`.
- *Stage 2, numeric* (`decide_numeric`): a symbol-free input is evaluated at
  machine precision and checked against a cancellation threshold, then up an
  MPFR ladder (53 → 200 → 500 → 1000 bits) while the verdict is ambiguous.
  Once the operand scale passes `ZT_HUGE_SCALE_BITS` (2^40), the machine
  rejection is not trusted and `decide_numeric_huge_scale` reads the trend
  of the residual across the ladder: a plateau means nonzero, shrinking to
  the top rung means zero.
- *Stage 3, Schwartz–Zippel* (`decide_schwartz_zippel_core`): each free
  symbol is replaced by a random real point with |v| ≥ 1, or by a point
  drawn from a per-symbol `SampleSpec` (real, integer or complex channels,
  with bounds) when assumptions apply. `ZT_SCREEN_SAMPLES` = 24 screening
  points are followed by `ZT_CONFIRM_SAMPLES` = 4 confirmations. IEEE
  overflow and numeric-function declines
  (`residue_is_numeric_decline`: `PolyLog[2, z]` returning itself for a
  large argument) are re-drawn on smaller magnitude shells, not
  treated as undecidable. An inequality that couples two or more symbols
  (`x > y`) is enforced by rejection sampling (`assignment_conforms`, up to
  `ZT_CONFORM_MAX_TRIES` = 128 tries), so every test point lies inside the
  region. A measure-zero fact (`a == b`, `Or`) cannot be sampled into, so
  under such a fact a `False` is downgraded to unknown.

Under assumptions, only `True` from the exact stages is trusted: an
unconditional `False` can be wrong on an integer lattice, as with
`E^(2 Pi I k) - 1`. Constrained sampling alone settles `False`.

**Result mapping.** `False` is returned only for a definite `False`. Both
`True` and unknown give `True`, the documented "assume zero when
uncertain" behaviour, and no `PossibleZeroQ::ztest1` message is emitted.
The attribute set is just `Protected`.

**Complexity / limits.** The cost is a few symbolic normalisations plus at
most about 28 sampled evaluations, with MPFR only on ambiguous rungs.
Without assumptions the free symbols are sampled on the real line. A
function that is not analytic in its argument (`Conjugate`, `Re`, `Im`)
therefore looks like the identity there: `PossibleZeroQ[Conjugate[x] - x]`
returns `True`. Supply an assumption such as `Element[x, Complexes]`, or
expand with `ComplexExpand`, before relying on that verdict.
