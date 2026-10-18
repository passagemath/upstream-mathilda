---
references:
  - "Joel S. Cohen, *Computer Algebra and Symbolic Computation: Mathematical Methods* (A K Peters, 2003)."
  - "NIST Digital Library of Mathematical Functions (DLMF), §5.11, §6.6, §6.12, §7.12 — asymptotic and logarithmic expansions used by the special-function emitters."
source: src/calculus/series.c
---
**Algorithm.** `Series` computes a truncated power series by recursive *series
algebra*. Its only attribute is `Protected`. The arguments evaluate normally,
because the stray `HoldAll` was removed in v0.236. `builtin_series` first checks
that every expansion variable is a symbol (`series_spec_vars_ok`; otherwise
`Series::ivar`). Inexact inputs are rationalised and the result is numericalised
afterwards (`internal_rationalize_then_numericalize`). A `List` first argument
is threaded. An `Assumptions -> assm` option is pulled out of the trailing specs
by its LHS symbol, so it is never mistaken for a leading-term spec `x -> x0`.
`series_effective_assumptions` turns it into an `AssumeCtx`, falling back to the
ambient `$Assumptions` or `Assuming[...]`. With several specs the leftmost
variable is expanded first. Each coefficient is then re-expanded in the remaining
variables, and the assumptions are forwarded.

`do_series_single` does the work for one variable. A body free of `x` is returned
verbatim. `WhittakerM`/`WhittakerW` and `MeijerG`/`BesselK` at `0` or `Infinity`
are rewritten to closed or Slater forms, which the generic path can expand. An
expansion at `-Infinity` is computed *by reflection*: the function expands
`f /. x -> -x` at `+Infinity` and flips the odd-exponent coefficients. At
`+Infinity`, a bank of dedicated emitters (`try_series_*_at_infinity`) supplies
asymptotic expansions that have no Laurent series. The essential factor is kept
symbolic. Covered heads: `ExpIntegralEi`, `ExpIntegralE[1, x]`, `Erf`/`Erfc`/`Erfi`,
`LogGamma` (Stirling), `PolyGamma`, `Zeta` (a truncated Dirichlet head),
`SinIntegral`/`CosIntegral`/`Sinh-`/`CoshIntegral`, `FresnelS`/`FresnelC`, the Airy
functions, `BesselJ`/`Y`/`I`/`K`, `HermiteH` and `ProductLog`. A scaled argument
`c x` with the sign of `c` known from the assumptions selects the branch
(`ArcTan[c x] -> Sign[c] Pi/2 - ...`) or the saturating constant (`Erf`, `Tanh`,
`Fresnel*`, `Erfc`). `ArcCosh`/`ArcSinh` of `c x` with `c < 0` reflect to the
principal branch. At `0` the logarithmic expansions of `ExpIntegralEi`,
`LogIntegral` (a generalized series in `1/Log[x]`), `CosIntegral`, `CoshIntegral`,
`BesselK`/`BesselY` and `LerchPhi` are emitted directly. When the assumptions
prove `x < 0`, these emitters use the `Log[-x]` branch.

On the generic path a rational function in `x` is first `Apart`-decomposed, and a
symbolic `x^alpha` prefactor is peeled off at `0` or `Infinity`. At `Infinity`
the engine then substitutes `x -> 1/u`. The body is expanded to a padded internal
order (`order + 12` at a numeric point, `+ 2` at a symbolic one, so that symbolic
convolution cannot explode) and truncated back to the requested order. The core
routine `series_expand` combines child series with `so_add`, `so_mul`, `so_inv`,
`so_div`, `so_pow_int` and `so_pow_expr`. The binomial `(a + b x^m)^alpha` has a
monomial fast path. Elementary and inverse heads compose known Taylor kernels
(`kernel_coefs` → `so_compose_scalar_kernel`) with the inner series. Reciprocal
heads are rewritten first (`rewrite_reciprocal_head`). The inverse trig and
hyperbolic heads have MMA-faithful branch-point handlers: Puiseux for
`ArcSin`/`ArcCos`/`ArcSinh`/`ArcCosh`, and an explicit `Log[x - x0]` term for
`ArcTan`/`ArcCot`/`ArcTanh`/`ArcCoth`. Both carry the
`(-1)^Floor[(Pi/2 - Arg[x - x0])/(2 Pi)]` discriminator. A head that cannot be
expanded otherwise falls back to `series_taylor_via_D`, which computes
`D^k f(x0)/k!` for up to 20 iterations and gives up when a derivative is infinite.

Assumptions (v0.368) act in four places:

- They give the sign of `x` for the log-branch emitters.
- They resolve the non-analytic heads `Abs`, `Sign`, `UnitStep`, `Sqrt[g^2]` and
  `Conjugate` by the sign of their argument at `x0`, decided through the
  Reduce/CAD bridge.
- They flip the even-order-root sign under `x < 0` (`Sqrt[1 - Cos[x]] -> -x/Sqrt[2]`).
- They clean up the final coefficients with `apply_assumption_rules`, for
  example `Sqrt[a^2] -> a` and `Sqrt[a] -> I Sqrt[-a]` for `a < 0`. No
  `Simplify` search is involved.

**Data structures.** The internal `SeriesObj` holds `x`, `x0`, an owned array of
coefficient `Expr*`, the leading exponent numerator `nmin`, the exclusive O-term
numerator `order`, and a common exponent denominator `den`. Coefficient `i`
multiplies `(x - x0)^((nmin+i)/den)`, so Laurent and Puiseux series share one
representation. `so_rescale`/`so_align_den` reconcile denominators before
arithmetic. `so_to_expr` emits `SeriesData[x, x0, {coefs}, nmin, nmax, den]`. At
infinity the variable is `1/x`. A symbolic prefactor (`E^x`, `x^a`) or an
additive growth head (`(x - 1/2) Log[x] - x + ...`) rides outside the
`SeriesData` at the expression level. Coefficient arithmetic promotes to
GMP-backed rationals on overflow.

**Complexity / limits.** Kernel-path heads are exact and fast. Convolution is
`O(N^2)` in the padded order. The D-path is capped at 20 derivatives and gives
up at singular points. For example, `Series[Gamma[x], {x, 0, n}]` stays
unevaluated, as does `Gamma` at any of its poles, and `Series[Gamma[x], {x,
Infinity, n}]` also stays unevaluated because `Gamma` has no asymptotic series
there. The general non-terminating `WhittakerM`/`WhittakerW` asymptotics are
deferred. A Puiseux base at a genuine branch point is handled only where an
explicit handler exists. This is a symbolic routine with no packed, NDArray or
`Compile[]` surface.
