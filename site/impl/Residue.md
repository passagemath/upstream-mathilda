---
source: src/calculus/residue.c
---
**Algorithm.** `builtin_residue` handles only the two-argument form
`Residue[f, {z, z0}]` (fewer arguments emit `Residue::argm`; a non-symbol `z`
or a malformed spec leaves the call unevaluated) and hands it to
`residue_compute`, which reads the coefficient of `(z - z0)^-1` out of a Laurent
expansion. The integrand is first classified by `residue_is_rational_in` (is
`Together[f]` a ratio of polynomials in `z`?). For a rational `f`,
`residue_simple_pole` tries the fast path `P(z0)/Q'(z0)` (`P/Q = Together[f]`,
`Q(z0) == 0` and `Q'(z0) != 0` decided by `PossibleZeroQ`), which never builds a
series. Otherwise `residue_shift_form` substitutes `z -> z0 + w` (a private
symbol ``Residue`$w``) and `Expand`s the shifted numerator and denominator
*separately*, so `Q`'s constant term collapses to a literal `0` and an algebraic
`z0` keeps one radical spelling. Transcendental integrands (`Cot`, `Zeta`,
unknown `f[z]/z^n`) are expanded directly about `z0`, so the series engine's own
Laurent knowledge applies.

In both cases the coefficient comes from `residue_extract`, which is
*self-validating* since v0.268. An order-0 probe (`residue_coeff_once`) learns
the pole order `m` and short-circuits the analytic (residue `0`) and branch-point
(`den > 1`, unevaluated) cases. `residue_coeff_resolved` then expands from order
`m - 1`, raising the order by 3 while the O-term still hides exponent `-1`.
Successive orders are compared (`residue_coeffs_agree`, a `PossibleZeroQ` of the
difference), and a coefficient is accepted only when two consecutive orders
agree. A padded `Series` can return a `-1` coefficient that is present but
inaccurate at a higher-order pole: the cofactor is truncated, so product-rule
cross terms drop out, and this happened at symbolic double poles such as
`E^(I z)/(z^2 + a^2)^2` at `I a`. The agreement test catches that case. A
second guard covers a *dropped pole*. When the direct expansion reports `0` but
`residue_has_pole_at` finds that `Denominator[Together[f]]` vanishes at `z0`,
the function re-expands `f /. z -> z0 + w` about `w = 0`. There `z^p` becomes the
analytic binomial `(z0 + w)^p`, which fixes the silent `0` previously returned
for `z^(1/3)/(1 + z^2)^2` at `z = I`.

**Data structures.** `Expr` trees throughout. The helpers `residue_eval1`/`residue_eval2`
build and evaluate calls to `Together`, `Numerator`, `Denominator`, `Expand`,
`D`, `ReplaceAll` and `Series`. The expansion is a `SeriesData[z, z0, {coefs},
nmin, nmax, den]`, and the residue is the coefficient at index `-1 - nmin`.

**Complexity / limits.** The simple-pole path costs one differentiation and two
substitutions. The series path costs a few Laurent expansions. The resolve loop
gives up past order 256, and the agreement loop is bounded at order 64 / 12
iterations, returning the last estimate. The residue is defined only for an
ordinary Laurent expansion, so a Puiseux expansion (a branch point, as in
`1/Sqrt[z]` at `0`) is left unevaluated. So is any point where `Series` itself has
no expansion. The essential singularity `E^(1/z)` at `0` and the poles of `Gamma`
at the non-positive integers fall in this class (`Series[Gamma[z], {z, 0, n}]`
is itself unevaluated). For those cases, use the numerical `NResidue`.
