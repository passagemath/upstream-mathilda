# Residue

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Residue[f, {z, z0}]`**

gives the residue of f at the isolated singularity z = z0 -- the coefficient of (z - z0)^-1 in the Laurent expansion of f.

<details>
<summary>Notes</summary>

Computed by power-series expansion, so a residue is found only where f admits a Laurent series at z0. Returns unevaluated at branch points (fractional-power expansions) and when no series can be produced. See NResidue for a numerical alternative that also handles essential singularities.

</details>

## Examples (19)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (9)

```mathematica
In[1]:= Residue[1/z, {z, 0}]
Out[1]= 1

In[2]:= Residue[1/z^2, {z, 0}]
Out[2]= 0

In[3]:= Residue[1/Sin[z]^5, {z, 0}]
Out[3]= 3/8
```

Order-2 pole

```mathematica
In[4]:= Residue[(z + 1)/(z^2 (z - 2)), {z, 0}]
Out[4]= -3/4
```

Complex pole

```mathematica
In[5]:= Residue[1/(z^2 + 1), {z, I}]
Out[5]= -1/2*I
```

Algebraic pole

```mathematica
In[6]:= Residue[x^3/(x^4 - 2), {x, 2^(1/4)}]
Out[6]= 1/4
```

Unknown numerator

```mathematica
In[7]:= Residue[f[z]/z^5, {z, 0}]
Out[7]= 1/24 Derivative[4][f][0]
```

```mathematica
In[8]:= Residue[Zeta[z]/(z - 1)^10, {z, 1}]
Out[8]= -1/362880 StieltjesGamma[9]
```

Branch point

```mathematica
In[9]:= Residue[1/Sqrt[z], {z, 0}]
Out[9]= Residue[1/Sqrt[z], {z, 0}]
```

### Applications (10)

A simple pole at the origin

```mathematica
In[10]:= Residue[1/z, {z, 0}]
Out[10]= 1
```

Simple pole at a complex location, via P(z0)/Q'(z0)

```mathematica
In[11]:= Residue[1/(z^2 + 1), {z, I}]
Out[11]= -1/2*I
```

An algebraic pole location, kept in one radical spelling

```mathematica
In[12]:= Residue[1/(z^4 + 1), {z, Exp[I Pi/4]}]
Out[12]= -1/4 (-1)^(1/4)
```

Fifth-order pole of a transcendental integrand

```mathematica
In[13]:= Residue[1/Sin[z]^5, {z, 0}]
Out[13]= 3/8
```

The summation kernel Pi Cot[Pi z] at a pole of the rational factor

```mathematica
In[14]:= Residue[Pi Cot[Pi z]/(z^2 + 1), {z, I}]
Out[14]= -1/2 Pi Coth[Pi]
```

Residue 1 at the only pole of Zeta

```mathematica
In[15]:= Residue[Zeta[s], {s, 1}]
Out[15]= 1
```

Symbolic double pole, self-validated across expansion orders

```mathematica
In[16]:= Residue[Exp[I z]/(z^2 + a^2)^2, {z, I a}]
Out[16]= ((-1/4*I) E^(-a))/a^3 + ((-1/4*I) E^(-a))/a^2
```

Fractional power times a double pole, via the dropped-pole guard

```mathematica
In[17]:= Residue[z^(1/3)/(1 + z^2)^2, {z, I}]
Out[17]= (-1/6*I) I^(1/3)
```

A branch point has no residue, so the call stays unevaluated

```mathematica
In[18]:= Residue[1/Sqrt[z], {z, 0}]
Out[18]= Residue[1/Sqrt[z], {z, 0}]
```

Limitation: Series has no expansion of Gamma at its poles, so this stays unevaluated; NResidue gives -0.166667

```mathematica
In[19]:= Residue[Gamma[z], {z, -3}]
Out[19]= Residue[Gamma[z], {z, -3}]
```

## Algorithm

residue.c -- Residue[expr, {z, z0}], the symbolic residue.

The residue of f at an isolated singularity z = z0 is the coefficient of (z - z0)^-1 in the Laurent expansion of f. We obtain it directly from the series engine: expand f to order (z - z0)^0 (which always spans the -1 term, however deep the pole), then read the coefficient at exponent -1 out of the resulting SeriesData[z, z0, {coefs}, nmin, nmax, den].

A residue is well defined only for an ordinary Laurent expansion (den == 1). A fractional-power (Puiseux) expansion, den > 1, signals a branch point, where the residue is undefined -- we leave the call unevaluated, matching Mathematica (e.g. Residue[1/Sqrt[z], {z, 0}]).

```text
Algebraic pole locations.  The series engine decides whether z0 is a pole by
```

evaluating the denominator there and testing it against zero; but for a pole whose location is a SUM of radicals (e.g. z0 = -2 + Sqrt[3], a root of 1 + 4 z + z^2), Denominator(z0) is an expression like 1 + 4 (-2 + Sqrt[3]) + (-2 + Sqrt[3])^2 that does not auto-simplify to 0, so

```text
the pole is missed and the residue wrongly comes out 0.  We defeat this by
```

expanding about z0 EXPLICITLY: substitute z -> z0 + w, then Expand the denominator of the result -- polynomial expansion collapses the radical arithmetic (Sqrt[3]^2 -> 3, ...) so the vanishing constant term becomes a

```text
literal 0 and the w-factor of the pole is exposed.  Reading the (z-z0)^-1
```

coefficient is then a plain Series-at-0 of the expanded form.

## Implementation notes

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

**Attributes:** `Protected`.

## References

**See also:** [NResidue](../../numerical-calculus/NResidue/), [Together](../../algebra/Together/), [Zeta](../../special-functions/Zeta/), [Series](../../power-series/Series/)

- Source: [`src/calculus/residue.c`](https://github.com/stblake/mathilda/blob/main/src/calculus/residue.c)
- Specification: [`docs/spec/builtins/calculus.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/calculus.md)
- Tests: [`tests/test_residue.c`](https://github.com/stblake/mathilda/blob/main/tests/test_residue.c)

## Notes & additional examples

### Notes

- `Residue[f, {z, z0}]` is the coefficient of `(z - z0)^-1` in the Laurent
  expansion of `f` at `z0`. An analytic point gives `0`. Poles of any order, at
  numeric or symbolic centres, are handled.
- For higher-order poles the coefficient is accepted only after two consecutive
  expansion orders agree. This closes the silent wrong answers at symbolic
  double and triple poles that existed before v0.268. Those wrong residues also
  fed `Integrate[..., Method -> "Residue"]`.
- A branch point (a fractional-power expansion) is left unevaluated, matching
  Mathematica. An essential singularity (`E^(1/z)` at `0`) and the poles of
  `Gamma` are also left unevaluated, because `Series` has no expansion there.
  Mathematica returns `1` and `(-1)^n/n!` for these. Use `NResidue` for a
  numerical value.
- Only the two-argument form is supported. Attributes: `Protected`.
- See also: `NResidue`, `Series`, `SeriesCoefficient`, `Integrate` (contour /
  residue-theorem definite integration).
