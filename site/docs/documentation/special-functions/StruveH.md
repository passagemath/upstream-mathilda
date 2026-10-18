# StruveH

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`StruveH[nu, z]`**

gives the Struve function H\_nu(z), the particular solution of the inhomogeneous Bessel equation with forcing (z/2)^(nu+1).

**`StruveH[nu, 0] = 0 for Re nu > -1. Numericizes to arbitrary precision via`**

<details>
<summary>Notes</summary>

H\_nu(z) = (z/2)^(nu+1) (2/(Sqrt\[Pi\] Gamma\[nu+3/2\])) 1F2(1; 3/2, nu+3/2; -z^2/4). FLINT/arb (the regularized 1F2, rigorous, all nu and complex z), and at machine precision via a fast real double kernel (convergent series / Y\_nu+K\_nu asymptotic) that also serves the packed/NDArray and Compile\[\] paths. Listable.

</details>

## Examples (13)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= StruveH[0, 1.5]
Out[1]= 0.736723

In[2]:= N[StruveH[1, 13/10], 30]
Out[2]= 0.3201231459764202395239403961771

In[3]:= StruveH[1, {0.5, 1., 2.}]
Out[3]= {0.0521737, 0.198457, 0.646764}
```

### Applications (10)

No closed form is substituted for exact symbolic input

```mathematica
In[4]:= StruveH[1/2, x]
Out[4]= StruveH[1/2, x]
```

Machine value from the double series kernel

```mathematica
In[5]:= StruveH[0, 1.5]
Out[5]= 0.736723
```

The Poisson-type integral (2/π)∫ sin(z cos θ) dθ gives the same number

```mathematica
In[6]:= NIntegrate[2/Pi Sin[1.5 Cos[t]], {t, 0, Pi/2}]
Out[6]= 0.736723
```

Agreement with the elementary half-integer form, to rounding

```mathematica
In[7]:= StruveH[1/2, 2.] - Sqrt[2/(Pi 2.)] (1 - Cos[2.])
Out[7]= 1.11022e-16
```

Arbitrary precision through the regularized 1F2 in arb

```mathematica
In[8]:= N[StruveH[1, 13/10], 40]
Out[8]= 0.32012314597642023952394039617728160520309
```

The arb path at ν = -3/2, checked against H_{-3/2} = -J_{3/2}

```mathematica
In[9]:= N[StruveH[-3/2, 2], 20] - N[-BesselJ[3/2, 2], 20]
Out[9]= 0.0
```

Large argument: the Y_ν + K_ν asymptotic regime

```mathematica
In[10]:= StruveH[0, 30.]
Out[10]= -0.0960984
```

Independent check of the previous line

```mathematica
In[11]:= NIntegrate[2/Pi Sin[30 Cos[t]], {t, 0, Pi/2}]
Out[11]= -0.0960984
```

Listable; a packed list runs the same kernel

```mathematica
In[12]:= StruveH[1, {0.5, 1., 2.}]
Out[12]= {0.0521737, 0.198457, 0.646764}
```

Complex z needs N: the machine kernel is real-only

```mathematica
In[13]:= N[StruveH[0, 1 + I], 20]
Out[13]= 0.766322794092971087217 + 0.484304188659985084847*I
```

## Algorithm

Mathilda -- the Struve function H_nu(z).

```text
  StruveH[nu, z]   H_nu(z) = (z/2)^(nu+1) (2/(Sqrt[Pi] Gamma[nu+3/2]))
                             1F2(1; 3/2, nu+3/2; -z^2/4)
```

Evaluation:

```text
  StruveH[nu, 0] = 0            exact, for a real nu > -1
  numeric nu and inexact z      via the 1F2 representation (reusing the tested
                                HypergeometricPFQ numeric evaluator)
  everything else               stays symbolic (return NULL)
```

The 1F2 identity is also the reduction used by the half-line Mellin transform (reduce_to_hypergeometric in integrate_ramanujan.c), so Integrate[x^(s-1) StruveH[0, a x], {x,0,Infinity}] closes through rec_pfq.

Attributes: Listable, NumericFunction, Protected.

## Implementation notes

**Algorithm.** `builtin_struveh` first returns the exact `0` when z is the integer 0 and ν is real with ν > −1. It then numericizes only when some argument is inexact and both are numeric quantities. The output precision is the smallest inexact bit count among the arguments.

- **Above 53 bits**, `flint_num_struveh` evaluates `(z/2)^(ν+1) · acb_hypgeom_pfq(1; 3/2, ν+3/2; −z²/4, regularized)` in arb. The regularized 1F2 absorbs the `1/Γ(ν+3/2)` pole, so every complex ν works. The starting working precision includes a cancellation allowance of about 1.45·|z| + 8 bits. The `nb_should_retry` / `nb_emit_bits` ladder doubles the precision until the ball radius certifies the requested digits.
- **At machine precision** with real ν and z, `sf_machine_struve_h` runs two regimes. Each returns a value together with a predicted relative error:
  - `struve_series`, the convergent 1F2 series. Its error estimate is the measured cancellation `max|term|/|sum| · ε`.
  - `struve_asymp`, which computes `Y_ν + K_ν` from two optimally truncated divergent series. Its error estimate is the size of the smallest omitted term.

  The kernel keeps the regime with the smaller predicted error and declines when even that one exceeds `STRUVE_ACCEPT = 1e-7`. A declined real z is retried through an expression template built on `HypergeometricPFQ`.

**Data structures.** The machine path works only on scalar doubles and keeps nothing between calls. The arb path uses `acb_t` balls for ν, z and the result, plus `acb_struct` parameter arrays for the 1F2. In `src/ndkernels.c`, `NDK_BIN2(StruveH, sf_machine_struve_h)` and `REG_B` register the same kernel for packed lists and `NDArray`, and `Compile[]` lowers it at scalar and rank-1 shapes. Every surface calls the same function, so the results are bit-identical.

**Complexity / limits.** The series costs O(|z|) terms, with a cap of 8000. The asymptotic regime stops at its least term, at most 400 terms. The two error curves cross near z ≈ 12–24 for small ν, and there the double result is only good to about 1e-9 absolute. There is no complex machine kernel, so a machine complex z stays unevaluated; an exact complex z under `N` works through arb. There is no `D` rule. One defect: at a non-positive half-integer ν (where Γ(ν+3/2) has a pole), the machine fallback template returns an unevaluated `0.0 HypergeometricPFQ[...]` instead of a number. The arb path gets these values right.

- Exact value: `StruveH[nu, 0] = 0` for a real ν > −1.
- **Arbitrary precision** (`N[StruveH[ν, z], p]` with an exact or MPFR argument):
  the rigorous FLINT/arb kernel `H_ν(z) = (z/2)^(ν+1) · ₁F̃₂(1; 3/2, ν+3/2; −z²/4)`
  (the *regularized* ₁F₂, so it is entire in ν — no pole to dodge — and valid for
  all ν over the whole complex plane). The working precision carries a
  cancellation allowance (~1.44·|z| bits) so a large |z| is resolved without
  blind doubling; the ball ladder renders only correct digits. Example:
  `N[StruveH[1, 13/10], 50] = 0.32012314597642023952394039617728160520311127053816`
  (matches mpmath).
- **Machine precision**: a fast real `double` kernel — the convergent ₁F₂-form
  series for small/moderate |z|, the `Y_ν + K_ν` asymptotic (DLMF 11.6.1 +
  10.17.4) for large |z| — shared bit-for-bit with the packed/NDArray kernel and
  `Compile[]` (`StruveH[0, 1.5] = 0.736723`, matching the integral representation
  (2/π)∫₀^(π/2) sin(z cos θ) dθ). Accurate to the double floor everywhere; in the
  narrow transition band z ∈ [~17, ~24] (small ν) both elementary series are
  intrinsically ill-conditioned and the result degrades to ~1e-9 absolute (as it
  does for every double implementation, scipy included) — the arb path above
  gives full precision there. A real argument the kernel declines falls back to
  the ₁F₂ template; a complex machine argument stays unevaluated (use `N` on an
  exact argument: `N[StruveH[0, 1 + I], 20]` goes through arb).
- Mellin transform: `∫₀^∞ x^(s-1) StruveH[0, a x] dx =
  π 2^(s-1) a^(-s)/(Cos(πs/2) Γ(1−s/2)²)`, closed by `Method -> "Mellin"` through
  the ₁F₂ reduction.

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## References

**See also:** [N](../../arithmetic/N/)

- DLMF §11.2.1 — the power series H_ν(z) = (z/2)^(ν+1) Σ (−1)^k (z/2)^(2k)/(Γ(k+3/2) Γ(k+ν+3/2)).
- DLMF §11.6.1 — the large-argument asymptotic expansion of K_ν = H_ν − Y_ν; DLMF §10.17.4 — the Hankel asymptotic expansion of Y_ν.
- Source: [`src/special_functions/struveh.c`](https://github.com/stblake/mathilda/blob/main/src/special_functions/struveh.c)
- Specification: [`docs/spec/builtins/special-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/special-functions.md)
- Tests: [`tests/test_integrate_ramanujan.c`](https://github.com/stblake/mathilda/blob/main/tests/test_integrate_ramanujan.c)
- Tests: [`tests/test_struveh.c`](https://github.com/stblake/mathilda/blob/main/tests/test_struveh.c)

## Notes & additional examples

### Notes

- `Listable`, `NumericFunction`, `Protected`. Exact input is never rewritten: `StruveH[1/2, x]` stays as it is.
- The only exact value is `StruveH[nu, 0] = 0` for a real ν > −1.
- Machine-precision real input runs `sf_machine_struve_h`: the convergent 1F2 series, or for large z the `Y_ν + K_ν` asymptotic, whichever predicts the smaller error. The same kernel backs packed lists, `NDArray` and `Compile[]` (scalar and rank-1), so all of them agree bit for bit. Between about z = 17 and z = 24 (small ν) the double answer is good to about 1e-9 absolute; `N[..., p]` gives full precision there.
- With precision above machine (`N[StruveH[ν, z], p]`, or an MPFR argument) the value comes from FLINT/arb as `(z/2)^(ν+1) · 1F2~(1; 3/2, ν+3/2; −z²/4)`, using the regularized 1F2, so it holds for any complex ν and z.
- A machine complex argument (`StruveH[0, 1. + I]`) stays unevaluated. Use `N` on exact input to get a complex value.
- No derivative rule exists yet: `D[StruveH[n, x], x]` returns an inert `Derivative[0, 1][StruveH][n, x]`, so `Series` at a regular point yields only `Derivative` coefficients.
- See also: `BesselJ`, `BesselY`, `HypergeometricPFQ`.
