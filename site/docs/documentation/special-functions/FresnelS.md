# FresnelS

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`FresnelS[z]`**

gives the Fresnel integral S(z) = Integral\_0^z Sin\[Pi t^2/2\] dt.

**`FresnelS[+-Infinity] = +-1/2, FresnelS[+-I Infinity] = -+I/2.`**

<details>
<summary>Notes</summary>

An entire, odd function with no branch cuts. FresnelS\[0\] = 0, Real and complex inputs evaluate numerically at machine or arbitrary (MPFR) precision; D\[FresnelS\[z\], z\] = Sin\[Pi z^2/2\]. Listable.

</details>

## Examples (13)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= FresnelC[1.8]
Out[1]= 0.333633

In[2]:= FresnelS[1.8]
Out[2]= 0.450939

In[3]:= FresnelC[2.5 + I]
Out[3]= 116.648 - 105.229*I

In[4]:= FresnelS[-I Infinity]
Out[4]= 1/2*I
```

### Applications (9)

Note the sign flip on the imaginary directions

```mathematica
In[5]:= {FresnelS[0], FresnelS[Infinity], FresnelS[-I Infinity]}
Out[5]= {0, 1/2, 1/2*I}
```

Machine real on the shared Fresnel kernel

```mathematica
In[6]:= FresnelS[1.8]
Out[6]= 0.450939
```

The exactly reduced phase keeps a huge argument accurate

```mathematica
In[7]:= FresnelS[123456.789]
Out[7]= 0.500002
```

Machine complex

```mathematica
In[8]:= FresnelS[2.5 + I]
Out[8]= 105.729 + 116.148*I
```

The tiny imaginary component is accurate relative to itself

```mathematica
In[9]:= FresnelS[1.*^10 + 1.*^-10 I]
Out[9]= 0.5 - 3.25467e-30*I
```

Arbitrary precision

```mathematica
In[10]:= N[FresnelS[2], 30]
Out[10]= 0.3434156783636982421953008159582
```

The packed list and the scalar calls agree

```mathematica
In[11]:= FresnelS[{0.5, 1., 1.5}] == Map[FresnelS, {0.5, 1., 1.5}]
Out[11]= True
```

The Maclaurin series begins with Pi x^3/6

```mathematica
In[12]:= Series[FresnelS[x], {x, 0, 7}]
Out[12]= 1/6 Pi x^3 - 1/336 Pi^3 x^7 + O[x]^8
```

The integrand

```mathematica
In[13]:= D[FresnelS[x], x]
Out[13]= Sin[1/2 Pi x^2]
```

## Algorithm

Mathilda -- the Fresnel integrals (Pi/2-normalized).

```text
  FresnelC[z] = Int_0^z Cos[Pi t^2 / 2] dt
  FresnelS[z] = Int_0^z Sin[Pi t^2 / 2] dt
```

Both are entire and odd, with no branch cuts. FresnelC and FresnelS share one numeric kernel: the pair (C, S) is computed together and each builtin returns its component. Evaluation is layered so each kind of argument takes the cheapest route:

```text
  exact special values     ->  0, +-1/2, +-I/2, Indeterminate
  machine / arbitrary real ->  MPFR series (small |x|) or asymptotic (large)
  complex (any precision)  ->  the paired A/B ncpx series with guard bits
  everything else          ->  stays symbolic (return NULL)
```

Convergent Maclaurin series (valid for all z):

```text
  C(z) = Sum_{m>=0} (-1)^m (Pi/2)^(2m)   z^(4m+1) / ((2m)!   (4m+1)),
  S(z) = Sum_{m>=0} (-1)^m (Pi/2)^(2m+1) z^(4m+3) / ((2m+1)! (4m+3)).
```

Equivalently A(z) = C(z) + i S(z) = Sum_{k>=0} (i Pi/2)^k z^(2k+1)/(k!(2k+1)) and B(z) = C(z) - i S(z) is the same with i -> -i; then C = (A+B)/2 and S = (A-B)/(2i). The real path sums C and S as two real series directly; the complex path sums A and B together in one ncpx loop.

The partial sums reach magnitude ~e^((Pi/2)|z|^2) before the O(1)-sized answer emerges, so the MPFR paths add ~(Pi/2)|z|^2/ln2 guard bits to absorb that cancellation exactly. For large real |x| the convergent series is infeasible; there we use the asymptotic expansion (DLMF 7.12)

```text
  C(x) = 1/2 + f(x) sin(Pi x^2/2) - g(x) cos(Pi x^2/2),
  S(x) = 1/2 - f(x) cos(Pi x^2/2) - g(x) sin(Pi x^2/2),
    f(x) ~ (1/(Pi x))   Sum_j (-1)^j (4j-1)!! / (Pi x^2)^(2j),
    g(x) ~ (1/(Pi^2 x^3)) Sum_j (-1)^j (4j+1)!! / (Pi x^2)^(2j),
```

summed to the smallest term (optimal truncation). The asymptotic constant 1/2 is the value only in a sector around the real axis (Stokes phenomenon), so it is used for real inputs only; complex inputs always use the convergent series (correct everywhere, merely costlier for large |z|).

Attributes: Listable, NumericFunction, Protected.

## Implementation notes

**Algorithm.** `builtin_fresnels` returns `S(z) = ∫₀ᶻ sin(π t²/2) dt`, using the π/2 convention. `S` is entire and odd, and it shares every kernel with `FresnelC`: the pair `(C, S)` is computed together and this builtin keeps the `S` component.

A machine call first goes through `fresnel_machine`, which uses the `sf_machine_scalar_args` gate.

**Real argument** (`sf_fresnel_s_r`):

- `|x| <= 1`: `S = x³ Σ s_k x^{4k}`, with `x³` formed to within a quarter ulp.
- `1 < |x| < 8`: `S = 1/2 - f cos(πx²/2) - g sin(πx²/2)`, where the auxiliary functions `f` and `g` come from 0.5-wide tables.
- `|x| >= 8`: the asymptotic series for `f` and `g`.
- Throughout, the phase `πx²/2` is reduced exactly: `x²` is held as a Dekker double-double and reduced mod 4.

**Complex argument** (`sf_fresnel_s_c`):

- `|z| <= 1.3`: the Maclaurin series.
- Next to an axis, and for huge `|z|` with `π|xy| < 8`: the scaled axis Taylor series. The imaginary axis is handled through `S(iz) = -iS(z)`.
- Elsewhere: `S = [(1+i)E₁ - (1-i)E₂]/(4i)`, where `E₁` and `E₂` are erf values of the exactly rotated arguments.

**Declined and MPFR calls** go to `fresnel_one_arg`:

- **Exact values:** `0`; `±Infinity → ±1/2`; `±I Infinity → ∓I/2` (note the sign flip); `ComplexInfinity` and `Indeterminate → Indeterminate`.
- **Odd fold** for symbolic negative arguments.
- **MPFR real:** the guarded series, or the DLMF 7.12 asymptotic form at large `|x|`.
- **MPFR complex:** the paired `A`/`B` series in `ncpx`, with `S = (A - B)/(2i)`.

**Data structures.** `Expr`; `ncpx` for MPFR complex values; double-double helpers plus generated tables on the machine side. The ND kernel is `NDKU_FresnelS = { sf_fresnel_s_c, ndk_FresnelS_r }`, registered with `REG_U`, and covers real and complex buffers. `Compile[]` lowers `FresnelS` through it for real and complex scalars and for rank-1 arrays. All surfaces call the same helper, so their results match bit for bit.

**Complexity / limits.** The machine kernels run in `O(1)`. The real kernel stays within 2.5 ulp. Away from the zeros the complex kernel stays within 3.5 ulp, and near an axis each component is accurate relative to its own size (`FresnelS[1.*^10 + 1.*^-10 I]` has an imaginary part of `-3.25*10^-30`). The MPFR cost grows with `|z|²` and with the requested precision. Exact non-special arguments stay unevaluated. `D[FresnelS[x], x] = Sin[π x²/2]`. Attributes: `Listable`, `NumericFunction`, `Protected`.

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## References

**See also:** [FresnelC](../../special-functions/FresnelC/), [NDArray](../../linear-algebra/NDArray/), [Piecewise](../../control-flow/Piecewise/)

- DLMF §7.2(iii) — the Fresnel integral S(z) = Int_0^z sin(π t^2/2) dt.
- DLMF §7.3.22, §7.5.3–4 — Fresnel integrals through erf, and through the auxiliary functions f, g.
- DLMF §7.12 — the asymptotic expansion of the Fresnel integrals.
- T. J. Dekker, *A floating-point technique for extending the available precision*, Numer. Math. 18 (1971) 224–242.
- Source: [`src/special_functions/fresnel.c`](https://github.com/stblake/mathilda/blob/main/src/special_functions/fresnel.c)
- Specification: [`docs/spec/builtins/special-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/special-functions.md)
- Tests: [`tests/test_compile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile.c)
- Tests: [`tests/test_fresnels.c`](https://github.com/stblake/mathilda/blob/main/tests/test_fresnels.c)
- Tests: [`tests/test_numeric_stress.c`](https://github.com/stblake/mathilda/blob/main/tests/test_numeric_stress.c)
- Tests: [`tests/test_series_assumptions.c`](https://github.com/stblake/mathilda/blob/main/tests/test_series_assumptions.c)

## Notes & additional examples

### Notes

- `FresnelS[z] = ∫₀ᶻ Sin[Pi t²/2] dt`, entire and odd; `FresnelS[±Infinity] = ±1/2` but `FresnelS[±I Infinity] = ∓I/2`.
- `FresnelS` and `FresnelC` share every kernel: the pair (C, S) is computed together. Machine arguments run on the `sf_machine_error.c` kernels, shared bit for bit with packed lists, `NDArray` and `Compile[]`; arbitrary precision uses the guarded MPFR series or the real-axis asymptotic expansion.
- Attributes: `Listable`, `NumericFunction`, `Protected`. See also `FresnelC`, `Erf`, `Erfi`.
