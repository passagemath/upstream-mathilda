# FresnelC

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`FresnelC[z]`**

gives the Fresnel integral C(z) = Integral\_0^z Cos\[Pi t^2/2\] dt.

**`FresnelC[+-Infinity] = +-1/2, FresnelC[+-I Infinity] = +-I/2.`**

<details>
<summary>Notes</summary>

An entire, odd function with no branch cuts. FresnelC\[0\] = 0, Real and complex inputs evaluate numerically at machine or arbitrary (MPFR) precision; D\[FresnelC\[z\], z\] = Cos\[Pi z^2/2\]. Listable.

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

Odd and entire, with limits 1/2 and I/2

```mathematica
In[5]:= {FresnelC[0], FresnelC[Infinity], FresnelC[I Infinity]}
Out[5]= {0, 1/2, 1/2*I}
```

Machine real: auxiliary f, g tables with an exactly reduced phase

```mathematica
In[6]:= FresnelC[1.8]
Out[6]= 0.333633
```

A huge argument is still accurate because the phase Pi x^2/2 is reduced exactly

```mathematica
In[7]:= FresnelC[123456.789]
Out[7]= 0.499999
```

Machine complex through the exactly rotated erf kernel

```mathematica
In[8]:= FresnelC[2.5 + I]
Out[8]= 116.648 - 105.229*I
```

The scaled axis series handles a huge argument just off the real axis

```mathematica
In[9]:= FresnelC[1.*^10 + 1.*^-10 I]
Out[9]= 0.5 + 3.67608e-10*I
```

Arbitrary precision through the guarded MPFR series

```mathematica
In[10]:= N[FresnelC[2], 50]
Out[10]= 0.488253406075340754500223503357261037688367154509216
```

Threads element-wise over a packed real list

```mathematica
In[11]:= FresnelC[{0.5, 1.0, 1.5}]
Out[11]= {0.492344, 0.779893, 0.445261}
```

The integrand, with the chain rule

```mathematica
In[12]:= D[FresnelC[x^2], x]
Out[12]= 2 x Cos[1/2 Pi x^4]
```

The odd-symmetry fold

```mathematica
In[13]:= FresnelC[-x]
Out[13]= -FresnelC[x]
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

**Algorithm.** `builtin_fresnelc` returns `C(z) = ∫₀ᶻ cos(π t²/2) dt` (Mathematica's π/2 convention). `C` is entire and odd. `FresnelS` shares every kernel: the pair `(C, S)` is computed together and each builtin keeps its own component.

A machine call first goes to `fresnel_machine`, which uses the `sf_machine_scalar_args` gate and the kernels in `sf_machine_error.c`.

**Real argument** (`sf_fresnel_c_r`):

- `|x| <= 1`: the Maclaurin series, with `x³` formed exactly.
- `1 < |x| < 8`: `C = 1/2 + f sin(πx²/2) - g cos(πx²/2)`, with the auxiliary functions `f` and `g` read from 0.5-wide polynomial tables (the Cephes `fresnl` structure).
- `|x| >= 8`: the nine-term asymptotic series for `f` and `g`.
- `|x| > 2^54`: exactly `1/2`.

The phase `πx²/2` is reduced exactly: `x²` is split as a Dekker double-double and reduced mod 4. A plain libm `sin(π x x/2)` would carry about `x²` ulps of phase error, so this keeps even `FresnelC[123456.789]` accurate to within a few ulp.

**Complex argument** (`sf_fresnel_c_c`):

- `|z| <= 1.3`: the Maclaurin series in `z⁴`.
- Next to an axis (`π|xy| < 1/2`), and also for huge `|z|` with `π|xy| < 8`: a scaled Taylor series about the real axis. The imaginary axis is mapped onto it with `C(iz) = iC(z)`. This keeps the small component relatively accurate and stops `FresnelC[1.*^10 + 1.*^-10 I]` from overflowing.
- Elsewhere: `C ± iS = ((1 ± i)/2) erf((√π/2)(1 ∓ i) z)` through the complex erf core. Every quantity that gets exponentiated (`X² = (π/4)(x+y)²`, `XY`, `Y² - X² = ∓πxy`) is built from the exact `x` and `y` rather than from a rounded rotated argument, which avoids about `π|z|²` ulps of error.

**Declined and MPFR calls** go to `fresnel_one_arg`:

- **Exact special values:** `0`, `±Infinity → ±1/2`, `±I Infinity → ±I/2`, and `ComplexInfinity` or `Indeterminate → Indeterminate`.
- **Odd fold:** `FresnelC[-x] = -FresnelC[x]`.
- **MPFR real:** the convergent series with about `(π/2)x²/ln2` guard bits, or the DLMF 7.12 asymptotic expansion at optimal truncation for large `|x|`. The asymptotic form is used on real inputs only, because the constant `1/2` changes across Stokes lines.
- **MPFR complex:** the paired `A = C + iS`, `B = C - iS` series in the `ncpx` toolkit, which is valid everywhere.

**Data structures.** `Expr`; `ncpx` (`numeric_complex.h`) for MPFR complex values. The machine side uses double-double helpers and the tables `FRF*_T` and `FRG*_T`. The ND kernel is `NDKU_FresnelC = { sf_fresnel_c_c, ndk_FresnelC_r }`, registered with `REG_U`, and covers real and complex buffers. `Compile[]` lowers `FresnelC` through it for real and complex scalars and for rank-1 arrays. Because every surface uses the same helper, scalar, packed, `NDArray` and `Compile[]` results match bit for bit.

**Complexity / limits.** The machine kernels are `O(1)`: real error is at most 2.5 ulp (median 0.3), and complex error is at most 3.5 ulp away from the zeros. The MPFR cost grows with `|z|²` and with the requested precision. Exact non-special arguments (`FresnelC[2]`) stay unevaluated. `D[FresnelC[x], x] = Cos[π x²/2]`. Attributes: `Listable`, `NumericFunction`, `Protected`.

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## References

**See also:** [FresnelS](../../special-functions/FresnelS/), [NDArray](../../linear-algebra/NDArray/), [Piecewise](../../control-flow/Piecewise/)

- DLMF §7.2(iii) — the Fresnel integral C(z) = Int_0^z cos(π t^2/2) dt.
- DLMF §7.3.22, §7.5.3–4 — Fresnel integrals through erf, and through the auxiliary functions f, g.
- DLMF §7.12 — the asymptotic expansion of the Fresnel integrals.
- T. J. Dekker, *A floating-point technique for extending the available precision*, Numer. Math. 18 (1971) 224–242.
- Source: [`src/special_functions/fresnel.c`](https://github.com/stblake/mathilda/blob/main/src/special_functions/fresnel.c)
- Specification: [`docs/spec/builtins/special-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/special-functions.md)
- Tests: [`tests/test_compile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile.c)
- Tests: [`tests/test_fresnelc.c`](https://github.com/stblake/mathilda/blob/main/tests/test_fresnelc.c)
- Tests: [`tests/test_numeric_stress.c`](https://github.com/stblake/mathilda/blob/main/tests/test_numeric_stress.c)
- Tests: [`tests/test_sf_machine_error.c`](https://github.com/stblake/mathilda/blob/main/tests/test_sf_machine_error.c)

## Notes & additional examples

### Notes

- `FresnelC[z] = ∫₀ᶻ Cos[Pi t²/2] dt` (Mathematica's π/2 convention), entire and odd; `FresnelC[±Infinity] = ±1/2`, `FresnelC[±I Infinity] = ±I/2`.
- Machine real and complex arguments use the `sf_machine_error.c` kernels — Maclaurin series, auxiliary-function tables with exact phase reduction, the axis Taylor series, and the complex erf core fed exactly rotated arguments — shared bit for bit with packed lists, `NDArray` and `Compile[]`.
- Arbitrary precision uses the convergent series with `(π/2)|z|²/ln2` guard bits, or (real arguments only) the DLMF 7.12 asymptotic expansion.
- Attributes: `Listable`, `NumericFunction`, `Protected`. See also `FresnelS`, `Erf`, `Erfi`.
