# Erfi

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Erfi[z]`**

gives the imaginary error function erfi(z) = -I Erf\[I z\] = (2/Sqrt\[Pi\]) Integral\_0^z e^(t^2) dt.

**`Erfi[0] = 0, Erfi[Infinity] = Infinity, Erfi[I Infinity] = I. An entire`**

<details>
<summary>Notes</summary>

function, odd in z. Real and complex inputs evaluate numerically at machine or arbitrary (MPFR) precision; D\[Erfi\[z\], z\] = (2/Sqrt\[Pi\]) E^(z^2). Listable.

</details>

## Examples (12)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= Erfi[I Infinity]
Out[1]= I

In[2]:= Erfi[2.5]
Out[2]= 130.396

In[3]:= Erfi[1.5 - I]
Out[3]= -0.70136 - 1.84683*I
```

### Applications (9)

The imaginary-axis limit is finite, and the function is odd

```mathematica
In[4]:= {Erfi[0], Erfi[I Infinity], Erfi[-x]}
Out[4]= {0, I, -Erfi[x]}
```

Machine real: e^(x^2) times a tabulated Dawson factor

```mathematica
In[5]:= Erfi[2.5]
Out[5]= 130.396
```

Still a machine double near the top of the range

```mathematica
In[6]:= Erfi[26.]
Out[6]= 8.31464e+291
```

Past the double range the MPFR series returns an extended-exponent number

```mathematica
In[7]:= Erfi[30.]
Out[7]= 1.379047205635398e+389
```

Machine complex through the complex erf kernel

```mathematica
In[8]:= Erfi[1.5 - I]
Out[8]= -0.70136 - 1.84683*I
```

The all-positive Maclaurin series at arbitrary precision

```mathematica
In[9]:= N[Erfi[1/2], 50]
Out[9]= 0.614952094696510980839681185623641393051345617895404
```

Another arbitrary-precision value

```mathematica
In[10]:= N[Erfi[1], 30]
Out[10]= 1.650425758797542876025337729561
```

Like Erf's series without the alternating signs

```mathematica
In[11]:= Series[Erfi[x], {x, 0, 7}]
Out[11]= (2 x)/Sqrt[Pi] + (2 x^3)/(3 Sqrt[Pi]) + x^5/(5 Sqrt[Pi]) + x^7/(21 Sqrt[Pi]) + O[x]^8
```

The derivative has the positive exponent

```mathematica
In[12]:= D[Erfi[z], z]
Out[12]= (2 E^z^2)/Sqrt[Pi]
```

## Algorithm

Mathilda -- the imaginary error function.

```text
  Erfi[z]   imaginary error function   erfi(z) = erf(i z)/i
                                               = (2/sqrt(pi)) Int_0^z e^t^2 dt
```

erfi is an entire function (no branch cuts) and odd in z. Evaluation is layered so each kind of argument takes the cheapest route:

```text
  exact special values   ->  0, +-Infinity, +-I (the imaginary-axis limits
                             are FINITE: erfi(i y) -> i as y -> +Infinity)
  symbolic odd argument   ->  Erfi[-x] = -Erfi[x]
  machine / arbitrary real -> the all-positive Maclaurin series
                             erfi(x) = (2/sqrt(pi)) Sum x^(2n+1)/(n!(2n+1)),
                             evaluated in MPFR. For real x every term shares
                             x's sign, so the partial sums climb monotonically
                             to the result -- no cancellation, no e^x^2
                             prefactor, pure-real arithmetic.
  complex (any precision)  -> erfi(z) = -i erf(i z), reusing the
                             cancellation-aware erf series (DLMF 7.6.2) in
                             MPFR with |z|^2/ln2 guard bits, so even
                             machine-precision complex results carry full
                             accuracy. Double-precision series is the
                             USE_MPFR=0 fallback.
  everything else         ->  stays symbolic (return NULL)
```

There is no libm erfi and no mpfr_erfi, so both numeric kernels are hand-rolled here.

Attributes: Listable, NumericFunction, Protected.

## Implementation notes

**Algorithm.** `builtin_erfi` handles `Erfi[z] = -i erf(iz)`, which is entire and odd. Neither libm nor MPFR provides `erfi`, so Mathilda implements every numeric path itself. A machine call first passes through the `sf_machine_scalar_args` gate and the shared kernels in `sf_machine_error.c`:

- **Real** (`sf_erfi_r`):
  - For `|x| < 1/2`, the all-positive Maclaurin polynomial.
  - For `1/2 <= |x| < 10`, `erfi(x) = e^{x²}·g(x)`. Here `g = e^{-x²} erfi(x) = (2/√π)·Dawson(x)` comes from 0.5-wide polynomial tables (`DAW1_T`, `DAW2_T`), and `e^{x²}` comes from the exact Dekker square.
  - For `|x| >= 10`, a 14-term asymptotic series for `g`.
  - When `e^{x²}` alone would overflow but `erfi` still fits (`x² ∈ (709, 713.6)`), the kernel scales by `e^{-64}` first.

  The error is ≤ 2.6 ulp; the old summed series reached 30 ulp at `|x| ≈ 5`. The kernel declines `|x| >= 26.72`.
- **Complex** (`sf_erfi_c`): `-i erf(iz)`, routed through the complex `erf` kernel (Maclaurin, axis Taylor series, Faddeeva table, Salzer's series).

Anything declined, or any MPFR call, goes to `erfi_one_arg`:

- **Exact special values:** `0`, `±Infinity → ±Infinity`, and the finite imaginary-axis limits `±I Infinity → ±I`.
- **Odd fold:** `Erfi[-x] = -Erfi[x]`.
- **Real at any precision:** the all-positive Maclaurin series `(2/√π) Σ x^{2n+1}/(n!(2n+1))` in MPFR. Nothing cancels, so a flat 64-bit guard is enough.
- **Complex at any precision:** `-i erf(iz)` through the guarded erf series in `ecx`.

Results pass through `numeric_promote_result_if_degenerate`, so `Erfi[30.]` comes back as an extended-exponent number, not `Infinity`.

**Data structures.** `Expr`; `ecx` (a pair of `mpfr_t`) and plain MPFR scalars. The ND kernel is `NDKU_Erfi = { sf_erfi_c, ndk_Erfi_r }`, registered with `REG_U`, covering real and complex packed arrays and `NDArray`s. `Compile[]` lowers `Erfi` through the kernel for real and complex scalars and for rank-1 arrays. Every surface shares one helper, so the answers agree bit for bit.

**Complexity / limits.** The machine kernels are `O(1)`. The MPFR real series needs `O(x²)` terms; the complex series' term count and guard bits grow with `|z|²` and with the precision. Exact non-special arguments (`Erfi[2]`) stay unevaluated. Attributes: `Listable`, `NumericFunction`, `Protected`.

- Exact special values: `Erfi[0] = 0`, `Erfi[Infinity] = Infinity`,
  `Erfi[-Infinity] = -Infinity`. The imaginary-axis limits are **finite** (unlike
  `Erf`): `Erfi[I Infinity] = I`, `Erfi[-I Infinity] = -I`, since
  erfi(i y) = -i·erf(-y) → i as y → ∞. `ComplexInfinity` and `Indeterminate`
  pass through.
- Odd symmetry for symbolic arguments: `Erfi[-x] = -Erfi[x]`,
  `Erfi[-2 x] = -Erfi[2 x]`.
- Numeric evaluation (there is no libm/MPFR `erfi`, so the kernels are
  hand-rolled):
  - **Machine precision** → the shared kernels of `sf_machine_error.c`. Real:
    |x| < 1/2 the positive Maclaurin series; otherwise
    `erfi(x) = e^(x²)·g(x)` with g = e^(−x²) erfi(x) = (2/√π)·Dawson(x) from
    0.5-wide tables (x < 10) or its 14-term asymptotic series, e^(x²) from the
    exact square (scaled by e^(−64) where e^(x²) alone would overflow,
    x² ∈ (709, 713.6)); ≤ 2.6 ulp (the old summed series reached 30 ulps at
    |x| ≈ 5, and cost ~80 ns). Complex: `−i erf(iz)` through the complex erf
    kernel. `|x| > 26.72` declines to MPFR.
  - Real, arbitrary precision → the all-positive Maclaurin series
    erfi(x) = (2/√π) Σ x^(2n+1)/(n!(2n+1)). For real x every term shares x's
    sign, so the partial sums climb monotonically to the result — no
    cancellation, no e^(x²) prefactor — evaluated in MPFR with output precision
    tracking the input, e.g. `Erfi[2.5] = 130.396`, `Erfi[0.5] = 0.614952`,
    `N[Erfi[1/2], 50] = 0.61495209469651098083968118562364139305134561789540`.
  - **Complex, arbitrary precision** → erfi(z) = -i erf(i z),
    reusing the cancellation-aware erf series (DLMF 7.6.2) in MPFR with
    `|z|²/ln2` guard bits, e.g. `Erfi[1.5 - I] = -0.70136 - 1.84683 I`,
    `N[Erfi[1/2 + I], 30] = 0.187973467223383313628263810077 + 0.950709728318957173804611826379 I`.
    A `double complex` series is the fallback for `USE_MPFR=0` builds.
- Derivative: `D[Erfi[z], z] = (2/Sqrt[Pi]) E^(z^2)` (positive exponent, vs
  `Erf`'s E^(−z²); chain rule applies), so the origin Taylor series follows from
  the generic `D`-based fallback, e.g. `Series[Erfi[x], {x, 0, 7}]` begins
  `2/Sqrt[Pi] x + 2/(3 Sqrt[Pi]) x^3 + …`.
- All other arguments (symbolic `Erfi[x]`, exact `Erfi[2]`) stay unevaluated.

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## References

**See also:** [Erf](../../special-functions/Erf/), [D](../../calculus/D/)

- DLMF §7.2, §7.6.2 — the imaginary error function erfi(z) = erf(iz)/i = (2/sqrt(π)) Int_0^z e^{t^2} dt.
- DLMF §7.2(ii), §7.12 — Dawson's integral and the asymptotic expansion behind the large-x tail.
- H. E. Salzer, *Formulas for calculating the error function of a complex variable*, Math. Tables Aids Comput. 5 (1951) 67–70.
- Source: [`src/special_functions/erfi.c`](https://github.com/stblake/mathilda/blob/main/src/special_functions/erfi.c)
- Specification: [`docs/spec/builtins/special-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/special-functions.md)
- Tests: [`tests/test_cherry_ei.c`](https://github.com/stblake/mathilda/blob/main/tests/test_cherry_ei.c)
- Tests: [`tests/test_compile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile.c)
- Tests: [`tests/test_erfi.c`](https://github.com/stblake/mathilda/blob/main/tests/test_erfi.c)
- Tests: [`tests/test_interval.c`](https://github.com/stblake/mathilda/blob/main/tests/test_interval.c)

## Notes & additional examples

### Notes

- `Erfi[z] = -I Erf[I z] = (2/Sqrt[Pi]) ∫₀ᶻ e^(t²) dt` is entire and odd; `Erfi[±I Infinity] = ±I`.
- There is no libm or MPFR `erfi`: machine reals use a Dawson-factor table times `e^(x²)` (exact square), machine complex the complex erf kernel, both shared bit for bit with packed lists, `NDArray` and `Compile[]`; arbitrary precision sums the all-positive Maclaurin series (real) or the guarded erf series (complex).
- Attributes: `Listable`, `NumericFunction`, `Protected`. See also `Erf`, `Erfc`, `FresnelC`.
