# Erf

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Erf[z]`**

gives the error function erf(z) = (2/Sqrt\[Pi\]) Integral\_0^z e^(-t^2) dt.

**`Erf[z0, z1]`**

gives the generalized error function erf(z1) - erf(z0).

**`Erf[0] = 0, Erf[Infinity] = 1, Erf[-Infinity] = -1. An entire function,`**

<details>
<summary>Notes</summary>

odd in z. Real and complex inputs evaluate numerically at machine or arbitrary (MPFR) precision; D\[Erf\[z\], z\] = (2/Sqrt\[Pi\]) E^(-z^2). Listable.

</details>

## Examples (13)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= Erf[0.95]
Out[1]= 0.820891

In[2]:= Erf[1.*^-300 + 2. I]
Out[2]= 6.16074e-299 + 18.5648*I

In[3]:= Erf[1.5, 2]
Out[3]= 0.0292171
```

### Applications (10)

Exact special values and the odd-symmetry fold

```mathematica
In[4]:= {Erf[0], Erf[Infinity], Erf[-z]}
Out[4]= {0, 1, -Erf[z]}
```

Machine real on the shared erfcx-table kernel

```mathematica
In[5]:= Erf[0.95]
Out[5]= 0.820891
```

Next to the imaginary axis the tiny real part is still accurate to every digit

```mathematica
In[6]:= Erf[1.*^-300 + 2. I]
Out[6]= 6.16074e-299 + 18.5648*I
```

Arbitrary-precision real through mpfr_erf

```mathematica
In[7]:= N[Erf[1], 40]
Out[7]= 0.84270079294971486934122063508260925929605
```

Arbitrary-precision complex through the guarded Maclaurin series

```mathematica
In[8]:= N[Erf[1 + I], 20]
Out[8]= 1.31615128169794764488 + 0.190453469237834686284*I
```

The generalized form erf(z1) - erf(z0)

```mathematica
In[9]:= Erf[1.5, 2]
Out[9]= 0.0292171
```

The packed list runs the same kernel as the scalar call

```mathematica
In[10]:= Erf[{0.5, 1.5, 2.5}] == Table[Erf[x], {x, {0.5, 1.5, 2.5}}]
Out[10]= True
```

The odd Maclaurin series

```mathematica
In[11]:= Series[Erf[x], {x, 0, 7}]
Out[11]= (2 x)/Sqrt[Pi] - (2 x^3)/(3 Sqrt[Pi]) + x^5/(5 Sqrt[Pi]) - x^7/(21 Sqrt[Pi]) + O[x]^8
```

The chain rule on the Gaussian derivative

```mathematica
In[12]:= D[Erf[x^2], x]
Out[12]= (4 x E^(-x^4))/Sqrt[Pi]
```

An exact non-special argument stays unevaluated

```mathematica
In[13]:= Erf[2]
Out[13]= Erf[2]
```

## Algorithm

Mathilda -- the error function.

```text
  Erf[z]       error function       erf(z)  = (2/sqrt(pi)) Int_0^z e^-t^2 dt
  Erf[z0, z1]  generalized error    erf(z1) - erf(z0)
```

erf is an entire function (no branch cuts) and odd in z. Evaluation is layered so each kind of argument takes the cheapest route:

```text
  exact special values   ->  0, +-1, DirectedInfinity[+-I], ...
  symbolic odd argument  ->  Erf[-x] = -Erf[x]
  machine real           ->  libm   erf
  arbitrary real         ->  MPFR   mpfr_erf
  complex (any precision) ->  the cancellation-aware Maclaurin series
                              (DLMF 7.6.2) evaluated in MPFR with guard
                              bits, so even machine-precision complex
                              results carry full accuracy. A double-complex
                              series is the fallback for USE_MPFR=0 builds.
  everything else        ->  stays symbolic (return NULL)

The series  erf(z) = (2/sqrt(pi)) e^-z^2 Sum_{n>=0} t_n,
  t_0 = z,  t_n = t_{n-1} (2 z^2)/(2n+1),
```

is convergent for every z. For real z all terms share a sign (no cancellation); for complex z the partial sums can reach magnitude ~e^|z|^2 before the e^-z^2 prefactor brings them back, so the MPFR path adds |z|^2/ln2 guard bits to absorb that cancellation exactly.

Attributes: Listable, NumericFunction, Protected.

## Implementation notes

**Algorithm.** `builtin_erf` handles `Erf[z]` and `Erf[z0, z1] = erf(z1) - erf(z0)`. A one-argument call first passes through the shared machine-scalar gate, `sf_machine_scalar_args`. It accepts machine `Real`, exact, or `Complex` leaves, needs at least one inexact leaf and no MPFR, and runs before any special-value check. The kernels live in `sf_machine_error.c`:

- **Real** (`sf_erf_r`). For `|x| < 1/2`, a 13-term Maclaurin polynomial. Otherwise `erf = 1 - e^{-x²} erfcx(x)`, with `erfcx` from 0.5-wide polynomial tables, and for `6 <= x <= 27.4` as a polynomial in `1/x²`. `e^{-x²}` is taken from the exact Dekker square, so there is no `x²·ulp` exponent error. The error is ≤ 1.6 ulp.
- **Complex** (`sf_erf_c`). Odd and conjugate symmetry reduce the call to the first quadrant, then `erf_quadrant` picks one of four methods:
  - `|z| <= 0.8`: the Maclaurin series.
  - `|Re z|, |Im z| < 4` with `XY < 1/4`: the Taylor series about the nearer axis (a Hermite recurrence), so the small component keeps full relative accuracy.
  - `|Re z|, |Im z| < 4` otherwise: `1 - e^{-z²} w(iz)`, with the Faddeeva function `w` from a 0.5×0.5 grid of degree-20 Taylor expansions, and with `e^{-z²}` built from an exact exponent and an exactly reduced phase.
  - Elsewhere: Salzer's series (A&S 7.1.29), an exponentially convergent trapezoidal sum whose relative error is `2e^{-4π²}`.

  It declines `|Im z| > 26`, where `e^{y²}` overflows.

A declined call, and any MPFR call, goes to `erf_one_arg`:

- **Exact special values:** `0`, `±Infinity → ±1`, `±I Infinity → DirectedInfinity[±I]`, plus `ComplexInfinity` and `Indeterminate`.
- **Odd fold:** a symbolic negative-leading argument folds by oddness (`Erf[-x] = -Erf[x]`).
- **Arbitrary-precision real:** `mpfr_erf`.
- **Arbitrary-precision complex:** the Maclaurin series `erf(z) = (2/√π) e^{-z²} Σ t_n` with `t_n = t_{n-1}·2z²/(2n+1)`, summed in the file-local `ecx` toolkit with `|z|²/ln2` guard bits to absorb the `~e^{|z|²}` cancellation.

`Erf[z0, z1]` commits only when both erf calls reduce to something concrete. An `Interval` goes to `interval_apply_function`.

**Data structures.** `Expr`; `ecx` (a pair of `mpfr_t`) for the MPFR complex series. The machine kernels use only `<math.h>` (no `<complex.h>`) with double-double helpers; the generated tables are in `sf_machine_error_tables.h` (`tools/sf_fit/gen_tables.py`). The ND kernel is `NDKU_Erf = { sf_erf_c, ndk_Erf_r }`, registered with `REG_U`, covering real and complex packed arrays and `NDArray`s. `Compile[]` lowers `Erf` with the dedicated opcodes `OP_ERF_R` and `OP_VERF_R` for reals, and through the kernel for complex values, at scalar and rank-1 shapes. Every surface calls the same helper, so scalar, `List`, packed, `NDArray` and `Compile[]` agree bit for bit.

**Complexity / limits.** The machine paths are `O(1)`. Salzer's series needs about `2|y| + 13` terms. Against mpmath the complex error has a median of 0.4 ulp and is ≤ 3.5 ulp away from the zeros. The MPFR series' term count and guard bits grow with `|z|²` and with the precision. Known gap: when `|Im z| > 26` and the result exceeds the double range, the machine kernel declines and the MPFR fallback answers `Overflow[]`. `Erfc` and `Erfi` return extended-exponent values in the same region. Exact non-special arguments (`Erf[2]`) stay unevaluated. Attributes: `Listable`, `NumericFunction`, `Protected`.

- Exact special values: `Erf[0] = 0`, `Erf[Infinity] = 1`,
  `Erf[-Infinity] = -1`, `Erf[I Infinity] = DirectedInfinity[I]`,
  `Erf[-I Infinity] = DirectedInfinity[-I]`, plus `ComplexInfinity` and
  `Indeterminate` pass through.
- Odd symmetry for symbolic arguments: `Erf[-x] = -Erf[x]`,
  `Erf[-2 x] = -Erf[2 x]`.
- Numeric evaluation:
  - **Machine precision** (real and complex) → the raw-double kernels of
    `src/special_functions/sf_machine_error.c`, shared bit-for-bit by the
    scalar builtin, the packed/`NDArray` element-wise map and `Compile[]`
    (scalar and rank-1, real and complex). The machine gate runs before any
    special-value ladder.
    - Real: |x| < 1/2 the Maclaurin series; otherwise `erf = 1 − e^(−x²)·erfcx(x)`
      with erfcx from 0.5-wide polynomial tables (6 ≤ x ≤ 27.4: x·erfcx as a
      polynomial in 1/x²), e^(−x²) from the exact (Dekker) square — ≤ 1.6 ulp,
      slightly faster than libm. E.g. `Erf[0.95] = 0.820891`, `Erf[1.5] = 0.966105`.
    - Complex: |z| ≤ 0.8 the Maclaurin series (19 terms, Estrin); next to an
      axis (2|xy| < 1/2) the Taylor series about that axis (Hermite recurrence),
      so the SMALL component keeps full relative accuracy
      (`Erf[1.*^-300 + 2. I] = 6.16074*10^-299 + 18.5648 I`, not `1.1*10^-16 + …`);
      |Re z|, |Im z| < 4 `erf = 1 − e^(−z²) w(iz)` with the Faddeeva function w
      from a 0.5×0.5 grid of degree-20 Taylor expansions (generated at 60 digits)
      and e^(−z²) from an exact exponent and exactly reduced phase; elsewhere
      Salzer's series (A&S 7.1.29, an exponentially convergent trapezoidal sum,
      error 2e^(−4π²) relative), with `cosh(nY)` generated outward from the
      Gaussian peak and compensated sums. Odd/conjugate symmetry fold the plane.
      Measured against mpmath: median 0.4 ulp, ≤ 3.5 ulp away from the complex
      zeros (absolute error scaled by |z erf′(z)| ≤ 2.7 ulp everywhere).
      `|Im z| > 26` (e^(y²) overflow) declines to the MPFR path.
  - Arbitrary precision (MPFR) real → `mpfr_erf`, output precision tracking the
    input, e.g.
    `N[Erf[3/2], 50] = 0.96610514647531072706697626164594785868141047925764`
    and `Erf[0.95`100]`.
  - **Complex, arbitrary precision** → the cancellation-aware
    Maclaurin series erf(z) = (2/√π) e^(−z²) Σ 2ⁿ z^(2n+1)/(1·3···(2n+1))
    (DLMF 7.6.2), evaluated in MPFR with `|z|²/ln2` guard bits, e.g.
    `Erf[1.5 - I] = 1.0784 + 0.0279637 I`,
    `N[Erf[1/2 + I], 30] = 1.20484755831421800270211268210 + 1.02440088160844588172486045441 I`.
    A `double complex` series is the fallback for `USE_MPFR=0` builds.
- The two-argument form `Erf[z0, z1]` reduces to `erf(z1) − erf(z0)` only when
  both reduce to something concrete, e.g. `Erf[1.5, 2] = 0.0292171`,
  `Erf[-Infinity, Infinity] = 2`; exact/symbolic pairs such as `Erf[2, 3]` and
  `Erf[a, b]` stay unevaluated.
- Derivative: `D[Erf[z], z] = (2/Sqrt[Pi]) E^(−z²)` (chain rule applies), so
  the origin Taylor series follows from the generic `D`-based fallback, e.g.
  `Series[Erf[x], {x, 0, 5}]` begins `2/Sqrt[Pi] x − 2/(3 Sqrt[Pi]) x^3 + …`.
- All other arguments (symbolic `Erf[x]`, exact `Erf[2]`) stay unevaluated.

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## References

**See also:** [NDArray](../../linear-algebra/NDArray/), [D](../../calculus/D/)

- DLMF §7.2, §7.6.2 — the error function and its Maclaurin series erf(z) = (2/sqrt(π)) e^{-z^2} Sum_n t_n.
- H. E. Salzer, *Formulas for calculating the error function of a complex variable*, Math. Tables Aids Comput. 5 (1951) 67–70 (A&S 7.1.29).
- T. J. Dekker, *A floating-point technique for extending the available precision*, Numer. Math. 18 (1971) 224–242.
- Source: [`src/special_functions/erf.c`](https://github.com/stblake/mathilda/blob/main/src/special_functions/erf.c)
- Specification: [`docs/spec/builtins/special-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/special-functions.md)
- Tests: [`tests/test_compile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile.c)
- Tests: [`tests/test_compiledfunction.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compiledfunction.c)
- Tests: [`tests/test_dsolve.c`](https://github.com/stblake/mathilda/blob/main/tests/test_dsolve.c)
- Tests: [`tests/test_erf.c`](https://github.com/stblake/mathilda/blob/main/tests/test_erf.c)

## Notes & additional examples

### Notes

- `Erf[z] = (2/Sqrt[Pi]) ∫₀ᶻ e^(-t²) dt` is entire and odd; `Erf[z0, z1] = erf(z1) - erf(z0)` reduces only when both values are concrete.
- Machine real and complex arguments use the `sf_machine_error.c` kernels (erfcx tables; Maclaurin, axis Taylor series, a Faddeeva-function table and Salzer's series in the complex plane), shared bit for bit by the scalar call, packed lists, `NDArray` and `Compile[]`.
- Arbitrary precision uses `mpfr_erf` (real) and a guarded DLMF 7.6.2 series (complex).
- Known gap: a machine complex argument with `|Im z| > 26` whose value leaves the double range answers `Overflow[]` (e.g. `Erf[3. + 30. I]`), where `Erfc` and `Erfi` return extended-exponent numbers.
- Attributes: `Listable`, `NumericFunction`, `Protected`. See also `Erfc`, `Erfi`, `InverseErf`, `FresnelC`, `FresnelS`.
