# Erfc

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Erfc[z]`**

gives the complementary error function erfc(z) = 1 - erf(z).

**`Erfc[0] = 1, Erfc[Infinity] = 0, Erfc[-Infinity] = 2. An entire`**

**`D[Erfc[z], z] = -(2/Sqrt[Pi]) E^(-z^2). Listable.`**

<details>
<summary>Notes</summary>

function. Real inputs evaluate via libm/MPFR erfc (cancellation-free); complex inputs via 1 - erf(z) at machine or arbitrary (MPFR) precision.

</details>

## Examples (12)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= Erfc[10.]
Out[1]= 2.08849e-45

In[2]:= Erfc[1.5 - I]
Out[2]= -0.0783992 - 0.0279637*I

In[3]:= Erfc[-x]
Out[3]= Erfc[-x]
```

### Applications (9)

Exact special values

```mathematica
In[4]:= {Erfc[0], Erfc[Infinity], Erfc[-Infinity]}
Out[4]= {1, 0, 2}
```

Computed as e^(-x^2) erfcx(x), with no cancellation against 1

```mathematica
In[5]:= Erfc[10.]
Out[5]= 2.08849e-45
```

The naive complement loses every digit

```mathematica
In[6]:= 1 - Erf[10.]
Out[6]= 0.0
```

The tail is rounded once into the subnormal range

```mathematica
In[7]:= Erfc[27.]
Out[7]= 5.23705e-319
```

Arbitrary-precision real through mpfr_erfc

```mathematica
In[8]:= N[Erfc[3/2], 50]
Out[8]= 0.0338948535246892729330237383540521413185895207423632
```

Machine complex on the shared kernel

```mathematica
In[9]:= Erfc[1.5 - I]
Out[9]= -0.0783992 - 0.0279637*I
```

Arbitrary-precision complex

```mathematica
In[10]:= N[Erfc[1 + I], 25]
Out[10]= -0.31615128169794764488027107 - 0.19045346923783468628410886*I
```

The Maclaurin series of 1 - erf

```mathematica
In[11]:= Series[Erfc[x], {x, 0, 5}]
Out[11]= 1 - (2 x)/Sqrt[Pi] + (2 x^3)/(3 Sqrt[Pi]) - x^5/(5 Sqrt[Pi]) + O[x]^6
```

No symmetry fold: erfc(-x) = 2 - erfc(x) is not applied

```mathematica
In[12]:= Erfc[-x]
Out[12]= Erfc[-x]
```

## Algorithm

Mathilda -- the complementary error function.

```text
  Erfc[z]   complementary error function   erfc(z) = 1 - erf(z)
```

erfc is an entire function (no branch cuts). It is the complement of erf; unlike erf it has no symmetry that simplifies Erfc[-x] (erfc(-x) = 2 - erfc(x), which Mathilda leaves unexpanded). Evaluation is layered so each kind of argument takes the cheapest, most accurate route:

```text
  exact special values    ->  1, 0, 2, DirectedInfinity[-+I], ...
  machine real            ->  libm   erfc
  arbitrary real          ->  MPFR   mpfr_erfc   (cancellation-free even for
                              large positive z, where 1 - erf(z) would lose
                              all significance)
  complex (any precision) ->  1 - erf(z), with erf(z) from the
                              cancellation-aware Maclaurin series (DLMF
                              7.6.2) evaluated in MPFR with guard bits; the
                              complement is formed at working precision
                              before rounding, so even machine-precision
                              complex results carry full accuracy. A
                              double-complex series is the USE_MPFR=0
                              fallback.
  everything else         ->  stays symbolic (return NULL)

The erf series  erf(z) = (2/sqrt(pi)) e^-z^2 Sum_{n>=0} t_n,
  t_0 = z,  t_n = t_{n-1} (2 z^2)/(2n+1),
```

is convergent for every z. For complex z the partial sums can reach magnitude ~e^|z|^2 before the e^-z^2 prefactor brings them back, so the MPFR path adds |z|^2/ln2 guard bits to absorb that cancellation exactly.

Attributes: Listable, NumericFunction, Protected.

## Implementation notes

**Algorithm.** `builtin_erfc` handles `Erfc[z] = 1 - erf(z)`. It is entire and, unlike `Erf`, has no symmetry fold: `Erfc[-x]` stays as it is. A machine call first passes through the `sf_machine_scalar_args` gate and the shared kernels in `sf_machine_error.c`:

- **Real** (`sf_erfc_r`). For `x >= 1/2`, `erfc = e^{-x²} erfcx(x)`, with `erfcx` from the same tables `Erf` uses and `e^{-x²}` from the exact Dekker square, so nothing cancels. The error is ≤ 3.1 ulp, against 5.1 ulp measured for the macOS libm. Past `x ≈ 26.5` the value drops into the subnormal range in a single rounding (`Erfc[27.] ≈ 5.2*10^-319`), and it reaches 0 at `x = 27.4`. For `x < 0` the kernel returns `2 - erfc(-x)`, and for `|x| < 1/2` it returns `1 - erf`.
- **Complex** (`sf_erfc_c`), for `Re z >= 0` outside a small disc:
  - Next to the real axis: `erfc(X) - D`, using the axis Taylor series, so each component keeps its own relative accuracy.
  - Next to the imaginary axis: `1 - erf`.
  - Inside the `|Re|, |Im| < 4` box: `e^{-z²} w(iz)`, from the Faddeeva table.
  - Beyond the box: `erfc(X) - (erf(z) - erf(X))`, taken from Salzer's series, so a small `erfc` does not cancel against 1.

  The kernel declines when `Re z > 26` or `|Im z| > 26`. For `Re z < 0` it uses `1 - erf(z)`, which lies between 1 and 2, so nothing cancels.

Anything else goes to `erfc_one_arg`:

- **Exact special values:** `0 → 1`, `Infinity → 0`, `-Infinity → 2`, `±I Infinity → DirectedInfinity[∓I]`, plus `ComplexInfinity` and `Indeterminate`.
- **Arbitrary-precision real:** `mpfr_erfc`, which does not cancel even for large `z`.
- **Arbitrary-precision complex:** `erf` from the guarded DLMF 7.6.2 series in `ecx`, with the complement formed at working precision.

The result passes through `numeric_promote_result_if_degenerate`, so a value that leaves the double range comes back as an extended-exponent number. An `Interval` goes to `interval_apply_function`.

**Data structures.** `Expr`; `ecx` (a pair of `mpfr_t`). The ND kernel is `NDKU_Erfc = { sf_erfc_c, ndk_Erfc_r }`, registered with `REG_U`, covering real and complex buffers. `Compile[]` uses `OP_ERFC_R` and `OP_VERFC_R` for reals and the kernel for complex values, at scalar and rank-1 shapes. Every surface shares one helper, so the answers agree bit for bit.

**Complexity / limits.** The machine paths are `O(1)`. The MPFR complex series costs more as `|z|²` and the precision grow. Exact non-special arguments (`Erfc[2]`) and symbolic arguments stay unevaluated, with no `2 - Erfc[x]` rewrite. Attributes: `Listable`, `NumericFunction`, `Protected`.

- Exact special values: `Erfc[0] = 1`, `Erfc[Infinity] = 0`,
  `Erfc[-Infinity] = 2`, `Erfc[I Infinity] = DirectedInfinity[-I]`,
  `Erfc[-I Infinity] = DirectedInfinity[I]` (negated relative to `Erf`), plus
  `ComplexInfinity` and `Indeterminate` pass through.
- Numeric evaluation:
  - **Machine precision** (real and complex) → the shared kernels of
    `sf_machine_error.c` (see `Erf`): real `erfc = e^(−x²)·erfcx(x)` for x ≥ 1/2
    (no cancellation; ≤ 3.1 ulp where the macOS libm measured 5.1), `2 − erfc(−x)`
    for x < 0; complex `erfc = e^(−z²) w(iz)` for Re z ≥ 0 inside the table box,
    `erfc(X) − (erf(z) − erf(X))` from Salzer's series beyond it (Re z > 26
    declines), `1 − erf(z)` for Re z < 0. E.g. `Erfc[0.95] = 0.179109`,
    `Erfc[1.5] = 0.0338949`, `Erfc[1.5 - I] = -0.0783992 - 0.0279637 I`.
  - Arbitrary precision (MPFR) real → `mpfr_erfc`, which is cancellation-free
    even for large positive z (where `1 − erf(z)` would lose all significance),
    output precision tracking the input, e.g.
    `N[Erfc[3/2], 50] = 0.033894853524689272933023738354052141318589520742363`.
  - **Complex, arbitrary precision** → `1 − erf(z)`, with erf(z)
    from the cancellation-aware DLMF 7.6.2 series evaluated in MPFR; the
    complement is formed at working precision (with `|z|²/ln2` guard bits) before
    rounding. A `double complex` series is the fallback for `USE_MPFR=0` builds.
- Derivative: `D[Erfc[z], z] = -(2/Sqrt[Pi]) E^(−z²)` (chain rule applies), so
  the origin Taylor series follows from the generic `D`-based fallback, e.g.
  `Series[Erfc[x], {x, 0, 3}]` begins `1 − 2/Sqrt[Pi] x + …`.
- All other arguments (symbolic `Erfc[x]`, exact `Erfc[2]`) stay unevaluated.

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## References

**See also:** [Erf](../../special-functions/Erf/), [D](../../calculus/D/)

- DLMF §7.2, §7.6.2 — the complementary error function erfc(z) = 1 - erf(z).
- H. E. Salzer, *Formulas for calculating the error function of a complex variable*, Math. Tables Aids Comput. 5 (1951) 67–70 (A&S 7.1.29).
- Source: [`src/special_functions/erfc.c`](https://github.com/stblake/mathilda/blob/main/src/special_functions/erfc.c)
- Specification: [`docs/spec/builtins/special-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/special-functions.md)
- Tests: [`tests/test_erfc.c`](https://github.com/stblake/mathilda/blob/main/tests/test_erfc.c)
- Tests: [`tests/test_fullsimplify.c`](https://github.com/stblake/mathilda/blob/main/tests/test_fullsimplify.c)
- Tests: [`tests/test_gruntz.c`](https://github.com/stblake/mathilda/blob/main/tests/test_gruntz.c)
- Tests: [`tests/test_gruntz_stress.c`](https://github.com/stblake/mathilda/blob/main/tests/test_gruntz_stress.c)

## Notes & additional examples

### Notes

- `Erfc[z] = 1 - Erf[z]` is entire; unlike `Erf` it is not odd, and `Erfc[-x]` is left as written.
- Machine arguments use the `sf_machine_error.c` kernels shared bit for bit with packed lists, `NDArray` and `Compile[]`: `e^(-x²) erfcx(x)` on the real line (accurate far into the tail), and `e^(-z²) w(iz)` or Salzer's series in the complex plane.
- Arbitrary precision uses `mpfr_erfc` (real) and the guarded complex erf series with the complement formed at working precision.
- Attributes: `Listable`, `NumericFunction`, `Protected`. See also `Erf`, `Erfi`, `InverseErfc`.
