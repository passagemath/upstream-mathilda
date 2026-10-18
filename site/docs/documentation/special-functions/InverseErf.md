# InverseErf

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`InverseErf[s]`**

gives the inverse error function: the z solving s = Erf\[z\].

**`InverseErf[z0, s]`**

gives the inverse of the generalized error function Erf\[z0, z\].

**`InverseErf[0] = 0, InverseErf[1] = Infinity, InverseErf[-1] = -Infinity.`**

<details>
<summary>Notes</summary>

Odd in s. Numerical values are given only for real s in \[-1, 1\], at machine or arbitrary (MPFR) precision; D\[InverseErf\[z\], z\] = (Sqrt\[Pi\]/2) E^(InverseErf\[z\]^2). Listable.

</details>

## Examples (12)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= InverseErf[0.6]
Out[1]= 0.595116

In[2]:= InverseErf[0.4, 0.2]
Out[2]= 0.631776

In[3]:= InverseErf[1]
Out[3]= Infinity
```

### Applications (9)

Exact values at the centre and the endpoints

```mathematica
In[4]:= {InverseErf[0], InverseErf[1], InverseErf[-1]}
Out[4]= {0, Infinity, -Infinity}
```

Machine real on the piecewise-polynomial kernel

```mathematica
In[5]:= InverseErf[0.6]
Out[5]= 0.595116
```

Listable, element by element

```mathematica
In[6]:= InverseErf[1/{2., 3., 4., 5.}]
Out[6]= {0.476936, 0.30457, 0.225312, 0.179143}
```

Deep in the tail, and the round trip back through Erf

```mathematica
In[7]:= {InverseErf[0.999999], Erf[InverseErf[0.999999]]}
Out[7]= {3.45891, 0.999999}
```

Arbitrary precision by Newton's method with precision doubling

```mathematica
In[8]:= N[InverseErf[33/100], 50]
Out[8]= 0.301332146133705826128502718158394773965824282828535
```

The 95% quantile of the standard normal distribution

```mathematica
In[9]:= N[Sqrt[2] InverseErf[2 (95/100) - 1], 30]
Out[9]= 1.644853626951472714863848907992
```

The generalized form: the z with Erf[0.4, z] = 0.2

```mathematica
In[10]:= InverseErf[0.4, 0.2]
Out[10]= 0.631776
```

The Taylor series through the derivative rule

```mathematica
In[11]:= Series[InverseErf[x], {x, 0, 7}]
Out[11]= 1/2 Sqrt[Pi] x + 1/24 Pi^(3/2) x^3 + 7/960 Pi^(5/2) x^5 + 127/80640 Pi^(7/2) x^7 + O[x]^8
```

Outside [-1, 1], and complex arguments, stay unevaluated

```mathematica
In[12]:= {InverseErf[2], InverseErf[0.5 + 0.1 I]}
Out[12]= {InverseErf[2], InverseErf[0.5 + 0.1*I]}
```

## Algorithm

Mathilda -- the inverse error function.

```text
  InverseErf[s]       inverse of erf: the z solving s = erf(z)
  InverseErf[z0, s]   inverse of the generalized error function:
                      the z solving s = Erf[z0, z] = erf(z) - erf(z0)
```

Per Mathematica, explicit numerical values are produced only for *real* s in [-1, 1]; complex and out-of-domain (|s| > 1) inputs stay symbolic. Evaluation is therefore much simpler than Erf -- no complex machinery:

```text
  exact special values   ->  InverseErf[0]=0, InverseErf[1]=Infinity,
                             InverseErf[-1]=-Infinity
  symbolic odd argument  ->  InverseErf[-x] = -InverseErf[x]
  machine real (|s|<1)   ->  Winitzki seed + Newton polish on libm erf
  arbitrary real (|s|<1) ->  Newton with precision doubling on mpfr_erf
  everything else        ->  stays symbolic (return NULL)
```

Neither C99 libm nor MPFR ships an inverse-erf, so the kernel is Newton's method on f(z) = erf(z) - s, with f'(z) = (2/sqrt(pi)) e^{-z^2}, i.e.

```text
  z <- z - (erf(z) - s) (sqrt(pi)/2) e^{z^2}.
```

Newton converges quadratically; the MPFR path doubles the working precision each step so the final full-precision erf dominates the cost.

The derivative D[InverseErf[z], z] = (sqrt(pi)/2) e^{InverseErf[z]^2} lives in src/calculus/deriv.c; Series follows from it via the generic Taylor-via-D fallback.

Attributes: Listable, NumericFunction, Protected.

## Implementation notes

**Algorithm.** `builtin_inverf` solves `s = erf(z)`. Following Mathematica, it returns explicit values only for real `s` in `[-1, 1]`; complex `s` and `|s| > 1` stay unevaluated. A one-argument machine-real call (passing the `sf_machine_scalar_args` gate, not complex) goes first to `sf_inverse_erf_r` in `sf_machine_error.c`. With `p = |s|` and its exact complement `q = 1 - p`:

- For `p < 0.875`, a degree-14 polynomial for `erfinv(p)/p` on 1/16-wide pieces, with no logarithm.
- Otherwise the kernel switches to Giles' variable `w = -log(q(2 - q))`, in which `erfinv(p)/p` is analytic in the strip `|Im w| < π`. For `w < 6` it uses unit-width degree-11 pieces. Past that it takes a `1e-7` seed polynomial in `1/√w`, then one Halley step on `erfc(z) = q`, written with the relative residual and `erfcx` so that neither `e^{±z²}` is ever formed.

No `erf` is evaluated on the main path. The kernel declines `±1`, `|s| > 1` and a subnormal `q`.

Declined calls go to `inverf_one_arg`:

- **Exact values:** `0 → 0`, `±1 → ±Infinity` (also for `±1.`).
- **Odd fold:** `InverseErf[-x] = -InverseErf[x]`.
- **Arbitrary-precision real** `|s| < 1`: Newton's method on `mpfr_erf`, `z ← z - (erf(z) - s)(√π/2) e^{z²}`, doubling the working precision at each step. The seed comes from `inverf_double`, which is the polynomial kernel above, falling back to Winitzki's closed form plus libm-`erf` Newton steps where the kernel declines.

The two-argument form `InverseErf[z0, s]` rewrites to `InverseErf[s + Erf[z0]]`, and stays in two-argument form if `Erf[z0]` does not reduce. The derivative `(√π/2) e^{InverseErf[z]²}` is defined in `calculus/deriv.c`, and `Series` builds on it through the generic Taylor-via-`D` route.

**Data structures.** `Expr`; `mpfr_t` for the Newton loop; the generated polynomial tables `INVP_T`, `INVW_T` and `INVT_T` in `sf_machine_error_tables.h`. The ND kernel is real-only, `NDKU_InverseErf = { NULL, ndk_InverseErf_r }`, registered with `REG_U`, so it is packed-aware. `Compile[]` lowers `InverseErf` at real scalar and rank-1 shapes through that kernel; a `_Complex` argument does not compile, since the function has no complex values. Scalar, packed, `NDArray` and `Compile[]` share `sf_inverse_erf_r` and agree bit for bit.

**Complexity / limits.** The machine kernel is a few polynomial evaluations, with error ≤ 1.5 ulp. The old kernel was Winitzki's closed form plus up to four Newton steps, each costing an `erf` and an `exp`. The MPFR cost is dominated by the final full-precision `mpfr_erf`. Exact rationals (`InverseErf[1/2]`), complex inputs and `|s| > 1` stay unevaluated. Attributes: `Listable`, `NumericFunction`, `Protected`.

- Exact special values: `InverseErf[0] = 0`, `InverseErf[1] = Infinity`,
  `InverseErf[-1] = -Infinity` (likewise `InverseErf[1.] = Infinity`), plus
  `Indeterminate` passes through.
- Odd symmetry for symbolic arguments: `InverseErf[-x] = -InverseErf[x]`.
- Numeric evaluation (neither C99 libm nor MPFR ships an inverse-erf, so the
  kernel is Newton's iteration on `f(z) = erf(z) − s`,
  `z ← z − (erf(z) − s)(Sqrt[Pi]/2) e^{z²}`):
  - Machine-precision real → `sf_inverse_erf_r` (`sf_machine_error.c`, shared
    with the packed/`NDArray` map and `Compile[]`): |s| < 0.875 a degree-14
    polynomial for erfinv(s)/s on 1/16-wide pieces (no logarithm); beyond, in
    Giles' variable w = −log(1 − s²) (erfinv(s)/s is analytic in the strip
    |Im w| < π), unit-width degree-11 pieces for w < 6 and, past that, a 1e-7 seed
    in 1/√w plus one Halley step on erfc written with the relative residual. No
    erf evaluation on the main path; ≤ 1.5 ulp (the old Winitzki + up-to-four
    Newton steps cost an erf and an exp each). E.g. `InverseErf[0.6] = 0.595116`,
    `InverseErf[1/{2.,3.,4.,5.}] = {0.476936, 0.30457, 0.225312, 0.179143}`.
  - Arbitrary precision (MPFR) real → Newton with precision doubling on
    `mpfr_erf`, output precision tracking the input, e.g.
    `N[InverseErf[33/100], 50] = 0.30133214613370582612850271815839477396582428282853`.
- Two-argument form: `InverseErf[0.4, 0.2] = 0.631776`; if `Erf[z0]` does not
  reduce the call stays in two-argument form, while a reducible `z0` collapses to
  the one-argument inverse, e.g. `InverseErf[0, 1.3] = InverseErf[1.3]`.
- Derivative: `D[InverseErf[z], z] = (Sqrt[Pi]/2) E^(InverseErf[z]²)` (chain rule
  applies), so the origin Taylor series follows from the generic `D`-based
  fallback, e.g. `Series[InverseErf[x], {x, 0, 8}] = (Sqrt[Pi]/2) x +
  (Pi^(3/2)/24) x³ + (7 Pi^(5/2)/960) x⁵ + (127 Pi^(7/2)/80640) x⁷ + O[x]^9`.
- All other arguments (symbolic `InverseErf[x]`, exact `InverseErf[1/2]`,
  out-of-domain `InverseErf[2]`) stay unevaluated.

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## References

**See also:** [NDArray](../../linear-algebra/NDArray/), [D](../../calculus/D/)

- DLMF §7.17 — inverse error functions.
- M. Giles, *Approximating the erfinv function*, in GPU Computing Gems Jade Edition (Morgan Kaufmann, 2011) — the w = -log(1 - p^2) variable.
- S. Winitzki, *A handy approximation for the error function and its inverse* (2008) — the closed-form fallback seed.
- Source: [`src/special_functions/inverf.c`](https://github.com/stblake/mathilda/blob/main/src/special_functions/inverf.c)
- Specification: [`docs/spec/builtins/special-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/special-functions.md)
- Tests: [`tests/test_compile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile.c)
- Tests: [`tests/test_interval.c`](https://github.com/stblake/mathilda/blob/main/tests/test_interval.c)
- Tests: [`tests/test_inverf.c`](https://github.com/stblake/mathilda/blob/main/tests/test_inverf.c)
- Tests: [`tests/test_numeric_stress.c`](https://github.com/stblake/mathilda/blob/main/tests/test_numeric_stress.c)

## Notes & additional examples

### Notes

- `InverseErf[s]` solves `s = Erf[z]`; as in Mathematica it gives explicit values only for real `s` in `[-1, 1]`, and `InverseErf[z0, s]` inverts the generalized `Erf[z0, z]`.
- Machine reals use `sf_inverse_erf_r` (polynomials in `s` and in Giles' variable `w = -log(1 - s²)`, one Halley step in the far tail), shared bit for bit with packed lists, `NDArray` and `Compile[]`; arbitrary precision uses Newton's method on `mpfr_erf`.
- `D[InverseErf[z], z] = (Sqrt[Pi]/2) E^(InverseErf[z]^2)`. Exact rationals such as `InverseErf[1/2]` stay unevaluated.
- Attributes: `Listable`, `NumericFunction`, `Protected`. See also `InverseErfc`, `Erf`, `Quantile`.
