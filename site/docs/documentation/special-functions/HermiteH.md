# HermiteH

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`HermiteH[n, x]`**

gives the physicists' Hermite polynomial H\_n(x).

<details>
<summary>Notes</summary>

A non-negative integer n yields the explicit polynomial (exact for an exact x, and via a stable recurrence for an inexact x -- the recurrence avoids the cancellation that corrupts the polynomial form for large n, e.g. HermiteH\[100, N\[5,16\]\]); a non-integer, complex, or negative order with an inexact argument evaluates numerically at machine or arbitrary (MPFR) precision, real or complex, through the confluent-hypergeometric form 2^n Sqrt\[Pi\] (Hypergeometric1F1\[-n/2, 1/2, x^2\]/Gamma\[(1-n)/2\] - 2 x Hypergeometric1F1\[(1-n)/2, 3/2, x^2\]/Gamma\[-n/2\]). The special value H\_n(0) = 2^n Sqrt\[Pi\]/Gamma\[(1-n)/2\]. D\[HermiteH\[n,x\],x\] = 2 n HermiteH\[n-1,x\]; Series (including at Infinity) and SeriesCoefficient are supported. Listable.

</details>

## Examples (16)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= HermiteH[10, x]
Out[1]= -30240 + 302400 x^2 - 403200 x^4 + 161280 x^6 - 23040 x^8 + 1024 x^10

In[2]:= Table[HermiteH[10, x], {x, 1, 5}]
Out[2]= {8224, 200416, -3093984, 125984224, 3275529760}

In[3]:= HermiteH[100, N[5, 16]]
Out[3]= 8.1931937130591502e+98

In[4]:= D[HermiteH[n, z], z]
Out[4]= 2 n HermiteH[-1 + n, z]
```

### Applications (12)

The physicists' polynomial with integer coefficients

```mathematica
In[5]:= HermiteH[6, x]
Out[5]= -120 + 720 x^2 - 480 x^4 + 64 x^6
```

The value at the origin for a symbolic order

```mathematica
In[6]:= HermiteH[n, 0]
Out[6]= (Sqrt[Pi] 2^n)/Gamma[1/2 (1 - n)]
```

The Hermite equation y'' - 2x y' + 2n y = 0

```mathematica
In[7]:= Simplify[D[HermiteH[5, x], {x, 2}] - 2 x D[HermiteH[5, x], x] + 10 HermiteH[5, x]]
Out[7]= 0
```

The norm 2^n n! Sqrt[Pi]

```mathematica
In[8]:= Integrate[HermiteH[3, x]^2 Exp[-x^2], {x, -Infinity, Infinity}]
Out[8]= 48 Sqrt[Pi]
```

Distinct degrees are orthogonal

```mathematica
In[9]:= Integrate[HermiteH[3, x] HermiteH[4, x] Exp[-x^2], {x, -Infinity, Infinity}]
Out[9]= 0
```

A degree-40 value at machine precision

```mathematica
In[10]:= HermiteH[40, 2.5]
Out[10]= -7.60498e+30
```

The same value at 20 digits

```mathematica
In[11]:= N[HermiteH[40, 5/2], 20]
Out[11]= -7.60498442897757420469e+30
```

A fractional order through the 1F1 form

```mathematica
In[12]:= HermiteH[1/2, 0.8]
Out[12]= 1.35032
```

A negative order: H_{-1}(x) = Sqrt[Pi] e^(x^2) Erfc[x]/2

```mathematica
In[13]:= {HermiteH[-1, 0.5], Sqrt[Pi] Exp[0.25] Erfc[0.5]/2}
Out[13]= {0.545641, 0.545641}
```

A packed Range input runs on the buffer kernel

```mathematica
In[14]:= HermiteH[4, Range[0., 2., 0.25]]
Out[14]= {12.0, 9.0625, 1.0, -9.9375, -20.0, -23.9375, -15.0, 15.0625, 76.0}
```

The algebraic asymptotic series at infinity

```mathematica
In[15]:= Series[HermiteH[n, x], {x, Infinity, 2}]
Out[15]= (2 x)^n (1 + 1/4 n (1 - n) (1/x)^2 + O[1/x]^4)
```

A fractional order with a symbolic argument stays unevaluated

```mathematica
In[16]:= HermiteH[1/2, x]
Out[16]= HermiteH[1/2, x]
```

## Algorithm

Mathilda -- HermiteH, the physicists' Hermite polynomial H_n(x).

```text
  HermiteH[n, x]

The Hermite polynomials satisfy  y'' - 2x y' + 2n y = 0  and are orthogonal on
```

(-Inf, Inf) with weight e^(-x^2). Evaluation is layered so each argument shape takes the cheapest exact or numeric route (mirroring the ChebyshevU design in chebyshev_u.c, whose integer-coefficient, two-argument substrate HermiteH shares via chebyshev_common.{c,h}):

```text
  non-negative integer order n
    symbolic / exact-number x  ->  the explicit degree-n polynomial in x,
                                   with exact integer coefficients, from
                                   FLINT's fmpz_poly_hermite_h (or, without
                                   FLINT, the three-term recurrence); an exact
                                   numeric x then evaluates to an exact number.
    inexact x                  ->  the stable three-term recurrence evaluated
                                   in machine/MPFR complex arithmetic, NOT by
                                   substituting into the polynomial (whose
                                   cancellation ruins accuracy for large n,
                                   e.g. HermiteH[100, N[5,16]]).
  non-integer / complex / negative order, some inexact argument
                               ->  H_n(z) = 2^n Sqrt[Pi] (
                                       1F1(-n/2; 1/2; z^2) / Gamma[(1-n)/2]
                                     - 2z 1F1((1-n)/2; 3/2; z^2) / Gamma[-n/2]),
                                   built as an expression over
                                   Hypergeometric1F1 so it inherits its
                                   machine + MPFR, real + complex numerics.
  special value
    H_0(x) = 1,  H_n(0) = 2^n Sqrt[Pi] / Gamma[(1-n)/2].
  everything else              ->  stays symbolic (return NULL).
```

Attributes: Listable, NumericFunction, Protected.

The machine-precision double / double-complex recurrence lives in sf_machine.c (sf_machine_hermite_h[_c]) and is shared bit-for-bit with the NDArray element-wise kernel and the Compile[] VM (ndkernels.c). The scalar machine path here dispatches to the same helpers.

Series / SeriesCoefficient: the derivative rule D[HermiteH[n,z],z] = 2 n HermiteH[n-1,z] lives in calculus/deriv.c and, with H_n(0) above, gives the Taylor expansion at a finite point for free via series.c's Taylor-via-D; the asymptotic expansion at Infinity and the symbolic-index SeriesCoefficient are bespoke kernels in calculus/series.c.

Deferred (left symbolic): negative-integer order with an exact argument (no Erfc closed form; a numeric argument is handled by the 1F1 form above), and bespoke tight Interval enclosures for non-integer order (the generic NumericFunction + D[] interval threading in eval.c handles those best-effort; an integer order threads exactly through the explicit polynomial).

## Implementation notes

**Algorithm.** `builtin_hermite_h` calls `hermite_h_two_arg`. (1) A non-integer order at `x = 0` returns `2^n Sqrt[Pi]/Gamma[(1-n)/2]` (`hermite_value_at_zero`). (2) A non-negative integer order with an inexact `x`, at machine precision with degree ≤ 1000, goes to `sf_machine_hermite_h` / `sf_machine_hermite_h_c`, the real and complex raw-double recurrences `H_k = 2x H_{k-1} - 2(k-1) H_{k-2}`. Otherwise it runs `hermite_numeric_recurrence` in MPFR. (3) An exact non-negative integer order with an exact or symbolic `x` builds the polynomial from FLINT's `fmpz_poly_hermite_h` (the integer recurrence when FLINT is absent). (4) A non-integer, complex or negative order with an inexact argument goes to `hermite_numeric_hyp`, which builds `2^n Sqrt[Pi] (1F1(-n/2; 1/2; z²)/Γ((1-n)/2) - 2z 1F1((1-n)/2; 3/2; z²)/Γ(-n/2))` (`hermite_hyp_form`). It numericalises that form at the output precision plus `hermite_hyp_guard_bits`: about `2·Re(z²)·log₂e + 64` bits, capped at 8192, which absorbs the cancellation between two terms that each grow like `e^{Re z²}`. It then rounds back to the input precision.

**Data structures.** `mpz_t` coefficient arrays hold the exact polynomial, `ncpx` MPFR pairs carry the recurrence, and `NumericSpec` targets (`NUMERIC_MODE_MPFR` then `NUMERIC_MODE_MPFR_CAP`) handle the guarded 1F1 evaluation. The binary ND kernel `ndk_HermiteH_c` (`src/ndkernels.c`) calls the same `sf_machine_hermite_h*` helpers, so packed and `NDArray` buffers, real or complex, stay on the buffer. `Compile[]` lowers the head at scalar and array shapes, and every surface agrees bit for bit.

**Complexity / limits.** The recurrence costs `O(n)` operations. The explicit polynomial is capped at degree 2000 (`HERMITE_POLY_CAP`). A non-integer order with an exact or symbolic argument stays unevaluated, and so does a negative integer order with an exact argument (there is no `Erfc` closed form). The derivative `2 n H_{n-1}` is in `src/calculus/deriv.c`. The asymptotic series at `Infinity` and the symbolic-index `SeriesCoefficient` are bespoke kernels in `src/calculus/series.c`; the asymptotic series omits the subdominant `e^{x²}` branch.

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## References

**See also:** [SeriesCoefficient](../../power-series/SeriesCoefficient/), [Piecewise](../../control-flow/Piecewise/)

- DLMF §18 — classical orthogonal polynomials (§18.9 recurrence, §18.3 orthogonality); §13.2 for the Kummer-function form of H_ν.
- Source: [`src/special_functions/hermite_h.c`](https://github.com/stblake/mathilda/blob/main/src/special_functions/hermite_h.c)
- Specification: [`docs/spec/builtins/special-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/special-functions.md)
- Tests: [`tests/test_hermite_h.c`](https://github.com/stblake/mathilda/blob/main/tests/test_hermite_h.c)
- Tests: [`tests/test_possiblezeroq_stress.c`](https://github.com/stblake/mathilda/blob/main/tests/test_possiblezeroq_stress.c)

## Notes & additional examples

### Notes

- `HermiteH[n, x]` is the physicists' Hermite polynomial, orthogonal on `(-∞, ∞)` with weight `e^(-x^2)`. For the probabilists' `He_n(x)`, use `2^(-n/2) HermiteH[n, x/Sqrt[2]]`.
- A non-integer, complex or negative order with an inexact argument goes through a two-term `Hypergeometric1F1` form, evaluated with extra working precision that scales with `Re(x^2)` to cover the cancellation.
- Numeric surfaces: the scalar, packed/`NDArray` (real and complex) and `Compile[]` paths share one `double` recurrence for integer orders up to 1000.
- The asymptotic series at `Infinity` omits the subdominant `e^(x^2)` branch.
- Attributes: `Listable`, `NumericFunction`, `Protected`.
- See also: `LaguerreL`, `Hypergeometric1F1`, `Erfc`.
