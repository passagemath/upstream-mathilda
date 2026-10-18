# ChebyshevU

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`ChebyshevU[n, x]`**

gives the Chebyshev polynomial of the second kind U\_n(x).

**`Cos[n Pi/2]`**

and U\_n(1) = n + 1. D\[ChebyshevU\[n,x\],x\] is supported.

<details>
<summary>Notes</summary>

Integer n yields the explicit polynomial (exact for an exact x, and via a stable recurrence for an inexact x); a non-integer or complex order with an inexact argument evaluates numerically at machine or arbitrary (MPFR) precision, real or complex. Special values include U\_n(0) = Listable.

</details>

## Examples (14)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= ChebyshevU[7, x]
Out[1]= -8 x + 80 x^3 - 192 x^5 + 128 x^7

In[2]:= Table[ChebyshevU[n, Pi/2], {n, 0, 3}]
Out[2]= {1, Pi, -1 + Pi^2, -2 Pi + Pi^3}

In[3]:= ChebyshevU[3, 2. + I]
Out[3]= 8.0 + 84.0*I
```

### Applications (11)

The explicit degree-6 polynomial

```mathematica
In[4]:= ChebyshevU[6, x]
Out[4]= -1 + 24 x^2 - 80 x^4 + 64 x^6
```

Negative orders reflect: U_{-1} = 0, U_{-m} = -U_{m-2}

```mathematica
In[5]:= {ChebyshevU[-1, x], ChebyshevU[-4, x]}
Out[5]= {0, 1 - 4 x^2}
```

The defining identity U_n(cos t) sin t = sin((n+1) t)

```mathematica
In[6]:= Simplify[TrigExpand[ChebyshevU[4, Cos[t]] Sin[t] - Sin[5 t]]]
Out[6]= 0
```

The bridge T_n = (U_n - U_{n-2})/2

```mathematica
In[7]:= Expand[ChebyshevT[5, x] - (ChebyshevU[5, x] - ChebyshevU[3, x])/2]
Out[7]= 0
```

A high degree at machine precision via the stable recurrence

```mathematica
In[8]:= {ChebyshevU[60, 0.7], Sin[61 ArcCos[0.7]]/Sin[ArcCos[0.7]]}
Out[8]= {-1.37881, -1.37881}
```

The same value at 25 digits

```mathematica
In[9]:= N[ChebyshevU[60, 7/10], 25]
Out[9]= -1.3788057165596781581949098
```

A non-integer order

```mathematica
In[10]:= {ChebyshevU[0.5, 0.3], Sin[1.5 ArcCos[0.3]]/Sin[ArcCos[0.3]]}
Out[10]= {0.992278, 0.992278}
```

A machine complex argument uses the complex recurrence

```mathematica
In[11]:= ChebyshevU[3, 2. + I]
Out[11]= 8.0 + 84.0*I
```

The orthogonality norm with weight Sqrt[1 - x^2]

```mathematica
In[12]:= Integrate[ChebyshevU[3, x]^2 Sqrt[1 - x^2], {x, -1, 1}]
Out[12]= 1/2 Pi
```

The derivative for a symbolic order

```mathematica
In[13]:= D[ChebyshevU[n, x], x]
Out[13]= (n x ChebyshevU[n, x] + (-1 - n) ChebyshevU[-1 + n, x])/(-1 + x^2)
```

The compiled array result is bit-identical to the packed one

```mathematica
In[14]:= cu = Compile[{{v, _Real, 1}}, ChebyshevU[9, v]]; cu[Range[-1., 1., 0.1]] === ChebyshevU[9, Range[-1., 1., 0.1]]
Out[14]= True
```

## Algorithm

Mathilda -- ChebyshevU, the Chebyshev polynomial of the second kind U_n(x).

```text
  ChebyshevU[n, x]
```

Evaluation is layered so each argument shape takes the cheapest exact or numeric route (mirroring the LegendreP design in legendre.c):

```text
  exact integer order n
    symbolic / exact-number x  ->  the explicit degree-|n'| polynomial in x,
                                   with exact integer coefficients, from
                                   FLINT's fmpz_poly_chebyshev_u (or, without
                                   FLINT, the three-term recurrence); an exact
                                   numeric x then evaluates to an exact number.
    inexact x                  ->  the stable three-term recurrence evaluated
                                   in machine/MPFR complex arithmetic, NOT by
                                   substituting into the polynomial (whose
                                   cancellation ruins accuracy for large n).
  non-integer / complex order, some inexact argument
                               ->  U_v(x) = (v+1) 2F1(-v, v+2; 3/2; (1-x)/2),
                                   built as an expression over Hypergeometric2F1
                                   so it inherits its machine + MPFR, real +
                                   complex numerics (and analytic continuation).
  special values
    U_0(x) = 1,  U_n(0) = Cos[n Pi/2],  U_n(1) = n + 1.
  everything else              ->  stays symbolic (return NULL).
```

Negative integer order reflects to a non-negative degree:

```text
  U_{-1}(x) = 0,  and  U_{-m}(x) = -U_{m-2}(x)  for m >= 2.
```

Attributes: Listable, NumericFunction, Protected.

Deferred (left symbolic): symbolic Series / SeriesCoefficient, D[] of a symbolic order (the D[ChebyshevU[n,x],x] recurrence lives in calculus/deriv.c), half-integer closed forms, and bespoke tight Interval enclosures (the generic NumericFunction + D[] interval threading in eval.c handles those best-effort).

The argument-shape predicates, the explicit-polynomial assembler, the MPFR result builder, the output-precision rule, and U_n(0) = Cos[n Pi/2] are shared with ChebyshevT via chebyshev_common.{c,h}; only the recurrence seed, the negative-order reflection, and the non-integer-order closed forms live here.

## Implementation notes

**Algorithm.** `builtin_chebyshev_u` dispatches through `chebyshev_u_two_arg`, which mirrors `ChebyshevT` and shares its substrate in `chebyshev_common.c`. (1) A non-integer order at `x = 0` or `x = 1` returns `Cos[n Pi/2]` or `n + 1`. (2) An integer order with an inexact `x` is reflected with `U_{-1} = 0` and `U_{-m} = -U_{m-2}` (`cheb_reflect` returns a degree and a sign). At machine precision with degree ≤ 1000 it goes to `sf_machine_chebyshev_u` / `sf_machine_chebyshev_u_c`, the raw-double real and complex recurrences `U_0 = 1`, `U_1 = 2x`, `U_k = 2x U_{k-1} - U_{k-2}`. Otherwise it runs the same recurrence in MPFR. (3) An exact integer order with an exact or symbolic `x` builds the explicit polynomial from FLINT's `fmpz_poly_chebyshev_u` (the integer recurrence when FLINT is absent). (4) A non-integer or complex order with an inexact argument is built as `(v+1) Hypergeometric2F1[-v, v+2, 3/2, (1-x)/2]` when `|(1-x)/2| < 1`, and as `Sin[(v+1) ArcCos[x]]/Sin[ArcCos[x]]` elsewhere.

**Data structures.** These match `ChebyshevT`: `mpz_t` coefficient arrays for the exact polynomial, MPFR or `ncpx` values for the recurrence, and plain `double`s on the machine path. The binary ND kernel `ndk_ChebyshevU_c` (`src/ndkernels.c`) calls the same helpers and has a complex arm, so a packed complex buffer is not materialised. `Compile[]` lowers the head at scalar and array shapes through that registered kernel.

**Complexity / limits.** The numeric recurrence costs `O(n)` operations. The caps are the same as for `ChebyshevT`: degree 2000 for the explicit polynomial and |n| ≤ 2·10⁶ for the recurrence. An exact fractional order with a symbolic argument stays unevaluated, and symbolic series and half-integer closed forms are not implemented. The derivative `((-1-n) U_{n-1} + n x U_n)/(x² - 1)` is produced by `src/calculus/deriv.c`.

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## References

**See also:** [Series](../../power-series/Series/), [SeriesCoefficient](../../power-series/SeriesCoefficient/)

- DLMF §18 — classical orthogonal polynomials (§18.9 recurrences; U_n(cos θ) = sin((n+1)θ)/sin θ and the 2F1 representation).
- Source: [`src/special_functions/chebyshev_u.c`](https://github.com/stblake/mathilda/blob/main/src/special_functions/chebyshev_u.c)
- Specification: [`docs/spec/builtins/special-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/special-functions.md)
- Tests: [`tests/test_chebyshev_t.c`](https://github.com/stblake/mathilda/blob/main/tests/test_chebyshev_t.c)
- Tests: [`tests/test_chebyshev_u.c`](https://github.com/stblake/mathilda/blob/main/tests/test_chebyshev_u.c)

## Notes & additional examples

### Notes

- `ChebyshevU[n, x]` is `U_n(x)`, with `U_n(cos θ) sin θ = sin((n+1)θ)`; it is orthogonal on `[-1, 1]` with weight `Sqrt[1 - x^2]`. Negative orders reflect as in Mathematica.
- Special values hold for a symbolic order: `ChebyshevU[n, 1] = n + 1` and `ChebyshevU[n, 0] = Cos[n Pi/2]`.
- Numeric surfaces: the scalar, packed/`NDArray` (including complex buffers) and `Compile[]` paths share one `double` recurrence. A fractional order goes through `Hypergeometric2F1` or the `ArcCos` closed form.
- Attributes: `Listable`, `NumericFunction`, `Protected`.
- See also: `ChebyshevT`, `GegenbauerC`, `LegendreP`.
