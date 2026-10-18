# ChebyshevT

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`ChebyshevT[n, x]`**

gives the Chebyshev polynomial of the first kind T\_n(x).

**`Cos[n Pi/2]`**

and T\_n(1) = 1. D\[ChebyshevT\[n,x\],x\] = n ChebyshevU\[n-1,x\].

<details>
<summary>Notes</summary>

Integer n yields the explicit polynomial (exact for an exact x, and via a stable recurrence for an inexact x); a non-integer or complex order with an inexact argument evaluates numerically at machine or arbitrary (MPFR) precision, real or complex. Special values include T\_n(0) = Listable.

</details>

## Examples (15)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= ChebyshevT[7, x]
Out[1]= -7 x + 56 x^3 - 112 x^5 + 64 x^7

In[2]:= Table[ChebyshevT[n, Pi/2], {n, 0, 3}]
Out[2]= {1, 1/2 Pi, -1 + 1/2 Pi^2, -3/2 Pi + 1/2 Pi^3}

In[3]:= {ChebyshevT[60, 0.7], Cos[60 ArcCos[0.7]]}
Out[3]= {-0.825313, -0.825313}
```

### Applications (12)

The explicit degree-8 polynomial with integer coefficients

```mathematica
In[4]:= ChebyshevT[8, x]
Out[4]= 1 - 32 x^2 + 160 x^4 - 256 x^6 + 128 x^8
```

The order index is even: T_{-n} = T_n

```mathematica
In[5]:= ChebyshevT[-3, x]
Out[5]= -3 x + 4 x^3
```

The defining identity T_n(cos t) = cos(n t)

```mathematica
In[6]:= Simplify[TrigExpand[ChebyshevT[5, Cos[t]] - Cos[5 t]]]
Out[6]= 0
```

The composition law T_m(T_n(x)) = T_{mn}(x)

```mathematica
In[7]:= Expand[ChebyshevT[3, ChebyshevT[4, x]] - ChebyshevT[12, x]]
Out[7]= 0
```

A high degree at machine precision via the stable recurrence

```mathematica
In[8]:= {ChebyshevT[60, 0.7], Cos[60 ArcCos[0.7]]}
Out[8]= {-0.825313, -0.825313}
```

The same value to 30 digits through the MPFR recurrence

```mathematica
In[9]:= N[ChebyshevT[60, 7/10], 30]
Out[9]= -0.8253130244838482248764016705343
```

A non-integer order, through 2F1 / Cos[v ArcCos[x]]

```mathematica
In[10]:= {ChebyshevT[1/3, 0.5], Cos[1/3 ArcCos[0.5]]}
Out[10]= {0.939693, 0.939693}
```

Outside [-1, 1] a fractional order grows like Cosh[v ArcCosh[x]]

```mathematica
In[11]:= ChebyshevT[2.5, 3.]
Out[11]= 41.0122
```

The orthogonality norm, pi/2 for n >= 1

```mathematica
In[12]:= Integrate[ChebyshevT[4, x]^2/Sqrt[1 - x^2], {x, -1, 1}]
Out[12]= 1/2 Pi
```

A packed Range input stays on the buffer kernel

```mathematica
In[13]:= ChebyshevT[5, Range[0., 1., 0.25]]
Out[13]= {0.0, 0.953125, 0.5, -0.890625, 1.0}
```

The head lowers in Compile at array shape

```mathematica
In[14]:= CompileDiagnostics[{{v, _Real, 1}}, ChebyshevT[7, v]]
Out[14]= {"Compiled" -> True, "ResultType" -> "Array", "Instructions" -> 14, "CommonSubexpressions" -> 0, "InstructionsUnoptimized" -> 15}
```

An exact fractional order with a symbolic argument stays unevaluated

```mathematica
In[15]:= ChebyshevT[1/2, x]
Out[15]= ChebyshevT[1/2, x]
```

## Algorithm

Mathilda -- ChebyshevT, the Chebyshev polynomial of the first kind T_n(x).

```text
  ChebyshevT[n, x]
```

Evaluation mirrors ChebyshevU (chebyshev_u.c): each argument shape takes the cheapest exact or numeric route.

```text
  exact integer order n
    symbolic / exact-number x  ->  the explicit degree-|n| polynomial in x,
                                   with exact integer coefficients, from
                                   FLINT's fmpz_poly_chebyshev_t (or, without
                                   FLINT, the three-term recurrence); an exact
                                   numeric x then evaluates to an exact number.
    inexact x                  ->  the stable three-term recurrence evaluated
                                   in machine/MPFR complex arithmetic, NOT by
                                   substituting into the polynomial (whose
                                   cancellation ruins accuracy for large n --
                                   e.g. ChebyshevT[60, N[2/3,16]]).
  non-integer / complex order, some inexact argument
                               ->  T_v(x) = 2F1(-v, v; 1/2; (1-x)/2), built as
                                   an expression over Hypergeometric2F1 so it
                                   inherits its machine + MPFR, real + complex
                                   numerics; for |(1-x)/2| >= 1 the whole-plane
                                   form T_v(x) = Cos[v ArcCos[x]] takes over.
  special values
    T_0(x) = 1,  T_n(0) = Cos[n Pi/2],  T_n(1) = 1.
  everything else              ->  stays symbolic (return NULL).
```

Unlike U, T is even in its order index -- T_{-n}(x) = T_n(x) -- so a negative order reflects to |n| with no sign and no special zero case (T_{-1}(x) = x).

The recurrence seed is the only arithmetic difference from U: T_1 = x where U_1 = 2x; the step T_k = 2x T_{k-1} - T_{k-2} is identical.

Attributes: Listable, NumericFunction, Protected.

Deferred (left symbolic), matching ChebyshevU: symbolic Series / SeriesCoefficient, asymptotic expansions, half-integer closed forms, exact Interval[{-1,1}] at infinite order, and the Cosh[..]-style symbolic simplification of a fully-exact complex-order call (e.g. ChebyshevT[2-7I, 5+3I] stays symbolic). D[ChebyshevT[n,x],x] = n U_{n-1}(x) lives in calculus/deriv.c.

## Implementation notes

**Algorithm.** `builtin_chebyshev_t` checks arity and calls `chebyshev_t_two_arg`, which picks one route per argument shape. (1) A non-integer order at `x = 0` or `x = 1` returns the special values `Cos[n Pi/2]` (`cheb_value_at_zero`) and `1`. (2) An integer order (exact, or an inexact integer-valued one via `cheb_inexact_int_order`) with an inexact `x` is reflected with `T_{-n} = T_n` (`cheb_reflect`). At machine precision with degree ≤ 1000 it goes to the raw-double recurrence `sf_machine_chebyshev_t` (real) or `sf_machine_chebyshev_t_c` (complex, re/im pairs, no `<complex.h>`). Otherwise, or when the helper declines on overflow, it runs `cheb_numeric_recurrence`, the three-term recurrence `T_k = 2x T_{k-1} - T_{k-2}` in MPFR, real or `ncpx` complex, with 32 guard bits. It never substitutes into the expanded polynomial, which loses all its digits to cancellation at high degree. (3) An exact integer order with an exact or symbolic `x` builds the explicit polynomial from FLINT's `fmpz_poly_chebyshev_t` (the integer recurrence when FLINT is absent) and assembles it with `cheb_poly_from_coeffs`. (4) A non-integer or complex order with some inexact argument (`cheb_numeric_general`) is built as the expression `Hypergeometric2F1[-v, v, 1/2, (1-x)/2]` when `|(1-x)/2| < 1`, and as `Cos[v ArcCos[x]]` otherwise, so it inherits those heads' machine and MPFR numerics and Mathematica's branch cut on `(-∞, -1)`.

**Data structures.** Exact coefficients live in an `mpz_t` array of length `deg+1`, filled from an `fmpz_poly`. Numeric values are MPFR pairs (`ncpx`), or plain `double`s on the machine path. The result precision comes from `cheb_out_prec`, the smallest precision among the inexact inputs. The ND kernel `ndk_ChebyshevT_c` (`src/ndkernels.c`, registered with `REG_B`) wraps the same two `sf_machine_*` helpers. Packed and `NDArray` buffers, real or complex, therefore stay on the buffer, and `Compile[]` lowers the head at both scalar and rank-1 shapes. Scalar, packed, `NDArray` and compiled results agree bit for bit.

**Complexity / limits.** The numeric recurrence costs `O(n)` operations at the working precision. The explicit polynomial is capped at degree 2000 (`CHEB_POLY_CAP`) and the numeric recurrence at |n| ≤ 2·10⁶ (`CHEB_NUMERIC_DEG_CAP`); beyond either cap the call stays unevaluated. An exact non-integer order with an exact or symbolic argument also stays unevaluated, and so does a fully exact complex order. Symbolic `Series` and the half-integer closed forms are not implemented. `D[ChebyshevT[n, x], x] = n ChebyshevU[n-1, x]` lives in `src/calculus/deriv.c`.

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## References

**See also:** [Series](../../power-series/Series/), [SeriesCoefficient](../../power-series/SeriesCoefficient/)

- DLMF §18 — classical orthogonal polynomials (§18.9 recurrences, §18.3 and §18.5.11 for T_n(cos θ) = cos nθ and the 2F1 form).
- Source: [`src/special_functions/chebyshev_t.c`](https://github.com/stblake/mathilda/blob/main/src/special_functions/chebyshev_t.c)
- Specification: [`docs/spec/builtins/special-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/special-functions.md)
- Tests: [`tests/test_chebyshev_t.c`](https://github.com/stblake/mathilda/blob/main/tests/test_chebyshev_t.c)
- Tests: [`tests/test_possiblezeroq_stress.c`](https://github.com/stblake/mathilda/blob/main/tests/test_possiblezeroq_stress.c)

## Notes & additional examples

### Notes

- `ChebyshevT[n, x]` is `T_n(x)`, defined by `T_n(cos θ) = cos(nθ)` and orthogonal on `[-1, 1]` with weight `1/Sqrt[1 - x^2]`. It agrees with Mathematica, including the branch cut on `(-∞, -1)` for non-integer orders.
- An integer order with an exact argument gives the exact polynomial; with an inexact argument it uses the stable recurrence, never the expanded polynomial. A fractional order is evaluated numerically only when some argument is inexact.
- Numeric surfaces: the machine `double` recurrence (real and complex) is shared by the scalar, packed/`NDArray` and `Compile[]` paths, so all of them agree bit for bit. An arbitrary-precision argument, a degree above 1000, or an overflow uses the MPFR recurrence.
- Attributes: `Listable`, `NumericFunction`, `Protected`.
- See also: `ChebyshevU`, `LegendreP`, `Hypergeometric2F1`.
