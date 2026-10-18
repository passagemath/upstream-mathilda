# LaguerreL

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`LaguerreL[n, x]`**

gives the Laguerre polynomial L\_n(x).

**`LaguerreL[n, a, x]`**

gives the generalized Laguerre polynomial L\_n^a(x)

**`Hypergeometric1F1[-n, a+1, x]. Special values include L_n(0) =`**

**`Binomial[n+a, n]. D[LaguerreL[n,x],x] = -LaguerreL[n-1,1,x] and`**

**`D[LaguerreL[n,a,x],x] = -LaguerreL[n-1,a+1,x]. Listable.`**

<details>
<summary>Notes</summary>

(with LaguerreL\[n, x\] == LaguerreL\[n, 0, x\]). A non-negative integer n yields the explicit polynomial (exact for an exact argument, and via a stable recurrence for an inexact one -- the recurrence avoids the cancellation that corrupts the polynomial form for large n); a non-integer or complex order with an inexact argument evaluates numerically at machine or arbitrary (MPFR) precision, real or complex, through the confluent-hypergeometric form Binomial\[n+a,n\]

</details>

## Examples (14)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= LaguerreL[7, x]
Out[1]= 1 - 7 x + 21/2 x^2 - 35/6 x^3 + 35/24 x^4 - 7/40 x^5 + 7/720 x^6 - 1/5040 x^7

In[2]:= Table[LaguerreL[n, a, x], {n, 0, 3}, {a, 0, 1}]
Out[2]= {{1, 1}, {1 - x, 2 - x}, {1 - 2 x + 1/2 x^2, 3 - 3 x + 1/2 x^2}, {1 - 3 x + 3/2 x^2 - 1/6 x^3, 4 - 6 x + 2 x^2 - 1/6 x^3}}

In[3]:= LaguerreL[100, 0, N[5, 20]]
Out[3]= 1.45552716253288126069
```

### Applications (11)

The ordinary Laguerre polynomial with rational coefficients

```mathematica
In[4]:= LaguerreL[4, x]
Out[4]= 1 - 4 x + 3 x^2 - 2/3 x^3 + 1/24 x^4
```

The generalized polynomial with a symbolic parameter

```mathematica
In[5]:= LaguerreL[3, a, x]
Out[5]= 1 + 11/6 a + a^2 + 1/6 a^3 - 3 x - 5/2 a x - 1/2 a^2 x + 3/2 x^2 + 1/2 a x^2 - 1/6 x^3
```

An exact value with a rational parameter

```mathematica
In[6]:= LaguerreL[3, 1/2, 2]
Out[6]= -43/48
```

The value at the origin is a binomial coefficient

```mathematica
In[7]:= LaguerreL[n, a, 0]
Out[7]= Binomial[a + n, n]
```

The Laguerre equation is satisfied identically

```mathematica
In[8]:= Simplify[x D[LaguerreL[4, a, x], {x, 2}] + (a + 1 - x) D[LaguerreL[4, a, x], x] + 4 LaguerreL[4, a, x]]
Out[8]= 0
```

The norm Gamma[n + a + 1]/n! = 20

```mathematica
In[9]:= Integrate[LaguerreL[3, 2, x]^2 x^2 Exp[-x], {x, 0, Infinity}]
Out[9]= 20
```

Distinct degrees are orthogonal

```mathematica
In[10]:= Integrate[LaguerreL[2, x] LaguerreL[3, x] Exp[-x], {x, 0, Infinity}]
Out[10]= 0
```

A degree-30 value at machine and 25-digit precision

```mathematica
In[11]:= {LaguerreL[30, 2, 10.5], N[LaguerreL[30, 2, 21/2], 25]}
Out[11]= {65.1827, 65.182670822327225144593153}
```

A fractional order through Binomial[n + a, n] Hypergeometric1F1[-n, a + 1, x]

```mathematica
In[12]:= LaguerreL[0.5, 1.2]
Out[12]= 0.28796
```

The generalized form over a packed Range input

```mathematica
In[13]:= LaguerreL[3, 0.5, Range[0., 2., 0.25]]
Out[13]= {2.1875, 1.20052, 0.416667, -0.179688, -0.604167, -0.872396, -1.0, -1.0026, -0.895833}
```

The three-argument form does not yet lower at array shape

```mathematica
In[14]:= CompileDiagnostics[{{v, _Real, 1}}, LaguerreL[5, 0.5, v]]
Out[14]= {"Compiled" -> False, "Reason" -> "no machine lowering for this head at these argument types", "Subexpression" -> "LaguerreL[5, 0.5, v]"}
```

## Algorithm

Mathilda -- LaguerreL, the (generalized) Laguerre polynomial.

```text
  LaguerreL[n, x]      L_n(x)       == LaguerreL[n, 0, x]
  LaguerreL[n, a, x]   L_n^(a)(x)   generalized / associated Laguerre
```

The polynomials are orthogonal with weight x^a e^-x and satisfy the ODE

```text
  x y'' + (a + 1 - x) y' + n y = 0.
```

They have the closed forms

```text
  L_n^(a)(x) = Sum_{i=0}^n (-1)^i Binomial[n+a, n-i] x^i / i!
             = Binomial[n+a, n] Hypergeometric1F1[-n, a+1, x]
```

(a = 0 recovers the ordinary L_n).

Evaluation is layered so each argument shape takes the cheapest exact or numeric route (mirroring the ChebyshevU / LegendreP design):

```text
  exact non-negative integer order n
    some inexact argument  ->  the stable three-term recurrence
                                 (k+1) L_{k+1} = (2k+1+a-x) L_k - (k+a) L_{k-1}
                               evaluated in machine / MPFR complex arithmetic,
                               NOT by substituting into the explicit polynomial
                               (whose cancellation ruins accuracy for large n,
                                e.g. LaguerreL[100,0,N[5,20]]).
    symbolic / exact x, numeric a  ->  the explicit degree-n polynomial in x
                               with exact rational coefficients, from the same
                               recurrence on coefficient vectors.
    symbolic parameter a    ->  the explicit polynomial in a and x, built from
                               the sum form and expanded.
  non-integer / complex order, some inexact argument
                           ->  the confluent-hypergeometric closed form, built
                               as an expression over Hypergeometric1F1 so it
                               inherits its machine + MPFR, real + complex
                               numerics. LaguerreL is entire in x, so a single
                               expression covers the whole plane.
  special values
    L_0(x) = 1,  L_n^(a)(0) = Binomial[n+a, n]  (= 1 when a = 0).
  everything else          ->  stays symbolic (return NULL).
```

Attributes: Listable, NumericFunction, Protected.

The machine-precision double / double-complex recurrence lives in sf_machine.c (sf_machine_laguerre_l[_assoc][_c]) and is shared bit-for-bit with the NDArray element-wise kernel and the Compile[] VM (ndkernels.c). The scalar machine path here dispatches to the same helpers.

Deferred (left symbolic, matching ChebyshevU): symbolic Series / SeriesCoefficient closed forms, the Series expansion at infinity, D[] of a symbolic order (the D[LaguerreL[...],x] rules live in calculus/deriv.c), and bespoke tight Interval enclosures (the generic NumericFunction + D[] interval threading in eval.c handles those best-effort). Negative integer order with an inexact argument also stays symbolic.

## Implementation notes

**Algorithm.** `builtin_laguerre_l` accepts `LaguerreL[n, x]` (with `a = 0`) and `LaguerreL[n, a, x]`. An `NDArray` or packed `x` with scalar `n` and `a` goes to `laguerre_l_ndarray`, which threads the machine helper over the buffer and delists to the scalar path if any element declines. Otherwise `laguerre_l_eval` picks a route by shape. (1) An integer order with an inexact argument, at machine precision with `n ≤ 1000`, calls `sf_machine_laguerre_l_assoc` (or its complex twin), the forward recurrence `(k+1) L_{k+1} = (2k+1+a-x) L_k - (k+a) L_{k-1}` in `double`s. Otherwise `lag_numeric_recurrence` runs the same recurrence in MPFR. (2) An exact order with an exact or symbolic `x` and a numeric `a` gives the explicit polynomial with rational coefficients from that recurrence on `mpq_t` coefficient vectors (`lag_poly_numeric_a`). A symbolic `a` gives the expanded sum `Σ (-1)^i Binomial[n+a, n-i] x^i/i!` (`lag_poly_symbolic_a`). (3) A non-integer or complex order with an inexact argument is built as `Binomial[n+a, n] Hypergeometric1F1[-n, a+1, x]` (`lag_numeric_hyp`), which is valid on the whole plane because `L_n` is entire in `x`. (4) The special values `L_0 = 1` and `L_n^{(a)}(0) = Binomial[n+a, n]` hold for a symbolic order.

**Data structures.** `mpq_t` arrays hold the coefficients for numeric `a`, and `Expr` sums hold them for symbolic `a`. The recurrence uses MPFR or `ncpx` values, and plain `double`s on the machine path. The two-argument form has the binary ND kernel `ndk_LaguerreL_c`, real and complex. The three-argument form has the n-ary kernel `ndk_LaguerreL_n` (`REG_N`), which `Compile[]` uses for the scalar shape.

**Complexity / limits.** The recurrence costs `O(n)` operations. The explicit polynomial is capped at `n ≤ 2000` (`LAG_POLY_CAP`), or `n ≤ 200` when `a` is symbolic (`LAG_SYM_POLY_CAP`), and the numeric recurrence at `2·10⁶`. `Compile[]` does not lower the three-argument form at array shape (`CompileDiagnostics` reports `Compiled -> False`), although the interpreter keeps it on the buffer. An exact fractional order with a symbolic argument stays unevaluated, and symbolic series and the expansion at infinity are not implemented. The derivative `-L_{n-1}^{(a+1)}` is in `src/calculus/deriv.c`.

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## References

**See also:** [NDArray](../../linear-algebra/NDArray/), [Series](../../power-series/Series/), [SeriesCoefficient](../../power-series/SeriesCoefficient/)

- DLMF §18 — classical orthogonal polynomials (§18.9.13 recurrence, §18.5.12 explicit sum); §13.6.19 for L_n^(α) = Binomial[n+α, n] M(-n, α+1, x).
- Source: [`src/special_functions/laguerre_l.c`](https://github.com/stblake/mathilda/blob/main/src/special_functions/laguerre_l.c)
- Specification: [`docs/spec/builtins/special-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/special-functions.md)
- Tests: [`tests/test_laguerre_l.c`](https://github.com/stblake/mathilda/blob/main/tests/test_laguerre_l.c)

## Notes & additional examples

### Notes

- `LaguerreL[n, x]` equals `LaguerreL[n, 0, x]`. `LaguerreL[n, a, x]` is the generalized (associated) polynomial `L_n^(a)(x)`, orthogonal on `(0, ∞)` with weight `x^a e^(-x)`.
- A symbolic parameter `a` gives the expanded polynomial in `a` and `x` (up to degree 200). A numeric `a` gives exact rational coefficients (up to degree 2000).
- A fractional or complex order is evaluated as `Binomial[n + a, n] Hypergeometric1F1[-n, a + 1, x]` when some argument is inexact.
- Numeric surfaces: both forms stay on the buffer for packed and `NDArray` input. `Compile[]` lowers `LaguerreL[n, x]` at scalar and array shapes, but lowers `LaguerreL[n, a, x]` only at scalar shape.
- Attributes: `Listable`, `NumericFunction`, `Protected`.
- See also: `HermiteH`, `Hypergeometric1F1`, `Binomial`.
