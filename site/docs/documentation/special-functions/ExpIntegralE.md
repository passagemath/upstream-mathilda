# ExpIntegralE

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`ExpIntegralE[n, z]`**

gives the generalized exponential integral E\_n(z) = Integral\_1^Infinity E^(-z t)/t^n dt.

**`ExpIntegralE[0, z] = E^-z/z, ExpIntegralE[n, 0] = 1/(n-1) for n > 1,`**

**`ExpIntegralE[1, z] = Gamma[0, z]. Integer order with an inexact real z > 0`**

<details>
<summary>Notes</summary>

evaluates numerically; D\[ExpIntegralE\[n, z\], z\] = -ExpIntegralE\[n-1, z\]. Listable.

</details>

## Examples (14)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= ExpIntegralE[0, z]
Out[1]= E^(-z)/z

In[2]:= {ExpIntegralE[2, 3/2], N[ExpIntegralE[2, 3/2]]}
Out[2]= {ExpIntegralE[2, 3/2], 0.0731008}

In[3]:= ExpIntegralE[3, -1.]
Out[3]= 1.77072 - 1.5708*I

In[4]:= D[ExpIntegralE[n, z], z]
Out[4]= -ExpIntegralE[-1 + n, z]
```

### Applications (10)

Order 0 is elementary

```mathematica
In[5]:= ExpIntegralE[0, z]
Out[5]= E^(-z)/z
```

The value at the origin is 1/(n - 1) for n >= 2

```mathematica
In[6]:= ExpIntegralE[3, 0]
Out[6]= 1/2
```

An exact argument stays symbolic until N

```mathematica
In[7]:= {ExpIntegralE[2, 3/2], N[ExpIntegralE[2, 3/2]]}
Out[7]= {ExpIntegralE[2, 3/2], 0.0731008}
```

The order-1 case is the incomplete gamma Gamma[0, z]

```mathematica
In[8]:= ExpIntegralE[1, 2.] == Gamma[0, 2.]
Out[8]= True
```

A negative real argument lands on the upper side of the cut

```mathematica
In[9]:= ExpIntegralE[3, -1.]
Out[9]= 1.77072 - 1.5708*I
```

A machine complex argument

```mathematica
In[10]:= ExpIntegralE[1, 2. + I]
Out[10]= 0.00938816 - 0.044463*I
```

The derivative lowers the order

```mathematica
In[11]:= D[ExpIntegralE[n, z], z]
Out[11]= -ExpIntegralE[-1 + n, z]
```

The Mellin transform Gamma[s]/(s + n - 1)

```mathematica
In[12]:= Integrate[x^(s - 1) ExpIntegralE[2, x], {x, 0, Infinity}, Assumptions -> s > 0]
Out[12]= Gamma[s]/(1 + s)
```

The compiled array result matches the packed one bit for bit

```mathematica
In[13]:= ce = Compile[{{v, _Real, 1}}, ExpIntegralE[3, v]]; ce[Range[0.5, 5., 0.5]] === ExpIntegralE[3, Range[0.5, 5., 0.5]]
Out[13]= True
```

A non-integer order stays unevaluated; E_n(z) = z^(n-1) Gamma[1 - n, z] gives the value

```mathematica
In[14]:= {ExpIntegralE[5/2, 1.2], 1.2^(3/2) Gamma[-3/2, 1.2]}
Out[14]= {ExpIntegralE[5/2, 1.2], 0.0958257}
```

## Algorithm

Mathilda -- the generalized exponential integral E_n.

```text
  ExpIntegralE[n, z]   E_n(z) = Integrate[E^(-z t)/t^n, {t, 1, Infinity}]
```

Layered evaluation (each argument kind takes the cheapest correct route):

```text
  E_0(z) = e^(-z)/z                      exact, any z
  E_n(0) = 1/(n-1),  integer n >= 2      exact
  integer n >= 0, machine real x > 0     sf_machine_expn (sf_machine_expint.c):
                                         recurrence/series (x <= 1) / continued
                                         fraction (x > 1)
  integer n >= 0, machine complex z      sf_machine_expn_c (and real x < 0,
  (or real x < 0)                        whose value is complex)
  everything else                        stays symbolic (return NULL)
```

The derivative rule d/dz E_n(z) = -E_(n-1)(z) (with the base E_0 above) lives in deriv.c; the half-line Mellin transform of E_n is produced by the operational-calculus fallback in integrate_ramanujan.c via that chain.

Attributes: Listable, NumericFunction, Protected.

## Implementation notes

**Algorithm.** `builtin_expintegrale` works in layers. (1) The **machine-scalar gate** (`sf_machine_scalar_args`) sends a machine-precision call with a real integer order to the `double` kernels in `sf_machine_expint.c`. For real `x > 0`, `sf_machine_expn` handles each case separately. `n = 0` is `e^{-x}/x`, and `n = 1` is the Cody–Thacher `E1` kernel. For `x ≤ 1` and `n ≤ 24` it runs the upward recurrence `E_{k+1} = (e^{-x} - x E_k)/k` from `E1`, which measured about 2 ulp. For `x ≤ 1` and larger `n` it sums the power series (DLMF 8.19.8) with the `ψ(n) - log x` term. For `x > 1` it evaluates the continued fraction 8.19.17 backward from the fixed depth `105/x + 12`, which avoids the error that the Lentz and Wallis forms accumulate. A complex `z`, or a real `x < 0` (where the value lies on the upper side of the cut), goes to `sf_machine_expn_c`. That routine uses the complex continued fraction when `|z| + Re z > 1`, with depth `240/s + 8` capped at 400, and the series elsewhere. (2) Exact rules: `E_0(z) = E^(-z)/z`, and `E_n(0) = 1/(n-1)` for integer `n ≥ 2`. (3) Everything else stays symbolic.

**Data structures.** The kernels use plain `double`s, with complex values spelled out as (re, im) pairs (`ExC`) and no `<complex.h>`. The binary ND kernel `ndk_ExpIntegralE_c` (`src/ndkernels.c`, `REG_B`) calls the same two functions. Packed and `NDArray` buffers therefore stay on the buffer, and `Compile[]` lowers the head at scalar and array shapes. Scalar, packed, `NDArray` and compiled results agree bit for bit. Measured accuracy against mpmath has a median of 0.4 ulp and a maximum of 4.3 ulp on the real line.

**Complexity / limits.** Each evaluation takes `O(1)` to `O(100)` floating-point steps. There is **no arbitrary-precision path**, so `N[ExpIntegralE[2, 3/2], 30]` stays unevaluated. A non-integer order, or an order above 10⁶, is left unevaluated even at machine precision, and so is `ExpIntegralE[1, 0]`. A machine-real underflow returns `0.`. The derivative `-E_{n-1}` lives in `src/calculus/deriv.c`, and the Mellin transform is handled by `Integrate`'s operational-calculus layer.

- Exact values: `ExpIntegralE[0, z] = E^(−z)/z` (the derivative base case),
  `ExpIntegralE[n, 0] = 1/(n−1)` for integer n ≥ 2.
- Numeric: an integer order n ≥ 1 with an inexact real z > 0 evaluates by the
  standard continued-fraction (z > 1) / power-series (0 < z ≤ 1) split
  (`ExpIntegralE[1, 2.] = 0.0489005`, `ExpIntegralE[2, 1.5] = 0.0731008`).
  Exact arguments stay symbolic; `N` forces evaluation.
- **Machine precision** (integer order n ≥ 0): double kernels shared bit-for-bit
  with the packed/NDArray kernels and `Compile[]`
  (`src/special_functions/sf_machine_expint.c`; the scalar builtin calls the same
  function, so scalar == List == packed == NDArray == Compile). A decline (a pole,
  a non-integer order) falls through unchanged. Real x > 0: the upward recurrence E_(k+1) =
  (e^−x − x E_k)/k from E1 for x ≤ 1 (the power series for n > 24), and for x > 1
  the continued fraction (DLMF 8.19.17) evaluated backward from depth 105/x + 12
  (forward Lentz and Wallis forms lose 1e-14 to accumulated rounding). Machine
  **complex** z — previously unevaluated — and a real x < 0 (complex-valued, the
  upper side of the cut, e.g. `ExpIntegralE[2, -1.] = 0.823164 - 3.14159 I`) use the
  complex continued fraction for |z| + Re z > 1 and the series elsewhere. A
  machine-real underflow returns 0., as before (there is no arbitrary-precision path).
  `ExpIntegralE[n, packed]` now has a binary NDArray kernel (it threaded element by
  element before). Accuracy vs mpmath: median 0.4 ulp, max 4.3 ulp (real); max 3.8
  ulp (complex).
- Derivative: `D[ExpIntegralE[n, z], z] = -ExpIntegralE[n-1, z]`.
- Mellin transform: `∫₀^∞ x^(s-1) ExpIntegralE[n, a x] dx = a^(-s) Γ(s)/(s+n-1)`
  (`Re s > 0`), closed by the operational-calculus layer of `Method -> "Mellin"`.

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## References

**See also:** [N](../../arithmetic/N/)

- DLMF §8.19 — the generalized exponential integral (8.19.8 power series, 8.19.12 recurrence, 8.19.17 continued fraction).
- W. H. Press et al., *Numerical Recipes*, §6.3 (expint) — the series / continued-fraction split for E_n.
- S. L. Moshier, *Methods and Programs for Mathematical Functions* (1989), the Cephes `expn` routine.
- Source: [`src/special_functions/expintegrale.c`](https://github.com/stblake/mathilda/blob/main/src/special_functions/expintegrale.c)
- Specification: [`docs/spec/builtins/special-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/special-functions.md)
- Tests: [`tests/test_expint_machine.c`](https://github.com/stblake/mathilda/blob/main/tests/test_expint_machine.c)
- Tests: [`tests/test_integrate_ramanujan.c`](https://github.com/stblake/mathilda/blob/main/tests/test_integrate_ramanujan.c)
- Tests: [`tests/test_series_assumptions.c`](https://github.com/stblake/mathilda/blob/main/tests/test_series_assumptions.c)

## Notes & additional examples

### Notes

- `ExpIntegralE[n, z] = ∫₁^∞ e^(-z t)/t^n dt`, so `ExpIntegralE[1, z] = Gamma[0, z]`. In general `E_n(z) = z^(n-1) Gamma[1-n, z]`.
- Only an integer order `n ≥ 0` is evaluated numerically, and only at machine precision. Mathematica also handles non-integer orders and arbitrary precision; here `ExpIntegralE[5/2, 1.2]` and `N[ExpIntegralE[2, 3/2], 30]` stay unevaluated, and `Gamma[1 - n, z]` gives those values.
- A machine-real underflow returns `0.`; for example `ExpIntegralE[2, 800.]` gives `0.0`, because the true value is near 10⁻³⁵¹.
- Numeric surfaces: one `double` kernel (real and complex) is shared by the scalar, packed/`NDArray` and `Compile[]` paths.
- Attributes: `Listable`, `NumericFunction`, `Protected`.
- See also: `ExpIntegralEi`, `Gamma`, `LogIntegral`.
