# CoshIntegral

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`CoshIntegral[z]`**

gives the hyperbolic cosine integral Chi(z) = EulerGamma + Log\[z\] + Integral\_0^z (Cosh\[t\] - 1)/t dt.

**`CoshIntegral[0] = -Infinity, CoshIntegral[Infinity] = Infinity,`**

**`CoshIntegral[+-I Infinity] = +-I Pi/2.`**

<details>
<summary>Notes</summary>

Has a logarithmic singularity at 0 and a branch cut on (-Infinity, 0\]. Real and complex inputs evaluate numerically at machine or arbitrary (MPFR) precision; D\[CoshIntegral\[z\], z\] = Cosh\[z\]/z. Listable.

</details>

## Examples (13)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= CoshIntegral[2.8]
Out[1]= 4.33122

In[2]:= CoshIntegral[-2.]
Out[2]= 2.45267 + 3.14159*I

In[3]:= CoshIntegral[3.5 + I]
Out[3]= 5.36668 + 4.28423*I
```

### Applications (10)

The logarithmic singularity and the limits

```mathematica
In[4]:= {CoshIntegral[0], CoshIntegral[Infinity], CoshIntegral[I Infinity]}
Out[4]= {-Infinity, Infinity, (1/2*I) Pi}
```

Machine real

```mathematica
In[5]:= CoshIntegral[2.8]
Out[5]= 4.33122
```

Beside the zero x0 = 0.52382..., still relatively accurate

```mathematica
In[6]:= CoshIntegral[0.5238]
Out[6]= -4.91387e-05
```

A negative real is on the cut: the value from above, Chi(2) + I Pi

```mathematica
In[7]:= CoshIntegral[-2.]
Out[7]= 2.45267 + 3.14159*I
```

On the imaginary axis, CosIntegral[3] + I Pi/2

```mathematica
In[8]:= CoshIntegral[3. I]
Out[8]= 0.11963 + 1.5708*I
```

Machine complex

```mathematica
In[9]:= CoshIntegral[3.5 + I]
Out[9]= 5.36668 + 4.28423*I
```

Arbitrary precision

```mathematica
In[10]:= N[CoshIntegral[2], 50]
Out[10]= 2.45266692264691452190613264749949287660178068872851
```

Threads element-wise over a packed real list

```mathematica
In[11]:= CoshIntegral[{1.0, 2.0, 3.0}]
Out[11]= {0.837867, 2.45267, 4.96039}
```

The logarithmic Taylor series

```mathematica
In[12]:= Series[CoshIntegral[x], {x, 0, 4}]
Out[12]= EulerGamma + Log[x] + 1/4 x^2 + 1/96 x^4 + O[x]^5
```

The derivative is Cosh[x]/x

```mathematica
In[13]:= D[CoshIntegral[x], x]
Out[13]= Cosh[x]/x
```

## Algorithm

Mathilda -- the hyperbolic cosine integral

```text
  Chi(z) = EulerGamma + Log[z] + Int_0^z (Cosh[t] - 1)/t dt.

  CoshIntegral[z]   Chi(z)
```

Chi has a logarithmic singularity at z = 0 and a branch cut running along the

```text
negative real axis (-Infinity, 0].  It is the imaginary-axis sibling of the
cosine integral: Chi(z) = Ci(i z) - i Pi/2.  Evaluation is layered so each
```

argument takes the cheapest route:

```text
  exact special values     ->  -Infinity, Infinity, +-I Pi/2, Indeterminate
  machine real             ->  MPFR series/asymptotic at 53 bits
  arbitrary real           ->  the same, at the input precision
  complex (any precision)  ->  the ncpx series/asymptotic with guard bits
  everything else          ->  stays symbolic (return NULL)
```

The convergent Maclaurin series (the trig series with the alternating sign removed) is

```text
  Chi(z) = EulerGamma + Log(z) + Sum_{k>=1} z^(2k) / (2k (2k)!),
```

valid on the principal branch for all z != 0 -- the principal Log(z) supplies the correct +-i Pi jump across the cut, so no folding is needed for the

```text
convergent path.  Every series term is positive, so on the real axis there is
NO catastrophic cancellation and only a small fixed guard is needed.  For large
|z| the convergent series is infeasible; there we use the asymptotic expansion

  Chi(z) ~ sinh(z) F(z) + cosh(z) G(z)   (+ Stokes constant, see below),
    F(z) ~ Sum (2k)!   / z^(2k+1) = 1/z + 2!/z^3 + ...,
    G(z) ~ Sum (2k+1)! / z^(2k+2) = 1/z^2 + 3!/z^4 + ...,
```

(the same F, G SinhIntegral computes; only the combination differs -- Shi uses

```text
cosh F + sinh G).  The bare part sinh F + cosh G is even and real on the
positive real axis (no constant).  Off the real axis, matching Mathematica's
```

Series[CoshIntegral[z],{z,Infinity,k}], a piecewise Stokes term restores the principal branch:

```text
  Chi(z) = B2(z) + i K,   B2(z) = sinh(z) F(z) + cosh(z) G(z),
    K = -Pi/2                       for Im(z) < 0,
    K = Pi sgn+(Re z) - Pi/2        for Im(z) > 0   (sgn+(0) = +1, from above).
```

(On the imaginary axis Im(z) > 0, Re(z) = 0 gives K = +Pi/2, so

```text
Chi(+- i Infinity) = +- i Pi/2.)  A negative real x is handled by the real
```

path as the from-above branch value Complex[Chi(|x|), Pi], matching the log jump; and Chi(-Infinity) -> Infinity (the real part dominates).

Machine-real results that overflow a C double (e.g. CoshIntegral[10.^6] ~ 1.5*10^434288) are emitted as a 53-bit MPFR real, whose exponent range easily holds them, rather than as an inf-valued double.

Attributes: Listable, NumericFunction, Protected.

## Implementation notes

**Algorithm.** `builtin_coshintegral` evaluates `Chi(z)`, the imaginary-axis sibling of `Ci` (`Chi(z) = Ci(iz) - iπ/2`). It has a logarithmic singularity at 0 and a branch cut along `(-∞, 0]`.

A machine call first passes the `sf_machine_scalar_args` gate and goes to the kernels in `sf_machine_expint.c`:

- **Real** (`sf_machine_chi`, Cephes `shichi` structure). For `0 < x <= 4` it uses `γ + log x + x²·Ch(x²)`; for larger `x`, `(Ei(x) - E1(x))/2`. It declines `x > 720` (overflow). For `x < 0` the scalar path asks `sf_machine_chi_c` for the complex from-above value `Chi(|x|) + iπ`.
- **Complex** (`sf_machine_chi_c`). It uses `Ci(iz) + (log z - log iz)`. The constant `-iπ/2`, or `3iπ/2` past `arg z = π/2`, is inserted exactly rather than obtained by subtraction.

A declined call, or an MPFR call, goes to `chi_one_arg`:

- **Exact values:** `Chi[0] = -Infinity`, `±Infinity → Infinity`, `±I Infinity → ±I π/2`, and `ComplexInfinity` or `Indeterminate → Indeterminate`.
- **MPFR series, moderate `|z|`.** The convergent series `γ + Log z + Σ z^{2k}/(2k (2k)!)`. All its terms are positive on the real axis, so a fixed 64-bit guard is enough. The principal `Log` supplies the jump across the cut.
- **MPFR asymptotic, large `|z|`.** `sinh(z) F(z) + cosh(z) G(z)` plus a piecewise Stokes constant.

A machine result that overflows a double comes back as a 53-bit MPFR real.

**Data structures.** `Expr`; `ncpx` (`numeric_complex.h`) on the MPFR side; refitted tables in `sf_machine_expint_coef.h`. The ND kernel is `NDKU_CoshIntegral = { coshintegral_machine_complex, ndk_CoshIntegral_r }`, registered with `REG_U`. Its real arm declines `x <= 0` and hands that buffer to the `List` path, which returns the complex values. `Compile[]` lowers `CoshIntegral` through the kernel for real and complex scalars and for rank-1 arrays. All surfaces call the same helpers and agree bit for bit.

**Complexity / limits.** The machine kernels are `O(1)`. Against mpmath the real kernel has median error 0.45 ulp and maximum 1.6 ulp, scaled at the zero `0.5238…`; the complex kernel's maximum is 2.6 ulp. The cost of the MPFR path grows with `|z|` and with the precision. Exact non-special arguments, and `CoshIntegral[-x]`, stay unevaluated. Attributes: `Listable`, `NumericFunction`, `Protected`.

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## References

**See also:** [CosIntegral](../../special-functions/CosIntegral/), [Log](../../elementary-functions/Log/)

- DLMF §6.2.16 — the hyperbolic cosine integral Chi(z) = γ + Log(z) + Int_0^z (cosh t - 1)/t dt.
- S. L. Moshier, *Methods and Programs for Mathematical Functions* (Ellis Horwood, 1989) — the Cephes `shichi` structure of the real kernel.
- W. J. Cody and H. C. Thacher, *Chebyshev approximations for the exponential integral Ei(x)*, Math. Comp. 23 (1969) 289–303.
- Source: [`src/special_functions/coshintegral.c`](https://github.com/stblake/mathilda/blob/main/src/special_functions/coshintegral.c)
- Specification: [`docs/spec/builtins/special-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/special-functions.md)
- Tests: [`tests/test_compile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile.c)
- Tests: [`tests/test_coshintegral.c`](https://github.com/stblake/mathilda/blob/main/tests/test_coshintegral.c)
- Tests: [`tests/test_expint_machine.c`](https://github.com/stblake/mathilda/blob/main/tests/test_expint_machine.c)
- Tests: [`tests/test_numeric_stress.c`](https://github.com/stblake/mathilda/blob/main/tests/test_numeric_stress.c)

## Notes & additional examples

### Notes

- `CoshIntegral[z] = EulerGamma + Log[z] + ∫₀ᶻ (Cosh[t] - 1)/t dt` has a logarithmic singularity at 0 and a branch cut on `(-∞, 0]`; it is the imaginary-axis sibling of `CosIntegral` (`Chi(z) = Ci(I z) - I Pi/2`).
- Machine arguments use the `sf_machine_expint.c` kernels (Cephes `shichi` on the real line, `(Ei - E1)/2` beyond 4; `Ci(iz)` with an exactly placed constant in the complex plane), shared bit for bit with packed lists, `NDArray` and `Compile[]`; overflow and arbitrary precision go to MPFR.
- Attributes: `Listable`, `NumericFunction`, `Protected`. See also `SinhIntegral`, `CosIntegral`, `ExpIntegralEi`.
