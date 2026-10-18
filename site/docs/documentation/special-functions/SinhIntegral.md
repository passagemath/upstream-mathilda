# SinhIntegral

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`SinhIntegral[z]`**

gives the hyperbolic sine integral Shi(z) = Integral\_0^z Sinh\[t\]/t dt.

**`SinhIntegral[+-Infinity] = +-Infinity, SinhIntegral[+-I Infinity] = +-I Pi/2.`**

<details>
<summary>Notes</summary>

An entire, odd function with no branch cuts. SinhIntegral\[0\] = 0, Real and complex inputs evaluate numerically at machine or arbitrary (MPFR) precision; D\[SinhIntegral\[z\], z\] = Sinh\[z\]/z. Listable.

</details>

## Examples (12)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= SinhIntegral[2.8]
Out[1]= 4.34808

In[2]:= SinhIntegral[2.5 + I]
Out[2]= 2.84649 + 2.22177*I

In[3]:= SinhIntegral[-2 x]
Out[3]= -SinhIntegral[2 x]
```

### Applications (9)

Exact special values

```mathematica
In[4]:= {SinhIntegral[0], SinhIntegral[Infinity], SinhIntegral[I Infinity]}
Out[4]= {0, Infinity, (1/2*I) Pi}
```

Odd in its argument

```mathematica
In[5]:= SinhIntegral[-2 x]
Out[5]= -SinhIntegral[2 x]
```

Machine real

```mathematica
In[6]:= SinhIntegral[2.8]
Out[6]= 4.34808
```

On the imaginary axis, I SinIntegral[3]

```mathematica
In[7]:= SinhIntegral[3. I]
Out[7]= 0.0 + 1.84865*I
```

Machine complex: an exact rotation of the Si kernel

```mathematica
In[8]:= SinhIntegral[2.5 + I]
Out[8]= 2.84649 + 2.22177*I
```

Past the double range the result is an extended-exponent number, not Infinity

```mathematica
In[9]:= SinhIntegral[10.^6]
Out[9]= 1.516609215011776e+434288
```

Arbitrary precision

```mathematica
In[10]:= N[SinhIntegral[2], 50]
Out[10]= 2.50156743335497564147337248272754239891627287369153
```

All-positive Taylor coefficients

```mathematica
In[11]:= Series[SinhIntegral[x], {x, 0, 7}]
Out[11]= x + 1/18 x^3 + 1/600 x^5 + 1/35280 x^7 + O[x]^8
```

The integral that defines it

```mathematica
In[12]:= Integrate[Sinh[a x]/x, x]
Out[12]= SinhIntegral[a x]
```

## Algorithm

```text
 Mathilda -- the hyperbolic sine integral  Shi(z) = Int_0^z Sinh[t]/t dt.

  SinhIntegral[z]   Shi(z)
```

Shi is entire and odd, with no branch cuts. It is the imaginary-axis sibling of the sine integral: Shi(z) = -i Si(i z). Evaluation is layered so each kind of argument takes the cheapest route:

```text
  exact special values     ->  0, +-Infinity, +-I Pi/2, Indeterminate
  machine real             ->  MPFR series/asymptotic at 53 bits
  arbitrary real           ->  the same, at the input precision
  complex (any precision)  ->  the ncpx series/asymptotic with guard bits
  everything else          ->  stays symbolic (return NULL)
```

The convergent Maclaurin series (the trig series with the alternating sign removed) is

```text
  Shi(z) = Sum_{k>=0} z^(2k+1) / ((2k+1) (2k+1)!),
```

valid for all z. Unlike Si, every term is positive, so on the real axis there is NO catastrophic cancellation -- the partial sums climb monotonically to the O(e^|z|)-sized answer and only a small fixed guard is needed. (Complex inputs still cancel like Si, so the complex paths keep the ~|z|/ln2 guard.) For large

```text
|z| the convergent series is infeasible; there we use the asymptotic expansion

  Shi(z) ~ cosh(z) F(z) + sinh(z) G(z)   (+ Stokes constant, see below),
    F(z) ~ Sum (2k)!   / z^(2k+1) = 1/z + 2!/z^3 + ...,
    G(z) ~ Sum (2k+1)! / z^(2k+2) = 1/z^2 + 3!/z^4 + ...,
```

summed to the smallest term (optimal truncation). F, G are Si's f, g without the (-1)^k. The odd-symmetry reduction Shi(-z) = -Shi(z) folds negative real / left-half-plane inputs onto Re >= 0, keeping the asymptotic within its valid sector. On the real axis the bare form above is exact (real, no constant); off it, matching Mathematica's Series[SinhIntegral[z],{z,Infinity,k}], a Stokes constant restores the analytic value:

```text
  Shi(z) = B(z) + i (Pi/2) sign(Im z),   B(z) = cosh(z) F(z) + sinh(z) G(z).
```

(B is odd; on the imaginary axis B ~ O(1/z) vanishes and the i Pi/2 sign(Im z) is all that survives, giving Shi(+- i Infinity) = +- i Pi/2.)

Machine-real results that overflow a C double (e.g. SinhIntegral[10.^6] ~ 1.5*10^434288) are emitted as a 53-bit MPFR real, whose exponent range easily holds them, rather than as an inf-valued double.

Attributes: Listable, NumericFunction, Protected.

## Implementation notes

**Algorithm.** `builtin_sinhintegral` evaluates `Shi(z) = ∫₀ᶻ sinh t/t dt`, which is entire and odd. It is the imaginary-axis sibling of `Si`: `Shi(z) = -i Si(iz)`. A machine call first passes the `sf_machine_scalar_args` gate and goes to the kernels in `sf_machine_expint.c`:

- **Real** (`sf_machine_shi`, the Cephes `shichi` structure). For `|x| <= 4` it uses `x·Sh(x²)` with positive coefficients. Beyond that it uses `(Ei(x) + E1(x))/2`, taking `Ei` and `E1` from the Cody–Thacher kernels. It declines `x > 720`, where the value overflows.
- **Complex** (`sf_machine_shi_c`). It uses `-i Si(iz)`, an exact rotation of the `SinIntegral` kernel.

A declined call, or an MPFR call, goes to `shi_one_arg`:

- **Exact values:** `0`, `±Infinity`, `±I Infinity → ±I π/2`, and `ComplexInfinity` or `Indeterminate → Indeterminate`.
- **Odd fold.** A negative-leading argument folds out by oddness.
- **MPFR, Maclaurin series.** On the real axis every term is positive, so only a small fixed guard is needed; off the axis it keeps `~|z|/ln2` guard bits.
- **MPFR, large `|z|`.** The asymptotic `cosh(z) F(z) + sinh(z) G(z)`, plus the Stokes constant `i(π/2) sign(Im z)` off the real axis.

When a machine-real result overflows a double (`SinhIntegral[10.^6]`), it comes back as a 53-bit MPFR real with an extended exponent instead of `Infinity`.

**Data structures.** `Expr`; `mpfr_t` and `ncpx` on the MPFR side; `(re, im)` doubles and the refitted tables in `sf_machine_expint_coef.h` on the machine side. The ND kernel is `NDKU_SinhIntegral = { sinhintegral_machine_complex, ndk_SinhIntegral_r }`, registered with `REG_U`, so it is packed-aware. `Compile[]` lowers `SinhIntegral` through that kernel for real and complex scalars and for rank-1 arrays. Every surface uses the same helpers, so the results agree bit for bit.

**Complexity / limits.** The machine kernels are `O(1)`. Against mpmath the real kernel has a median error of 0.4 ulp and a maximum of 2.3 ulp; the complex kernel's maximum is 2.4 ulp. The MPFR series costs grow with `|z|` and with the precision. Exact non-special arguments (`SinhIntegral[2]`) stay unevaluated. `D[SinhIntegral[z], z] = Sinh[z]/z`. Attributes: `Listable`, `NumericFunction`, `Protected`.

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## References

**See also:** [SinIntegral](../../special-functions/SinIntegral/)

- DLMF §6.2 — the hyperbolic sine integral Shi.
- S. L. Moshier, *Methods and Programs for Mathematical Functions* (Ellis Horwood, 1989) — the Cephes `shichi` structure of the real kernel.
- W. J. Cody and H. C. Thacher, *Chebyshev approximations for the exponential integral Ei(x)*, Math. Comp. 23 (1969) 289–303.
- Source: [`src/special_functions/sinhintegral.c`](https://github.com/stblake/mathilda/blob/main/src/special_functions/sinhintegral.c)
- Specification: [`docs/spec/builtins/special-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/special-functions.md)
- Tests: [`tests/test_compile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile.c)
- Tests: [`tests/test_expint_machine.c`](https://github.com/stblake/mathilda/blob/main/tests/test_expint_machine.c)
- Tests: [`tests/test_numeric_stress.c`](https://github.com/stblake/mathilda/blob/main/tests/test_numeric_stress.c)
- Tests: [`tests/test_sinhintegral.c`](https://github.com/stblake/mathilda/blob/main/tests/test_sinhintegral.c)

## Notes & additional examples

### Notes

- `SinhIntegral[z] = ∫₀ᶻ Sinh[t]/t dt` is entire and odd, the imaginary-axis sibling of `SinIntegral` (`Shi(z) = -I Si(I z)`); `SinhIntegral[±I Infinity] = ±I Pi/2`.
- Machine arguments use the `sf_machine_expint.c` kernels (Cephes `shichi` on the real line, `(Ei + E1)/2` beyond 4), shared bit for bit with packed lists, `NDArray` and `Compile[]`; overflow past `x ≈ 720` and arbitrary precision go to MPFR.
- Attributes: `Listable`, `NumericFunction`, `Protected`. See also `CoshIntegral`, `SinIntegral`, `ExpIntegralEi`.
