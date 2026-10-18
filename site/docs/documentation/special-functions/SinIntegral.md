# SinIntegral

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`SinIntegral[z]`**

gives the sine integral Si(z) = Integral\_0^z Sin\[t\]/t dt.

**`SinIntegral[+-Infinity] = +-Pi/2, SinIntegral[+-I Infinity] = +-I Infinity.`**

<details>
<summary>Notes</summary>

An entire, odd function with no branch cuts. SinIntegral\[0\] = 0, Real and complex inputs evaluate numerically at machine or arbitrary (MPFR) precision; D\[SinIntegral\[z\], z\] = Sinc\[z\]. Listable.

</details>

## Examples (13)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= SinIntegral[2.8]
Out[1]= 1.8321

In[2]:= SinIntegral[2.5 + I]
Out[2]= 1.99549 + 0.222995*I

In[3]:= SinIntegral[-x]
Out[3]= -SinIntegral[x]
```

### Applications (10)

Limits along the axes

```mathematica
In[4]:= SinIntegral[{-Infinity, Infinity, -I Infinity, I Infinity}]
Out[4]= {-1/2 Pi, 1/2 Pi, -I Infinity, I Infinity}
```

Machine real

```mathematica
In[5]:= SinIntegral[2.8]
Out[5]= 1.8321
```

Approaches Pi/2 for large x

```mathematica
In[6]:= SinIntegral[1.*^6]
Out[6]= 1.5708
```

On the imaginary axis, I SinhIntegral[3]

```mathematica
In[7]:= SinIntegral[3. I]
Out[7]= 0.0 + 4.97344*I
```

Machine complex

```mathematica
In[8]:= SinIntegral[2.5 + I]
Out[8]= 1.99549 + 0.222995*I
```

Arbitrary precision

```mathematica
In[9]:= N[SinIntegral[2], 50]
Out[9]= 1.6054129768026948485767201481985889408485834223285
```

The odd-symmetry fold

```mathematica
In[10]:= SinIntegral[-x]
Out[10]= -SinIntegral[x]
```

The Taylor series at the origin

```mathematica
In[11]:= Series[SinIntegral[x], {x, 0, 7}]
Out[11]= x - 1/18 x^3 + 1/600 x^5 - 1/35280 x^7 + O[x]^8
```

The asymptotic expansion

```mathematica
In[12]:= Normal[Series[SinIntegral[x], {x, Infinity, 3}]]
Out[12]= 1/2 Pi - Sin[x]/x^2 + Cos[x] (-1/x + 2/x^3)
```

The derivative is Sinc

```mathematica
In[13]:= D[SinIntegral[x], x]
Out[13]= Sinc[x]
```

## Algorithm

```text
 Mathilda -- the sine integral  Si(z) = Int_0^z Sin[t]/t dt.

  SinIntegral[z]   Si(z)
```

Si is entire and odd, with no branch cuts. Evaluation is layered so each kind of argument takes the cheapest route:

```text
  exact special values     ->  0, +-Pi/2, +-I Infinity, Indeterminate
  machine real / complex   ->  the double kernels in sf_machine_expint.c
                               (shared with NDArray and Compile[])
  arbitrary real           ->  MPFR series/asymptotic at the input precision
  arbitrary complex        ->  the ncpx series/asymptotic with guard bits
  everything else          ->  stays symbolic (return NULL)
```

The convergent Maclaurin series (DLMF 6.6.5) is

```text
  Si(z) = Sum_{k>=0} (-1)^k z^(2k+1) / ((2k+1) (2k+1)!),
```

valid for all z. Its partial sums can reach magnitude ~e^|z| before the O(1)-sized answer emerges, so the MPFR paths add ~|z|/ln2 guard bits to absorb that cancellation exactly. For large |z| the convergent series is infeasible; there we use the asymptotic expansion (DLMF 6.12.3)

```text
  Si(z) = Pi/2 - cos(z) f(z) - sin(z) g(z),
    f(z) ~ Sum (-1)^k (2k)!   / z^(2k+1) = 1/z - 2!/z^3 + ...,
    g(z) ~ Sum (-1)^k (2k+1)! / z^(2k+2) = 1/z^2 - 3!/z^4 + ...,
```

summed to the smallest term (optimal truncation). Both odd-symmetry reductions Si(-z) = -Si(z) fold negative real / left-half-plane inputs onto the right half-plane, keeping the asymptotic within its valid sector.

Attributes: Listable, NumericFunction, Protected.

## Implementation notes

**Algorithm.** `builtin_sinintegral` evaluates `Si(z) = ∫₀ᶻ sin t/t dt`, an entire and odd function. A machine call first goes through the `sf_machine_scalar_args` gate and the kernels in `sf_machine_expint.c`:

- **Real** (`sf_machine_si`, Cephes `sici` structure, coefficients refitted from mpmath). For `|x| <= 4` it uses `x·S(x²)`, an economised Maclaurin series. Beyond that it uses `π/2 - f cos x - g sin x`, with the auxiliary functions written as `x f = 1 - w·RF(w)` and `x² g = 1 - w·RG(w)`, where `w = 1/x²`, on `[4, 8]` and `[8, ∞)`.
- **Complex** (`sf_machine_si_c`). It uses the Maclaurin series by Horner where its terms do not cancel (`|z| - |Im z| <= 2` or `|z| <= 2.5`). Otherwise it uses `π/2 + (E1(iz) - E1(-iz))/(2i)`, with `E1` from the continued fraction evaluated backward. For `|z| > 45` it uses the `f`/`g` asymptotic series.

If the kernel declines, or the call is MPFR, `si_one_arg` handles it:

- **Exact values:** `0`, `±Infinity → ±π/2`, `±I Infinity → ±I Infinity`, and `ComplexInfinity` or `Indeterminate → Indeterminate`.
- **Odd fold:** a negative-leading `Times` folds by oddness.
- **MPFR, moderate `|z|`:** the convergent Maclaurin series (DLMF 6.6.5) with about `|z|/ln2` guard bits, since the partial sums reach about `e^|z|`.
- **MPFR, large `|z|`:** the asymptotic expansion `π/2 - cos(z) f(z) - sin(z) g(z)`, summed to its smallest term.

Complex values use the `ncpx` toolkit.

**Data structures.** `Expr`; `mpfr_t` and `ncpx` on the MPFR side; refitted tables in `sf_machine_expint_coef.h` on the machine side. The ND kernel is `NDKU_SinIntegral = { sinintegral_machine_complex, ndk_SinIntegral_r }`, registered with `REG_U`, so it is packed-aware. `Compile[]` lowers `SinIntegral` through that kernel for real and complex scalars and for rank-1 arrays. All five surfaces (scalar, `List`, packed, `NDArray`, `Compile[]`) share the helpers and agree bit for bit.

**Complexity / limits.** The machine kernels are `O(1)`. Against mpmath, the real kernel has median 0.3 ulp and max 1.8 ulp error, and the complex kernel has max 1.5 ulp scaled by `|z f′|` near its zeros. The MPFR convergent regime costs more as `|z|` and the precision grow; the asymptotic regime uses optimal truncation. Exact non-special arguments (`SinIntegral[2]`) stay unevaluated. Attributes: `Listable`, `NumericFunction`, `Protected`.

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## References

- DLMF §6.2 — the sine integral Si.
- DLMF §6.6.5 — convergent Maclaurin series; §6.9.1 — the E1 continued fraction; §6.12.3 — asymptotic expansion.
- S. L. Moshier, *Methods and Programs for Mathematical Functions* (Ellis Horwood, 1989) — the Cephes `sici` structure of the real kernel.
- Source: [`src/special_functions/sinintegral.c`](https://github.com/stblake/mathilda/blob/main/src/special_functions/sinintegral.c)
- Specification: [`docs/spec/builtins/special-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/special-functions.md)
- Tests: [`tests/test_compile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile.c)
- Tests: [`tests/test_expint_machine.c`](https://github.com/stblake/mathilda/blob/main/tests/test_expint_machine.c)
- Tests: [`tests/test_integrate_divergent.c`](https://github.com/stblake/mathilda/blob/main/tests/test_integrate_divergent.c)
- Tests: [`tests/test_integrate_ramanujan.c`](https://github.com/stblake/mathilda/blob/main/tests/test_integrate_ramanujan.c)

## Notes & additional examples

### Notes

- `SinIntegral[z] = ∫₀ᶻ Sin[t]/t dt` is entire and odd, with `SinIntegral[±Infinity] = ±Pi/2`.
- Machine arguments use the `sf_machine_expint.c` kernels (Cephes `sici` structure on the real line; series, backward E1 continued fraction, or asymptotic series in the complex plane), shared bit for bit with packed lists, `NDArray` and `Compile[]`; arbitrary precision uses the guarded Maclaurin series or the asymptotic expansion.
- Attributes: `Listable`, `NumericFunction`, `Protected`. See also `CosIntegral`, `SinhIntegral`, `Sinc`, `ExpIntegralEi`.
