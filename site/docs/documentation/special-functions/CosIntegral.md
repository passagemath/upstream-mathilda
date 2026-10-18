# CosIntegral

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`CosIntegral[z]`**

gives the cosine integral Ci(z) = -Integral\_z^Infinity Cos\[t\]/t dt.

**`CosIntegral[0] = -Infinity, CosIntegral[Infinity] = 0,`**

**`CosIntegral[-Infinity] = I Pi, CosIntegral[+-I Infinity] = Infinity.`**

<details>
<summary>Notes</summary>

Has a logarithmic singularity at 0 and a branch cut on (-Infinity, 0\]. Real and complex inputs evaluate numerically at machine or arbitrary (MPFR) precision; D\[CosIntegral\[z\], z\] = Cos\[z\]/z. Listable.

</details>

## Examples (13)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= CosIntegral[2.8]
Out[1]= 0.186488

In[2]:= CosIntegral[-2.]
Out[2]= 0.422981 + 3.14159*I

In[3]:= CosIntegral[3. I]
Out[3]= 4.96039 + 1.5708*I
```

### Applications (10)

The limits along the four axis directions

```mathematica
In[4]:= CosIntegral[{-Infinity, Infinity, -I Infinity, I Infinity}]
Out[4]= {I Pi, 0, Infinity, Infinity}
```

Machine real

```mathematica
In[5]:= CosIntegral[2.8]
Out[5]= 0.186488
```

Inside [2, 4] the real kernel uses the factored form (x - x0) Q(x) through the first zero x0 = 3.3842...

```mathematica
In[6]:= CosIntegral[3.5]
Out[6]= -0.0321285
```

A negative real lies on the cut: the value from above, Ci(2) + I Pi

```mathematica
In[7]:= CosIntegral[-2.]
Out[7]= 0.422981 + 3.14159*I
```

On the imaginary axis, CoshIntegral[3] + I Pi/2

```mathematica
In[8]:= CosIntegral[3. I]
Out[8]= 4.96039 + 1.5708*I
```

Machine complex

```mathematica
In[9]:= CosIntegral[2.5 + I]
Out[9]= 0.331465 - 0.388237*I
```

Arbitrary precision

```mathematica
In[10]:= N[CosIntegral[2], 50]
Out[10]= 0.422980828774864995698565153198255894135737756306189
```

The logarithmic Taylor series at the origin

```mathematica
In[11]:= Series[CosIntegral[x], {x, 0, 6}]
Out[11]= EulerGamma + Log[x] - 1/4 x^2 + 1/96 x^4 - 1/4320 x^6 + O[x]^7
```

The trig-prefactored asymptotic expansion

```mathematica
In[12]:= Normal[Series[CosIntegral[x], {x, Infinity, 3}]]
Out[12]= -Cos[x]/x^2 + Sin[x] (1/x - 2/x^3)
```

No symmetry fold

```mathematica
In[13]:= CosIntegral[-x]
Out[13]= CosIntegral[-x]
```

## Algorithm

```text
 Mathilda -- the cosine integral  Ci(z) = -Int_z^Inf Cos[t]/t dt.

  CosIntegral[z]   Ci(z)
```

Ci has a logarithmic singularity at z = 0 and a branch cut running along the negative real axis (-Infinity, 0]. It is NOT entire and NOT odd, so -- unlike SinIntegral -- there is no odd-symmetry fold; instead negative / left-half-plane inputs are handled through the principal Log branch and an explicit reflection. Evaluation is layered so each argument takes the cheapest route:

```text
  exact special values     ->  -Infinity, 0, I Pi, Infinity, Indeterminate
  machine real             ->  MPFR series/asymptotic at 53 bits
  arbitrary real           ->  the same, at the input precision
  complex (any precision)  ->  the ncpx series/asymptotic with guard bits
  everything else          ->  stays symbolic (return NULL)
```

The convergent Maclaurin series (DLMF 6.6.5) is

```text
  Ci(z) = EulerGamma + Log(z) + Sum_{k>=1} (-1)^k z^(2k) / (2k (2k)!),
```

valid on the principal branch for all z != 0 -- the principal Log(z) already supplies the correct +-i Pi jump across the cut, so no folding is needed for the convergent path. Its partial sums can reach magnitude ~e^|z| before the O(1)-sized answer emerges, so the MPFR paths add ~|z|/ln2 guard bits to absorb that cancellation exactly. For large |z| the convergent series is infeasible; there we use the asymptotic expansion (DLMF 6.12.4)

```text
  Ci(z) = sin(z) f(z) - cos(z) g(z)   (+ Stokes constant, see below),
    f(z) ~ Sum (-1)^k (2k)!   / z^(2k+1) = 1/z - 2!/z^3 + ...,
    g(z) ~ Sum (-1)^k (2k+1)! / z^(2k+2) = 1/z^2 - 3!/z^4 + ...,
```

(the same f, g SinIntegral computes; only the combination differs -- Si uses Pi/2 - cos f - sin g). The bare sin f - cos g asymptotic is valid in the open right half plane; on/left of the imaginary axis a piecewise-constant Stokes term restores the principal branch:

```text
  Ci(z) = sin(z) f(z) - cos(z) g(z) + i C,
    C =  0                for Re(z) > 0,
    C = +-Pi/2            for Re(z) = 0 (sign of Im z),
    C = +-Pi              for Re(z) < 0 (sign of Im z; +Pi from above on the cut).
```

(Equivalently, for Re(z) < 0 we fold w = -z into the right half plane, use the reflection Ci(z) = Ci(-z) + i Pi sign(Im z), and evaluate Ci(-z) by the clean right-half asymptotic.) For a negative real x the result is the from-above limit Complex[Ci(|x|), Pi], matching Ci(-Infinity) = I Pi.

Attributes: Listable, NumericFunction, Protected.

## Implementation notes

**Algorithm.** `builtin_cosintegral` evaluates `Ci(z) = -∫_z^∞ cos t/t dt`. The function has a logarithmic singularity at 0 and a branch cut along `(-∞, 0]`, and it is not odd. A machine call first passes the `sf_machine_scalar_args` gate and runs the kernels in `sf_machine_expint.c`.

**Real kernel** (`sf_machine_ci`), which follows the Cephes `sici` structure with coefficients refitted from mpmath:

- `0 < x <= 2`: `γ + log x + x²·C(x²)`.
- `2 < x <= 4`: `(x - x₀)·Q(x)` through the first zero `x₀ = 3.3842…`. The series form would cancel 15-fold by `x = 4`.
- `x > 4`: `f sin x - g cos x`, where the auxiliary functions use `x f = 1 - w·RF(w)` and `x² g = 1 - w·RG(w)` with `w = 1/x²`, rational on `[4, 8]` and `[8, ∞)`.
- `x < 0`: the real kernel declines, and the scalar path asks `sf_machine_ci_c` for the complex from-above value `Ci(|x|) + iπ`.

**Complex kernel** (`sf_machine_ci_c`):

- The Maclaurin series by Horner, but only where it does not cancel (`|z| <= 1` or `|z| - |Im z| <= 1`).
- Otherwise `Ci = -(E1(iz) + E1(-iz))/2`, with `E1` from the backward continued fraction.
- The reflection `Ci(z) = Ci(-z) ± iπ` in the left half-plane.
- The asymptotic series for `|z| > 45`.

A declined call goes to `ci_one_arg`:

- **Exact values:** `Ci[0] = -Infinity`, `Infinity → 0`, `-Infinity → I π`, `±I Infinity → Infinity`, and `ComplexInfinity`/`Indeterminate → Indeterminate`.
- **MPFR, moderate `|z|`:** the convergent series `γ + Log z + Σ (-1)^k z^{2k}/(2k (2k)!)` with `~|z|/ln2` guard bits. The principal `Log` supplies the `±iπ` jump across the cut.
- **MPFR, large `|z|`:** the asymptotic `sin(z) f(z) - cos(z) g(z)`, plus a piecewise Stokes constant (`0`, `±π/2` or `±π`) that restores the principal branch.

Complex MPFR values use the shared `ncpx` toolkit.

**Data structures.** `Expr`; `ncpx` (`numeric_complex.h`); refitted coefficient tables in `sf_machine_expint_coef.h`. The ND kernel is `NDKU_CosIntegral = { cosintegral_machine_complex, ndk_CosIntegral_r }`, registered with `REG_U`. Its real arm declines `x <= 0`, which hands the buffer to the `List` path, and that path returns the complex values. `Compile[]` lowers `CosIntegral` through the kernel for real and complex scalars and for rank-1 arrays. The scalar call, `List`, packed arrays, `NDArray` and `Compile[]` all share these helpers and agree bit for bit.

**Complexity / limits.** The machine kernels are `O(1)`. The real kernel's median error is 0.45 ulp, and its maximum is 3.5 ulp when scaled by `|x Ci′(x)|` (the relative error grows without bound only at the zeros). The complex kernel's maximum is 1.8 ulp, scaled the same way. The cost of the MPFR series grows with `|z|` and the precision. Exact non-special arguments, and `CosIntegral[-x]`, stay unevaluated. Attributes: `Listable`, `NumericFunction`, `Protected`.

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## References

**See also:** [Log](../../elementary-functions/Log/), [SinIntegral](../../special-functions/SinIntegral/)

- DLMF §6.2.11, §6.6.5 — the cosine integral Ci(z) and its convergent Maclaurin series.
- DLMF §6.2.17–6.2.18, §6.12 — the auxiliary functions f, g and the asymptotic expansion Ci(z) ~ sin(z) f(z) - cos(z) g(z).
- S. L. Moshier, *Methods and Programs for Mathematical Functions* (Ellis Horwood, 1989) — the Cephes `sici` structure of the real kernel.
- Source: [`src/special_functions/cosintegral.c`](https://github.com/stblake/mathilda/blob/main/src/special_functions/cosintegral.c)
- Specification: [`docs/spec/builtins/special-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/special-functions.md)
- Tests: [`tests/test_compile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile.c)
- Tests: [`tests/test_cosintegral.c`](https://github.com/stblake/mathilda/blob/main/tests/test_cosintegral.c)
- Tests: [`tests/test_expint_machine.c`](https://github.com/stblake/mathilda/blob/main/tests/test_expint_machine.c)
- Tests: [`tests/test_integrate_ramanujan.c`](https://github.com/stblake/mathilda/blob/main/tests/test_integrate_ramanujan.c)

## Notes & additional examples

### Notes

- `CosIntegral[z] = -∫_z^∞ Cos[t]/t dt` has a logarithmic singularity at 0 and a branch cut on `(-∞, 0]`; a negative real gives the from-above value `Ci(|x|) + I Pi`.
- Machine arguments use the `sf_machine_expint.c` kernels (Cephes `sici` structure on the real line with a factored form through the first zero; series or E1 continued fraction in the complex plane), shared bit for bit with packed lists, `NDArray` and `Compile[]`; arbitrary precision uses the guarded series or the Stokes-corrected asymptotic expansion.
- `D[CosIntegral[z], z] = Cos[z]/z`. Attributes: `Listable`, `NumericFunction`, `Protected`. See also `SinIntegral`, `CoshIntegral`, `ExpIntegralEi`.
