---
source: src/special_functions/sinintegral.c
references:
  - "DLMF §6.2 — the sine integral Si."
  - "DLMF §6.6.5 — convergent Maclaurin series; §6.9.1 — the E1 continued fraction; §6.12.3 — asymptotic expansion."
  - "S. L. Moshier, *Methods and Programs for Mathematical Functions* (Ellis Horwood, 1989) — the Cephes `sici` structure of the real kernel."
---
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
