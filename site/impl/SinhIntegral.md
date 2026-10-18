---
source: src/special_functions/sinhintegral.c
references:
  - "DLMF §6.2 — the hyperbolic sine integral Shi."
  - "S. L. Moshier, *Methods and Programs for Mathematical Functions* (Ellis Horwood, 1989) — the Cephes `shichi` structure of the real kernel."
  - "W. J. Cody and H. C. Thacher, *Chebyshev approximations for the exponential integral Ei(x)*, Math. Comp. 23 (1969) 289–303."
---
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
