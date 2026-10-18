---
source: src/special_functions/coshintegral.c
references:
  - "DLMF §6.2.16 — the hyperbolic cosine integral Chi(z) = γ + Log(z) + Int_0^z (cosh t - 1)/t dt."
  - "S. L. Moshier, *Methods and Programs for Mathematical Functions* (Ellis Horwood, 1989) — the Cephes `shichi` structure of the real kernel."
  - "W. J. Cody and H. C. Thacher, *Chebyshev approximations for the exponential integral Ei(x)*, Math. Comp. 23 (1969) 289–303."
---
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
