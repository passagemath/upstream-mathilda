---
source: src/special_functions/cosintegral.c
references:
  - "DLMF §6.2.11, §6.6.5 — the cosine integral Ci(z) and its convergent Maclaurin series."
  - "DLMF §6.2.17–6.2.18, §6.12 — the auxiliary functions f, g and the asymptotic expansion Ci(z) ~ sin(z) f(z) - cos(z) g(z)."
  - "S. L. Moshier, *Methods and Programs for Mathematical Functions* (Ellis Horwood, 1989) — the Cephes `sici` structure of the real kernel."
---
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
