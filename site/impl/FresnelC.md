---
source: src/special_functions/fresnel.c
references:
  - "DLMF §7.2(iii) — the Fresnel integral C(z) = Int_0^z cos(π t^2/2) dt."
  - "DLMF §7.3.22, §7.5.3–4 — Fresnel integrals through erf, and through the auxiliary functions f, g."
  - "DLMF §7.12 — the asymptotic expansion of the Fresnel integrals."
  - "T. J. Dekker, *A floating-point technique for extending the available precision*, Numer. Math. 18 (1971) 224–242."
---
**Algorithm.** `builtin_fresnelc` returns `C(z) = ∫₀ᶻ cos(π t²/2) dt` (Mathematica's π/2 convention). `C` is entire and odd. `FresnelS` shares every kernel: the pair `(C, S)` is computed together and each builtin keeps its own component.

A machine call first goes to `fresnel_machine`, which uses the `sf_machine_scalar_args` gate and the kernels in `sf_machine_error.c`.

**Real argument** (`sf_fresnel_c_r`):

- `|x| <= 1`: the Maclaurin series, with `x³` formed exactly.
- `1 < |x| < 8`: `C = 1/2 + f sin(πx²/2) - g cos(πx²/2)`, with the auxiliary functions `f` and `g` read from 0.5-wide polynomial tables (the Cephes `fresnl` structure).
- `|x| >= 8`: the nine-term asymptotic series for `f` and `g`.
- `|x| > 2^54`: exactly `1/2`.

The phase `πx²/2` is reduced exactly: `x²` is split as a Dekker double-double and reduced mod 4. A plain libm `sin(π x x/2)` would carry about `x²` ulps of phase error, so this keeps even `FresnelC[123456.789]` accurate to within a few ulp.

**Complex argument** (`sf_fresnel_c_c`):

- `|z| <= 1.3`: the Maclaurin series in `z⁴`.
- Next to an axis (`π|xy| < 1/2`), and also for huge `|z|` with `π|xy| < 8`: a scaled Taylor series about the real axis. The imaginary axis is mapped onto it with `C(iz) = iC(z)`. This keeps the small component relatively accurate and stops `FresnelC[1.*^10 + 1.*^-10 I]` from overflowing.
- Elsewhere: `C ± iS = ((1 ± i)/2) erf((√π/2)(1 ∓ i) z)` through the complex erf core. Every quantity that gets exponentiated (`X² = (π/4)(x+y)²`, `XY`, `Y² - X² = ∓πxy`) is built from the exact `x` and `y` rather than from a rounded rotated argument, which avoids about `π|z|²` ulps of error.

**Declined and MPFR calls** go to `fresnel_one_arg`:

- **Exact special values:** `0`, `±Infinity → ±1/2`, `±I Infinity → ±I/2`, and `ComplexInfinity` or `Indeterminate → Indeterminate`.
- **Odd fold:** `FresnelC[-x] = -FresnelC[x]`.
- **MPFR real:** the convergent series with about `(π/2)x²/ln2` guard bits, or the DLMF 7.12 asymptotic expansion at optimal truncation for large `|x|`. The asymptotic form is used on real inputs only, because the constant `1/2` changes across Stokes lines.
- **MPFR complex:** the paired `A = C + iS`, `B = C - iS` series in the `ncpx` toolkit, which is valid everywhere.

**Data structures.** `Expr`; `ncpx` (`numeric_complex.h`) for MPFR complex values. The machine side uses double-double helpers and the tables `FRF*_T` and `FRG*_T`. The ND kernel is `NDKU_FresnelC = { sf_fresnel_c_c, ndk_FresnelC_r }`, registered with `REG_U`, and covers real and complex buffers. `Compile[]` lowers `FresnelC` through it for real and complex scalars and for rank-1 arrays. Because every surface uses the same helper, scalar, packed, `NDArray` and `Compile[]` results match bit for bit.

**Complexity / limits.** The machine kernels are `O(1)`: real error is at most 2.5 ulp (median 0.3), and complex error is at most 3.5 ulp away from the zeros. The MPFR cost grows with `|z|²` and with the requested precision. Exact non-special arguments (`FresnelC[2]`) stay unevaluated. `D[FresnelC[x], x] = Cos[π x²/2]`. Attributes: `Listable`, `NumericFunction`, `Protected`.
