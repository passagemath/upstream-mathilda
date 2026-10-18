---
source: src/special_functions/fresnel.c
references:
  - "DLMF §7.2(iii) — the Fresnel integral S(z) = Int_0^z sin(π t^2/2) dt."
  - "DLMF §7.3.22, §7.5.3–4 — Fresnel integrals through erf, and through the auxiliary functions f, g."
  - "DLMF §7.12 — the asymptotic expansion of the Fresnel integrals."
  - "T. J. Dekker, *A floating-point technique for extending the available precision*, Numer. Math. 18 (1971) 224–242."
---
**Algorithm.** `builtin_fresnels` returns `S(z) = ∫₀ᶻ sin(π t²/2) dt`, using the π/2 convention. `S` is entire and odd, and it shares every kernel with `FresnelC`: the pair `(C, S)` is computed together and this builtin keeps the `S` component.

A machine call first goes through `fresnel_machine`, which uses the `sf_machine_scalar_args` gate.

**Real argument** (`sf_fresnel_s_r`):

- `|x| <= 1`: `S = x³ Σ s_k x^{4k}`, with `x³` formed to within a quarter ulp.
- `1 < |x| < 8`: `S = 1/2 - f cos(πx²/2) - g sin(πx²/2)`, where the auxiliary functions `f` and `g` come from 0.5-wide tables.
- `|x| >= 8`: the asymptotic series for `f` and `g`.
- Throughout, the phase `πx²/2` is reduced exactly: `x²` is held as a Dekker double-double and reduced mod 4.

**Complex argument** (`sf_fresnel_s_c`):

- `|z| <= 1.3`: the Maclaurin series.
- Next to an axis, and for huge `|z|` with `π|xy| < 8`: the scaled axis Taylor series. The imaginary axis is handled through `S(iz) = -iS(z)`.
- Elsewhere: `S = [(1+i)E₁ - (1-i)E₂]/(4i)`, where `E₁` and `E₂` are erf values of the exactly rotated arguments.

**Declined and MPFR calls** go to `fresnel_one_arg`:

- **Exact values:** `0`; `±Infinity → ±1/2`; `±I Infinity → ∓I/2` (note the sign flip); `ComplexInfinity` and `Indeterminate → Indeterminate`.
- **Odd fold** for symbolic negative arguments.
- **MPFR real:** the guarded series, or the DLMF 7.12 asymptotic form at large `|x|`.
- **MPFR complex:** the paired `A`/`B` series in `ncpx`, with `S = (A - B)/(2i)`.

**Data structures.** `Expr`; `ncpx` for MPFR complex values; double-double helpers plus generated tables on the machine side. The ND kernel is `NDKU_FresnelS = { sf_fresnel_s_c, ndk_FresnelS_r }`, registered with `REG_U`, and covers real and complex buffers. `Compile[]` lowers `FresnelS` through it for real and complex scalars and for rank-1 arrays. All surfaces call the same helper, so their results match bit for bit.

**Complexity / limits.** The machine kernels run in `O(1)`. The real kernel stays within 2.5 ulp. Away from the zeros the complex kernel stays within 3.5 ulp, and near an axis each component is accurate relative to its own size (`FresnelS[1.*^10 + 1.*^-10 I]` has an imaginary part of `-3.25*10^-30`). The MPFR cost grows with `|z|²` and with the requested precision. Exact non-special arguments stay unevaluated. `D[FresnelS[x], x] = Sin[π x²/2]`. Attributes: `Listable`, `NumericFunction`, `Protected`.
