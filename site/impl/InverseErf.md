---
source: src/special_functions/inverf.c
references:
  - "DLMF §7.17 — inverse error functions."
  - "M. Giles, *Approximating the erfinv function*, in GPU Computing Gems Jade Edition (Morgan Kaufmann, 2011) — the w = -log(1 - p^2) variable."
  - "S. Winitzki, *A handy approximation for the error function and its inverse* (2008) — the closed-form fallback seed."
---
**Algorithm.** `builtin_inverf` solves `s = erf(z)`. Following Mathematica, it returns explicit values only for real `s` in `[-1, 1]`; complex `s` and `|s| > 1` stay unevaluated. A one-argument machine-real call (passing the `sf_machine_scalar_args` gate, not complex) goes first to `sf_inverse_erf_r` in `sf_machine_error.c`. With `p = |s|` and its exact complement `q = 1 - p`:

- For `p < 0.875`, a degree-14 polynomial for `erfinv(p)/p` on 1/16-wide pieces, with no logarithm.
- Otherwise the kernel switches to Giles' variable `w = -log(q(2 - q))`, in which `erfinv(p)/p` is analytic in the strip `|Im w| < π`. For `w < 6` it uses unit-width degree-11 pieces. Past that it takes a `1e-7` seed polynomial in `1/√w`, then one Halley step on `erfc(z) = q`, written with the relative residual and `erfcx` so that neither `e^{±z²}` is ever formed.

No `erf` is evaluated on the main path. The kernel declines `±1`, `|s| > 1` and a subnormal `q`.

Declined calls go to `inverf_one_arg`:

- **Exact values:** `0 → 0`, `±1 → ±Infinity` (also for `±1.`).
- **Odd fold:** `InverseErf[-x] = -InverseErf[x]`.
- **Arbitrary-precision real** `|s| < 1`: Newton's method on `mpfr_erf`, `z ← z - (erf(z) - s)(√π/2) e^{z²}`, doubling the working precision at each step. The seed comes from `inverf_double`, which is the polynomial kernel above, falling back to Winitzki's closed form plus libm-`erf` Newton steps where the kernel declines.

The two-argument form `InverseErf[z0, s]` rewrites to `InverseErf[s + Erf[z0]]`, and stays in two-argument form if `Erf[z0]` does not reduce. The derivative `(√π/2) e^{InverseErf[z]²}` is defined in `calculus/deriv.c`, and `Series` builds on it through the generic Taylor-via-`D` route.

**Data structures.** `Expr`; `mpfr_t` for the Newton loop; the generated polynomial tables `INVP_T`, `INVW_T` and `INVT_T` in `sf_machine_error_tables.h`. The ND kernel is real-only, `NDKU_InverseErf = { NULL, ndk_InverseErf_r }`, registered with `REG_U`, so it is packed-aware. `Compile[]` lowers `InverseErf` at real scalar and rank-1 shapes through that kernel; a `_Complex` argument does not compile, since the function has no complex values. Scalar, packed, `NDArray` and `Compile[]` share `sf_inverse_erf_r` and agree bit for bit.

**Complexity / limits.** The machine kernel is a few polynomial evaluations, with error ≤ 1.5 ulp. The old kernel was Winitzki's closed form plus up to four Newton steps, each costing an `erf` and an `exp`. The MPFR cost is dominated by the final full-precision `mpfr_erf`. Exact rationals (`InverseErf[1/2]`), complex inputs and `|s| > 1` stay unevaluated. Attributes: `Listable`, `NumericFunction`, `Protected`.
