---
source: src/special_functions/erfi.c
references:
  - "DLMF §7.2, §7.6.2 — the imaginary error function erfi(z) = erf(iz)/i = (2/sqrt(π)) Int_0^z e^{t^2} dt."
  - "DLMF §7.2(ii), §7.12 — Dawson's integral and the asymptotic expansion behind the large-x tail."
  - "H. E. Salzer, *Formulas for calculating the error function of a complex variable*, Math. Tables Aids Comput. 5 (1951) 67–70."
---
**Algorithm.** `builtin_erfi` handles `Erfi[z] = -i erf(iz)`, which is entire and odd. Neither libm nor MPFR provides `erfi`, so Mathilda implements every numeric path itself. A machine call first passes through the `sf_machine_scalar_args` gate and the shared kernels in `sf_machine_error.c`:

- **Real** (`sf_erfi_r`):
  - For `|x| < 1/2`, the all-positive Maclaurin polynomial.
  - For `1/2 <= |x| < 10`, `erfi(x) = e^{x²}·g(x)`. Here `g = e^{-x²} erfi(x) = (2/√π)·Dawson(x)` comes from 0.5-wide polynomial tables (`DAW1_T`, `DAW2_T`), and `e^{x²}` comes from the exact Dekker square.
  - For `|x| >= 10`, a 14-term asymptotic series for `g`.
  - When `e^{x²}` alone would overflow but `erfi` still fits (`x² ∈ (709, 713.6)`), the kernel scales by `e^{-64}` first.

  The error is ≤ 2.6 ulp; the old summed series reached 30 ulp at `|x| ≈ 5`. The kernel declines `|x| >= 26.72`.
- **Complex** (`sf_erfi_c`): `-i erf(iz)`, routed through the complex `erf` kernel (Maclaurin, axis Taylor series, Faddeeva table, Salzer's series).

Anything declined, or any MPFR call, goes to `erfi_one_arg`:

- **Exact special values:** `0`, `±Infinity → ±Infinity`, and the finite imaginary-axis limits `±I Infinity → ±I`.
- **Odd fold:** `Erfi[-x] = -Erfi[x]`.
- **Real at any precision:** the all-positive Maclaurin series `(2/√π) Σ x^{2n+1}/(n!(2n+1))` in MPFR. Nothing cancels, so a flat 64-bit guard is enough.
- **Complex at any precision:** `-i erf(iz)` through the guarded erf series in `ecx`.

Results pass through `numeric_promote_result_if_degenerate`, so `Erfi[30.]` comes back as an extended-exponent number, not `Infinity`.

**Data structures.** `Expr`; `ecx` (a pair of `mpfr_t`) and plain MPFR scalars. The ND kernel is `NDKU_Erfi = { sf_erfi_c, ndk_Erfi_r }`, registered with `REG_U`, covering real and complex packed arrays and `NDArray`s. `Compile[]` lowers `Erfi` through the kernel for real and complex scalars and for rank-1 arrays. Every surface shares one helper, so the answers agree bit for bit.

**Complexity / limits.** The machine kernels are `O(1)`. The MPFR real series needs `O(x²)` terms; the complex series' term count and guard bits grow with `|z|²` and with the precision. Exact non-special arguments (`Erfi[2]`) stay unevaluated. Attributes: `Listable`, `NumericFunction`, `Protected`.
