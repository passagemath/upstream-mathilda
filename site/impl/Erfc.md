---
source: src/special_functions/erfc.c
references:
  - "DLMF §7.2, §7.6.2 — the complementary error function erfc(z) = 1 - erf(z)."
  - "H. E. Salzer, *Formulas for calculating the error function of a complex variable*, Math. Tables Aids Comput. 5 (1951) 67–70 (A&S 7.1.29)."
---
**Algorithm.** `builtin_erfc` handles `Erfc[z] = 1 - erf(z)`. It is entire and, unlike `Erf`, has no symmetry fold: `Erfc[-x]` stays as it is. A machine call first passes through the `sf_machine_scalar_args` gate and the shared kernels in `sf_machine_error.c`:

- **Real** (`sf_erfc_r`). For `x >= 1/2`, `erfc = e^{-x²} erfcx(x)`, with `erfcx` from the same tables `Erf` uses and `e^{-x²}` from the exact Dekker square, so nothing cancels. The error is ≤ 3.1 ulp, against 5.1 ulp measured for the macOS libm. Past `x ≈ 26.5` the value drops into the subnormal range in a single rounding (`Erfc[27.] ≈ 5.2*10^-319`), and it reaches 0 at `x = 27.4`. For `x < 0` the kernel returns `2 - erfc(-x)`, and for `|x| < 1/2` it returns `1 - erf`.
- **Complex** (`sf_erfc_c`), for `Re z >= 0` outside a small disc:
  - Next to the real axis: `erfc(X) - D`, using the axis Taylor series, so each component keeps its own relative accuracy.
  - Next to the imaginary axis: `1 - erf`.
  - Inside the `|Re|, |Im| < 4` box: `e^{-z²} w(iz)`, from the Faddeeva table.
  - Beyond the box: `erfc(X) - (erf(z) - erf(X))`, taken from Salzer's series, so a small `erfc` does not cancel against 1.

  The kernel declines when `Re z > 26` or `|Im z| > 26`. For `Re z < 0` it uses `1 - erf(z)`, which lies between 1 and 2, so nothing cancels.

Anything else goes to `erfc_one_arg`:

- **Exact special values:** `0 → 1`, `Infinity → 0`, `-Infinity → 2`, `±I Infinity → DirectedInfinity[∓I]`, plus `ComplexInfinity` and `Indeterminate`.
- **Arbitrary-precision real:** `mpfr_erfc`, which does not cancel even for large `z`.
- **Arbitrary-precision complex:** `erf` from the guarded DLMF 7.6.2 series in `ecx`, with the complement formed at working precision.

The result passes through `numeric_promote_result_if_degenerate`, so a value that leaves the double range comes back as an extended-exponent number. An `Interval` goes to `interval_apply_function`.

**Data structures.** `Expr`; `ecx` (a pair of `mpfr_t`). The ND kernel is `NDKU_Erfc = { sf_erfc_c, ndk_Erfc_r }`, registered with `REG_U`, covering real and complex buffers. `Compile[]` uses `OP_ERFC_R` and `OP_VERFC_R` for reals and the kernel for complex values, at scalar and rank-1 shapes. Every surface shares one helper, so the answers agree bit for bit.

**Complexity / limits.** The machine paths are `O(1)`. The MPFR complex series costs more as `|z|²` and the precision grow. Exact non-special arguments (`Erfc[2]`) and symbolic arguments stay unevaluated, with no `2 - Erfc[x]` rewrite. Attributes: `Listable`, `NumericFunction`, `Protected`.
