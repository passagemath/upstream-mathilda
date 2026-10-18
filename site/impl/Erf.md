---
source: src/special_functions/erf.c
references:
  - "DLMF §7.2, §7.6.2 — the error function and its Maclaurin series erf(z) = (2/sqrt(π)) e^{-z^2} Sum_n t_n."
  - "H. E. Salzer, *Formulas for calculating the error function of a complex variable*, Math. Tables Aids Comput. 5 (1951) 67–70 (A&S 7.1.29)."
  - "T. J. Dekker, *A floating-point technique for extending the available precision*, Numer. Math. 18 (1971) 224–242."
---
**Algorithm.** `builtin_erf` handles `Erf[z]` and `Erf[z0, z1] = erf(z1) - erf(z0)`. A one-argument call first passes through the shared machine-scalar gate, `sf_machine_scalar_args`. It accepts machine `Real`, exact, or `Complex` leaves, needs at least one inexact leaf and no MPFR, and runs before any special-value check. The kernels live in `sf_machine_error.c`:

- **Real** (`sf_erf_r`). For `|x| < 1/2`, a 13-term Maclaurin polynomial. Otherwise `erf = 1 - e^{-x²} erfcx(x)`, with `erfcx` from 0.5-wide polynomial tables, and for `6 <= x <= 27.4` as a polynomial in `1/x²`. `e^{-x²}` is taken from the exact Dekker square, so there is no `x²·ulp` exponent error. The error is ≤ 1.6 ulp.
- **Complex** (`sf_erf_c`). Odd and conjugate symmetry reduce the call to the first quadrant, then `erf_quadrant` picks one of four methods:
  - `|z| <= 0.8`: the Maclaurin series.
  - `|Re z|, |Im z| < 4` with `XY < 1/4`: the Taylor series about the nearer axis (a Hermite recurrence), so the small component keeps full relative accuracy.
  - `|Re z|, |Im z| < 4` otherwise: `1 - e^{-z²} w(iz)`, with the Faddeeva function `w` from a 0.5×0.5 grid of degree-20 Taylor expansions, and with `e^{-z²}` built from an exact exponent and an exactly reduced phase.
  - Elsewhere: Salzer's series (A&S 7.1.29), an exponentially convergent trapezoidal sum whose relative error is `2e^{-4π²}`.

  It declines `|Im z| > 26`, where `e^{y²}` overflows.

A declined call, and any MPFR call, goes to `erf_one_arg`:

- **Exact special values:** `0`, `±Infinity → ±1`, `±I Infinity → DirectedInfinity[±I]`, plus `ComplexInfinity` and `Indeterminate`.
- **Odd fold:** a symbolic negative-leading argument folds by oddness (`Erf[-x] = -Erf[x]`).
- **Arbitrary-precision real:** `mpfr_erf`.
- **Arbitrary-precision complex:** the Maclaurin series `erf(z) = (2/√π) e^{-z²} Σ t_n` with `t_n = t_{n-1}·2z²/(2n+1)`, summed in the file-local `ecx` toolkit with `|z|²/ln2` guard bits to absorb the `~e^{|z|²}` cancellation.

`Erf[z0, z1]` commits only when both erf calls reduce to something concrete. An `Interval` goes to `interval_apply_function`.

**Data structures.** `Expr`; `ecx` (a pair of `mpfr_t`) for the MPFR complex series. The machine kernels use only `<math.h>` (no `<complex.h>`) with double-double helpers; the generated tables are in `sf_machine_error_tables.h` (`tools/sf_fit/gen_tables.py`). The ND kernel is `NDKU_Erf = { sf_erf_c, ndk_Erf_r }`, registered with `REG_U`, covering real and complex packed arrays and `NDArray`s. `Compile[]` lowers `Erf` with the dedicated opcodes `OP_ERF_R` and `OP_VERF_R` for reals, and through the kernel for complex values, at scalar and rank-1 shapes. Every surface calls the same helper, so scalar, `List`, packed, `NDArray` and `Compile[]` agree bit for bit.

**Complexity / limits.** The machine paths are `O(1)`. Salzer's series needs about `2|y| + 13` terms. Against mpmath the complex error has a median of 0.4 ulp and is ≤ 3.5 ulp away from the zeros. The MPFR series' term count and guard bits grow with `|z|²` and with the precision. Known gap: when `|Im z| > 26` and the result exceeds the double range, the machine kernel declines and the MPFR fallback answers `Overflow[]`. `Erfc` and `Erfi` return extended-exponent values in the same region. Exact non-special arguments (`Erf[2]`) stay unevaluated. Attributes: `Listable`, `NumericFunction`, `Protected`.
