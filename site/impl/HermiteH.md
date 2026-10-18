---
source:
  - src/special_functions/hermite_h.c
  - src/special_functions/chebyshev_common.c
  - src/special_functions/sf_machine.c
references:
  - "DLMF §18 — classical orthogonal polynomials (§18.9 recurrence, §18.3 orthogonality); §13.2 for the Kummer-function form of H_ν."
---
**Algorithm.** `builtin_hermite_h` calls `hermite_h_two_arg`. (1) A non-integer order at `x = 0` returns `2^n Sqrt[Pi]/Gamma[(1-n)/2]` (`hermite_value_at_zero`). (2) A non-negative integer order with an inexact `x`, at machine precision with degree ≤ 1000, goes to `sf_machine_hermite_h` / `sf_machine_hermite_h_c`, the real and complex raw-double recurrences `H_k = 2x H_{k-1} - 2(k-1) H_{k-2}`. Otherwise it runs `hermite_numeric_recurrence` in MPFR. (3) An exact non-negative integer order with an exact or symbolic `x` builds the polynomial from FLINT's `fmpz_poly_hermite_h` (the integer recurrence when FLINT is absent). (4) A non-integer, complex or negative order with an inexact argument goes to `hermite_numeric_hyp`, which builds `2^n Sqrt[Pi] (1F1(-n/2; 1/2; z²)/Γ((1-n)/2) - 2z 1F1((1-n)/2; 3/2; z²)/Γ(-n/2))` (`hermite_hyp_form`). It numericalises that form at the output precision plus `hermite_hyp_guard_bits`: about `2·Re(z²)·log₂e + 64` bits, capped at 8192, which absorbs the cancellation between two terms that each grow like `e^{Re z²}`. It then rounds back to the input precision.

**Data structures.** `mpz_t` coefficient arrays hold the exact polynomial, `ncpx` MPFR pairs carry the recurrence, and `NumericSpec` targets (`NUMERIC_MODE_MPFR` then `NUMERIC_MODE_MPFR_CAP`) handle the guarded 1F1 evaluation. The binary ND kernel `ndk_HermiteH_c` (`src/ndkernels.c`) calls the same `sf_machine_hermite_h*` helpers, so packed and `NDArray` buffers, real or complex, stay on the buffer. `Compile[]` lowers the head at scalar and array shapes, and every surface agrees bit for bit.

**Complexity / limits.** The recurrence costs `O(n)` operations. The explicit polynomial is capped at degree 2000 (`HERMITE_POLY_CAP`). A non-integer order with an exact or symbolic argument stays unevaluated, and so does a negative integer order with an exact argument (there is no `Erfc` closed form). The derivative `2 n H_{n-1}` is in `src/calculus/deriv.c`. The asymptotic series at `Infinity` and the symbolic-index `SeriesCoefficient` are bespoke kernels in `src/calculus/series.c`; the asymptotic series omits the subdominant `e^{x²}` branch.
