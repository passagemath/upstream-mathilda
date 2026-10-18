---
source: src/special_functions/legendre.c
references:
  - "DLMF §14 — Legendre and related functions (§14.3 series, §14.7 integer degree, §14.10.3 the three-term recurrence)."
---
**Algorithm.** `builtin_legendre_p` handles `LegendreP[n, x]`,
`LegendreP[n, m, x]` and `LegendreP[n, m, a, x]`. For the two-argument form it
first tries the **machine gate** `legendre_machine`: when both arguments are
machine-scalar leaves with at least one machine `Real` (`sf_machine_scalar_args`),
the order is integer-valued and real, and `x` is not an exact `1`, it calls
`sf_machine_legendre_p_c` (`src/special_functions/sf_machine_plog.c`). That runs the
three-term recurrence `k P_k = (2k-1) x P_{k-1} - (k-1) P_{k-2}` in `double`
(real `x`, via `sf_machine_legendre_p` in `sf_machine.c`) or in re/im doubles
(complex `x`), maps negative orders with `P_{-n-1} = P_n`, and declines above
degree 1000. Otherwise `legendre_two_arg` does the exact work: an exact integer
order builds the explicit degree-`|n'|` polynomial from the same recurrence with
GMP `mpq` coefficients (`legendre_p_coeffs`, cap `LEG_POLY_CAP = 2000`);
`x == 1 -> 1` for any order; a non-integer order with an inexact argument sums the
Gauss series `P_n(x) = 2F1(-n, n+1; 1; (1-x)/2)` in the `ncpx` MPFR-complex toolkit
(`leg_numeric_p`, real or complex, machine or arbitrary precision; requires
`|(1-x)/2| < 1`). The associated forms (integer `n`, integer `m >= 0`): type 1 is
the Rodrigues derivative `(-1)^m (1-x^2)^(m/2) d^m/dx^m P_n(x)` (0 when
`m > |n'|`); types 2 and 3 are the regularized Gauss polynomial
`2F1Reg(-n, n+1, 1-m, (1-x)/2)` (`legendre_2f1reg_core`) times
`(1+x)^(m/2) (1-x)^(-m/2)` or `(1+x)^(m/2) (-1+x)^(-m/2)`. For a machine buffer
argument, `leg_poly_from_coeffs` evaluates these polynomials by Horner in `double`
instead of as a sum of monomials, each of which would cost a `pow()` per element.
`D` for symbolic order lives in `calculus/deriv.c`.

**Data structures.** `Expr`; GMP `mpq_t` coefficient arrays for the polynomials;
`ncpx` (`mpfr_t` re/im) for the non-integer-order series; plain doubles (no
`<complex.h>`) in the machine kernels. ND: binary kernel
`NDKB_LegendreP = { ndk_LegendreP_c, ... }`, which calls the same
`sf_machine_legendre_p_c`, registered `REG_B` (so `packed_aware`); it threads `n`
and `x` element-wise. Attributes: `Listable`, `NumericFunction`, `Protected`.

**Complexity / limits.** The machine recurrence is `O(n)` per point, up to degree
1000. The exact recurrence is `O(n^2)` in growing big rationals, capped at
`n = 2000`. The non-integer-order series costs `O(wp)` terms inside the
convergence disc. The Horner buffer path is `O(n)` fused multiply-adds per element
(degree 30 over 10^5 points about 10 ms). Left symbolic: non-integer or negative
`m`, symbolic `Series`, and non-integer order with `|(1-x)/2| >= 1`
(`LegendreP[0.5, 3.]`). `Compile[]` lowers `LegendreP[n, x]` at real and complex
scalar shapes and at rank-1 array shape (`Compiled -> True`), through the same
kernel. The associated three-argument form is not lowered.
