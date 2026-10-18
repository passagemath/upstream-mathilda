---
source: src/special_functions/gamma.c
references:
  - "DLMF §5.2, §5.5 — the Euler gamma function and the reflection formula."
  - "DLMF §5.11 — Stirling's series for ln Gamma."
  - "DLMF §8.2, §8.11 — the (upper) incomplete gamma function Gamma(a,z)."
  - "J. L. Spouge, *Computation of the gamma, digamma, and trigamma functions*, SIAM J. Numer. Anal. 31 (1994) 931–944."
  - "C. Lanczos, *A precision approximation of the gamma function*, SIAM J. Numer. Anal. Ser. B 1 (1964) 86–96."
  - "D. E. G. Hare, *Computing the principal branch of log-Gamma*, J. Algorithms 25 (1997) 221–236."
  - "T. J. Dekker, *A floating-point technique for extending the available precision*, Numer. Math. 18 (1971) 224–242."
---
**Algorithm.** `builtin_gamma` serves `Gamma[z]`, `Gamma[a, z]` (upper incomplete) and `Gamma[a, z0, z1] = Gamma[a, z0] - Gamma[a, z1]`; an `Interval` argument goes to `interval_apply_function`. `gamma_one_arg` tests the cheap numeric cases first, then the exact ones:

1. **Machine real** → libm `tgamma`. A pole (`NaN`) gives `ComplexInfinity`. An overflow (`Gamma[172.]` and beyond) is recomputed by `mpfr_gamma` at 53 bits, so the answer is a machine-precision number with an extended exponent and not `Infinity`.
2. **Complex with an MPFR part** → `gamma_mpfr_complex`: Spouge's approximation in the file-local complex-MPFR toolkit `gcx`. Its coefficients are computed at run time for the working precision (output precision + 64 guard bits), so the requested digits are correct.
3. **Machine complex** → `gamma_machine_complex`, which calls `sfg_gamma_c` in `sf_machine_gamma.c`. For `-5 <= Re z < 10`, `|Im z| < 19/4`, the exact recurrence moves `z` into the strip `1 <= Re < 5/2`, and lnΓ comes from a Taylor series about the nearest of 30 precomputed centres. Elsewhere the kernel shifts to `Re w >= 7` and sums Stirling's series. The `N` recurrence logarithms are taken as one log of the running product, with its windings counted (Hare's continuous branch). Γ is then `exp` of a double-double exponent, built with Dekker products, times the product correction. A real-valued result (`Im = 0`) comes back as a `Real`. Inputs the kernel declines (`Re z < -30`, overflow) fall back to the old `g = 7`, 9-term Lanczos series with reflection.
4. **Exact integer / half-integer** → `(z-1)!` through the `Factorial` machinery. This gives an exact integer, a `BigInt`, or a rational multiple of `Sqrt[Pi]`, and works for negative half-integers too. Non-positive integers are poles (`ComplexInfinity`).
5. **Symbolic infinities**: `Infinity → Infinity`, `-Infinity → Indeterminate`, `ComplexInfinity → ComplexInfinity`.
6. **Arbitrary-precision real** → `mpfr_gamma` at the input precision.

The incomplete form does the exact rewrites first: `Gamma[a, 0] = Gamma[a]` and `Gamma[a, Infinity] = 0`. A positive-integer `a <= 1000` with exact or symbolic `z` expands to `e^{-z} Sum_k ((n-1)!/k!) z^k`, with the coefficients built in GMP (`gamma_incomplete_int`). For a numeric real it calls `mpfr_gamma_inc` at the minimum input precision, keeping a 53-bit MPFR value when a double would overflow or underflow. For numeric complex it uses `gcx_inc_gamma`: the lower series `gamma(a,z) = z^a e^{-z} Sum z^n/(a)_{n+1}` when `Re z < Re a + 1`, otherwise a Lentz continued fraction, then multiplies by `e^{-z} z^a`.

**Data structures.** `Expr` trees; `gcx` (a pair of `mpfr_t`; alias-safe, explicit precision) for the complex and incomplete MPFR work; GMP for the factorial and closed-form coefficients. The machine kernel works on `(re, im)` doubles, without `<complex.h>`. Its tables (tile Taylor coefficients, Stirling `B_2k/(2k(2k-1))`) come from `tools/sf_accuracy/gen_gamma_tables.py` (mpmath), and `tools/sf_accuracy/gamma.py` gates their accuracy. The ND kernel is `NDKU_Gamma = { gamma_machine_complex, ndk_Gamma_r }`, registered with `REG_U`, so packed arrays and visible `NDArray`s are mapped directly: real elements use `tgamma` (a pole or overflow declines to the `List` path), complex elements use the same `gamma_machine_complex` as the scalar call. `Compile[]` lowers `Gamma` through that kernel for real scalars, complex scalars and rank-1 arrays. An integer argument uses the integer-closed opcode `OP_GAMMA_I`, so `Gamma[5]` stays the exact `24`. The scalar, packed, `NDArray` and `Compile[]` answers agree bit for bit.

**Complexity / limits.** The machine paths are `O(1)` per element. The complex kernel has a median error of 0.8 ulp and a maximum of about 4 ulp in the box `0.5 <= Re z <= 4`, `|Im z| <= 3` (the Lanczos series it replaced reached 50 ulp). The MPFR paths scale their term counts and guard bits with the precision. Exact complex arguments (`Gamma[3/2, I]`), exact non-integer `a` in `Gamma[a, z]`, and symbolic `z` stay unevaluated. The incomplete integer expansion is capped at `a <= 1000`. Attributes: `Listable`, `NumericFunction`, `Protected`.
