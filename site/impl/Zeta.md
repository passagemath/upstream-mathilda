---
source: src/special_functions/zeta.c
references:
  - "DLMF §25.2 — the Riemann zeta function; §25.4.1 — the functional equation; §25.11 — the Hurwitz zeta function."
  - "DLMF §25.11.5 / §25.2.9 — the Euler-Maclaurin summation formula used by the machine and complex-MPFR kernels."
---
**Algorithm.** `builtin_zeta` evaluates `Zeta[s]` (Riemann) and `Zeta[s, a]`
(Hurwitz) in a fixed dispatch order. (1) An `Interval` argument goes through
`interval_apply_function`. (2) **Machine-scalar gate**: when
`sf_machine_scalar_args` classifies the call as machine precision (some machine
`Real`, no MPFR leaf), the raw-`double` kernels of `sf_machine_zeta.c` answer —
`sf_zeta_r` / `sf_zeta_c` for one argument, `sf_hurwitz_r` / `sf_hurwitz_c` for
two (which only accept `Re a > 0`, where the symmetric and principal powers
coincide). `sf_zeta_r` reads tabulated ζ(2…64) and ζ(0…−60) at integers, sums the
Dirichlet series directly for `s ≥ 20`, uses Euler–Maclaurin for `s ≥ 1/2`, and
the functional equation (DLMF 25.4.1, with sin(πs/2) reduced exactly and the
rounding of `1 − s` corrected to first order in Γ and ζ) below 1/2; complex
`Re s < 1/2` declines. A decline (pole, domain, overflow) falls through. (3)
**Exact** one-argument closed forms: `1 → ComplexInfinity`, `0 → -1/2`, even
`2n > 0 →` rational·π²ⁿ from an exact Bernoulli cache, negative `−m →` rational
(zero at even `m`), odd positive integers stay symbolic, `Zeta[Infinity] = 1`.
(4) **Arbitrary precision** (an MPFR real, or a complex with an MPFR part) goes
to arb's `acb_dirichlet_zeta` via `flint_num_zeta`; then `mpfr_zeta` for any
remaining real `s`, and the complex-MPFR Euler–Maclaurin kernel
`zeta_hurwitz_cx` (with `a = 1`) for a complex `s` the gate declined. The
two-argument form reduces exactly at `a = 1`, `a = 1/2` (`(2^s − 1) Zeta[s]`)
and positive integer `a` (`Zeta[s] − Σ_{k<a} k⁻ˢ`); any other inexact input runs
`zeta_numeric` (complex-MPFR Euler–Maclaurin), which uses the **symmetric** power
`((a+k)²)^(−s/2)` so `Re a < 0` gives Mathematica's real value.

**Data structures.** `Expr`; an exact `mpq` Bernoulli cache for the closed forms;
the file-local `zcx` (`mpfr_t` re/im) toolkit for the Euler–Maclaurin kernel;
`double` tables (`ZM_ZPOS`, `ZM_ZNEG`, `ZM_B2K_FACT`) and TwoSum/double-double
helpers in the machine kernel; complex machine results share one `Complex` head
node (`sf_machine_complex_result`). ND: unary kernel `NDKU_Zeta` with both a real
arm (`sf_machine_zeta`) and a complex arm (`sf_zeta_c`), registered `REG_U`, so
packed and visible-`NDArray` arguments stay on the buffer. Attributes:
`Listable`, `NumericFunction`, `Protected`.

**Complexity / limits.** Machine values cost well under a microsecond (one table
read or a short Euler–Maclaurin sum, `w ≥ max(10, |s|/2 + 4)`, correction cut at
2⁻⁵⁶). Arbitrary-precision cost is arb's. Exact-integer cap
`ZETA_EXACT_INT_CAP = 10000`; Hurwitz-integer-`a` cap `100000`. Complex
`Re s < 1/2` at machine precision is computed on the MPFR path (slower but
accurate). `Compile[]` lowers `Zeta[s]` at scalar and rank-1 array shapes, real
and complex (`Compiled -> True`); the two-argument `Zeta[s, a]` does **not**
compile — use `HurwitzZeta`, which does (identical for `Re a > 0`).
