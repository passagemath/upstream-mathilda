---
source: src/special_functions/beta.c
references:
  - "DLMF §5.12 — the Euler beta function B(a,b) = Gamma(a) Gamma(b)/Gamma(a+b)."
  - "DLMF §5.15 — the digamma function psi, used for the first-order rounding correction of a+b."
  - "DLMF §8.17 — the incomplete beta function B_z(a,b) and its relation to the Gauss hypergeometric function."
---
**Algorithm.** `builtin_beta` checks the machine case before anything else. When `sf_machine_scalar_args` accepts a two-argument call (machine `Real`, exact, or `Complex` leaves, with at least one inexact leaf and no MPFR), it computes the value directly in double:

- `sfg_beta` (real) takes `tgamma(a)/tgamma(a+b)·tgamma(b)`. Because `a + b` is rarely representable, it undoes that rounding to first order through `ψ(a+b)`, using `Γ(s+e) ≈ Γ(s)(1 + ψ(s) e)`.
- `sfg_beta_c` (complex) combines the three Gammas under one exponential, summing their double-double exponents from the `sf_machine_gamma.c` tiles or Stirling series, with the same ψ correction.

A decline falls through to the layered path below: a pole, `|a|`, `|b|` or `|a+b| > 170`, or underflow.

The layered path builds on the `Gamma` machinery. `beta_lattice_case` first sorts out the gamma poles on the integer lattice: a pole that survives gives `ComplexInfinity`, a cancelling pair gives its finite limit (`Beta[-2, 1] = -1/2`), and a pole only in the denominator gives `0`. Purely symbolic arguments stay unevaluated. If one argument is a positive integer `n`, the call collapses to `(n-1)!/Pochhammer[b, n]`, which is exact for rational `b` (`Beta[3, 1/3] = 27/14`). Otherwise the call is rebuilt as `Gamma[a] Gamma[b]/Gamma[a+b]` and evaluated. That keeps half-integers exact (`Beta[1/2, 1/2] = Pi`), keeps other rationals as a Gamma quotient, and sends MPFR and complex arguments to `Gamma`'s arbitrary-precision paths.

The incomplete form uses `Beta[z, a, b] = z^a/a · 2F1(a, 1-b; a+1; z)` through `Hypergeometric2F1`. The series terminates to a closed form when `b` is a positive integer, even for symbolic `z` and `a`. Otherwise it reduces only when every argument is numeric. Two special values are handled directly: `Beta[0, a, b] = 0` and `Beta[1, a, b] = Beta[a, b]`. The four-argument form `Beta[z0, z1, a, b]` evaluates as `Beta[z1, a, b] - Beta[z0, a, b]` when every argument is numeric.

**Data structures.** `Expr` trees driven through `eval_and_free`. The machine kernel works on `(re, im)` doubles in `sf_machine_gamma.c` and shares that file's Gamma and digamma tables. The ND kernel is the binary `NDKB_Beta`, registered with `REG_B`; its `ndk_Beta_c` calls `sfg_beta_c` for real and complex buffers, including the case where both arguments are arrays. So packed arrays and `NDArray`s are handled directly, and the scalar, packed, `NDArray` and `Compile[]` answers agree bit for bit. `Compile[]` lowers `Beta` through that kernel at real and complex scalar shapes, and at rank-1 shapes when one argument is an array and the other a scalar (`Beta[v, 2.5]`). An array in both arguments is not yet lowered.

**Complexity / limits.** The machine kernel is `O(1)`. Measured error: real median 0.8 ulp, max 4.1 ulp on `[0.3, 5]²`; complex median 1.3 ulp, max 6.3 ulp on `[0.5, 3] × ±2i`. The exact path costs what `Gamma`/`Pochhammer` cost, and the incomplete form costs one `2F1` evaluation. Symbolic two-argument calls, and the four-argument form with any symbolic argument, stay unevaluated; symbolic derivatives of the four-argument form are handled in `deriv.c`. Attributes: `Listable`, `NumericFunction`, `Protected`.
