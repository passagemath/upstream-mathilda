---
source: src/special_functions/expintegralei.c
references:
  - "DLMF §6.2.5, §6.6.2 — the exponential integral Ei(z) and its convergent series γ + Log(z) + Sum z^k/(k k!)."
  - "DLMF §6.9.1 — the continued fraction for E1; DLMF §6.12.2 — the asymptotic expansion Ei(z) ~ (e^z/z) Sum k!/z^k."
  - "W. J. Cody and H. C. Thacher, *Rational Chebyshev approximations for the exponential integral E1(x)*, Math. Comp. 22 (1968) 641–649."
  - "W. J. Cody and H. C. Thacher, *Chebyshev approximations for the exponential integral Ei(x)*, Math. Comp. 23 (1969) 289–303."
---
**Algorithm.** `builtin_expintegralei` evaluates `Ei(z)`, which has a branch cut along `(-∞, 0)`; on the cut it returns the principal value. A machine call first passes the `sf_machine_scalar_args` gate and goes to the double kernels in `sf_machine_expint.c`.

**Real machine kernel** (`sf_machine_ei`), following Cody & Thacher's scheme with coefficients refitted from mpmath by `tools/sf_accuracy/fit_expint.py`:

- `0 < x <= 16`: `γ + log x + x·R(x)`, where `R` is split into six pieces, each a polynomial with positive coefficients (so Horner never cancels).
- `16 < x <= 64`: `e^x/x·(1 + u·H(u))` with `u = 1/x`, using rationals on two octaves.
- `x > 64`: the asymptotic series.
- `x < 0`: `-E1(-x)`, where `E1` is a polynomial for `x <= 0.6` and `e^{-x}/x·(1 - u·J(u))` beyond.

`e^x` is formed as `e^{x/2}·e^{x/2}`, so it does not overflow early.

**Complex machine kernel** (`sf_machine_ei_c`):

- The power series by Horner wherever its terms do not cancel (`|z| - Re z <= 2`).
- Otherwise `-E1(-z) ± iπ`, with `E1` from the continued fraction (DLMF 6.9.1). The fraction is evaluated backward from a fixed depth `n = 240/s + 8`, where `s = |w| + Re w`; forward Lentz evaluation measured 55 ulp.
- The asymptotic series for `|z| > 45`.

A decline (the pole at 0, overflow, underflow) falls through to `ei_one_arg`:

- **Exact special values:** `Ei[0] = -Infinity`, `Infinity → Infinity`, `-Infinity → 0`, `±I Infinity → ±I π`, and `ComplexInfinity`/`Indeterminate → Indeterminate`.
- **Real, `x > 0`, any precision:** `mpfr_eint`, which is correctly rounded. This also gives the extended-exponent answer for `ExpIntegralEi[800.]`.
- **Real, `x < 0`:** the series `γ + ln|x| + Σ x^k/(k·k!)` with `|x|/ln2` guard bits, or the asymptotic form when `|x|` is large.
- **Complex:** the same series with the principal `Log z` in the `ecx` toolkit, using `(|z| + |Re z|)/ln2` guard bits, or for large `|z|` the asymptotic form plus the `iπ sign(Im z)` jump.

An explicit low-precision request (an MPFR leaf below 53 bits) gets an MPFR result, so the requested digits print. Every result passes through `numeric_promote_result_if_degenerate`.

**Data structures.** `Expr` trees and the `ecx` toolkit (pairs of `mpfr_t`) for the MPFR paths. The machine side keeps coefficient tables in `sf_machine_expint_coef.h` and does complex arithmetic on explicit `(re, im)` doubles. The ND kernel is `NDKU_ExpIntegralEi = { expintegralei_machine_complex, ndk_ExpIntegralEi_r }`, registered with `REG_U`, so it is packed-aware. Its real arm declines the pole, and its complex arm is `sf_machine_ei_c`. `Compile[]` lowers `ExpIntegralEi` through this kernel for real and complex scalars and for rank-1 arrays. The scalar call, packed arrays, `NDArray` and `Compile[]` all use the same helpers and agree bit for bit.

**Complexity / limits.** The machine kernels are `O(1)`. Against mpmath they show a median error of 0.4 ulp (real) and 0.6 ulp (complex), with a maximum of 2.6 ulp away from the zero `x₀ = 0.3725`. The MPFR series costs about `2|z|` terms plus guard bits until the asymptotic regime takes over. Exact non-zero arguments (`ExpIntegralEi[2]`) stay unevaluated. Attributes: `Listable`, `NumericFunction`, `Protected`.
