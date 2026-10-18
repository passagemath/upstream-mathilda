---
source:
  - src/special_functions/struveh.c
  - src/special_functions/sf_machine.c
  - src/flint_num_bridge.c
  - src/ndkernels.c
references:
  - "DLMF §11.2.1 — the power series H_ν(z) = (z/2)^(ν+1) Σ (−1)^k (z/2)^(2k)/(Γ(k+3/2) Γ(k+ν+3/2))."
  - "DLMF §11.6.1 — the large-argument asymptotic expansion of K_ν = H_ν − Y_ν; DLMF §10.17.4 — the Hankel asymptotic expansion of Y_ν."
---
**Algorithm.** `builtin_struveh` first returns the exact `0` when z is the integer 0 and ν is real with ν > −1. It then numericizes only when some argument is inexact and both are numeric quantities. The output precision is the smallest inexact bit count among the arguments.

- **Above 53 bits**, `flint_num_struveh` evaluates `(z/2)^(ν+1) · acb_hypgeom_pfq(1; 3/2, ν+3/2; −z²/4, regularized)` in arb. The regularized 1F2 absorbs the `1/Γ(ν+3/2)` pole, so every complex ν works. The starting working precision includes a cancellation allowance of about 1.45·|z| + 8 bits. The `nb_should_retry` / `nb_emit_bits` ladder doubles the precision until the ball radius certifies the requested digits.
- **At machine precision** with real ν and z, `sf_machine_struve_h` runs two regimes. Each returns a value together with a predicted relative error:
  - `struve_series`, the convergent 1F2 series. Its error estimate is the measured cancellation `max|term|/|sum| · ε`.
  - `struve_asymp`, which computes `Y_ν + K_ν` from two optimally truncated divergent series. Its error estimate is the size of the smallest omitted term.

  The kernel keeps the regime with the smaller predicted error and declines when even that one exceeds `STRUVE_ACCEPT = 1e-7`. A declined real z is retried through an expression template built on `HypergeometricPFQ`.

**Data structures.** The machine path works only on scalar doubles and keeps nothing between calls. The arb path uses `acb_t` balls for ν, z and the result, plus `acb_struct` parameter arrays for the 1F2. In `src/ndkernels.c`, `NDK_BIN2(StruveH, sf_machine_struve_h)` and `REG_B` register the same kernel for packed lists and `NDArray`, and `Compile[]` lowers it at scalar and rank-1 shapes. Every surface calls the same function, so the results are bit-identical.

**Complexity / limits.** The series costs O(|z|) terms, with a cap of 8000. The asymptotic regime stops at its least term, at most 400 terms. The two error curves cross near z ≈ 12–24 for small ν, and there the double result is only good to about 1e-9 absolute. There is no complex machine kernel, so a machine complex z stays unevaluated; an exact complex z under `N` works through arb. There is no `D` rule. One defect: at a non-positive half-integer ν (where Γ(ν+3/2) has a pole), the machine fallback template returns an unevaluated `0.0 HypergeometricPFQ[...]` instead of a number. The arb path gets these values right.
