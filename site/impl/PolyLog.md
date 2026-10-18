---
source: src/special_functions/polylog.c
references:
  - "DLMF §25.12 — the polylogarithm and Jonquiere's function (§25.12.11/12 for the zeta expansion, §25.12.4 and §25.12.13 for the inversion formulas)."
  - "G. 't Hooft and M. Veltman, Scalar one-loop integrals, Nucl. Phys. B153 (1979) 365-401 — the Bernoulli series for Li_2 in u = -log(1-z)."
  - "R. E. Crandall, Note on fast polylogarithm computation (2006) — the dilogarithm series and the integer-order inversion formula."
---
**Algorithm.** `builtin_polylog` evaluates `PolyLog[n, z] = Li_n(z)`. (1)
**Machine-scalar gate** (any `z` except the exact `0`): `sf_machine_scalar_args`
routes a machine-precision call to the `double` kernels of `sf_machine_zeta.c` —
`sf_machine_polylog` (= `sf_polylog_r`) for a real order and real `z ∈ [−1, 1]`
(every real `z` for the dilogarithm), `sf_polylog_c` otherwise. `PolyLog[2, z]`
uses `sf_dilog_r` / `sf_dilog_c`: the Bernoulli series in `u = −log(1−z)` after
the maps `z → 1/z` (`|z| > 1`) and `z → 1−z` (`Re z > 1/2`), so `|u| ≤ 1.25`.
Other real orders: closed forms at `s = 0, 1`, the power series for `|z| ≤ 1/2`,
the duplication formula for real `x ∈ [−1, −1/2)` with `s ≥ 1`, the Jonquière
expansion in `μ = log z` for `|μ| < 3.6`, and for integer `s ≥ 2` and `|z| > 1`
the Bernoulli-polynomial inversion formula. A real `z > 1` answers from below the
cut (`Im z = −0`); a real `z < −1` keeps only the (real) real part. A complex
order declines. (2) **Exact closed forms**: `Li_n(0) = 0`; `Li_1(z) = −Log[1−z]`;
`Li_0(z) = z/(1−z)`; `Li_{−m}(z)` the Eulerian-number rational function (cap
`POLYLOG_NEGINT_CAP = 400`); for integer `n ≥ 2`, `Li_n(1) = Zeta[n]`,
`Li_n(−1) = (2^(1−n)−1) Zeta[n]`, `Li_2(1/2)`, `Li_3(1/2)`, and `Li_2(±I) =
−π²/48 ± I Catalan`. (3) **Numeric** (`polylog_numeric`, at least one inexact
operand): real order with real `−1 < z < 1` uses arb's `acb_polylog`
(`flint_num_polylog`) above machine precision, else the real-MPFR series
`pl_real_series`; everything else runs the complex-MPFR kernel — power series
for `|z| ≤ 1/2`, the Jonquière/zeta expansion (`ζ(s−k)` and `Γ(1−s)`, reflected
through the functional equation in the left half-plane) for `|ln z| < 2π − 0.05`,
and the `Li_2` inversion `pl_li2_inversion` for `|z| > 1`; otherwise symbolic.
`PolyLog[n, p, z]` (Nielsen) is accepted and left symbolic.

**Data structures.** `Expr`; the local `pcx` (`mpfr_t` re/im) toolkit; real paths
use `mpfr_zeta`/`mpfr_gamma`, complex-order paths reuse the `Zeta`/`Gamma`
builtins. The machine kernel carries `static const` tables of ζ at the integers
(`ZM_ZPOS`, `ZM_ZNEG`), Bernoulli coefficients (`ZM_DILOG_C`, `ZM_BN`) and π²/6
as a two-double constant. ND: a hand-written binary kernel `NDKB_PolyLog` with a
complex arm (`ndk_PolyLog_c` → `sf_polylog_c`), registered `REG_B`, so packed
and visible-`NDArray` arguments, real or complex, stay on the buffer.
Attributes: `Listable`, `NumericFunction`, `Protected`.

**Complexity / limits.** Machine values ≈ 0.15 µs real, ≈ 0.4 µs complex; a
10⁶-element `PolyLog[2, v]` takes a few ms. Each MPFR regime is gated by
`|z|`/`|ln z|`; outside them (e.g. non-integer order with large `|z|`, or an
order the machine kernel declines) the call stays symbolic. `Compile[]` lowers
`PolyLog[n, z]` at scalar and rank-1 array shapes, real and complex
(`Compiled -> True`).
