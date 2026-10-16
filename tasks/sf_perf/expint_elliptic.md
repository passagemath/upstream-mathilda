# EllipticK / EllipticE (complete) — machine kernels (sf-perf/expint-ell)

## What changed
- New `src/special_functions/sf_machine_ellk.{c,h}`: one machine kernel per
  (function, real | complex) for the complete integrals `K(m)`, `E(m)`. The
  scalar builtins, the packed/NDArray element-wise path and `Compile[]` all call
  it, so scalar == List == packed == NDArray == Compile bit for bit (tested).
- **Real `m < 1`**: Cephes' form (Moshier, `ellpk`/`ellpe`; A&S 17.3.34/17.3.36)
  `K = P(t) − log(t) Q(t)`, `E = P(t) − t log(t) Q(t)`, `t = 1 − m`, degree-10
  polynomials re-fitted by `tools/sf_accuracy/fit_ellk.py` (linear relative-error
  minimax, Lawson IRLS, mpmath 50 digits, points down to `t = 1e-300`; fit error
  ≤ 1.2e-18). `m < 0` by DLMF 19.7.5 (same polynomials at `t' = 1/t`).
  `elliptic_machine_k` / `elliptic_machine_e_complete` (used by the incomplete
  integrals' quasi-period too) now delegate here; Carlson `R_F`/`R_D`/`R_J`
  remain for `EllipticF`, `EllipticE[φ, m]`, `EllipticPi`.
- **Complex `m`**: complex AGM in raw double re/im (no `<complex.h>`),
  `K = π/(2 M(1, √(1−m)))`, `E = K(1 − Σ 2ⁿ⁻¹cₙ²)` (DLMF 19.8.5–6), right-choice
  square roots (Cox 1984) → principal branch, cut `[1, ∞)`, value on the cut =
  limit from below (`EllipticK[2.] = 1.31103 − 1.31103 I`, matching Mathematica,
  mpmath and the previous Arb answer). Final `π/(2a)` by Smith's division.
- **Wiring**: `EllipticK[m]` / `EllipticE[m]` take the shared gate
  `sf_machine_scalar_args` right after the pole check; a machine complex `m` and
  a machine real `m > 1` return a **machine** `Complex` (previously Arb returned
  53-bit `EXPR_MPFR` — `Precision[EllipticK[0.6 + 0.5 I]]` was 15.95, now
  `MachinePrecision`). Declines (pole `K(1)`, `|m| > 1e300`, non-finite) fall
  through unchanged. `NDKU_EllipticK/E` gained complex arms, so complex NDArrays
  and `Compile[{{z, _Complex}}, EllipticK[z]]` run on the kernel.

## Accuracy (tools/sf_accuracy/expint.py, 2000 pts/regime, vs mpmath 30 digits)
| head | MR med/max | C med/max | MRW (−50..1) max | CW (±20) med/max |
|---|---|---|---|---|
| EllipticK | 0.36 / 1.73 | 0.62 / 3.42 | 2.42 | 0.72 / 3.55 |
| EllipticE | 0.32 / 1.27 | 0.78 / 3.89 | 2.34 | 1.05 / 4.72 |

## Speed (net µs per call, run_shootout.py, same session, Python = scipy (MR) / mpmath (C))
| head | MR before → after | Python | C before → after | Python |
|---|---|---|---|---|
| EllipticK | 0.193 → 0.123 | 0.148 | 6.74 → 0.371 | 232 |
| EllipticE | 0.236 → 0.121 | 0.145 | 19.9 → 0.356 | 696 |

10⁶-element buffer: `EllipticK` 12.5 → 4.0 ms, `EllipticE` 18.5 → 3.3 ms
(SciPy 16.7 / 16.5 ms in the same run).
