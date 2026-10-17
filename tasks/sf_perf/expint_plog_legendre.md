# sf-perf / expint: ProductLog, LegendreP, LegendreQ — changelog prose

## ProductLog, LegendreP, LegendreQ: machine scalars on the shared double kernels

A machine-precision call (`ProductLog[x]`, `ProductLog[0, x]`, `LegendreP[n, x]`,
`LegendreQ[n, x]` with an integer-valued degree and a machine real or complex
argument) now goes straight to the raw-double kernels in the new
`src/special_functions/sf_machine_plog.c`. The NDArray element-wise kernels and
`Compile[]` call the same functions, so scalar == unpacked List == packed List ==
NDArray == Compile bit for bit, real and complex (tests/test_productlog_machine.c,
tests/test_legendre_machine.c). Before this change the scalar paths ran the
complex-MPFR Halley core (ProductLog) or built and evaluated an exact-coefficient
polynomial / closed-form expression (LegendreP/Q). A kernel that declines falls
through to those paths unchanged.

- **ProductLog, complex (new).** Principal branch W0 by Halley's iteration in
  re/im doubles (Corless, Gonnet, Hare, Jeffrey & Knuth 1996). The seed is the
  Padé form `z(3 + 6z + z²)/(3 + 9z + 5z²)` near the origin and the asymptotic
  `L1 - L2 + L2/L1` elsewhere. Every answer is verified twice: `w + Log w = Log z`
  (mod 2πi), and `w` must lie in W0's range (`|Im w| < π`, right of the curve
  `-v cot v + iv`). A real `x < -1/e` now returns its complex principal value
  through this kernel.
- **ProductLog near the branch point (real and complex).** For `|e z + 1| < 1/2`,
  Halley runs in `u = w + 1` on `g(u) - e(z + 1/e)`, with
  `g(u) = (u - 1)e^u + 1` summed from its Maclaurin series. `z + 1/e` uses `1/e`
  split into a double head and tail. The textbook residual `w e^w - z` subtracts
  two numbers close to `-1/e`, so near the branch point it is pure rounding.
  The split constant also puts the real/complex boundary in the right place: the
  double nearest `-1/e` is below it, so its W0 is `-1 + 8.22e-9 I`. The MPFR core
  always gave that answer; the old double kernel answered `-1.0`.
- **ProductLog hot path.** The exact `-1/E` and `-Pi/2` candidate expressions are
  now built only for a `Times[...]` argument. Before, they were built and compared
  on every call.
- **LegendreQ.** Uses Christoffel's form `Q_n = P_n Q_0 - W_{n-1}`, with `P` and
  `W` both computed by the three-term recurrence. The old recurrence on `Q` itself
  cancels in `Q_1 = x Q_0 - 1` and lost up to 9 ulp near |x| ≈ 0.8.
  `Q_0 = (Log[1+x] - Log[1-x])/2` is evaluated without cancellation (Kahan's
  log1p form for |x| < 1/2, a difference of exact-argument logs otherwise). It
  uses the closed form's branch, so a real |x| > 1 gives the same complex value
  as before (`LegendreQ[2, 1.5] = 0.063567 - 4.51604 I`). `x = ±1` declines
  (Indeterminate, unchanged).
- **LegendreP.** Complex recurrence added. Negative degree uses `P_{-n-1} = P_n`.
  A machine-integer degree (`3.`) also takes the kernel. `LegendreP[n, 1]` with an
  exact 1 stays exact.
- **ND/Compile.** `NDKU_ProductLog` gains a complex arm. `NDKB_LegendreP/Q` accept
  complex arguments; a real element whose value leaves the real axis (LegendreQ
  at a real |x| > 1) declines to the List path. `CompileDiagnostics` reports
  `Compiled -> True` for real scalar, real rank-1 and complex shapes of all three.

### Accuracy (mpmath 30 digits, 2000 points per regime; ulp, rel = relative, scl = scaled by max(|f|, |x f'|) near zeros)

| head | regime | rel med | rel max | scl med | scl max |
|---|---|---|---|---|---|
| ProductLog | MR [-0.3, 10] | 0.28 | 2.77 | 0.28 | 1.38 |
| ProductLog | C  | 0.34 | 1.60 | 0.34 | 1.60 |
| ProductLog | MRW [-1/e, 1e6] | 0.23 | 0.61 | 0.23 | 0.61 |
| ProductLog | CW ±30 | 0.30 | 1.16 | 0.30 | 1.16 |
| LegendreP[3,x] | MR | 0.54 | 318 (at the zero ±0.7746) | 0.28 | 1.66 |
| LegendreP[3,x] | C | 0.58 | 334 (at the zero) | 0.28 | 1.72 |
| LegendreP[3,x] | CW ±10 | 0.67 | 3.26 | 0.23 | 1.00 |
| LegendreQ[2,x] | MR | 0.31 | 4.44 | 0.24 | 1.24 |
| LegendreQ[2,x] | C | 0.38 | 2.46 | 0.30 | 2.02 |
| LegendreQ[2,x] | CW ±5 | 0.62 | 3.05 | 0.31 | 1.33 |

### Net latency (µs, run_shootout.py; Python figures from the same run)

| head | MR before → after | Python MR | C before → after | Python C | best other (baseline) |
|---|---|---|---|---|---|
| ProductLog | 27.7 → 0.30 | 2.48 | 51.7 → 0.75 | 2.32 | MR 1.84, C 1.77 |
| LegendreP[3,x] | 7.22 → 0.135 | 1.27 | 10.1 → 0.285 | 1.44 | MR 0.818, C 1.17 |
| LegendreQ[2,x] | 12.5 → 0.156 | 3558 (mpmath) | 23.3 → 0.363 | 8159 (mpmath) | MR 4.56, C (winning) |

Vectorised 10^6 buffer (ms): ProductLog 21.7 vs SciPy 300; LegendreP 4.0 vs 7.7;
LegendreQ 6.0 (no SciPy equivalent).
