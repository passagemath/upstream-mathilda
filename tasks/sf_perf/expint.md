# sf-perf / expint — changelog prose (for the lead to merge)

Family "expint": SinIntegral, CosIntegral, SinhIntegral, CoshIntegral,
ExpIntegralEi, ExpIntegralE, LogIntegral, EllipticK, EllipticE (complete),
ProductLog, LegendreP, LegendreQ. Three parts: this file covers the
exponential-integral family; `expint_elliptic.md` and `expint_plog_legendre.md`
(same directory) cover the other two sub-branches, merged into
`sf-perf/expint`.

## Exponential-integral family: one double kernel per head, every surface

A machine-precision `SinIntegral`, `CosIntegral`, `SinhIntegral`, `CoshIntegral`,
`ExpIntegralEi`, `LogIntegral` or `ExpIntegralE[n, .]` call, real or complex, now
goes straight to the raw-double kernels in the new
`src/special_functions/sf_machine_expint.c`. The shared gate
`sf_machine_scalar_args` sits at the top of each builtin, ahead of the
special-value ladder. The NDArray element-wise kernels and `Compile[]` call the
same functions, so scalar == unpacked List == packed List == NDArray == Compile
bit for bit, real and complex (`tests/test_expint_machine.c`).

Before this change:
- The scalar builtins computed in MPFR (or ncpx complex MPFR pairs) and rounded.
- `LogIntegral` built `ExpIntegralEi[Log[z]]` and evaluated it.
- The existing double kernels in `sf_machine.c` were used only by NDArray and
  Compile. They were not accurate enough to wire in: up to 10 ulp for Si, 648 ulp
  for Ci near its zero, 24 ulp for E1 at x ≈ 1, and 54 ulp for li(x < 1).
- The complex `*_machine_complex` series declined wherever the terms cancelled,
  and were only about 1e-3 × eps accurate elsewhere.

A kernel that declines falls through to the MPFR path unchanged:
- poles: Ei(0), li(1), Ci(0);
- overflow and underflow: Ei(±800) and Shi(800) keep their extended-exponent
  MPFR answers;
- non-finite input.

### Real kernels

All coefficients are refitted from mpmath at 50 digits by
`tools/sf_accuracy/fit_expint.py` (near-minimax relative-error fits:
Sanathanan–Koerner + Lawson, `tools/sf_accuracy/fitlib.py`). The structures
follow the classical double designs; no external tables are copied.

- **Si, Ci** (structure of Cephes `sici`):
  - Si = x·S(x²) for |x| ≤ 4.
  - Ci = γ + log x + x²·C(x²) for x ≤ 2.
  - On [2, 4], Ci = (x − x₀)·Q(x) through its first zero x₀ = 3.38418. The series
    form cancels 15-fold by x = 4 and measured 5 ulp there.
  - For x > 4, the auxiliary f, g of DLMF 6.2.17 are written x f = 1 − w RF(w) and
    x² g = 1 − w RG(w), with w = 1/x². RF and RG are Stieltjes-type rationals on
    [4, 8] and [8, ∞), with one-signed coefficients.
- **Shi, Chi** (structure of Cephes `shichi`):
  - Shi = x·Sh(x²) and Chi = γ + log x + x²·Ch(x²) for x ≤ 4. The polynomials are
    fitted in v = x²/16, so every coefficient is positive and Horner never cancels.
  - For x > 4, (Ei(x) ± E1(x))/2.
- **Ei, E1** (Cody & Thacher, Math. Comp. 1968/1969 structure):
  - For 0 < x ≤ 16, Ei = γ + log x + x·R(x), with R on six pieces. Each piece is a
    positive-coefficient polynomial fitted from its LEFT endpoint. Rationals in
    1/x were tried and rejected there: they evaluate with 1–2e-15 error because
    x e^−x Ei(x) is not a Stieltjes function.
  - For 16 < x ≤ 64, Ei = e^x/x·(1 + u H(u)).
  - Beyond 64, the asymptotic series.
  - E1 = x·R(−x) − γ − log x for x ≤ 0.6, then e^−x/x·(1 − u J(u)) rationals.
- **li(x) = Ei(log x).** The rounding of log x costs about log x ulps (9 at
  x = 1e4) through Ei's e^y, so it is recovered as dy = (x − e^y)/x and added back
  as (x/y)·dy. In the large-argument forms, x itself replaces e^(log x).
- **E_n(x)**, integer n ≥ 0:
  - For x ≤ 1, the upward recurrence E_(k+1) = (e^−x − x E_k)/k from E1. The power
    series' first three terms cancel 5–11-fold (4.75 ulp).
  - For x > 1, the continued fraction (DLMF 8.19.17) evaluated BACKWARD from depth
    105/x + 12, calibrated against mpmath for 1 ≤ n ≤ 50. Forward Lentz
    accumulates one rounding per level. Forward Wallis (no divisions) cancels:
    8e-15 at x = 2.4.

### Complex kernels (raw re/im doubles, no `<complex.h>`)

- **Power series by Horner, with tabulated coefficients.** The term recurrence
  measured 5 ulp; Horner measures 2. The series is used only where its terms do
  not cancel:
  - |z| − Re z ≤ 2 for Ei;
  - |z| − |Im z| ≤ 2 or |z| ≤ 2.5 for Si, Shi, Chi;
  - ≤ 1 for Ci, whose log cancels against the sum near the real axis (measured
    7 ulp with the Si thresholds).
- **Elsewhere, E1 by its continued fraction (DLMF 6.9.1), evaluated backward from
  depth 240/s + 8, s = |w| + Re w.** Forward Lentz measured 55 ulp at depth 60.
  The other functions are derived from E1:
  - Si = π/2 + (E1(iz) − E1(−iz))/(2i) and Ci = −(E1(iz) + E1(−iz))/2
    (DLMF 6.5.5–6);
  - Ei(z) = −E1(−z) ± iπ.
- **|z| > 45:** the asymptotic expansions.
- **Rotations and compositions:**
  - Shi(z) = −i Si(iz).
  - Chi(z) = Ci(iz) + (log z − log iz), with the constant taken exactly.
  - li(z) = Ei(log z), with the same log-rounding correction as the real kernel.
- **ExpIntegralE** with a machine **complex** argument now evaluates; it was
  unevaluated before. A real x < 0 is complex-valued, on the upper side of the
  cut: `ExpIntegralE[2, -1.] = 0.823164 - 3.14159 I`.
  - `ExpIntegralE[n, packed]` gets a binary NDArray kernel; it used to thread
    element by element.
  - `CompileDiagnostics` reports `Compiled -> True` for the `_Real`,
    `{_Real, 1}` and `_Complex` shapes.

### Accuracy

Measured with `tools/sf_accuracy/expint.py`: mpmath at 30 digits, 2000 points
per regime. MR/C are the benchmark domains; MRW/CW are wider bands. "scl" is
the error in ulps of max(|f|, |x f′|), which differs from the relative error
only near a zero.

| head | MR rel med / max (scl max) | C rel med / max (scl max) | MRW scl max | CW scl max |
|---|---|---|---|---|
| SinIntegral | 0.30 / 1.73 (1.73) | 0.46 / 24.9 (1.54)¹ | 1.79 | 1.38 |
| CosIntegral | 0.45 / 256¹ (3.50) | 0.69 / 5.45¹ (1.82) | 4.24 | 0.87 |
| SinhIntegral | 0.43 / 2.34 (1.21) | 0.48 / 2.40 (1.95) | 1.14 | 1.32 |
| CoshIntegral | 0.45 / 47¹ (1.64) | 0.51 / 2.55 (1.42) | 1.36 | 1.85 |
| ExpIntegralEi | 0.41 / 61¹ (1.60) | 0.56 / 2.53 (1.72) | 1.70 | 1.23 |
| ExpIntegralE[2,·] | 0.43 / 4.26 (3.02) | 0.68 / 3.82 (2.05) | 3.46 | 2.25 |
| LogIntegral | 0.42 / 1.89 (1.86) | 0.54 / 2.71 (2.71) | 2.67 | 3.33 |

¹ The relative maximum sits at a zero of the function (Ci at 3.38/12.6, Chi at
0.5238, Ei at 0.3725, Si's complex zeros), where the relative error is
unbounded. The scaled column is the meaningful one there.

### Net latency

Net µs per call from `run_shootout.py`. The before figures and the best-other
engine come from the baseline JSON; Python is from the same run.

| head | MR before → after | Python | best other | C before → after | Python | best other |
|---|---|---|---|---|---|---|
| SinIntegral | 7.34 → 0.158 | 1.51 | 1.23 | 30.9 → 1.21 | 145 | 25.7 |
| CosIntegral | 9.99 → 0.172 | 1.53 | 1.26 | 39.3 → 1.88 | 159 | 25.9 |
| SinhIntegral | 5.12 → 0.143 | 1.68 | 1.29 | 25.3 → 0.334 | 261 | 38.5 |
| CoshIntegral | 7.73 → 0.148 | 1.62 | 1.31 | 35.7 → 0.367 | 254 | 42.1 |
| ExpIntegralEi | 14.3 → 0.147 | 0.308 | 0.260 | 64.1 → 0.451 | 109 | 18.6 |
| ExpIntegralE[2,·] | 0.246 → 0.315 | 1.41 | 1.19 | (n/a, was unevaluated) | | |
| LogIntegral | 12.4 → 0.164 | 0.369 | 0.321 | 54.7 → 0.477 | 130 | 9.38 |

10⁶-element buffer, Mathilda (the packed NDArray path) vs SciPy:

| head | before | after | SciPy |
|---|---|---|---|
| SinIntegral | 76.5 | 5.94 | 32.5 |
| CosIntegral | 73.0 | 8.14 | 32.8 |
| SinhIntegral | 33.2 | 4.84 | 123 |
| CoshIntegral | 32.1 | 5.33 | 120 |
| ExpIntegralEi | 11.8 | 5.31 | 177 |
| ExpIntegralE[2,·] | 557 | 22.3 | 255 |
| LogIntegral | 12.2 | 6.95 | — |

All figures are in ms. ExpIntegralE was threaded element by element before.
Single-threaded, its kernel costs about 0.22 µs per element, the continued
fraction near x = 1 being the cost; SciPy's `expn` costs 0.25 µs.

### Notes

- `ExpIntegralE[2, x]` scalar is 0.07 µs slower than before: 0.315 µs against
  0.246. The old continued fraction stopped at a 2e-15 tolerance. It is still
  four times faster than the best other engine.
- `tests/test_compile.c`'s "expint cancellation gate" used to require the
  complex kernels to DECLINE at 40 + 30i, where the old series lost every digit.
  They now answer accurately there, through the asymptotic expansion. The test
  now requires an answer, and checks it against mpmath.
- Gates: `make check-c99`, `check-messages`, `check-packed-aware`,
  `check-nd-surfaces`, `check-array-exactness` and `check-compile-coverage` all
  pass.
