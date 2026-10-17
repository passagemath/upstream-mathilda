# Special-function perf campaign — hypergeometric family

Heads: `Hypergeometric0F1`, `Hypergeometric1F1`, `Hypergeometric2F1`,
`HypergeometricPFQ`, `HypergeometricU`, `WhittakerM`, `WhittakerW`.
Branch `sf-perf/hyper`. No version bump or tag (the lead does those at merge).

## What changed

- **One machine kernel per head, on every surface.** New
  `src/special_functions/sf_machine_hyper.{c,h}` holds the real `double` and raw
  re/im complex kernels. The real 0F1/1F1/2F1/pFq definitions moved there from
  `sf_machine.c` and keep their names. The scalar builtins call these kernels
  for every machine-precision call, exact Integer/Rational parameters included,
  *before* the exact rewrites:
  - 0F1/1F1/2F1 no longer rebuild a `HypergeometricPFQ` expression and
    re-evaluate it.
  - `HypergeometricU` no longer evaluates `b - a - 1` on every call.
  - Whittaker no longer evaluates an expression template for a rational `m`.

  The kernels are bit-identical across five surfaces (scalar, unpacked List,
  packed List, NDArray, `Compile[]`), real and complex. This is asserted in
  `tests/test_sf_machine_hyper.c`.
- **Buffers.** The PFQ, U and Whittaker ND paths now:
  - accept exact parameters;
  - accept complex buffers;
  - widen to complex for points on the cut.

  A constant-parameter vector driver (`sf_machine_pfq_vec`) tabulates the 2F1
  series coefficients once per buffer, so the per-element loop has no division.
  When the parameters make the series one-signed and non-growing
  (`ab <= c`, `a + b <= c + 1`), it sums with Fast2Sum on two chains and runs
  four elements in lock-step using the scalar loop's exact arithmetic. The
  0F1 ND kernel gained a complex arm.
- **Compile.** New opcode `OP_KERNN_C`, the complex arm of the n-ary kernel. A
  complex operand of 1F1/2F1/PFQ/U/Whittaker now compiles
  (`CompileDiagnostics[{{x,_Complex}}, …]` gives `Compiled -> True`). The
  rank-1 array shape of the n-ary heads is still the pre-existing BASELINE gap.
- **Attributes.** `Hypergeometric0F1`, `1F1` and `2F1` are now `Listable`, as in
  Mathematica.
- **Semantics.**
  - `Hypergeometric2F1` is continued analytically to the whole plane. These
    used to stay unevaluated: real `z < -1`, complex `|z| >= 1`, and real
    `z > 1`.
  - On the cut `z > 1`, a real argument returns the complex limit from below,
    the Mathematica and mpmath convention: `Hypergeometric2F1[1,2,3,2.5]` gives
    `-0.929749 - 1.00531 I`.
  - At real `z < 0`, `WhittakerM`, `WhittakerW` and `HypergeometricU` return
    the principal complex value from the machine kernel.

## Algorithms (all in double / raw re-im double, cited in the source)

- **Every series measures its own loss.** The series is summed in double with
  compensated (TwoSum) summation while tracking the peak term and its index.
  - If peak/|sum| > 2, or the peak lies past the 6th term, the same series is
    re-summed in double-double (Dekker / Hida–Li–Bailey).
  - Beyond 2^45 of cancellation the kernel declines to MPFR.
  - Compensation alone took `2F1(1,2;3;0.89)` (380 terms) from 13 ulp to 1.4 ulp.
- **0F1.**
  - z > 0: the series; for z > 900, the Hankel expansion of I_ν, with
    2√z formed in double-double.
  - z < 0, little cancellation: the series in double.
  - z < 0, otherwise: the double-double re-sum up to |z| ≈ 250, then Miller's
    backward recurrence for J_ν written in −z itself, normalised by the Neumann
    series (Gautschi 1967).
- **1F1.**
  - Closed forms: a = 0, a = b, and (1, 2) → `expm1(z)/z` (with a complex form
    that does not cancel).
  - Terminating polynomials, valid even past an unreached lower pole.
  - |z| > 20: DLMF 13.7.2 with both exponential sectors (Stokes-line `cos`
    factor), accepted only when the least term is below half an ulp.
  - Re z < 0: Kummer's transformation, with b − a carried in double-double.
- **2F1.**
  - Real z: series on [0, 0.9]; Pfaff on [−9, 0); Gauss's sum at z = 1.
  - Complex z: the cheapest of the direct and Pfaff series, then the 1−z and
    1/(1−z) connection formulas (A&S 15.3.6), including the degenerate
    logarithmic cases (A&S 15.3.10–12). Their digamma differences use their own
    recurrences and exact finite sums for integer offsets. A connection result
    is accepted only when its condition number κ ≤ 2.
  - Elsewhere: ODE analytic continuation by Taylor steps (Bühring 1987) from
    z0 = ½ ± 0.45i. This covers the e^{±iπ/3} band, ill-conditioned parameters,
    and the cut, where the limit from below is taken.
  - Prefactors are computed as pow(|w|, s)·cis(s·arg) rather than
    exp(s·log w).
- **U (real).**
  - Asymptotic series, now compensated and accepted only below half an ulp
    (it used to accept 1e-13).
  - Integer-b log series and non-integer-b connection, each now accepted only
    when its combination amplifies rounding ≤ 4×. These guards used to allow
    1e6 and 1e13: `HypergeometricU[5.5, 0.3, 3.]` was 1.5e-10 off, and
    `U(1.3, 5.5, 10)` 60 ulp off.
  - The moderate-z band in between, which used to fall to FLINT at 280 µs: the
    integral representation by double-exponential (exp-sinh) quadrature, step
    1/8 with 50–90 nodes, ≈4 µs, ≤ 2.3 ulp.
- **U (complex).**
  - Asymptotic series.
  - Integer-b DLMF 13.2.9.
  - Kummer connection via complex 1F1.
- **Whittaker M/W.**
  - Prefactor × 1F1/U.
  - On the real axis the phase is reduced exactly, so an integer m + ½ gives a
    real value.

## Net µs (shootout, MR/C; Mathilda vs Python same run; Mathematica from 2026-10-17 baseline)

| head | MR before | MR after | Python | best other | C before | C after | Python | best other |
|---|---|---|---|---|---|---|---|---|
| 0F1 | 7.22 | 0.83 | 1.44 | 1.55 (py) | 0.95 | 1.15 | 1.67 | 2.21 (py) |
| 1F1 | 1.29 | 0.14 | 1.61 | 1.79 (py) | 4.67 | 0.31 | 1.72 | 1.87 (py) |
| 2F1 | 2.89 | 0.31 | 2.16 | 2.36 (py) | 4.81 | 1.74 | 2.89 | 3.24 (py) |
| PFQ | 1.04 | 0.31 | 27.9 (mpmath) | 25.3 (py) | 1.10 | 0.66 | 63.3 | 58.7 (py) |
| U | 0.86 | 0.11 | 1.60 | 1.48 (py) | 0.97 | 0.33 | 36.6 | 8.32 (wl) |
| WhittakerM | 4.74 | 0.14 | 50.5 (mpmath) | 10.0 (wl) | 7.84 | 0.34 | 76.8 | 10.8 (wl) |
| WhittakerW | 3.32 | 0.15 | 53.6 (mpmath) | 12.3 (wl) | 6.33 | 0.38 | 91.9 | 12.6 (wl) |

10^6-buffer throughput (ms; SciPy same run):

| head | before | after | SciPy |
|---|---|---|---|
| 0F1 | 10.1* | 87.5 | 259 |
| 1F1 | 1262 | 11.1 | 27.0 |
| 2F1 | 2596 | 136 | 162 |
| PFQ | 1314 | 186 | — |
| U | 26.4 | 7.4 | 36.8 |
| WhittakerM | 5495 | 16.0 | — |
| WhittakerW | 3850 | 17.9 | — |

\* The old 0F1 buffer path summed the alternating series in double and
accepted up to 1e13 of cancellation, so it was fast and wrong at z < 0: the
scalar path re-summed those points in MPFR, and the buffer disagreed with it.
It is now accurate and bit-identical to the scalar, and still 3× faster than
SciPy.

## Accuracy (tools/sf_accuracy/hyper.py, mpmath 30 digits at the double parameters)

Benchmark domains, 2000 points each, max (median) ulp:

| head | MR | C |
|---|---|---|
| 0F1 | 0.85 (0.25) | 1.24 (0.29) |
| 1F1 | 1.13 (0.29) | 3.5 (0.64) |
| 2F1 | 2.74 (0.43) | 3.25 (0.50) |
| 2F2 | 0.87 (0.25) | 0.80 (0.27) |
| U | 0.50 (0.26) | 2.14 (0.48) |
| WhittakerM | 1.26 (0.29) | 2.52 (0.49) |
| WhittakerW | 1.26 (0.29) | 2.52 (0.49) |

The old binary, same sweep: 1F1 C max 9.2, 2F1 MR 17.7, 2F1 C 36.5 (with 10 of
400 points unevaluated).

Wider bands are documented as looser. The largest maxima are below:

| case | median | max | max scaled |
|---|---|---|---|
| 2F1(1.5, 2.25; 3.1; x), x ∈ [−30, 0) | 1.6 | 15 | 8 |
| 2F1(½, ½; 1; x), x < −9 | 2.4 | 6 | — |
| 0F1, z < −250 (Miller) | — | — | 6.7 |
| 2F1 cut | — | 10 | 5 |
| WhittakerM(0.3, 0.7, x) | 1.3 | 10 | 2.3 |

The 2F1(1.5, 2.25; 3.1) case is the connection formula's κ ≈ 2 times the
Gamma-ratio rounding. The 2F1(½, ½; 1) case is digamma(½) rounding. For
WhittakerM, the derived a = m − k + ½ is rounded.
