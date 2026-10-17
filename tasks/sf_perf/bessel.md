# Bessel / Airy family: machine-precision kernels (changelog prose)

## Machine-precision BesselJ / BesselY / BesselI / BesselK and the Airy functions

Every machine-number call of `BesselJ`, `BesselY`, `BesselI`, `BesselK` (with a
real order, machine or exact), `AiryAi`, `AiryBi`, `AiryAiPrime` and
`AiryBiPrime` — real or complex argument — is now answered by raw-double
kernels in the new `src/special_functions/sf_machine_bessel.c`, instead of the
complex-MPFR cores rounded to a double. One dispatch (`sf_bessel_nd` /
`sf_bessel_scalar`, `sf_airy_nd` / `sf_airy_scalar`) serves the scalar
builtins, the NDArray / packed element-wise kernels and `Compile[]`, so the
scalar, List, packed, `NDArray` and compiled answers agree bit for bit (real and
complex; `tests/test_sf_machine_bessel.c` checks all five). Arbitrary precision
is untouched: a kernel decline (order beyond 1000, overflow / underflow of the
double range, a complex order) falls through to the MPFR / arb path unchanged.

**Real axis.**
- Integer-order J and Y: libm `jn` / `yn` (≤ 3.6 ulp of the local magnitude
  over the benchmark domain and a wide band).
- I₀, I₁, K₀, K₁: Cephes-style regions with regenerated fits — the ascending
  series for small x (positive terms), piecewise Chebyshev fits of the
  exponentially scaled function in the middle, Hankel's expansion beyond 32;
  K₀/K₁ use the logarithmic series below x = 1/2. Kₙ by upward recurrence.
  ≤ 2.2 ulp (were 13 / 31 ulp in the old double kernels). Tables:
  `tools/sf_accuracy/gen_bessel01_tables.py`.
- Airy: piecewise Chebyshev polynomials on [−10, 10] (the functions on
  [−10, 2], the scaled Ai e^ζ, Bi e^−ζ on [2, 10]) and the Poincaré expansions
  beyond, with ζ = (2/3)x^{3/2} in double-double (a plain double ζ costs
  hundreds of ulps at large |x|). ≤ 3.7 ulp (the old kernel reached 270 ulp at
  x ≈ 2.3 and ~1000 ulp at x ≈ 100). Tables:
  `tools/sf_accuracy/gen_airy_tables.py`.

**Complex arguments and non-integer orders** (Amos's scheme, ACM TOMS 644):
J and Y through I and K of u = −iz in the right half plane, I/K continued to
the left half plane (DLMF 10.34), negative orders by the J/Y and I/K
reflections. K from Temme's series (|z| ≤ 1.25), from Taylor continuation of
the modified Bessel equation off tabulated nodes for 1.25 ≤ |z| < 21 (orders 0
and ±1/3 mod 1; `tools/sf_accuracy/gen_bessel_ktab.py`), Steed's CF2 for other
fractional orders, Hankel's expansion beyond; I from Miller's backward
recurrence (from the fractional part of the order, with a double-double 2/z),
normalised by whichever Gegenbauer / Neumann sum is better conditioned. Airy:
the Maclaurin series for |z| ≤ 1, else DLMF 9.6.1–9.6.3 via one core call at
order −1/3, with the continuation sector (Re ζ < 0) in a simplified form where
no two large terms cancel, and a double-double ζ.

**The `bessel.m` DownValues are now C.** The half-integer-order closed forms
(symbolic argument only) and the negative-integer-order reflections moved from
`src/internal/bessel.m` (deleted) into the builtins, run after the numeric
paths decline. As DownValues they were matched — six conditional patterns —
before the builtin on every numeric call. Symbolic results are unchanged
(`BesselJ[5/2, x]`, `BesselK[-3/2, x]`, `BesselJ[-3, x]`, …).
`tools/check_install.sh` now probes `HeavisideTheta[0]` (from
`distributions.m`) to detect a missing module tree.

**Fix.** `AiryAi[120.]` (and `AiryAiPrime` at large x) returned `0.`: the
53-bit evaluation flushed the value below the double range. It now promotes
and re-evaluates with an extended exponent: `2.1575×10^-382`.

### Accuracy (`tools/sf_accuracy/bessel.py`, 2000 points per head and regime, vs mpmath at 30 digits)

Ulps of the local magnitude (absolute error scaled by max(|f|, h|f'|), h the
local length scale; equals the relative error away from zeros).

| head | MR median / max | C median / max | MR wide | C wide |
|---|---|---|---|---|
| BesselJ[2, ·] | 0.42 / 3.6 | 0.92 / 7.5 | 0.05 / 3.7 | 0.97 / 7.9 |
| BesselY[1, ·] | 0.33 / 2.7 | 0.96 / 5.2 | 0.06 / 2.9 | 0.99 / 6.8 |
| BesselI[0, ·] | 0.38 / 2.3 | 0.80 / 4.1 | 0.52 / 2.7 | 1.18 / 12.2 |
| BesselK[0, ·] | 0.32 / 1.9 | 0.31 / 3.9 | 0.16 / 2.8 | 0.92 / 9.5 |
| AiryAi | 0.28 / 2.7 | 0.97 / 9.1 | 0.51 / 3.4 | 1.27 / 8.4 |
| AiryBi | 0.24 / 2.8 | 1.13 / 4.8 | 0.51 / 3.7 | 1.32 / 7.1 |
| AiryAiPrime | 0.25 / 2.3 | 1.04 / 6.7 | 0.48 / 3.4 | 1.23 / 6.8 |
| AiryBiPrime | 0.25 / 2.7 | 1.07 / 4.6 | 0.52 / 4.1 | 1.16 / 10.2 |

Domains: the benchmark's (`registry.py`); wide = J/Y x ∈ (10⁻³, 80), I/K up to
600, Airy [−60, 100]; complex |Re|, |Im| up to 15–30.

Non-integer real orders (random ν ∈ [−4.7, 9.3]), not in the benchmark: median
0.2–1.7, max ≤ ~18 ulp (the wider band). For comparison SciPy (Amos / Cephes)
on the same complex domains: J₂ max 21, Y₁ 26, K₀ 17, Airy 190–830 ulp.

### Speed (net µs per call, `run_shootout.py`, Mathilda vs SciPy in the same run)

| head | MR before → after | MR SciPy | C before → after | C SciPy | best other (baseline) MR / C |
|---|---|---|---|---|---|
| BesselJ | 14.4 → 0.21 | 1.83 | 88.1 → 0.64 | 1.74 | 2.10 / 1.91 |
| BesselY | 51.4 → 0.18 | 2.38 | 158 → 0.79 | 2.98 | 2.95 / 2.83 |
| BesselI | 30.5 → 0.14 | 1.56 | 51.0 → 0.47 | 1.63 | 1.43 / 1.50 |
| BesselK | 58.2 → 0.14 | 1.60 | 109 → 0.40 | 1.78 | 1.13 (Mathematica) / 1.65 |
| AiryAi | 88.8 → 0.12 | 1.91 | 130 → 0.69 | 5.47 | 0.234 (Mathematica) / 4.8 |
| AiryBi | 85.0 → 0.12 | 1.87 | 126 → 0.83 | 5.29 | 0.250 (Mathematica) / 4.7 |
| AiryAiPrime | 87.3 → 0.12 | 1.88 | 130 → 0.70 | 5.24 | 0.227 (Mathematica) / 4.8 |
| AiryBiPrime | 85.2 → 0.12 | 1.82 | 126 → 0.88 | 5.09 | 0.246 (Mathematica) / 4.7 |

10⁶-element machine-real buffer (ms): BesselJ 15.7 (SciPy 683), BesselY 9.1
(1037), BesselI 5.3 (382), BesselK 6.2 (344), Airy ×4 3.4–3.7 (137–159).
