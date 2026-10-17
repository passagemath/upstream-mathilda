# Special-function performance campaign: error-function family

Branch `sf-perf/error`. Heads: `Erf`, `Erfc`, `Erfi`, `FresnelC`, `FresnelS`,
`Sinc`, `InverseErf`, `InverseErfc`. Changelog prose for the lead to merge into
the weekly `docs/spec/changelog/` file.

## What changed

All machine-precision evaluation of the family now runs on raw-double kernels in
the new `src/special_functions/sf_machine_error.{c,h}` (plus the generated
`sf_machine_error_tables.h`). One helper per head and argument type serves every
surface: the scalar builtin (behind the shared `sf_machine_scalar_args` gate,
ahead of each special-value ladder), the packed / `NDArray` element-wise map
(`src/ndkernels.c`, now with complex arms) and `Compile[]` (scalar and rank-1,
real and complex: `Erf`/`Erfc` complex lower through the registered kernel, and
`OP_ERF_R`/`OP_ERFC_R`/`OP_VERF_R`/`OP_VERFC_R` call the shared kernels instead
of libm). Scalar == unpacked List == packed List == NDArray == Compile holds bit
for bit for every head, real and complex (`tests/test_sf_machine_error.c`).
Previously the scalar paths ran MPFR (real `Erfi`, `Fresnel*`, `Sinc`; every
complex call) and then rounded.

Algorithms (all in raw double re/im, no `<complex.h>`; sources cited in the code):

- **Complex erf / erfc / erfi.** |z| ≤ 0.8: the Maclaurin series (19 terms,
  Estrin). Next to an axis (2|xy| < 1/2): the Taylor series about that axis
  (erf^(k+1) = (2/√π)(−1)^k H_k e^(−t²), Hermite forward recurrence), which keeps
  the SMALL component accurate relative to itself (the w route below left it an
  absolute error of ~ε|e^(−z²)w|: `Erf[1.*^-300 + 2. I]` read `1.1*10^-16 + 18.56 I`
  in a first version; it is `6.16*10^-299 + 18.56 I`, as the MPFR path gave).
  |Re z|, |Im z| < 4: `erf = 1 − e^(−z²) w(iz)`, `erfc = e^(−z²) w(iz)`,
  with the Faddeeva function w from a 0.5×0.5 grid of degree-20 Taylor expansions
  (generated at 60 digits from w′ = −2ζw + 2i/√π; Poppe–Wijers / Zaghloul-style
  inner-region tabulation). Elsewhere: Salzer's series (A&S 7.1.29), re-derived as
  a trapezoidal sum with relative error 2e^(−4π²); its cosh(nY) factors are
  generated outward from the Gaussian peak with compensated sums. Every
  exponentiated quantity (x², y², 2xy) is formed as a Dekker double-double, which
  removes the |z|²-ulp loss of rounding the exponent. erfi(z) = −i erf(iz).
- **Real erf / erfc.** |x| < 1/2: Maclaurin; otherwise e^(−x²)·erfcx(x) from
  0.5-wide tables (6 ≤ x ≤ 27.4: x·erfcx as a polynomial in 1/x²). Faster than
  libm and more accurate (macOS libm erfc measured up to 5.1 ulp; now ≤ 3.1).
- **Real erfi.** e^(x²)·g(x), g = e^(−x²) erfi(x) (Dawson) from tables, 14-term
  asymptotic beyond x = 10, scaled where e^(x²) alone overflows. The old summed
  Maclaurin kernel (still used by nothing now) reached 30 ulps at |x| ≈ 5.
- **Fresnel.** (Complex next to an axis, π|xy| < 1/2: the Taylor series about
  the real axis, F(x+iy) = F(x) + φ(x)Σ q_n (iy)^{n+1}/(n+1)!, with C, S taking
  Re/Im of φq_n as real coefficients so the small component is exact-relative;
  C(iz) = iC(z), S(iz) = −iS(z) fold the imaginary axis.) Real: Maclaurin (|x| ≤ 1); auxiliary f, g tables (1 < |x| < 8,
  the Cephes `fresnl` structure with polynomial tables); asymptotic f, g beyond;
  the phase πx²/2 reduced exactly mod 4 from the double-double x², so large-x
  values are correctly rounded (`FresnelC[123456.789]`). Complex: Maclaurin in z⁴
  for |z| ≤ 1.3, else C ± iS = ((1±i)/2) erf((√π/2)(1∓i)z) through the erf
  kernels, fed the EXACT pieces X² = (π/4)(x+y)², XY = (π/4)(y²−x²),
  Y²−X² = ∓πxy instead of a rounded rotated argument (which costs π|z|² ulps).
- **Sinc.** Real sin(x)/x (unchanged value; scalar now off MPFR). Complex: even
  Maclaurin for |z| ≤ 1, else (sin x cosh y + i cos x sinh y)/z.
- **InverseErf / InverseErfc.** |p| < 0.875: degree-14 polynomials for erfinv(p)/p
  on 1/16-wide pieces (no logarithm); beyond, Giles' w = −log(1 − p²) (erfinv/p
  analytic in |Im w| < π) on unit pieces of degree 11 up to w = 6, then a 1e-7
  seed in 1/√w and one Halley step on erfc written with the relative residual.
  InverseErfc works from the exact complement min(q, 2 − q), so the tail is
  resolved down to the least normal double. Replaces Winitzki + up to four
  Newton steps (an erf and an exp each).
- **Bug fixed:** `InverseErfc` of a subnormal machine number answered
  `Indeterminate` (the Winitzki seed is log(0)); it now solves at 53 bits in MPFR
  (`InverseErfc[1.*^-320] = 27.0732`).
- Machine-complex inputs now give machine-complex results (`FresnelC[2. I] =
  0. + 0.488253 I`, as in Mathematica).

Tables are reproducible bit for bit: `benchmarks/.venv/bin/python -I
tools/sf_fit/gen_tables.py > src/special_functions/sf_machine_error_tables.h`
(fitlib: Chebyshev-node interpolation in a centred variable, constant term kept
as a double-double, every rounded piece re-measured in high precision).

## Accuracy (tools/sf_accuracy/error.py, 2000 points per band, mpmath 30 digits)

ulps of the result. "away" = relative max over points not near a zero
(|f| ≥ |z f′(z)|/4); "scaled" = |error| / ulp(max(|f|, |z f′(z)|)).

| head | regime / band | median | max (away) | max scaled | component median / p90 |
|---|---|---|---|---|---|
| Erf | MR [-3,3] | 0.28 | 1.61 | 1.61 | |
| Erf | C [-2,2]² | 0.42 | 3.45 | 2.70 | 0.93 / 2.71 |
| Erf | C near axes | 0.35 | 2.44 | 1.13 | 0.66 / 1.42 |
| Erfc | MR [-3,3] | 0.31 | 2.05 | 1.77 | |
| Erfc | C [-2,2]² | 0.55 | 3.47 | 1.98 | 1.01 / 3.26 |
| Erfc | C near axes | 0.48 | 2.22 | 1.06 | 0.67 / 1.38 |
| Erfi | MR [-3,3] | 0.46 | 2.17 | 1.86 | |
| Erfi | C [-2,2]² | 0.42 | 3.21 | 3.21 | 0.91 / 2.61 |
| FresnelC | MR [-3,3] / [-50,50] | 0.26 / 0.25 | 1.10 / 0.87 | 1.08 | |
| FresnelC | C [-2,2]² | 0.71 | 2.35 | 1.66 | 1.25 / 4.88 |
| FresnelC | C near axes | 0.28 | 1.26 | 1.19 | 0.54 / 1.07 |
| FresnelS | MR [-3,3] / [-50,50] | 0.32 / 0.25 | 2.41 / 1.98 | 1.14 | |
| FresnelS | C [-2,2]² | 0.83 | 3.49 | 1.90 | 1.37 / 4.72 |
| Sinc | MR [0.2,10] | 0.33 | 1.43 | 1.43 | |
| Sinc | C [0.5,5]×[-3,3] | 0.71 | 3.35 | 2.93 | 1.16 / 3.72 |
| InverseErf | MR [-0.95,0.95] / (-1,1) | 0.28 / 0.30 | 1.34 | 1.10 | |
| InverseErfc | MR [0.05,1.95] / 10^[-300,-1] | 0.31 / 0.25 | 1.99 / 0.88 | 1.94 | |

"near axes": one coordinate in [0.8, 4] (Fresnel [1.3, 6]), the other ±10^[-12,-1].
"component": each of Re, Im against its own ulp (a component crossing zero is
legitimately inexact relative to itself, hence median and p90 rather than max).
Wider bands (Erf/Erfc/Erfi C [-5,5]², Fresnel C [-4,4]², Sinc C [-20,20]²) are
within the same bounds (max away ≤ 3.3, scaled ≤ 2.5). The shootout's agreement
digits rose in every cell (e.g. Erf C 14.65 → 15.88, FresnelS MR 13.74 → 15.95).

## Speed (net µs, `run_shootout.py --regimes MR,C --engines mathilda,python`)

Same-run Python figures; "before" is the 2026-10-17 net baseline.

| head | reg | before | after | Python | best other (baseline) |
|---|---|---|---|---|---|
| Erf | MR | 0.102 | 0.103 | 0.124 | 0.125 |
| Erf | C | 64.6 | 0.264 | 0.248 | 0.243 |
| Erfc | MR | 0.103 | 0.107 | 0.123 | 0.124 |
| Erfc | C | 65.2 | 0.266 | 0.246 | 0.423 |
| Erfi | MR | 6.83 | 0.102 | 0.134 | 0.167 |
| Erfi | C | 65.3 | 0.264 | 0.248 | 0.259 |
| InverseErf | MR | 0.190 | 0.096 | 0.130 | 0.161 |
| InverseErfc | MR | 0.190 | 0.100 | 0.126 | 0.182 |
| FresnelC | MR | 11.9 | 0.121 | 1.29 | 1.23 |
| FresnelC | C | 119 | 0.305 | 1.93 | 1.93 |
| FresnelS | MR | 12.2 | 0.121 | 1.28 | 1.25 |
| FresnelS | C | 121 | 0.305 | 1.91 | 1.89 |
| Sinc | MR | 2.18 | 0.100 | 0.138 | 0.140 |
| Sinc | C | 5.03 | 0.234 | 0.158 | 0.191 |

10⁶ vectorised throughput (ms, Mathilda vs SciPy, same run): Erf 4.5 vs 14.2,
Erfc 3.2 vs 13.1, Erfi 3.1 vs 24.4, InverseErf 2.3 vs 18.7, InverseErfc 2.7 vs
15.3 (was tied 16.0 vs 15.9), FresnelC 5.3 vs 34.1, FresnelS 5.3 vs 34.2,
Sinc 2.7.

Kernel cost alone (C microbenchmark, 10⁶ random points per call site, one core,
same run; libm `erf` 19.6 ns as the yardstick): real erf 17.1, real erfc 18.4,
complex erf 59 (the first, Salzer-only version was ~3× that), complex erfc 62,
complex Fresnel C 103 (two erfs), complex Sinc 34, real Fresnel 33, real erfi
18.5, InverseErf 12.6, InverseErfc 15.6 ns.

## Still behind, and why

`Erf`/`Erfc`/`Erfi` C (by ~6%) and `Sinc` C (by ~48%): the kernels are 25–50 ns,
but a call that RETURNS a `Complex` pays ~115–170 ns more than one returning a
Real — the evaluator re-evaluates the fresh `Complex[r, i]` node (measured:
`Sinc[z]` 338 ns vs 178 ns when the same builtin returns a Real; `Re[z]`, complex
in / real out, 105 ns). Stamping the node as an evaluated GROUND fixed point in
`sf_machine_complex_result` (`r->last_evaluated_at = eval_clock_get() | GROUND`)
cut `Sinc[z]` from 338 to 224 ns in a trial build, which would put all four ahead
of SciPy. That is a shared-file / evaluator change (Phase 2), not made here.
