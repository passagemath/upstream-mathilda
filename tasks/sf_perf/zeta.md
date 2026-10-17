### Zeta family: machine-precision kernels for Zeta, HurwitzZeta, LerchPhi, PolyLog

Machine-precision calls to `Zeta[s]`, `Zeta[s, a]`, `HurwitzZeta[s, a]`, `LerchPhi[z, s, a]` and
`PolyLog[s, z]` used to go through MPFR, or through file-local complex-MPFR toolkits, and were then
rounded to a double. They now go through one set of raw-`double` kernels in the new
`src/special_functions/sf_machine_zeta.{c,h}`. These are the same helpers that the packed and
`NDArray` element-wise kernels and `Compile[]` call, so the five surfaces agree bit for bit:
scalar, unpacked List, packed List, `NDArray` and `Compile`. All four builtins put the shared
`sf_machine_scalar_args` gate ahead of their exact-value ladders. Anything a kernel declines (a
pole, a domain edge, cancellation, overflow) falls through to the arbitrary-precision path, which
is unchanged.

**Algorithms** (published, cited in the source):

- **Zeta and Hurwitz.**
  - Integer `s` returns a tabulated value: ζ(2…64) and ζ(0…−60).
  - For `s ≥ 20` the Dirichlet series is summed directly.
  - Otherwise Euler–Maclaurin (DLMF 25.2.9 / 25.11.5) is used. The direct sum runs to
    w ≥ max(10, |s|/2 + 4) and the Bernoulli correction is cut at 2⁻⁵⁶.
  - Each `k + a` is formed by an exact TwoSum, and its power is corrected to first order.
  - Real `s < ½` uses the functional equation:
    - sin(πs/2) is reduced exactly in `s`;
    - (2π)ˢ is formed from an FMA two-product;
    - the rounding of `1 − s` is corrected to first order in Γ and ζ. Without this correction
      `Zeta[-31.02]` was off by 100 ulp.
  - Complex `Re s < ½` declines. The Γ(1−s) factor cannot be formed in double: it measured
    18 ulp median, and the MPFR path stays.
  - A complex Hurwitz sum whose terms exceed 16× its value declines.
- **Dilogarithm.** The Bernoulli series Σ Bₙ uⁿ⁺¹/(n+1)! in u = −log(1−z) ('t Hooft–Veltman 1979;
  Crandall 2006). The maps z → 1/z and z → 1−z land every z where |u| ≤ 1.25, so at most ~14 terms
  are needed. π²/6 is carried as a two-double constant. A real x > 1 is taken from below the cut,
  as in Mathematica.
- **PolyLog, other real orders.**
  - s = 0 and s = 1: closed forms.
  - |z| ≤ ½: the power series, with a double-double running power and Neumaier summation.
  - |μ| < 3.6: the Jonquière expansion in μ = log z (DLMF 25.12.11/12), with a complex μ for
    negative x.
  - Real x < −½ with s ≥ 1: the duplication formula.
  - Integer s and larger |z|: the inversion formula with Bernoulli polynomials (DLMF 25.12.13).
  - Near an integer order the Jonquière pieces have cancelling poles. If the terms exceed 4× the
    sum (16× for integer s), the kernel declines.
  - This replaces the old real `sf_machine_polylog`. That kernel used the duplication formula for
    every negative x, which cancels catastrophically when s < 1: `PolyLog[-1.88, -0.987]` was off
    by about 4 × 10⁸ ulp on the packed path.
- **LerchPhi.**
  - The defining series, now also for complex z, s and a with `Re a > 0`.
  - zᵏ is carried in double-double, with an FMA residual per product. This brings the real kernel
    from 14 ulp max to 1.1.
  - TwoSum `k + a`, Neumaier summation, and a cancellation guard. The old real kernel returned
    garbage for s < 0 near z = −1: 17 575 ulp.

**Surfaces.**

- `ndkernels.c` gains complex arms for Zeta, HurwitzZeta, PolyLog and LerchPhi.
- `Compile[]` gains `OP_KERNN_C`, a complex n-ary kernel opcode. With it, `LerchPhi` with complex
  arguments compiles; until now the n-ary path was real only.
- `CompileDiagnostics` reports `Compiled -> True` for all four heads, real and complex, at scalar
  and rank-1 shapes.

**Net latency** (µs, `run_shootout.py`, the same run as Python; two runs):

| head | regime | before | after | Python (same run) | best other (baseline) |
|---|---|---|---|---|---|
| Zeta | MR | 74.7 | 0.31–0.37 | 1.6–1.9 (scipy) | 1.35 (scipy) |
| Zeta | C | 403 | 0.79–0.95 | 444–531 (mpmath) | 20.9 (Maxima) |
| HurwitzZeta | MR | 275 | 0.25–0.30 | 1.4–1.6 (scipy) | 1.19 (scipy) |
| HurwitzZeta | C | 709 | 0.97–1.17 | 2236–2470 (mpmath) | 316 (Mathematica) |
| LerchPhi | MR | 0.88 | 1.7–2.0 | 21 669–25 608 (mpmath) | 149 (Mathematica) |
| LerchPhi | C | 571 | 2.2–2.7 | 44 856–44 893 (mpmath) | 149 (Mathematica) |
| PolyLog[2, ·] | MR | 494 | 0.14–0.16 | 592–618 (mpmath) | 6.9 (Mathematica) |
| PolyLog[2, ·] | C | 6560 | 0.36–0.43 | 1542–1657 (mpmath) | 115 (Mathematica) |

LerchPhi MR is about 2× slower than before. That cost buys the accuracy fixes above, and the head
is still about 80× ahead of Mathematica.

**10⁶-element buffers:**

| call | Mathilda | SciPy |
|---|---|---|
| `Zeta[v]` | 29 ms | 277 ms (`zeta(v, 1)`) |
| `HurwitzZeta[2.5, a]` | 35 ms | 294 ms |
| `PolyLog[2, x]` | 4.8 ms | 27 ms (`spence(1 − x)`) |

**Accuracy** against mpmath (`tools/sf_accuracy/zeta.py`, 2000 points per row, through the real
builtins):

| row | median ulp | max ulp |
|---|---|---|
| Zeta MR [1.5, 9] | 0.30 | 1.03 |
| Zeta C [1.5, 5]×[−3, 3] | 0.26 | 1.46 |
| Zeta real [−40, ½] | 0.72 | 5.7 |
| HurwitzZeta MR | 0.35 | 1.43 |
| HurwitzZeta C bench | 0.78 | 4.7 (κ-scaled 0.88) |
| LerchPhi MR | 0.28 | 1.12 |
| LerchPhi C bench | 0.22 | 0.84 |
| Li₂ MR | 0.42 | 2.25 |
| Li₂ C bench | 0.39 | 3.94 |
| Li₂ C wide box ±20 | 0.39 | 2.6 |

In the wider bands the large errors follow the function's condition number:

- complex Hurwitz with large `Im s`: median 5, max 184 ulp, ≤ 3.9 κ;
- complex ζ near the critical line: max 82 ulp, ≤ 10 κ;
- PolyLog of a general real order: median 0.35, max 10 ulp in real, 29 in complex.
