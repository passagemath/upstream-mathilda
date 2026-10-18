# Mathar integral corpus: triage and extension plan

Source: R. J. Mathar, *Yet Another Table of Integrals*, arXiv:1207.5845, from its
TeX source. The corpus files' headers record the sha256 of that source.

Corpus:
- `tests/MatharDefiniteIntegralsCorpus.m`: 580 records.
- `tests/MatharIndefiniteIntegralsCorpus.m`: 135 records.
- `tests/MatharSkipped.md`: the remaining integral sites, each with its reason.

Runner: `tests/test_mathar_corpus.c`. Tooling: `tools/mathar/`.

## 1. Building the corpus

| Step | Tool | Result |
|---|---|---|
| Harvest every display block / inline fragment containing `\int` | `extract.py` | 1,007 sites |
| Translate TeX → InputForm (Pratt parser; Mathematica's TeX import cannot read these) | `tex2m.py` | 593 auto-translated |
| Numeric oracle (Mathematica): the paper's value must equal `NIntegrate` (no messages) at two parameter instantiations | `verify.wls` | 369 verified automatically |
| Hand repair of the rest, five parallel agents, every reading re-verified by the oracle | `overrides.json` | 574 entries reviewed; about 90 typos in the paper corrected |
| Emit the corpora | `instance.wls`, `build_corpus.py` | 715 records; 292 skips |

The 292 skips break down as follows. 108 are multiple integrals and 44 are recurrences or reduction formulas.
Another 23 are definitions and 21 involve generic functions (f, Z_ν). There are 15 unrepairable paper
errors, and 20 entries are numerically unverifiable (12 judged so by a hand-repair agent, 8 by the oracle).
The remaining 61 are miscellaneous; `tests/MatharSkipped.md` lists each one.

Oracle hardening, learned the hard way:
- Sum/Product dummy indices are not free symbols.
- Use `Sqrt[u^2]` for `Abs` in the derivative test.
- Draw integers for parameters declared integer.
- Use Landen forms of K′/E′ near 0.
- Reject an instantiation whose value is below 10⁻⁶. The relative test is vacuous there.

Each record runs twice:
- **S**: symbolic, with `Assumptions`.
- **N**: a numeric instance.

Gates: WRONG, NONNUMERIC and CRASH must be 0, and the CORRECT floors (S, N) ratchet in `tests/CMakeLists.txt`.

## 2. Bugs the corpus exposed (all fixed in v0.458)

| Bug | Where | Corpus symptom |
|---|---|---|
| `Series`/`Limit` of `ArcTan[k x]` at ∞ took +π/2 for a negative radical `k` | `series.c` | ∫₀^∞ dx/(1+2x²+x⁴/2) = 2.05 (true 0.850) |
| Bigint coefficients enrolled as polynomial variables | `poly.c` `is_number` | ∫ x/(625/4096+x⁴) dx = 0 |
| `FactorSquareFree` returned a product ≠ input (missing cofactor set to 1) | `facpoly_squarefree.inc` | ∫₀^∞ log²x/(a²+x²) = −π/(6a) |
| Symbolic residue families classified poles at one sample point | `integrate_residue.c` | 1/(1+2ax²+bx⁴) wrong for a² > b |
| ParallelMixedTower treated `g[Cos[t]]` as a constant, leaked its tower variable | `integrate.c` | antiderivative containing `r25` |
| ParallelMixedSpecial cascade accepted a half-evaluated package body | `integrate.c` | `{0, 0}` as a definite value |
| `Rational[INT64_MIN, d]` → SIGFPE | `arithmetic.c` | crash |
| Use-after-free in extension autodetect | `qafactor.c` | (debug build) |
| TimeConstrained unwind freed Block values still in use | `modular.c` | 8 timing-dependent segfaults |
| `E^(I Pi)` unevaluated; `PolyLog[2, I]` unevaluated | `power.c`, `polylog.c` | polylog values not closing |
| False `Integrate::idiv` at a radical endpoint | `integrate_newton_leibniz.c` | ∫₀^{2−√3} log x/(1+x²) "divergent" |

Each fix has a root-cause test in `tests/test_mathar_regressions.c`.

## 3. Results

| | baseline | + fixes | + PolyLogLadder (v0.458) |
|---|---|---|---|
| definite S | 49 | 49 | **80** |
| definite N | 124 | 124 | **166** |
| indefinite S | 63 | 63 | **66** |
| indefinite N | 79 | 79 | **82** |

v0.459 adds Log-of-trig: definite 93/191, indefinite unchanged 66/82.

WRONG, NONNUMERIC and CRASH are all 0. CMake floors: 78/163 and 65/80, a small margin for timeouts.

## 4. What remains: ranked by "Mathematica closes it"

`tools/mathar/mma_integrate.wls` runs Mathematica's own `Integrate` (60 s, parameters assumed positive) on
every record. Mathematica closes 414 of the 715. Of Mathilda's failures, 270 are integrals Mathematica
closes; those are the algorithmically reachable targets. The table counts symbolic-mode failures after
PolyLogLadder.

| Class | Records | Failing | …of which Mathematica closes |
|---|---|---|---|
| definite trig (rational/power/elliptic-type) | 99 | 80 | 42 |
| definite log | 89 | 71 | 41 |
| definite Bessel | 71 | 70 | 38 |
| definite log × trig | 57 | 55 | 30 |
| definite trig × algebraic | 23 | 21 | 16 |
| indefinite Bessel | 19 | 19 | 13 |
| definite exp | 26 | 23 | 10 |
| indefinite log | 22 | 10 | 9 |
| definite / indefinite other special (Gamma, Zeta, Erf, Ei, ...) | 37 | 35 | 16 |
| definite algebraic | 21 | 11 | 8 |
| definite elliptic | 38 | 38 | 5 |
| orthogonal polynomials (definite + indefinite) | 44 | 44 | 4 |

Elliptic and orthogonal-polynomial entries are mostly generic-order identities that Mathematica cannot
close either. They rank last.

## 5. Extension stages

### Done: PolyLogLadder (v0.458)

Module `src/internal/PolyLogLadder.m`; stage and `Method -> "PolyLogLadder"`. It handles two classes:
- `R(x) Log[a x + b]^k`, by partial fractions in the kernel variable and the polylogarithm ladder.
- A polynomial × trig/hyperbolic rational, rewritten in one exponential kernel, by the exponential ladder
  with ρ-derivatives for repeated roots.

Both are exact identities, accepted after a differentiation check. Supporting changes: `E^(I q Pi)`,
`PolyLog[2, ±I]`, and a real-form finishing pass in Newton–Leibniz.

### Next, in rank order

1. **Done in v0.459. Log of a trig rational: `P(x) Log[G(trig)]`** (log × trig, about 30 closable).
   - Write `G` in T = E^(i b x) and factor it: G = K T^e Π(T − ρᵢ)^{mᵢ}.
   - Then d/dx Log G = e i b + Σ mᵢ i b T/(T − ρᵢ), so Log G = C + e i b x + Σ mᵢ Log(1 − T/ρᵢ) on each
     interval where no factor crosses a branch cut.
   - Each `x^n Log(1 − T/ρ)` is the exponential ladder `x^n Li₁(T/ρ)`, which integrates by parts to
     `Σ x^{n−j} Li_{2+j}(T/ρ)/(i b)^{j+1}`.
   - The constant C must be fixed per branch interval. Evaluate it at an interior point and split the
     definite range at the zeros of the factors on the real axis.
   - Closes ∫ log sin, ∫ log cot, ∫ x log 2sin(πx/2), ∫ log(1 + a² + 2a cos t). The values come out as
     Catalan, Clausen, ζ.
2. **Bessel × Bessel × power on [0, ∞) by Mellin–Barnes / MeijerG** (38 definite).
   - Each kernel (J, Y, I, K, exp, power, (1+x²)^ν) is a MeijerG.
   - ∫₀^∞ G₁(ax) G₂(bx) x^{s−1} dx is a single MeijerG by the convolution theorem.
   - Reduce that G back to hypergeometric / Gamma / elementary forms.
   - `src/special_functions/meijer_g.c` already has the G-function. Needed: (a) a kernel → MeijerG table,
     (b) the convolution, (c) G → closed-form reduction for the common cases (Gamma ratios, ₂F₁ at
     special points).
3. **Indefinite Bessel × power** (13). These are the Lommel-type recurrences: ∫ x^{μ} Z_ν(x) by the
   recurrence on (μ, ν) down to the base cases ∫ x^{ν+1} Z_ν = x^{ν+1} Z_{ν+1}, plus Struve for the rest.
4. **Complete elliptic reduction on [0, π/2]**.
   - Applies to `sin^{2m} cos^{2n} (1 − k² sin²)^{±1/2}` (trig-only entries M8035, M8090, ...).
   - Reduce to K, E with the Byrd–Friedman recurrences for ∫ sin^{2m}/Δ and ∫ sin^{2m} Δ. Every term is a
     rational combination of K(k²) and E(k²).
   - These currently time out.
5. **Multi-kernel log products** (indefinite log, 9): `Log[x] Log[1+x]/x` and similar. They need Li₃ of
   cross-ratios or Nielsen S_{n,p}. Start with the pole-at-a-kernel-root case, which is the classical Li
   ladder by parts.
6. **Polylog constant reduction.** Values like `PolyLog[2, 1 ± I]` and `PolyLog[2, (1 ± I)/2]` should
   reduce to π, log 2, Catalan. This is cosmetic: the corpus already counts them correct. It needs a table
   of Li₂ at Gaussian rationals.

## 6. Open bugs found by the corpus, not yet fixed

- `fullsimplify_corpus_tests`: `FullSimplify[Log[Cosh[x] + Sinh[x]]]` gives `Log[E^x]`. This predates the
  campaign and also fails on HEAD.
- Expressions written as `E^(I q Pi)` with q not a half-integer stay symbolic (as in Mathematica). Whereas
  `Exp[I Pi/3]` → `(-1)^(1/3)` in Mathilda but stays `E^(I Pi/3)` in Mathematica: a cosmetic divergence.
- ParallelMixedSpecial emits `Set::wrsym` on `Plus`/`Missing` for √(1−k² sin²θ) integrands. This
  predates the campaign; package state is suspect.
