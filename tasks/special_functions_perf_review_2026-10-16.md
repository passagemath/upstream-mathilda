# Special-functions performance review (all heads) — 2026-10-16

Applied the Whittaker review method to all 61 special-function heads, at **machine**
and **MPFR** precision, vs SciPy / mpmath / wolframscript (Mathematica 14.0).
Machine axis driven by `tools/numeric_sweep.py --group special`; MPFR axis by a
50-digit spot-check (`$CLAUDE_JOB_DIR/tmp/mpfr_spot.{m,py}`, min/one pass).

## Machine precision (double) — vs SciPy

**Mathilda is faster than SciPy on the large majority**: Gamma 4.5×, Bessel
10–120×, Airy(+primes) 10×, ProductLog 15×, ExpIntegralEi 10×, LogIntegral 7×,
incomplete elliptics 8–9×, Beta 11×, Hyper0F1 32×, HurwitzZeta 4.9×, HermiteH 5.3×,
Laguerre 4.7×, Chebyshev 2.6–3.3×, Legendre 2.2–2.9×, Sinc/Fibonacci/Lucas ~6×,
HypergeometricU 2.4×, and more.

### Fixed this pass (each a tagged commit)
| Head | Before | After | Fix |
|---|---|---|---|
| PolyLog (v0.424) | 24× slower | **4.75×** | constant-ζ tables (no per-term `sf_machine_zeta`) + integer power in the defining series |
| FresnelC/S (v0.425) | 5.1× slower | **2.4×** | Lentz CF rewritten in real/imag pairs (no libgcc `__divdc3`) |
| LerchPhi (v0.427) | 10⁶ timed out | **1/1.1× (faster)** | scalar + buffer now use `sf_machine_lerchphi`; 628→2.5 µs/call |

### Residual machine gaps — documented, not fixed (accuracy-gated / architectural)
- **PolyLog 4.75×, Fresnel 2.4×**: SciPy's `spence`/`fresnl` are specialised
  minimax-rational approximations (no loop); Mathilda uses general series/CF. Full
  parity needs a dedicated real dilog and Fresnel auxiliary-rational path.
- **Zeta 2.58×**: `pow(n,−s)`×15 per element, `s` varies so the values cannot be
  cached across the buffer. The one cheap win (static `ln n` table + `exp(−s·ln n)`)
  trades ~1 ulp of accuracy — rejected (the kernel bit-tracks MPFR). Needs a tuned
  Euler–Maclaurin, not a constant-factor hack.
- **Hypergeometric2F1 1.62×**: the `pfq` ascending series; SciPy's `hyp2f1` is tuned.
- **InverseErf 1.23×, InverseErfc 1.31×, EllipticE 1.15×**: marginal, within
  measurement noise and the SIMD-prefactor effect seen on cheap kernels (NumPy
  vectorises `exp`/`pow` across the array; our per-element libm cannot).

### Coverage added (permanent probes in `numeric_sweep.py`)
16 kernel-backed heads that had no probe: HurwitzZeta, Hypergeometric0F1/PFQ/U,
Chebyshev T/U, HermiteH, LaguerreL, LegendreQ, AiryAiPrime/BiPrime, InverseErfc,
Sinc, Fibonacci, LucasL, LerchPhi. All agree with SciPy/numpy (0 value mismatches);
almost all are faster than the baseline.

## MPFR / arbitrary precision (50 digits) — CORRECTED + resolved (v0.428–v0.432)

**Correction.** This section first reported "wolframscript 3×–890× FASTER" — that
was a **measurement artifact**. The harness warmed up then timed the *same*
arguments, and Mathematica memoizes `N[expr, prec]`, so the timed pass measured
cache hits (~3µs), not computation. Re-measured **cache-safe** (distinct argument
per call, single timed pass; confirmed with `ClearSystemCache[]`), Mathematica's
real 50-digit times are 0.06–2.2 ms/call, and Mathilda was already faster on most
heads. Four heads were genuinely slower — closed by routing their transcendental
series through FLINT/arb (v0.428–v0.432).

**Final cache-safe standing (50 digits, ms/call; a broad 28-head sweep):** Mathilda
is **faster than Mathematica on every head measured.** Representative:

| head | Mathilda | wolfram | speedup | | head | Mathilda | wolfram | speedup |
|---|---|---|---|---|---|---|---|---|
| Gamma | 0.037 | 0.090 | 2.4× | | PolyLog | 0.027 | 0.747 | **28×** |
| Erf | 0.019 | 0.202 | 10.6× | | LerchPhi | 0.156 | 2.236 | 14× |
| Zeta | 0.098 | 0.435 | 4.4× | | Hyp2F1 | 0.422 | 0.723 | 1.7× |
| HurwitzZeta | 0.106 | 0.312 | 2.9× | | Hyp1F1 | 0.115 | 0.899 | 7.8× |
| PolyGamma | 0.102 | 0.341 | 3.3× | | HypergeometricU | 0.269 | 2.176 | 8.1× |
| EllipticK | 0.008 | 0.024 | 3.0× | | BesselJ | 0.037 | 0.101 | 2.7× |
| AiryAi | 0.014 | 0.065 | 4.6× | | FresnelC | 0.037 | 0.454 | 12× |
| SinIntegral | 0.015 | 0.393 | 27× | | Beta | 0.096 | 0.338 | 3.5× |

**What changed (arb routing via the `flint_num_*` bridge):** Zeta → `acb_dirichlet_zeta`
(v0.428); PolyLog real −1<z<1 → `acb_polylog` (v0.429, ~45×); AiryAi/Bi(+primes) →
`acb_hypgeom_airy` (v0.430, ~9×); LerchPhi |z|≤1 → `acb_dirichlet_lerch_phi` (v0.431,
~15×); HurwitzZeta → `acb_dirichlet_hurwitz` (v0.432, ~10×). Each keeps the machine
`double` path and branch-cut/pole/exact handling; accuracy exact vs 30 digits;
valgrind-clean (no Mathilda-code leaks). The feared central-`numeric.c` per-call
overhead rewrite was **not needed** — Mathematica carries comparable per-call cost.

**Lesson:** benchmark a CAS that memoizes with **distinct arguments per timed call**
(or `ClearSystemCache[]`); a warm-up + same-args timed pass measures cache hits.

**Residual (not speed):** `StruveH` does not numericalize at arbitrary precision
(`N[StruveH[1,1.3],50]` stays symbolic) — a coverage gap (FLINT has no Struve
kernel), tracked separately.

## Correctness fixes
- **LerchPhi (v0.426)**: a visible `NDArray` in any slot was left **unevaluated**
  (wrong answer per the packed contract). Now delists / threads the kernel; added to
  `pack.c` AWARE. (v0.427 then made that path fast.)
- **LegendreP mismatch** was a checksum artifact, not a bug: `LegendreP[2, u]` summed
  over the symmetric `u`-grid cancels to ~1e-6 (P₂ is orthogonal). `LegendreP` itself
  is correct to <6e-18 vs exact. Probe argument changed `u → p`.

## Feature gaps (not implemented — a separate task, not perf)
- **GammaRegularized / BetaRegularized**: not registered anywhere in `src/`; returned
  unevaluated. The two phantom sweep probes were removed. Implementing them (regularised
  incomplete Γ/B, `mpfr_gamma_inc` / `betainc` based) is a feature task.
- **EllipticPi incomplete third kind**: no SciPy one-liner, so no machine baseline
  (the complete case is covered via the Carlson identity).

## Scoped follow-ups (priority order)
1. **MPFR arbitrary-precision speed** (broad; the biggest open gap vs wolframscript).
2. Dedicated real **dilog** rational approx (PolyLog s=2 → beat SciPy's `spence`).
3. **Fresnel** auxiliary-rational path for `x > 1.5` (→ beat SciPy's `fresnl`).
4. **Zeta** tuned Euler–Maclaurin; **Hyper2F1** series tuning (accuracy-gated).

## Verification (per fix)
Accuracy vs 30-digit FLINT/MPFR (≤~2e-16 rel); scalar == NDArray == Compile
bit-identical where the scalar path shares the kernel; `check-c99`, `check-messages`,
`check-packed-aware`, `check-nd-surfaces`, `check-array-exactness` green; relevant
`tests/*_tests` pass. Tags v0.424–v0.427; the `numeric_sweep.py` coverage/harness
commit rides untagged (contributor-facing).
