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

## MPFR / arbitrary precision (50 digits) — HEADLINE FINDING

All values agree to 15+ digits across Mathilda / mpmath / wolframscript (correctness
is solid). **But wolframscript is 3×–890× FASTER than Mathilda at 50 digits on every
function tested:**

| 50-digit, ms/call | Mathilda | wolframscript | mpmath |
|---|---|---|---|
| Gamma | 0.035 | 0.003 | 0.028 |
| Erf | 0.019 | 0.003 | 0.031 |
| Zeta | 0.69 | 0.0035 | 0.53 |
| PolyLog | 4.10 | 0.0046 | 1.94 |
| LerchPhi | 2.90 | 0.010 | 47.9 |
| HypergeometricU | 0.28 | 0.012 | 4.50 |
| EllipticK | 0.009 | 0.003 | 0.029 |
| BesselJ | 0.040 | 0.004 | 0.024 |
| AiryAi | 0.21 | 0.003 | 0.13 |
| FresnelC | 0.037 | 0.009 | 0.073 |

Mathilda is competitive with / faster than **mpmath** on several (LerchPhi 16×,
HypergeometricU 16×, Erf/EllipticK/FresnelC 1.6–3×), but trails wolframscript
everywhere. Mathilda's per-call time tracks the function (Gamma 35µs vs PolyLog
4ms) while wolframscript's is a near-flat ~3–12µs — i.e. the cost is Mathilda's
per-function MPFR series/evaluation, not a fixed per-call overhead, and the
transcendental series heads (PolyLog, LerchPhi, Zeta, HypergeometricU, AiryAi) are
the worst. This **reverses** the narrower impression from the Whittaker review
(where wolframscript's own `WhittakerW` happened to be slow); on the bread-and-butter
heads wolframscript's arbitrary-precision kernels are dramatically faster.

**This is a broad gap across dozens of functions — out of scope for this pass.**
Recommended as a dedicated follow-up campaign (profile: per-call MPFR context/alloc
overhead vs series efficiency; candidates: cache MPFR scratch, raise FLINT/arb use
beyond elliptic+HypergeometricU, tune series term counts at low precision).

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
