---
title: Extreme definite-integral stress campaign (battery #2)
date_started: 2026-10-13
status: baseline-complete, awaiting scope
corpus: tests/ExtremeDefiniteIntegralsCorpus.m (to be promoted from scratchpad)
---

# Extreme / "out-of-sample" definite-integral campaign (50 cases)

Battery `extremeAssumedIntegralStressTest` — a second, harder 50-case battery:
Frullani integrals & advanced rational poles (1–10), fractional branch cuts /
Mellin kernels (11–20), essential singularities / oscillatory / hyperbolic
kernels (21–30), complex-contour deflections & Cauchy PVs (31–40),
polylogarithms / special removals / Gamma-Zeta bounds (41–50).

## FOUR BATTERY ORACLE ERRORS (verified by NIntegrate — do NOT re-fix Mathilda to them)
- **6**  `(e^{-ax²}-e^{-bx²})/x²` on [0,∞): battery `2√π(√b−√a)` is **2× too big**.
  True = **`√π(√b−√a)`** (NIntegrate a=2,b=3 → 0.56335, and dF/da = −½√(π/a) ⇒ F=√π(√b−√a)).
- **24** `x sin(ax)/cosh²(bx)` on [0,∞): battery `πa/(2b² sinh(aπ/2b))` is the formula for
  `∫cos(ax)sech²(bx)`. True = **`(π/(2b²))·csch(aπ/2b)·(aπ/2b·coth(aπ/2b) − 1)`**
  (= −d/da of the cos integral; matches NIntegrate to 30 digits at a=1,b=2).
- **36** `√x/(x²−a²)` PV on [0,∞): battery `−π/(2√a)` has the **wrong sign**.
  True = **`+π/(2√a)`** (ε=1e-6,hi=1e8 PV → +1.11052→+1.11072).
- **37** `x log x/(x²−a²)²` PV on [0,∞): a **double pole** (Hadamard finite part), oracle
  is regularization-convention-dependent; battery `−1/(4a²)=−0.0625` matches no standard
  reg (our 2A/ε-subtracted FP → −0.08664 ≈ −log2/8). **EXCLUDED from scoring** (informational).

## TWO DEGENERATE INSTANTIATIONS (fixed in corpus so an eval-to-0 can't score a false CORRECT)
- **12** `x^(a−1)/(1−x^b)` PV: a→3/2,b→3 ⇒ Tan[π/2]=∞ ⇒ expected 0. Changed to **a→1,b→3**.
- **19** `x^a/(1−x)²` PV: a→1/2 ⇒ Cot[π/2]=0 ⇒ expected 0. Changed to **a→1/3**.

## Baseline (v0.387, corrected/verified oracles, one process per case, 30s TMO)
**25 CORRECT / 0 WRONG / 1 FORM / 19 UNEVAL / 4 TLE / 1 EXCLUDED.**
Soundness bar already met (zero wrong). Work = converting UNEVAL/TLE/FORM → CORRECT.

- CORRECT (25): 1,3,4,5,7,8,9,10,11,13,14,15,16,17,18,20,22,23,39,40,42,43,45,48,49
- FORM (1): **21** `Exp[-I a x²]` whole-line → value-correct but wrapped as
  `2 ConditionalExpression[√π/(2√(I a)), I a>0]` (spurious condition, `√(I a)` unreduced).
- UNEVAL (19): 2,6,12,19,25,26,27,28,30,31,32,33,34,36,41,44,46,47,50
- TLE 30s (4): 24,29,35,38
- EXCLUDED (1): 37 (double-pole Hadamard FP, ambiguous oracle)

## Open cases grouped by machinery (24 scored targets)
- **Frullani** (2): `(cos ax−cos bx)/x → log(b/a)` — exp variant (1) & arctan (3) already pass;
  cos needs CosIntegral-limit handling. General Frullani recognizer.
- **Gaussian difference** (6): `(e^{-ax²}−e^{-bx²})/x² → √π(√b−√a)`.
- **Mellin PV** (12,19): `x^(a−1)/(1−x^b)` PV → `π/(b tan(aπ/b))`; `x^a/(1−x)²` PV → `πa cot(πa)`.
- **Hyperbolic kernels** (24,25,29,30): `∫₀^∞ {1,x,x²}·trig/hyperbolic` →
  π·csch/tanh/coth families (Mittag-Leffler / rectangular-contour residue sums). 23 passes.
- **Bessel Weber-Schafheitlin** (26,27,28): `J₀(ax){cos,sin}(bx)` discontinuous integrals →
  `1/√(a²−b²)` / `0` / `1/√(b²−a²)`.
- **Half-line rational PV, on-axis pole** (31,32,33,34): `{1,x,x²}/(x^n−a^n)` PV →
  indented-contour half-residue. The biggest PV gap.
- **Log/√ over (x²−a²) PV** (35,36): `log x/(x²−a²)→π²/(4a)`; `√x/(x²−a²)→+π/(2√a)`.
- **Whole-line log** (38,44): `log|x−a|/(x²+b²)→π log(a²+b²)/(2b)`;
  `log(x²+a²)/(x²+b²)→2π log(a+b)/b`.
- **Half-line log-rational / log-√** (41,46,47,50): `log³x/(1−x)→−π⁴/15`;
  `log(1+x²)/(1+x²)→π log2`; `log x/√(1−x²)→−(π log2)/2`; `log Γ(x)` on [0,1]→`log(2π)/2` (Raabe).
- **Fresnel whole-line form** (21): ConditionalExpression strip + `√(I a)` branch reduction.

## PROGRESS (live) — extreme 42/50 correct, 0 wrong; battery #1 44/50 (bonus)
| Version | Cases closed | Mechanism |
|---|---|---|
| v0.388 | 2 | Frullani cosine (Cesàro mean at ∞) + scaled-arg CosIntegral series |
| v0.389 | 24,25,29,30 | half-line hyperbolic-kernel table (residue_family_hyperbolic_halfline); fixed 2 hangs |
| v0.390 | 26,27,28 | Weber-Schafheitlin Bessel (rec_bessel_cossin in intrep) |
| v0.391 | 35 **+ battery#1 48,49** | a0_keyhole_log symbolic-param guard (un-hung the keyhole-log grind) |
| v0.392 | 31,32,33,34 | half-line rational PV, on-axis pole (residue_family_rational_halfline_pv → m=0 keyhole) |
| v0.393 | 36 | branch-power Mellin keyhole PV on-axis half-residue (mellin_core) |
| v0.394 | 44,46 | whole-line log-over-quadratic (residue_family_logquad) + even-fold |

Correct now (42): 1-5,7-18,20-34,36,39,40,42-46,48,49 (+21 value-correct, 37 excluded).
**Open (7 scored)**: 6 (Gaussian diff), 12,19 (Mellin PV), 38 (whole-line Log|x-a|, hangs in
diffunderint + piecewise grabs it), 41 (Log^3 x/(1-x) on [0,1]), 47 (Log x/Sqrt[1-x^2] on [0,1]),
50 (Raabe Log Gamma).
Diagnostic specs in hand (subagents): 41 (new integrate_ratlog01_try: (-1)^n n! HurwitzZeta[n+1,c];
also up-front Beta decline to stop Infinity::indet leak), 47 (extend integrate_beta_try with
(1-x^2)^e via t=x^2), 50 (new integrate_raabe_try: Log[Gamma[x]] on [t,t+1] → ½Log2π + t Log t - t),
38 (whole-line Log[Abs] = ½Log[(x-a)^2] via logquad on real double root; needs piecewise to decline +
diffunderint divergence gate). 6 (Mellin of exp-difference at s=-1: √π(√b-√a)). 12,19 (Mellin x^(a-1)/(1∓x^b) PV).

## Invariant
Zero wrong answers. Each substantive change: full corpus (no new WRONG) + touched unit suites
+ check-messages + check-c99 + valgrind new C; version bump + git tag + changelog + docs.
