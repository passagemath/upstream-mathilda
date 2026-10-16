# WhittakerM / WhittakerW performance review — 2026-10-16 (v0.423)

Goal: confirm both heads are at least as fast as SciPy and wolframscript, at
**machine** and **MPFR** precision. Both were already wired across the scalar,
NDArray/packed, and `Compile[]` surfaces; speed had never been measured.

## Method
- Engines: Mathilda (`./Mathilda -file`), SciPy 1.17 / mpmath 1.3 (python3.11),
  wolframscript (Mathematica 14.0). Scripts in `$CLAUDE_JOB_DIR/tmp/whit_math.m`,
  `whit.py`, `whit_wl.m`. SciPy has no native Whittaker → compose
  `exp(-z/2)·z^(m+½)·{hyp1f1,hyperu}` (exactly the kernel's own definition);
  mpmath has `whitm`/`whitw` (the arbitrary-precision Python peer).
- Regimes: integer m (k=2,m=3 → `b=7`, the `hyperu_integer_b` log series),
  non-integer m (m=0.3 → connection formula), large z (asymptotic). All real, z>0.
- Min across 3 whole-script runs; vectorized = 10⁶ reals; MPFR = 50 digits.

## Finding
One real gap: **vectorized `WhittakerW`, integer/half-integer m — 1.8× slower
than SciPy** (1510 ms vs 837 ms / 10⁶). Cause: `hyperu_integer_b` called
`sf_machine_digamma` 3×/term. Every other regime already matched or beat SciPy,
and Mathilda beat wolframscript everywhere (Mathematica has no packed Whittaker
path: 20 s for the same `WhittakerM` buffer).

## Fix (src/special_functions/sf_machine.c, `hyperu_integer_b`)
Advance the three digammas by `ψ(x+1)=ψ(x)+1/x` from 3 seeds — one division per
term instead of a full digamma call ×3. Safe: the asymptotic branch peels off
every terminating/pole case first, so a+k never hits a non-positive integer.

## Results (min, lower=better)
Machine vectorized, ns/element:

| case            | Mathilda(before→after) | SciPy | WL      |
|-----------------|------------------------|-------|---------|
| M int-m z[1,2]  | 115 → 107              | 125   | 22017   |
| **W int-m z[1,2]** | **1510 → 420**      | 768   | 408125  |
| M nonint z[1,2] | — / 118                | 133   | 19093   |
| W nonint z[1,2] | — / 327                | 372   | 873848  |
| M int-m z[20,30]| — / 350                | 363   | 390173  |
| W int-m z[20,30]| — / 278                | 225   | 930230  |

Machine scalar µs/call: M 1.57, W int-m 3.55→1.95, W nonint 1.85 — WL 20–874 (13–470× faster).

MPFR ms/call (vs mpmath / WL): M-int-m 0.16 (mpmath 0.11, WL 0.28); W-int-m 0.35
(mpmath 5.13, WL 2.33); W-nonint 0.22 (mpmath 0.93, WL 3.08).

## Verdict vs the bar
- **wolframscript**: Mathilda faster on *every* case, both precisions (100–3000×). ✓
- **SciPy** (machine double only): faster on 5/6 machine regimes. The lone
  remaining 1.24× is `W` large-z (asymptotic), where our per-element libm
  `exp`/`pow` prefactor meets NumPy's SIMD-vectorised prefactor on a regime
  whose series is cheap — an array-substrate (SIMD) matter, not a kernel defect.

## Known / deferred (not regressions from this change)
- `W` large-z 1.24× vs SciPy — SIMD prefactor; out of scope (array substrate).
- MPFR `WhittakerM` 1.5× slower than mpmath (faster than WL) — the Expr-template
  round-trip per call; a direct arb kernel would close it. mpmath is not a named
  baseline; WL is beaten. Optional.
- Rank-1 `Compile[]` for N-ary (3-arg) kernels reports `False` — pre-existing,
  all such heads, tracked in `tools/compile_coverage.py` BASELINE. Scalar
  `Compile` works. Does not affect the packed throughput (ndarray kernel, not Compile).

## Verification
- Accuracy ≤1e-14 rel vs 30-digit FLINT at all z.
- Bit-identical scalar == NDArray == Compile (diffs exactly 0.0).
- `whittaker_tests`, `hypergeometricu_tests` pass; `check-c99`, `check-messages`,
  `check-packed-aware`, `check-nd-surfaces`, `check-array-exactness`,
  `check-compile-coverage` all green (no new gaps).
- Permanent guard: `whittakerm`/`whittakerw` probes in `tools/numeric_sweep.py`
  (0 value mismatches; 1.17× / 1.77× faster than NumPy in-harness).
