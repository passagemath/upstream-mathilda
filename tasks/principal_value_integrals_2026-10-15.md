---
title: Principal-Value definite-integration campaign (50-case battery)
status: in-progress
started: 2026-10-15
corpus: tests/PrincipalValueIntegralsCorpus.m
runner: tests/test_principal_value_corpus.c  (principal_value_corpus_tests)
gate: WRONG <= ratcheting baseline (start 4 -> 0); CORRECT floor ratchets up
baseline_version: v0.406
---

# Principal-Value definite-integration campaign

User supplied 50 Cauchy principal-value definite integrals ("extremeComplexSingular
IntegralStressTest") to stress test Mathilda and guide *fully general, algorithmic*
improvements. Independent verification (analytic Mellin derivation + direct /
symmetric-exclusion quadrature) found the battery's oracle is ~40% unusable:
~16 expected values are wrong, 5 integrals diverge, 3 are convention-dependent
double-pole Hadamard finite parts, 2 are non-elementary (Si/Ci), 2 need an
infinite-pole-sum mechanism. Mathilda is ALREADY correct on ~20 of the 39 scorable
cases.

## Verdict table (TRUE value established independently; "oracle" = user's value)

Legend: ✅ Mathilda correct now · 🐞 Mathilda wrong now · ⬜ uneval/timeout now · ✎ oracle corrected · ⊘ excluded from scoring

| # | integrand (PV, domain) | TRUE value | oracle | Mathilda | notes |
|--|--|--|--|--|--|
|1|x^(a-1)/(x-1), 0∞|-π Cot[πa]|✓|✅||
|2|x^(a-1)/(1-x), 0∞|π Cot[πa]|✓|✅||
|3|x^(a-1)/(x²-1), 0∞|-(π/2)Cot[πa/2]|✓|✅||
|4|x^a/(x²-1), 0∞|+(π/2)Tan[πa/2]|✎ sign|✅|oracle had −|
|5|x^(a-1)/(x²-b²), 0∞|-(π/(2b^(2-a)))Cot[πa/2]|✓|🐞|spurious +Im (b≠1 bug)|
|6|x^a/(x²-b²), 0∞|+(π/2)b^(a-1)Tan[πa/2]|✎ sign|🐞|oracle sign + Im bug|
|7|x^(a-1)/(x⁴-1), 0∞|-(π/4)Cot[πa/4]|✎|✅|oracle had spurious Csc|
|8|x^a/(x⁴-1), 0∞|-(π/4)Cot[(a+1)π/4]|✎|✅|oracle spurious Csc|
|9|x^(a-1)/(x-b), 0∞|-π b^(a-1)Cot[πa]|✓|✅||
|10|x^(a-1)/(b-x), 0∞|π b^(a-1)Cot[πa]|✓|✅||
|11|Log x/(x-1), 0∞|DIVERGES (tail ln x/x)|✗|✅idiv|⊘ divergent|
|12|Log x/(x²-1), 0∞|π²/4|✓|✅||
|13|Log x/(x²-b²), 0∞|π²/(4b)|✓|✅||
|14|x Log x/(x²-1), 0∞|DIVERGES|✗|⬜|⊘ divergent|
|15|x Log x/(x²-b²), 0∞|DIVERGES|✗|⬜|⊘ divergent|
|16|Log²x/(x²-1), 0∞|0|✓|✅||
|17|Log²x/(x²-b²), 0∞|π²Log[b]/(2b)|✎|✅|oracle had 0|
|18|Log x/(x³-1), 0∞|4π²/27|✎|✅|oracle 2π²/(9√3)|
|19|x Log x/(x³-1), 0∞|4π²/27|✎|✅|oracle π²/(9√3)|
|20|Log x/(x⁴-1), 0∞|π²/8|✓|✅||
|21|Tan[ax]/x, 0∞|π/2|✓?|⬜|⊘ hard ∞-pole|
|22|Cot[ax]/x, 0∞|DIVERGES at 0 (~1/x²)|✗|⬜|⊘ divergent|
|23|Sin[ax]/(x²-b²), 0∞|(Cos[ab]Si[ab]-Sin[ab]Ci[ab])/b|✎|⬜|oracle πcos/(2b); non-elem|
|24|Cos[ax]/(x²-b²), 0∞|-π Sin[ab]/(2b)|✓|⬜|**Phase 1**|
|25|x Sin[ax]/(x²-b²), 0∞|π Cos[ab]/2|✓|⬜|**Phase 1**|
|26|x Cos[ax]/(x²-b²), 0∞|-(Cos[ab]Ci[ab]+Sin[ab]Si[ab])|✎|⬜|oracle -πsin/2; non-elem|
|27|Sin[ax]/(x(x²-b²)), 0∞|π(Cos[ab]-1)/(2b²)|✓|⬜|**Phase 1**|
|28|Cos[ax]/(x(x²-b²)), 0∞|DIVERGES at 0|✗|⬜|⊘ divergent|
|29|Tan[ax]/(x²-b²), 0∞|? (∞ poles)|?|⬜|⊘ hard|
|30|Cot[ax]/(x²-b²), 0∞|? (∞ poles)|?|⬜|⊘ hard|
|31|Exp[iax]/(x-b), ℝ|iπ e^{iab}|✓|⬜|**Phase 1**|
|32|Exp[-iax]/(x-b), ℝ|-iπ e^{-iab}|✓|⬜|**Phase 1**|
|33|Exp[iax]/(x²-b²), ℝ|-π Sin[ab]/b|✎|⬜|=2·(24); oracle iπe^{iab}/b|
|34|Exp[-iax]/(x²-b²), ℝ|-π Sin[ab]/b|✎|⬜|=2·(24)|
|35|x Exp[iax]/(x²-b²), ℝ|iπ Cos[ab]|✓|⬜TIMEOUT|=2i·(25); perf hang|
|36|x Exp[-iax]/(x²-b²), ℝ|-iπ Cos[ab]|✓|⬜TIMEOUT|=−2i·(25)|
|37|Sin[ax]/Sinh[bx], ℝ|(π/b)Tanh[πa/2b]|✓|✅||
|38|Cos[ax]/Sinh[bx], ℝ|0|✓|⬜||
|39|x Cos[ax]/Sinh[bx], ℝ|-π²/(2b²Cosh²[πa/2b])|✓?|⬜||
|40|Exp[-ax]/Sinh[bx], 0∞|(1/2b)(ψ[(a+b)/2b]-ψ[a/2b])|✓?|⬜||
|41|x^(a-1)/(x²-1)², 0∞|Hadamard FP (conv-dep)|n/a|⬜|⊘ double pole|
|42|x^a/(x²-1)², 0∞|Hadamard FP|n/a|⬜|⊘ double pole|
|43|x^(a-1)/(x³-1)², 0∞|Hadamard FP|n/a|⬜|⊘ double pole|
|44|x^(a-1)/(x²-2x cos b+1), 0∞|π Sin[(1-a)(π-b)]/(Sin[πa]Sin[b])|✎|✅*|correct but monster form|
|45|x^a/(x²+2x cos b+1), 0∞|π Sin[ab]/(Sin[πa]Sin[b])|✓|✅||
|46|Log x/(x²-2x cos b+1), 0∞|0 (x→1/x symmetry)|✗|✅~0|oracle (bπ-π²)/(2sin b)|
|47|Log x/(x²+2x cos b+1), 0∞|0|✗|✅~0|oracle bπ/(2sin b)|
|48|Log x·x^(a-1)/(x-1), 0∞|π²Csc²[πa]|✎ swap|⬜|oracle had Csc·Cot (that's 49)|
|49|Log x·x^(a-1)/(x+1), 0∞|-π²Csc[πa]Cot[πa]|✎ swap|🐞 ComplexInfinity|oracle had Csc² (that's 48)|
|50|Log²x·x^(a-1)/(x+1), 0∞|π³Csc³[πa](1+Cos²[πa])|✓|🐞 Indeterminate||

Excluded from scoring (11): 11,14,15 (divergent), 22,28 (divergent at 0), 21,29,30
(infinitely many poles — deferred), 41,42,43 (convention-dependent Hadamard FP).
Scored: 39.  Current: ~20 CORRECT, 4 WRONG (5,6,49,50), rest UNEVAL/TIMEOUT.

Generic instantiation: {a->37/100, b->8/10} (case 40: {a->13/10, b->1/2}, needs a>b>0).

## Phased plan (leverage-ordered general fixes)
- **Phase 0** — this corpus + runner + baseline.  [in progress]
- **Phase 1** — real-axis PV half-residue in the residue Fourier/rational families
  (`integrate_residue.c` fourier @738, fourier_sum @909, rational @615). Unlocks
  24,25,27,31,32,33,34,35,36 (and 23,26 if Si/Ci emitted). Fix 35/36 TIMEOUT.
- **Phase 2** — b≠1 keyhole-Mellin PV spurious imaginary (`mellin_core` @1922-1950).
  Fixes 5,6 (WRONG→CORRECT).
- **Phase 3** — (deferred/optional) order-2 axis pole / Hadamard FP — 41,42,43 excluded.
- **Phase 4** — Log×Mellin s=integer limit (`ratlogpow`): 48,49,50 (49,50 WRONG→CORRECT, 48 UNEVAL→CORRECT).
- **Phase 5** — simplify the quadratic-denominator Mellin monster — 44 (cosmetic).

## Progress log
| Version | Cases closed | Mechanism |
|--|--|--|
| v0.406 (baseline) | — | 20/39 correct, 4 wrong |
| v0.410 | 24,25,31,32,33,34,35,36 | real-axis PV half-residue in residue Fourier family (`integrate_residue.c`) |
| v0.410 | 5,6 | b≠1 keyhole-Mellin PV real projection (`mellin_core`) |
| v0.410 | 49,50 | Log-weighted Mellin `x^(c-1)Log^k/(1±x^d)` (`integrate_mellin_logpow_try`); keyhole-log declines symbolic exponent |

**Status after v0.410: 32/39 correct, 0 wrong** (corpus `principal_value_corpus_tests`
baseline ratcheted 4 → 0). Regression corpora green: extreme 48/0, singular 44/0,
mse 15/0, residue_tests pass.

Still open (UNEVALUATED, not wrong): 23, 26 (non-elementary Si/Ci targets — need
`SinIntegral`/`CosIntegral` emission); 27 (removable+PV mixed, declined); 38, 39, 40
(whole-line Sinh-denominator extensions); 48 (PV `Log·x^(a-1)/(x-1)` — the handler is
correct but its `prove_true` strip check intermittently fails after the PV residue
stage leaves global assumption state; a false *decline*, never wrong). 48's near-miss
and the residue-stage assumption-state leak are the natural next items.

**Deferred git step (concurrent tree):** main advanced v0.406→v0.409 during the work,
so the version bump to v0.410 is in the working tree but NOT committed/tagged — commit
on a branch and resolve any `version.h` collision against current main, then
`git tag v0.410`.

## Lessons
- The battery oracle was ~40% unusable; every "FAIL" had to be adjudicated
  independently (analytic Mellin + quadrature) before touching the CAS — most were
  Mathilda-correct against a wrong oracle. See `[[feedback_verify_reference_values_never_recall_them]]`.
- The residue series engine mishandles a **numeric-frequency kernel at a symbolic
  pole** (`Exp[I x]/(x²-b²)`: the negated pole's `e^{iz₀}` collapses to 1). Gated
  around in `residue_family_fourier` (`res_has_param`); a latent bug worth a real fix.
- `prove_true`/`Simplify` convergence checks are perturbed by global assumption state
  the residue PV stage leaves behind — a handler that runs after residue can see its
  strip proof flip. Prefer a self-contained proof or run before the state is set.
