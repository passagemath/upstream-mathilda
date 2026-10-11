# TODO: Integer-assumption Simplify campaign

Plan: `~/.claude/plans/here-s-some-problems-to-floating-milner.md`
Baseline: 51/100 on v0.346. Target ~90/100 (Phases 1–4).

## Phase 0 — harness (folded into per-phase; asserts active by default so each
## commit's test must be green → assert only cases that pass at that commit)
- Measurement throughout: `./Mathilda -file` heredoc over all 100 (non-aborting).
- `test_integer_assume_simplify.c` created in P1 with Cat A/B + declines + #29/#88
  non-identity asserts + #94; later phases APPEND their green cases.
- CMake: add_executable + include dirs + build_ecm dep + add_test.

## Phase 1 — root-of-unity / parity (Cat A, B, #94)  ✅ DONE v0.347 (51->64/100)
- [x] `simp_reduce_mod2` Boolean-ring primitive + `canon_neg1_pow` + Mod-absorption
- [x] New `Power[-1|I|-I|E^(IPi·)]` block in `assume_structural_rewrite`
- [x] Supersede even-only `prov_even` fold; added base-I/-I folding
- [x] De-restrict `trig_pi_shift_symbolic` (emit symbolic `(-1)^K`); prov_int(Mod)
- [x] test_integer_assume_simplify.c + CMake; declines; #94 Arg bug fixed
- [x] no regressions (logexp/invtrig/hyp/refine/element); valgrind clean; gates green

## Phase 2 — Cat C modular/rounding  ✅ DONE v0.348 (64->78/100)
- [x] (3a) confluence post-pass in builtin_simplify (fixpoint, gated to rounding heads)
- [x] (3b) integer-shift extraction, floor-division (Expand), Mod[P,2] reduce, Mod divisibility
- [x] tests + declines; no regressions (fixed logexp (ab)^c re-contraction); valgrind clean
- residual: #54 Mod[n,2]+Mod[n+1,2]=1 (parity complement), not asserted

## Phase 3 — Gamma/Factorial (Cat D)  ✅ DONE v0.349 (78->85/100)
- [x] Gamma/Pochhammer/Binomial->Factorial conversion in simp_factorial (ctx-gated)
      + Factorial->Gamma back-conversion; shape/gate widened to factorial-like
- [x] Factorial2[2n]=2^n n! in assume_structural_rewrite (even, prov_int-gated)
- [x] Infinity::indet gone (symbolic reduction fires before the pole sampler)
- [x] ctx-gating fixed the Gosper/Sum regression (found via expand_tests A/B)
- [x] tests + no regressions (gamma/pochhammer/fullsimplify/expand/sum); valgrind clean

## Phase 4 — LegendreP + GCD/LCM (Cat E-partial, F)  ✅ DONE v0.350 (85->90/100)
- [x] LegendreP[n,-1]=(-1)^n (assumption path, prov_int-gated)
- [x] GCD[a,a]=Abs, GCD[a,1]=1, LCM[a,1]=Abs, GCD[n,n+1]=1 (consecutive coprime)
- [x] post-pass gate widened to GCD/LCM/LegendreP; tests + soundness; valgrind clean

## CAMPAIGN COMPLETE: 51 -> 90/100 (target met). Remaining 10:
##   #29,#88 non-identities (Mathilda correct); #54 Mod parity-complement residual;
##   #76-81 ChebyshevT/U (unimplemented functions, separate feature); #87 binomial sum.
## Out of scope (per plan): ChebyshevT/U, binomial sum #87.

## Review (campaign complete)
Result: assumed-integer `Simplify` suite **51 → 90/100** across v0.347–v0.350,
four tagged commits. One elegant core primitive (`simp_reduce_mod2`, Boolean
ring mod 2) drives the root-of-unity work; the rest reuses existing engines
(`simp_factorial` ratio cancellation, `trig_pi_shift_symbolic`) by routing into
them. Each mechanism declines off the integer assumption (soundness tests in
`tests/test_integer_assume_simplify.c`). Zero regressions across 17 suites;
valgrind-clean; `check-c99`/`check-messages` green.

Key lessons (also in harness memory):
- A broad Gamma/Binomial→Factorial conversion in general `Simplify` silently
  corrupts the Sum engine's (unassumed) internal Simplify (Gosper) → gate it on
  an active integer-assumption `ctx`. Found via an `expand_tests` git-stash A/B.
- `Simplify` drops bare structural assumption rewrites that `Refine` surfaces →
  a confluence post-pass in `builtin_simplify` (gated to the relevant heads).

Not done (deliberate / out of scope): ChebyshevT/U (#76–81, unimplemented
functions), binomial sum #87 (Sum subsystem), #54 residual, #29/#88 non-identities.
Tags are LOCAL (not pushed) — push with `git push --follow-tags` when ready.
