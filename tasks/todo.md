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

## Phase 2 — Cat C modular/rounding
- [ ] (3a) root-cause + fix Simplify↔Refine structural-rewrite surfacing
- [ ] (3b) integer-shift extraction, floor-division, Mod[P,m] via reducer, Mod divisibility
- [ ] tests; build; bump+tag; docs

## Phase 3 — Gamma/Factorial (Cat D)
- [ ] Gamma↔Factorial normalizer routing through simp_factorial
- [ ] Fix zero_test integer-sampler Gamma-pole exclusion (Infinity::indet)
- [ ] tests; build; bump+tag; docs

## Phase 4 — LegendreP + GCD/LCM (Cat E-partial, F)
- [ ] LegendreP[n,-1], [n,0] special values
- [ ] GCD/LCM symbolic rules under integer assumptions
- [ ] tests; build; bump+tag; docs

## Out of scope: ChebyshevT/U (#76-81), binomial sum #87.
