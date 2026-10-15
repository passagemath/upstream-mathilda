# DSolve — next 500 of Nasser's corpus (§2.2.41–§2.2.45)

Plan: `~/.claude/plans/let-s-continue-our-implementation-sleepy-rain.md`
Start version: 0.411 (after M67 / §2.2.40 / problem 4000).

Per block: fetch → convert → indvar-audit → register gate+scoreboard → measure
baseline → diagnose dominant gap → general root-cause fix → re-measure + lower
baseline → regression A/B → track (STATUS/README/reports/changelog/PLAN) →
version bump + git tag.

## Block 1 — §2.2.41 (problems 4001–4100), milestone M69, target v0.412
- [x] Fetch `Ch2.S1.SS41.htm`
- [x] Convert → `DE_examples_2241.m` (label 2.2.41) — 100 scalar, 6 IVP
- [x] `make check-corpus-indvar` green + eyeball spot-check
- [x] Add `dsolve_corpus_2_2_41_tests` to CMakeLists (baseline 25) + STATUS/README blocks
- [x] Measure honest baseline → 71/100 (29 non-PASS); dominant gap = regular-singular 2nd-order linear
- [x] Land general fix: confluent-Whittaker recogniser in normal-form pre-pass → 75/100
- [x] No masked-wrong-answer (fully 1F1-verified; degenerate 2μ∈ℤ correctly declines)
- [x] Regression A/B on 7 prior sections: 0 regress, §2.2.20 −4 bonus; specialform unit-test exprs byte-identical HEAD
- [x] version bump 0.411→0.412; changelog + PLAN M69 + report written
- [x] gates: §2.2.41 ctest 25≤25 PASS; check-c99 rc=0; check-messages OK; dsolve_tests aborts on PRE-EXISTING exact-ODE test (byte-identical HEAD, not mine)
- [x] commit + tag v0.412

## Block 2 — §2.2.42 (4101–4200), M70
- [ ] (same pipeline)

## Block 3 — §2.2.43 (4201–4300), M71
- [ ] (same pipeline)

## Block 4 — §2.2.44 (4301–4400), M72
- [ ] (same pipeline)

## Block 5 — §2.2.45 (4401–4500), M73
- [ ] (same pipeline)

## Review
(to be filled at the end)
