# Assumption-Dependent Definite-Integral Stress Campaign (2026-10-13)

Plan: `~/.claude/plans/let-s-thoroughly-test-the-wild-harbor.md`
Findings: `tasks/assumed_definite_integrals_2026-10-13.md`
Baseline (v0.376, per-case isolated): ASSUM 35 PASS / 5 WRONG / 6 UNEVAL / 4 ABORTED.
40/50 values correct; 0 true wrong answers. Target: all 50 ASSUM-PASS, clean forms,
no uneval/abort/crash.

## Phase 0 — Diagnostic battery  [DONE]
- [x] per-case isolated harness (alldefs.m); both ASSUM + LITERAL metrics
- [x] categorise all 50; findings file written

## Theme 1 — Result cleanliness under Assumptions (Mathematica parity)
- [ ] 1a. driver post-simplify: integrate_definite -> Simplify[r, assumptions] (gated). Fixes 12-17,10,31,48,49
- [ ] 1b. ConditionalExpression discharge under entailed assumptions (8,9,18)
- [ ] 1c. nested-radical Sqrt[1+p/q] pull-out under positivity (8,9,37)
- [ ] 1d. Beta[x,y] -> Gamma expansion (18)
- [ ] 1e. trig periodicity Csc/Sec/Sin/Cos[Pi(k+a)], integer k (26)

## Theme 2 — Genuine integration gaps
- [x] 2a. piecewise handler integrate_piecewise.c (43-47)  [DONE v0.377]
      New Integrate`Piecewise (Abs/Sign/HeavisideTheta/UnitStep/Min/Max/Piecewise);
      breakpoint split, assumption-aware branch select, recurse+sum. Internal Solve
      Quiet-wrapped. ASSUM 35->40 PASS. Edge cases verified (bp-outside, 2-bp,
      product, undecidable declines). check-c99 + check-messages PASS. No regression:
      ramanujan/beta/intrep/newton_leibniz/symmetry suites clean; no existing integrate
      test uses branch heads. RELEASE (commit/tag v0.377) AWAITING USER GO.
- [ ] 2b. improper endpoint / ArcSin: 1/Sqrt[a^2-x^2] on [0,a] (27)
- [ ] 2c. sign reflection: Sin[a x]/x, a<0 (20)
- [ ] 2d. log-rational half-line: Log[x]/(x^2+a^2) (32)
- [ ] 2e. Poisson log on [0,Pi] (38,39)

## Theme 3 — Robustness
- [ ] characterise the TimeConstrained-accumulation segfault (likely pre-existing)

## Verification / release (per change)
- [ ] re-run full battery (both metrics, per-case) after each change — no regression
- [ ] tests/AssumedDefiniteIntegralsCorpus.m + runner (model test_crc_corpus.c)
- [ ] touched unit suites; check-messages; check-c99; valgrind new C handlers
- [ ] docs/spec + weekly changelog; version bump + tag per substantive commit
- [ ] update findings file + tasks/lessons.md on corrections
