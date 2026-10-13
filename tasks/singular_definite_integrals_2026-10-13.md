---
title: Singular / FTC-boundary definite-integral stress campaign
date_started: 2026-10-13
status: in-progress
plan: ~/.claude/plans/pasted-content-id-a30c-this-50-expressi-dazzling-sunset.md
---

# Singular / FTC-boundary definite-integral campaign (50 cases)

Battery: `assumedSingularIntegralStressTest` — 50 definite integrals at the
Newton–Leibniz boundaries: integrable endpoint singularities, in-domain
kinks (Abs/Sign/CubeRoot/Min/Max), **Cauchy principal values**, branch cuts in
the antiderivative (trig-rational over a period, Poisson kernels), and
essential/oscillatory endpoints. Scope (user): **FULL — close all 50.**

Corpus: `tests/SingularDefiniteIntegralsCorpus.m` (single source of truth).
Committed gate: `tests/test_singular_definite_corpus.c`
(`singular_definite_corpus_tests`), fork-per-case, WRONG baseline ratchets down.
Fast driver: `scratchpad/run_corpus.sh` (one process per case, `-file` mode).

## TWO BATTERY ORACLE ERRORS (verified by NIntegrate — do NOT re-fix Mathilda to them)
- **15** `Integrate[Abs[x^2-a^2],{x,-2a,2a}]` = **4a^3** (battery said 6a^3; NInt a=2 → 32).
  Mathilda already returns 4a^3 — correct.
- **41** `Integrate[Sin[a/x]/x^2,{x,a,Inf}]` = **(1-Cos[1])/a** (battery sign-flipped;
  NInt a=3 → +0.1532).

## Baseline (v0.381, corrected oracles, one process per case)
20 CORRECT, 4 WRONG (finite ≠ expected), 23 UNEVAL, 3 TIMEOUT.
- CORRECT incl. the value-correct-but-unreduced 1,2,3 (Beta form) and 33,34 (radical), + 15.
- **WRONG (soundness, Phase-1 targets): 14, 18, 24, 25** (finite wrong) **+ 12** (garbage with
  free x — diff-eval hangs → shows as TIMEOUT; same root bug).
- TIMEOUT: 12 (diff-eval), 48, 49 (`Log^n x/(a^2+x^2)` half-line hang).

## Architecture anchors (from exploration)
- Dispatch `integrate.c:1762`; definite `integrate_definite` `integrate.c:1233`; options
  parsed inline (`assumptions` `:1256`, `principal_value` by string `:1268`); driver post-Refine
  `:1422-1434`.
- Cascade `integrate.c:1316-1416`: a==b → dirac → piecewise(`:1325`) → residue(`:1333`) →
  symmetry(`:1353`) → beta/euler_2f1(`:1370`) → trigpower(`:1375`) → ramanujan(`:1382`) →
  intrep(`:1391`) → NewtonLeibniz(`:1397`) → sinpowmono/ratlogpow → diffunderint(`:1416`).
- Newton–Leibniz `integrate_newton_leibniz.c`: `_try_pv`(`:566`), `nl_collect_poles`(`:313`),
  `nl_eval_at`(`:405`)/`nl_limit`(`:388`), idiv `:656`, crosscheck `nl_crosscheck_ok`(`:454`),
  PV engine `nl_principal_value`(`:530`, finite-numeric-bounds gate `:533`).
  **No `assumptions` param** — symbolic bounds bail `NL_UNKNOWN`(`:237`).
- PV wired only into residue + N-L; whole-line PV absent
  (`integrate_residue.c:644-649`, `807-817`).
- Piecewise `integrate_piecewise.c` (`has_target:224`, `solve_roots:146`); trig-rational
  residue-on-unit-circle `integrate_residue.c:1099/1502`; Beta `integrate_beta.c:207` (returns
  `Beta[...]`); Fresnel/osc/Dirichlet `integrate_ramanujan.c` (`rec_trig:386`,`sinpowmono:2143`);
  log-rational complex poles via Mellin binomial `mellin_term:1786`.

## Gaps → phases (see plan)
- **P1 soundness**: √(g²)→drop-Abs (18); Tan/Cot through interior pole → finite complex −I Pi
  (24,25, even w/o PV); piecewise neg-power branch complex artifact (6,14) + symbolic-bound
  garbage (12); 48/49 hang.
- **P2 keystones**: assumption-aware Newton–Leibniz (unlocks 4,40 + symbolic-bound PV);
  PrincipalValue first-class across 21–30 incl. whole-line (30).
- **P3 recognizers**: Log[Sin]/Log[Cos] (9,10); Poisson log (36,37); reciprocal u=1/x (41–43);
  King-reflection (35); (a²−x²)^s (4,5,11,40,16); Tan²/Sin² (38,39).
- **P4 form**: Beta[k,1]→1/k (1,2,3); Sqrt[(1−a²)²] reduction (33,34).

## Progress (live)
| Version | Fix | Cases | Status |
|---|---|---|---|
| — | Phase 0: corpus + committed gate + driver | — | DONE (contributor-facing, no bump) |
| (pending bump) | Phase 1 soundness — eliminate all 5 wrong answers | 14,18,24,25,12 | DONE, verifying |

**Phase 1 soundness — WRONG 5 → 0** (corpus: CORRECT 20→24, UNEVAL 23→24, TIMEOUT 3→2).
Changes:
- `integrate_newton_leibniz.c`: (a) **radical branch-hazard net** — force the NIntegrate
  crosscheck when F carries a radical whose radicand vanishes inside (a,b)
  (`nl_find_radical_hazard`); (b) **trig-pole detection** — `Tan/Cot/Sec/Csc` poles are
  invisible to `Denominator[Together]` (atomic head; `Sin/Cos` re-canonicalizes to `Tan`),
  so collect `Cos[u]==0 / Sin[u]==0` zeros directly (`nl_collect_trig_poles`), with periodic
  `ConditionalExpression` family enumeration (`nl_enum_periodic_pole`); fixes Tan/Cot PV = 0
  (24,25) AND the non-PV `-I Pi` soundness (→ idiv); (c) **complex-value net** — a real-looking
  integrand with a provably-nonzero `Im` is a branch-crossing artifact → idiv.
- `integrate_piecewise.c`: (a) **`Sqrt[perfect square]→Abs` normalization** via
  `Simplify[Sqrt[g], Element[x,Reals]]` so `Sqrt[1-Cos[2x]]→√2 Abs[Sin x]` reaches the Abs
  machinery (18); (b) **periodic breakpoint enumeration** so `Abs[Sin x]` on `[0,2π]` splits at
  `π` (`enum_periodic_root`); (c) **reflect a nonpositive sub-interval** `x→-x` in
  `integrate_piece` so `(-x)^a` stays real instead of the complex `e^{iπa}x^a` branch (6,14).
- `integrate.c`: **dispatcher guard** — reject a real-axis definite result still containing the
  integration variable (leaked antiderivative), always wrong → unevaluated (12 → sound UNEVAL).
- Committed gate `CORPUS_WRONG_BASELINE` ratcheted 5 → 0 (assert-empty).

Known residual message leaks (value correct; Check[] cleanliness, to tidy):
`Integrate::idiv` under a successful Cot PV (25); `Integrate::nonelem` under case 12's decline.
Known PRE-EXISTING soundness hole (not a corpus case, not a regression): `Sqrt[2-2Cos x]` on
[0,2π] → 0 (Simplify can't crack the x/2 half-angle; radicand zero only at endpoints so the
N-L radical-hazard net does not fire).

## Phase 2 complete (v0.383 + v0.384) — corpus 36/50 correct, 0 wrong

**2a (v0.383) assumption-aware Newton–Leibniz** — `assumptions` threaded into pole
classification (`nl_classify_symbolic`, Refine-proved) and boundary eval (Limit + Refine);
unlocked 4, 40 AND the reciprocal `u=1/x` family 41, 42, 43 (via assumption-aware boundary
limits). Corpus 24→31.
**2b (v0.384) PrincipalValue with symbolic bounds** — `nl_principal_value` accepts symbolic
bounds (local odd-order probe + assumption-aware endpoints); numeric realness gate + literal-Re
fallback; `DiffUnderInt` skipped under PV; driver Refine guarded against INTRODUCING complexity
(it was stripping a realness-protecting `Re[CosIntegral[a]-CosIntegral[-a]]`). Unlocked 21, 22,
27, 28, 29. Corpus 31→36.
Regression: all 19 integration/limit/refine + CRC + intrat corpus suites green (0 new fails).

Remaining open (12 uneval + 2 timeout):
- **Phase 3 recognizers**: Log[Sin]/Log[Cos] (9,10); Poisson log (36,37); King-reflection (35);
  Tan²/Sin² over period with Solve::ifun (38,39); CubeRoot odd (16).
- **Symbolic |x|^-b** (12): reflected piece ∫₀^a x^(-b) still uneval under symbolic a.
- **Hard PV**: symbolic pole x=a (23,26); whole-line PV `Cos[x]/(x²-a²)` (30, residue
  indented contour).
- **48,49 residue hang** (pre-existing; Ramanujan solves instantly).

## Phase 1 release-readiness (v0.382) — AWAITING USER GO
Working tree ready to commit + tag `v0.382`:
- `src/calculus/integrate_newton_leibniz.c`, `integrate_piecewise.c`, `integrate.c`,
  `src/version.h` (0.381→0.382), `docs/spec/changelog/2026-10-12.md`,
  `docs/spec/builtins/calculus.md`, `tests/SingularDefiniteIntegralsCorpus.m`,
  `tests/test_singular_definite_corpus.c`, `tests/CMakeLists.txt`.
- Verification: corpus WRONG 0 (24 correct); 16 integration/limit/refine suites green;
  check-messages + check-c99 pass; valgrind shows **no leak in any new function** (reported
  bytes are Darwin ObjC/Accelerate false positives + the pre-existing startup baseline).
- NOT yet committed/tagged (per the await-user-go release pattern).

## Soundness invariant
Zero true wrong answers is the bar. Re-run the full corpus (committed gate + driver) after
each change; WRONG must only go down. Phase-2 core-FTC changes additionally gate on the full
integrate/limit/refine/zerotest suites + prior-campaign corpus A/B.
