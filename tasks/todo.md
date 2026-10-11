# Assumed-Hyperbolic Simplify Stress Corpus

Parallel to the trig campaign (v0.332–v0.338, 22→49/50). Build a 50-case
`Simplify[diff, assumptions] → 0` hyperbolic corpus + the C changes to prove it.
Plan: `/Users/user/.claude/plans/for-each-of-the-mighty-dove.md`.

## Phase 0 — construct + baseline (no commit)
- [ ] Write corpus `.m` to scratchpad ({id, asm, expr} + Do scorer)
- [ ] Numeric-verify all 50 identities in-region (reject wrong constructions)
- [ ] Run `./Mathilda` baseline; record pass/fail per case

## Phase A — M1 hyperbolic sign oracle (simp_assume.c)  [commit+tag]
- [ ] prov_pos: Sinh/Tanh/Coth/Csch > 0 when arg > 0; Sech > 0 for real arg
- [ ] prov_neg: Sinh/Tanh/Coth/Csch < 0 when arg < 0
- [ ] prov_nn: Sech >= 0 for real arg
- [ ] Re-score corpus (expect G2/G5 sqrt/abs cases)

## Phase B — M4 reciprocal-Pythag radicand prep (simp_logexp_abs.c)  [commit+tag]
- [ ] Add `1 - Sech^2 :> Tanh^2`, `-1 + Coth^2 :> Csch^2`
- [ ] Re-score (#17, #18, #45)

## Phase C — M5 inverse-of-forward (simp_assume_rewrite.c)  [commit+tag]
- [ ] ArcCsch[Csch]->e, ArcCoth[Coth]->e, ArcSech[Sech]->±e
- [ ] Re-score (#50)

## Phase D — G4 imaginary-period symbolic-integer shift  [commit+tag]
- [ ] New hyperbolic-period helper mirroring trig_pi_shift_symbolic
- [ ] Re-score (#34–#39; confirm #31–33/#40)

## Phase E — G3 residue (measure first)  [commit+tag if needed]
- [ ] Check #21–30 via existing roundtrip; #20 via TrigToExp
- [ ] Add M7 hyperbolic multiangle candidate only if gap remains

## Wrap-up
- [ ] New `tests/test_hyp_assume_simplify.c` + CMakeLists registration
- [ ] Soundness decline tests
- [ ] No trig regression (logexp/invtrig/trigexp_zero tests)
- [ ] Docs: changelog 2026-10-05.md + simplification.md
- [ ] version.h bumps + tags; make check-messages / check-c99 / full make
- [ ] Memory update; rebuild code-review graph

## Review (v0.346)

Built the assumed-hyperbolic Simplify stress corpus, sibling of the assumed-trig
campaign. Baseline **30/50 → 48/50**; all 50 numerically verified as true
identities first. New `tests/test_hyp_assume_simplify.c` (48 collapse asserts + 5
soundness declines), registered in CMakeLists.

Mechanisms (all sound, each with a decline test; verified off-identity behaviour):
- M1 sign oracle (`simp_assume.c`): `sign(Sinh/Tanh/Coth/Csch[u])=sign(u)`, `Sech>0`;
  + fixed a latent `prov_re` gap (a `Rational` literal was not recognised real, so
  `prov_re[x/2]` failed — closed #16). Unlocked Groups 2 & 5 sqrt/abs.
- M4 radicand prep (`simp_logexp_abs.c`): `1 ± Cosh[u] → 2{Cosh,Sinh}[u/2]²`.
- M5 inverse-of-forward (`simp_assume_rewrite.c`): ArcSech/ArcCsch/ArcCoth.
- Imaginary-period shift (`simp_assume_rewrite.c`): `H[t + m Pi I]`, symbolic int m;
  correctly rejects real-Pi shifts.

Residuals (true identities, documented, not hyperbolic-specific): #26, #29 — need
radical-gradient / inverse-hyperbolic constant zero-testing (as trig left #20).

Verification: all 4 suites green (logexp/invtrig/trigexp_zero/hyp_assume, 0 FAIL),
corpus 48/50, `make check-c99` (exit 0), `make check-messages` (OK), full build
clean. `$VersionNumber` 0.346.

