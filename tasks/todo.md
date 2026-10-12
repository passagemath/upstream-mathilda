# TODO: Integrate / PossibleZeroQ no spurious `Power::infy` (v0.362)

Plan: `~/.claude/plans/the-following-should-not-humble-hammock.md`

Bug: `Integrate[1/Sqrt[Coth[x]], x]` returns the correct antiderivative but leaks
two `Power::infy: Infinite expression 1/0 encountered.` lines. Not cosmetic —
`Check[Integrate[1/Sqrt[Coth[x]], x], MARKER]` returned `MARKER` (failure branch),
because the leaked diagnostic was *noted*. Root cause: `PossibleZeroQ`'s numeric
sampler (`src/zero_test.c`) did not mute its throwaway samples; close the whole
class of un-muted numeric probes in the Integrate pipeline.

## Fix
- [x] Tier 1 (root): `src/zero_test.c` — wrap `evaluate_rung` (per-sample numeric
      eval, the choke point for every PZQ sample) in `arith_warnings_mute_push/pop`
      via a single-exit wrapper (`evaluate_rung_core`).
- [x] Tier 2 (bespoke zero-tests): `integrate_goursat.c` `diff_back_ok`,
      `sample_clearly_nonzero`.
- [x] Tier 3 (bespoke sign/branch probes): `intsimp.c` `intsimp_numeric_sign` +
      `it_num_real` (covers the 6 `intrat.c` sign sites); goursat
      `involution_score` / `radical_sign_negative` / `numeric_sign_ratio` /
      `monomialize_root_of_unity`; `integrate_quadrad.c` `qr_sign_probe`;
      `integrate_jeffrey.c` `const_sign`.
- [x] Tier 4: regression `test_coth_no_infy_leak` in
      `tests/test_integrate_derivdivides.c` (asserts the integral closes AND
      `mth_msg_fired_count()` is unchanged — the direct Check-corruption guard).

## Verification
- [x] `make` clean (gcc-16, `-std=c99 -Wall -Wextra -Werror=...`), warning-free.
- [x] `Integrate[1/Sqrt[Coth[x]], x]` → 0 leaks (was 2); `Check[…, MARKER]` →
      antiderivative (was `MARKER`); `PossibleZeroQ[1/x]` → False, 0 leaks.
- [x] Answers byte-identical across the probe set (Coth, Tanh, Tan, Sqrt[(1+x)/x],
      Sqrt[x+Sqrt[x]], 1/(1-x^2), 1/(1+x^4), Sech, ArcTanh[x], Sqrt[Coth]).
- [x] No genuine diagnostic hidden: `1/0`, `1/(x-x)` still warn; `Check[1/0,
      FAILED]` still returns `FAILED`.
- [x] `make check-messages` (BASELINE empty), `make check-c99` (exit 0).
- [x] Affected + high-fanout suites: integrate_derivdivides (incl. new test),
      goursat, quadrad, jeffrey, zero_test, integrals, possiblezeroq_stress — pass.
- [ ] Full `tests/` suite run (in progress — zero_test.c is high-fanout).
- [x] Version bump `src/version.h` 0.361 → 0.362; changelog
      `docs/spec/changelog/2026-10-12.md`.
- [ ] Commit + tag `v0.362` + push `--follow-tags` (after full suite passes).
- [ ] Rebuild code-review-graph (SessionStart hook).

## Review
Root cause was the shared `PossibleZeroQ` sampler, not the Integrate call sites —
confirmed by a bare `PossibleZeroQ[…]` leaking and by Tier 1 alone dropping the
reported integral from 2 leaks to 0. The `evaluate_rung` mute is behaviour-
preserving (a pole sample is classified decisive-FALSE from its `Infinity` value,
never from message state), so no verdict changes. The remaining tiers close the
same latent class (bespoke `N[…]` probes that bypass PZQ) so other integrands can't
reproduce it. Architectural note: Goursat's `diff_back_ok` bypasses PZQ *by design*
(PZQ misfires/too-slow on its output); it was muted, not reconciled with PZQ —
reconciliation is a separate, riskier decision.
