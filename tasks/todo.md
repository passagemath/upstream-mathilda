# LaguerreL — todo

Implement `LaguerreL[n,x]` and generalized `LaguerreL[n,a,x]` to the ChebyshevU
bar (v0.373 → v0.374). Plan: `~/.claude/plans/pasted-content-id-8ac5-*.md`.

## Step 1 — machine helpers (sf_machine.{h,c})  [DONE]
- [x] `sf_machine_laguerre_l(n,x,*out)` + `_c` (a=0 recurrence)
- [x] `sf_machine_laguerre_l_assoc(n,a,x,*out)` + `_c` (generalized)

## Step 2 — module (src/special_functions/laguerre_l.{c,h})  [DONE]
- [x] shape predicates / mpq poly assembly / MPFR result helpers (lag_*)
- [x] explicit polynomial: mpq recurrence (numeric/absent a); Expr+Expand (symbolic a)
- [x] numeric dispatch: machine-gate → sf_machine; else MPFR ncpx recurrence
- [x] non-integer/complex order → Hypergeometric1F1 template
- [x] special values (L_0=1, L_n(0)=Binomial[n+a,n])
- [x] builtin_laguerre_l (2|3 args, builtin_arg_error), laguerre_l_init (attrs)
- [x] C99: fallbacks fine (no M_* consts used), int64 N/A

## Step 3 — ND/Compile kernels (ndkernels.c)  [DONE]
- [x] NDKB_LaguerreL (binary, hand-written complex arm) + REG_B
- [x] NDKN_LaguerreL (n-ary, nargs=3) + REG_N

## Step 4 — glue  [DONE]
- [x] SYM_LaguerreL (sym_names.h decl + sym_names.c def + intern)
- [x] core.c laguerre_l_init() + include
- [x] info.c docstring

## Step 5 — derivative (src/calculus/deriv.c)  [DONE]
- [x] n==2: -LaguerreL[N-1,1,Z]·dZ + inert Derivative[1,0]
- [x] n==3: -LaguerreL[N-1,A+1,Z]·dZ + inert Derivative terms

## Step 6 — docs/tests/version  [DONE]
- [x] tests/test_laguerre_l.c + CMakeLists.txt (incl. source in mathilda_common!)
- [x] docs/spec/builtins/special-functions.md (## LaguerreL)
- [x] docs/spec/changelog/2026-10-12.md (prepend ## LaguerreL ; v0.374)
- [x] src/version.h → 0.374 / "0.374"

## Verify
- [x] make -j (clean build, laguerre_l.o linked)
- [x] check-c99 PASS; check-messages PASS
- [x] behavioral parity — ALL spec examples match (incl. cancellation 1.45553 vs 1.45577)
- [x] multi-surface identity: scalar==List==NDArray==Compile bit-for-bit (delta 0.0);
      CompileDiagnostics True for 2-arg scalar+array, 3-arg scalar
- [x] check-array-exactness PASS; check-compile-coverage PASS (no new gaps)
- [~] check-packed-aware: fails only on PRE-EXISTING StandardForm (commit 4a66b3fe,
      not mine); LaguerreL NOT flagged
- [x] check-nd-surfaces PASS (no LaguerreL disagreement)
- [~] check-fastpath-sweep: fails on a PRE-EXISTING 27-head stale-baseline drift
      (graph/string/stats predicates, Xor, Refine, Predict, distances; OFF_BUFFER
      "last re-recorded 2026-08-03"); LaguerreL NOT flagged
- [x] laguerre_l_tests PASS (0 FAIL lines; fixed my 17/2→11/2 reference slip)
- [x] valgrind CLEAN (identical leak totals to baseline; no LaguerreL frames)
- [x] sibling tests PASS (chebyshev_u/t, legendre) — additive changes don't regress
- [x] rebuilt code-review graph (incremental, 13 files, no errors)

## Pre-existing red gates (NOT mine, out of scope)
- check-packed-aware: StandardForm (commit 4a66b3fe / v0.367, print.c)
- check-fastpath-sweep: 27-head stale OFF_BUFFER baseline (2026-08-03)

## Release (awaiting user go — outward-facing)
- [ ] commit (; v0.374), tag v0.374, push --follow-tags
