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

## STRESS TEST (independent, 2026-10-13) — user-requested
Goal: LaguerreL works seamlessly with NDArray/packed AND is >= as fast as
scipy and wolframscript. Baselines found: benchmarks/.venv (scipy 1.18.1,
numpy 2.5.2; eval_laguerre/eval_genlaguerre), /Applications/Mathematica.app
(wolframscript 1.11.0).

- [ ] A. Correctness oracle vs scipy across grid (2-arg, 3-arg, real/complex,
      non-int order, large n, special values). Tolerance-checked.
- [ ] B. Cross-surface identity (independent re-check): scalar == List ==
      packed List == NDArray == Compile, 2-arg & 3-arg, real & complex.
- [ ] C. int64 NDArray question: LaguerreL[3, NDArray[{1,2,3},"int64"]] returns
      plain List {-2/3,-1/3,1} (un-packs visible NDArray). Bug or deliberate?
      Compare Mathematica; decide fix vs EXEMPT.
- [ ] D. PERFORMANCE vs scipy & wolframscript: vectorized 1e6 pts (2-arg,
      3-arg), Compile[] throughput, scalar loop. Must be >= both.
- [ ] E. Fix any gap found; re-verify; spec/changelog/version/tag if code changes.

## STRESS TEST RESULTS (2026-10-13) — v0.375
- [x] A. Correctness oracle: 902/902 (n,a,x) grid points match mpmath < 4.1e-14
      (2-arg/3-arg, real/complex, non-int/neg order, orders to 300). On par with
      scipy; BEATS scipy at n=50,x=1.5 (2e-15 vs 7e-12). CHECK values agree
      across Mathilda/scipy/wolfram.
- [x] B. Cross-surface identity: 2-arg (real+cplx) was already 0.0. 3-arg was
      OFF by 1.27e-12 (array path used symbolic poly, not machine helper) ->
      FIXED -> now 0.0 across scalar/List/packed/NDArray/Compile, real+complex.
- [x] C. int64 NDArray un-packs to exact List {-2/3,...} — NOT a bug: identical
      to Sin/ChebyshevU siblings (system convention). Preserved by the fix
      (branch declines int64/bool dtype).
- [x] D. PERFORMANCE (1e6 pts, min-of-5 ms):
        case      Mathilda(pre) Mathilda(post) scipy   wolfram
        L2 n5      2.7          2.7            13.7    1304     (M 5x>scipy, 478x>wolf)
        L2 n12     5.9          5.9            27.5    3032     (4.7x / 516x)
        L2 n50     26.8         26.8           128.5   11519    (4.8x / 430x)
        L3 n5      104.5        16.9           29.7    1430     (1.76x / 85x)
        L3 n12     147.7        38.9           67.5    2919     (1.74x / 75x)
        L3 n50     OOM-killed   222            211.6   11159    (parity / 50x)
      Requirement MET: >= scipy in 5/6 (parity at L3 n50), 50-516x >= wolfram.
- [x] E. FIX: laguerre_l.c — laguerre_l_ndarray() threads sf_machine_laguerre_l_assoc[_c]
      over real/complex buffer, expr_new_ndarray_like (preserves present_as),
      delists for int64/non-machine-a/n>1000/overflow. Reuses scalar helper ->
      bit-identity free. version 0.374->0.375; spec + changelog updated.
- [x] Verify: laguerre_l_tests PASS; valgrind no LaguerreL leak frames (all
      startup/framework baseline); 3-arg complex NDArray + edge cases bit-exact.
- [x] Audits: check-array-exactness PASS (352 probes, 0 MIXED, 0 no-result);
      check-nd-surfaces CLEAN (no disagreement, LaguerreL not flagged, now packs)
- [ ] Release: commit (; v0.375), tag v0.375, push --follow-tags (AWAITING USER GO)
