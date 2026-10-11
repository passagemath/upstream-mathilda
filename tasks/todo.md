# TODO: Implement ChebyshevU (Chebyshev polynomial of the second kind)

Plan: `~/.claude/plans/pasted-content-id-4498-in-src-special-f-reactive-stream.md`
New builtin (no ChebyshevT exists). Template: `src/special_functions/legendre.c`.
Engine: FLINT `fmpz_poly_chebyshev_u` (user decision) + GMP-recurrence `#else`.

## Phase 1 — symbol + source + registration
- [x] `SYM_ChebyshevU` in sym_names.h / sym_names.c (3 sites)
- [x] `src/special_functions/chebyshev_u.{c,h}` — builtin + helpers + `chebyshev_u_init()`
- [x] `chebyshev_u_init()` wired into `src/core.c` (near `legendre_init();`)
- [x] attributes: Listable | NumericFunction | Protected

## Phase 2 — numeric + fast paths
- [x] `sf_machine_chebyshev_u` in sf_machine.{c,h} (forward recurrence, integer n gate)
- [x] `NDK_BIN2(ChebyshevU, ...)` + `REG_B(ChebyshevU, ...)` in ndkernels.c
- [x] MPFR/complex `ncpx` recurrence (integer order) + `cheb_make_result`
- [x] non-integer/complex order: 2F1 template; **verify gate fired** → added the
      ArcCos closed form `Sin[(n+1)ArcCos x]/Sin[ArcCos x]` (matches Mathematica
      branch cut; the algebraic power-form gave the wrong branch at x<-1)
- [x] inexact integer-valued order routes to the real recurrence (clean real)

## Phase 3 — symbolic surface
- [x] special values: U_0=1, U_n(0)=Cos[n Pi/2], U_n(1)=n+1
- [x] D[] rule in deriv.c (enables D + generic Interval threading)

## Phase 4 — docs / tests / release
- [x] docstring in info.c
- [x] `tests/test_chebyshev_u.c` (12 tests) + CMakeLists (source + add_executable + add_test)
- [x] docs/spec/builtins/special-functions.md (## ChebyshevU + intro mention)
- [x] docs/spec/changelog/2026-10-05.md new top entry (historical :20 line left intact)
- [x] version.h bump 0.350 -> 0.351

## Verification
- [x] clean build (gcc-16 -Wall -Wextra, 0 warnings) + `make check-c99`
- [x] all reference values match incl. cancellation ChebyshevU[60,N[2/3,16]]=1.15684 (not 0.)
- [x] fast paths: NDArray eval packs; CompileDiagnostics scalar/array -> Compiled -> True
- [x] `make check-packed-aware` ✓ `check-nd-surfaces` ✓ `check-messages` ✓ `check-c99` ✓
      (`check-compile-coverage` is PRE-EXISTING red on 21 unrelated Image/stats/Fit heads;
       ChebyshevU is NOT in the failure list and passes)
- [x] test suite passes (12/12); valgrind identical to control -> zero new leaks
- [x] graph rebuild; commit `; v0.351` + tag v0.351

## Review

Implemented `ChebyshevU[n, x]` as a new builtin (the request's premise that it
and `ChebyshevT` already existed was wrong — neither did; template was
`legendre.c`). Per user decision, exact integer-order polynomials come from
FLINT (`fmpz_poly_chebyshev_u`) with a GMP-recurrence fallback.

Key design decisions:
- **Three numeric routes**: stable recurrence for integer order + inexact x
  (fixes the degree-60 polynomial-cancellation case the request highlights);
  2F1 expression-template on its convergence disk; ArcCos closed form elsewhere.
  The ArcCos form was chosen over the algebraic power-form *because testing
  showed* the power-form takes the wrong branch on the cut (-∞,-1) — it gave
  -0.383 I where Mathematica (and ArcCos) give +0.597614 I for ChebyshevU[-1.5,-2.4].
- **Fast paths are generic**: one `REG_B` kernel registration delivered NDArray,
  packed, AND Compile[] with no edits to pack.c or compile_*.
- **Quality polish**: negative-order reflection folded into coefficients (prints
  `1 - 4 x^2`, not `-(-1 + 4 x^2)`); inexact integer-valued order routed to the
  real recurrence (clean `35.`, not `35. + 2e-15 I`).

All 10 request reference outputs reproduced exactly. No regressions.

**Flagged for maintainers (out of scope here):** `make check-compile-coverage`
is red on 21 pre-existing heads (Image*, Chop, Clip, Fit, PackedArrayQ, …) that
have a numeric fast path but no Compile lowering and aren't baselined. The gate
is not in CI. Unrelated to this change.
