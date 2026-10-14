# TODO: Windows (MSYS2/MinGW UCRT64) clean-compile

Context: passagemath PR #2966 — Mathilda v0.387 failed to build on native MSYS2/MinGW
because of unguarded POSIX time functions (regression). Make it compile+link cleanly on
UCRT64. Full plan: ~/.claude/plans/let-s-review-the-comment-partitioned-pixel.md

## Blockers (always-compiled)
- [x] 1. New `src/mtime.{c,h}` timing shim (mono / wall / sleep), POSIX + `_WIN32` impls
- [x] 1b. Replace 5 copy-pasted clock helpers with the shim (datetime, simp_util, integrate x2, core)
- [x] 2. `core.c` TimeConstrained: portable `TC_SETJMP`/`TC_LONGJMP` + `#ifndef _WIN32` on SIGPROF/itimer; cooperative deadline as `double` via shim
- [x] 3. `meminfo.c` getrusage → guard + Windows GetProcessMemoryInfo branch
- [x] 4. `files.c:86` lstat → `#define lstat stat` on `_WIN32`
- [x] 5. `ffi/mathilda_ffi.c:62` setenv → `_putenv` on `_WIN32`

## Build / CI / docs
- [x] makefile: `-lpsapi` on Windows LDFLAGS; MPFR-must-be-present note; `TARGET=Mathilda.exe` on Windows
- [x] `.github/workflows/build.yml`: MSYS2/UCRT64 build+smoke job (authoritative Windows gate)
- [x] `tools/check_c99_portability.py`: MinGW-missing symbol set, `_WIN32`-region-aware (assert-empty) + SPEC.md §10 note
- [x] `docs/building.md`: mpfr in pacman list; cooperative-only TimeConstrained + degradations note
- [x] `tests/CMakeLists.txt`: add `../src/mtime.c` to mathilda_common (explicit source list)

## Housekeeping
- [x] `src/version.h` → v0.406; changelog note under week-of-2026-10-12 (ISO Monday)
- [x] macOS sanity: clean `make`; `make check-c99` (incl. new MinGW scan) passes
- [x] behaviour verified unchanged: TimeConstrained (abort/failexpr/nested/Infinity), AbsoluteTiming, Timing, MemoryInUse, MaxMemoryUsed, FileExistsQ, Pause, SessionTime, DateList
- [x] unit suites: `time_constrained_tests` (incl. cooperative-backstop-without-sigprof) + `datetime_tests` pass

## Review

**Outcome.** All Tier-1 compile/link blockers for native Windows (MSYS2/MinGW
UCRT64) closed, plus the regression gate. macOS build is clean and behaviourally
identical; the MinGW surface is now either `_WIN32`-guarded or routed through the
new `src/mtime.c` seam.

**What shipped (v0.406):** `src/mtime.{c,h}` shim; `core.c` TimeConstrained split
into the cooperative deadline (portable, kept on Windows) and the `#ifndef _WIN32`
SIGPROF/itimer async layer, with a `TC_SETJMP`/`TC_LONGJMP` macro; `meminfo.c`
GetProcessMemoryInfo; `files.c` lstat alias; `ffi` `_putenv`; makefile `-lpsapi`
+ `Mathilda.exe` target; MinGW scan in `check_c99_portability.py`; MSYS2/UCRT64
CI job; docs (`building.md`, `SPEC.md §10`) + changelog.

**Verification gaps / notes.**
- The real UCRT64 compile is only provable in the new CI job (we build on macOS);
  that job is the authoritative gate.
- `TimeConstrained` on Windows is cooperative-only (documented degradation; cannot
  preempt a non-cooperative C kernel) — a deliberate scope decision, not a bug.
- Not done (out of scope, per approved plan): porting the C unit-test POSIX
  timeout harness (`alarm`/`fork`/`usleep`) so `ctest` runs on Windows.
- Did not rebuild the code-review graph (verification was via compile + targeted
  tests + the portability gate; a structural diff adds nothing here).
