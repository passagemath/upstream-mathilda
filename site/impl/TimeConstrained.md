---
source: src/core.c
---
**Algorithm.** `builtin_time_constrained` implements `TimeConstrained[expr, t]`
and `TimeConstrained[expr, t, failexpr]` under `HoldAll`. Only `t` is evaluated
up front. `tc_parse_time` accepts Integer, Real, BigInt, `Rational`, and
`Infinity`/`DirectedInfinity[1]`; anything else leaves the call unevaluated.
`Infinity` evaluates the body unconstrained. A non-positive or NaN budget aborts
immediately, returning `failexpr` or `$Aborted`, and a runaway budget is clamped
to 2147483 s.

For a finite budget two layers are armed:

1. **Async layer (POSIX only).** A `SIGPROF` handler (`tc_sigprof_handler`) and
   an `ITIMER_PROF` timer. The timer is read with `getitimer` first, and the
   armed value is `min(t, outer remaining)`, so nesting **clamps** rather than
   refunds.
2. **Cooperative wall-clock layer.** A deadline `tc_deadline_sec`, computed from
   `mth_mono_seconds` and min'd with any enclosing deadline. The evaluator checks
   it each rewrite step through `tc_check_deadline`, and heavy C kernels such as
   CAD poll the non-jumping `tc_deadline_passed` and decline cleanly. This layer
   is the backstop on WSL 1 and the only layer on native Windows.

The body runs in `tc_run_guarded`, which holds the `sigsetjmp` landing pad
(`TC_SETJMP`, plain `setjmp` on MinGW). The jump is not always immediate. While
a GMP or FLINT allocation is in flight (`tc_gmp_alloc_busy`), or inside a heavy
method's `tc_async_defer_push` region, the handler only flags the timeout. The
jump is then taken at the next safe point (`tc_alloc_safepoint`, or the
per-step check), and only on the main thread, never from a worker. On unwind,
`tc_run_guarded` resets the Orderless sort memo, drains `Block` frames, and
restores the async-defer count, the `Quiet` and `Solve::ifun` suppression depths,
and every counter registered with `tc_register_unwind_counter`. The caller
restores the recursion depth, disarms its timer, and reinstalls the outer timer
**charged** for the CPU time consumed (at least 1 µs, so an exhausted outer
budget aborts at its next step). It also restores the previous handler, the
`jmp_buf` and the deadline state. Finally it returns `failexpr`, evaluated
afterward with a fresh budget, or `$Aborted`.

**Data structures.** The state is a process-global `tc_jmp_env` plus
`sig_atomic_t` flags (`tc_timed_out`, `tc_gmp_alloc_busy`, `tc_async_deferred`)
and the cooperative deadline pair. Each call saves these on the C stack, so
nested calls compose. The unwind-counter registry is a fixed array of at most 32
`int*`. GMP and FLINT memory functions are replaced at startup by
`tc_install_alloc_guard` so that allocator calls are bracketed.

**Complexity / limits.** The overhead is a few system calls per call plus one
monotonic-clock read per rewrite step while a budget is active. On POSIX the
budget measures CPU time (`ITIMER_PROF`), not wall-clock time. On Windows it is
cooperative only: a single long C builtin that never returns to the evaluator or
polls the deadline runs to completion. On the timeout path the abandoned
partial tree is leaked, an inherent consequence of `siglongjmp`. Bounded
allocations made during a deferred jump leak the same way.
