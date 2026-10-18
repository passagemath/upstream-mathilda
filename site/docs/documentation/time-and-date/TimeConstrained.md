# TimeConstrained

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`TimeConstrained[expr, t]`**

evaluates expr, stopping after t seconds.

**`TimeConstrained[expr, t, failexpr]`**

returns failexpr if the time constraint is not met.

**`TimeConstrained[expr, Infinity]`**

imposes no time constraint.

<details>
<summary>Notes</summary>

TimeConstrained generates an interrupt to abort the evaluation of expr if the evaluation is not completed within the specified time. TimeConstrained evaluates failexpr only if the evaluation is aborted. TimeConstrained returns $Aborted if the evaluation is aborted and no failexpr is specified. TimeConstrained may give different results on different occasions within a single session, for example as a result of different conditions of internal system caches. TimeConstrained takes account only of CPU time spent inside the main Mathilda kernel process; it does not include additional threads or processes.

</details>

## Examples (14)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= TimeConstrained[2 + 2, 1, Print["never printed"]; -1]
Out[1]= 4

In[2]:= TimeConstrained[While[True], 0.2, "timed out"]
Out[2]= "timed out"

In[3]:= TimeConstrained[Print["not run"], 0]
Out[3]= $Aborted

In[4]:= TimeConstrained[{TimeConstrained[While[True], 0.1, "inner"], 7}, 5]
Out[4]= {"inner", 7}
```

### Applications (10)

A body that finishes inside its budget returns its value unchanged

```mathematica
In[5]:= TimeConstrained[Integrate[x^2 Exp[x], x], 10]
Out[5]= E^x (2 - 2 x + x^2)
```

An over-budget body is aborted, even inside an auto-compiled loop

```mathematica
In[6]:= TimeConstrained[Do[Sin[1.0], {10^9}], 0.5]
Out[6]= $Aborted
```

The three-argument form returns failexpr on timeout

```mathematica
In[7]:= TimeConstrained[While[True], 0.2, "timed out"]
Out[7]= "timed out"
```

Failexpr is held and evaluated only if the body is aborted

```mathematica
In[8]:= TimeConstrained[2 + 2, 1, Print["never printed"]; -1]
Out[8]= 4
```

Only the budget is evaluated up front, so it can be any expression

```mathematica
In[9]:= t = 0.1; TimeConstrained[While[True], 2 t, {"budget", 2 t}]
Out[9]= {"budget", 0.2}
```

An Infinity budget imposes no limit

```mathematica
In[10]:= TimeConstrained[FactorInteger[2^64 + 1], Infinity]
Out[10]= {{274177, 1}, {67280421310721, 1}}
```

A non-positive budget aborts at once, without evaluating the body

```mathematica
In[11]:= TimeConstrained[Print["not run"], 0]
Out[11]= $Aborted
```

Nesting clamps: the 30 s inner budget is cut to the outer 0.3 s

```mathematica
In[12]:= {t1, r} = AbsoluteTiming[TimeConstrained[TimeConstrained[While[True], 30, "inner"], 0.3, "outer"]]; {r, t1 < 1}
Out[12]= {"inner", True}
```

An inner timeout is local; the outer computation carries on

```mathematica
In[13]:= TimeConstrained[{TimeConstrained[While[True], 0.1, "inner"], 7}, 5]
Out[13]= {"inner", 7}
```

A budget that is not a number is left unevaluated

```mathematica
In[14]:= TimeConstrained[x^2, a]
Out[14]= TimeConstrained[x^2, a]
```

## Implementation notes

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

- `HoldAll`, `Protected`.
- Returns the value of `expr` if it completes within `t` seconds of CPU time.
- Returns `$Aborted` if the budget is exhausted and no `failexpr` is provided.
- Returns the (then-evaluated) value of `failexpr` if the budget is exhausted and the three-argument form is used. `failexpr` is **not** evaluated when the body completes in time.
- `TimeConstrained[expr, Infinity]` imposes no time constraint.
- Zero, negative, or non-numeric (NaN) time budgets abort immediately without evaluating `expr`; the abort path still produces `$Aborted` (or `failexpr` for the 3-argument form).
- The time budget is measured with `ITIMER_PROF` and counts only CPU time (user + kernel) consumed by the Mathilda kernel process. Wall-clock time spent in I/O, sleep, or other processes is not charged.
- A cooperative wall-clock deadline (`CLOCK_MONOTONIC` + budget) is also armed and checked once per evaluator rewrite step, as a portability backstop on hosts where `ITIMER_PROF` is unreliable (notably WSL 1, whose syscall-translation layer under-counts CPU time and delivers `SIGPROF` many seconds late). On real Linux / macOS the `SIGPROF` normally fires first and the cooperative check is a cheap no-op; on broken hosts the cooperative check enforces the deadline between rewrite steps. The only case that escapes both layers is a single long-running C builtin (e.g. `FactorInteger` on a huge composite) on a broken host -- it must wait for the late `SIGPROF`.
- May give different results on different occasions within a single session, for example as a result of different conditions of internal system caches.
- Nested `TimeConstrained` calls compose: each call saves and restores the previous `SIGPROF` handler, `ITIMER_PROF` state, and the cooperative-deadline state, so an inner abort does not disturb an outer time budget.
- A timeout also restores every C nesting counter a subsystem has registered with `tc_register_unwind_counter` (v0.451) — the unwind skips their decrements, and a counter left raised changes behaviour for the rest of the session (one aborted nested `DSolve` used to leave `DSolve`'s recursion depth at 2, disabling its top-level-only methods and its failure memo).
- **Nesting CLAMPS, it does not refund** (v0.253). An inner budget can only ever *shorten* the enclosing one: the timer actually armed is `min(inner request, outer remaining)`, and on exit the enclosing budget is reinstalled **charged** for the CPU time the inner call consumed. So `TimeConstrained[expr, t]` is an upper bound on the whole call, however many nested budgets `expr` enters. A shorter inner budget still wins, and a computation that fits inside its budget is unaffected. Before this, both layers lifted the caller's deadline merely by entering an inner scope — the `ITIMER_PROF` layer reinstalled the outer timer at the value it held when the inner call *started*, and the cooperative layer set an inner absolute deadline that could be later than the outer's — so `TimeConstrained[TimeConstrained[<long loop>, 30], 3]` ran the loop to completion in 24.5 s. Since `Integrate`, `Simplify`, `NIntegrate` and `DSolve` all bound their own sub-steps this way, no budget was an upper bound in practice.
- The abort is implemented by `siglongjmp`-ing out of the in-flight evaluator. Expression nodes allocated by the aborted computation are not reclaimed; this is the documented behaviour.

**Attributes:** `HoldAll`, `Protected`.

## References

**See also:** [HoldAll](../../expression-information/HoldAll/), [FactorInteger](../../number-theory/FactorInteger/), [DSolve](../../calculus/DSolve/), [Integrate](../../calculus/Integrate/), [Simplify](../../simplification/Simplify/), [NIntegrate](../../numerical-calculus/NIntegrate/)

- Source: [`src/core.c`](https://github.com/stblake/mathilda/blob/main/src/core.c)
- Specification: [`docs/spec/builtins/time-and-date.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/time-and-date.md)
- Tests: [`tests/test_dsolve.c`](https://github.com/stblake/mathilda/blob/main/tests/test_dsolve.c)
- Tests: [`tests/test_dsolve_m34_stress.c`](https://github.com/stblake/mathilda/blob/main/tests/test_dsolve_m34_stress.c)
- Tests: [`tests/test_dsolve_m62_stress.c`](https://github.com/stblake/mathilda/blob/main/tests/test_dsolve_m62_stress.c)
- Tests: [`tests/test_dsolve_m63_stress.c`](https://github.com/stblake/mathilda/blob/main/tests/test_dsolve_m63_stress.c)

## Notes & additional examples

### Notes

- `TimeConstrained[expr, t]` returns the value of `expr` if it finishes within `t` seconds, and `$Aborted` otherwise. `TimeConstrained[expr, t, failexpr]` returns `failexpr` instead, evaluated after the timeout with a fresh budget.
- Attributes `HoldAll`, `Protected`. Only `t` is evaluated up front. It may be an Integer, a Real, a big integer, a Rational or `Infinity`. Zero, negative and NaN budgets abort immediately. Budgets are capped near 2147483 s.
- On POSIX the budget is CPU time, measured with `ITIMER_PROF`/`SIGPROF`. This interrupts even a long-running C builtin or compiled loop. A cooperative wall-clock deadline is also polled at every evaluator step and by heavy kernels such as CAD, as a backstop on hosts with unreliable profiling timers (WSL 1). Native Windows has only the cooperative layer, so a single long C builtin there runs to completion before the deadline is noticed.
- Nesting clamps and never refunds (v0.253). An inner budget is `min(inner, outer remaining)`, and its CPU time is charged to the outer budget, so `TimeConstrained[expr, t]` bounds the whole call however many budgets `Integrate`, `Simplify`, `DSolve` or `NIntegrate` open inside it. When the clamped inner budget is what fires, the inner `failexpr` is the value you see.
- After a timeout the kernel restores the evaluator recursion depth, `Block` frames, message-suppression depth and every C counter registered with `tc_register_unwind_counter` (v0.451). An aborted call therefore cannot disable later behaviour in the session.
- The abort unwinds with `siglongjmp`, so expression nodes built by the abandoned computation are leaked by design. Results may differ between runs because the cut-off point depends on load and on cache state.
- See also: `AbsoluteTiming`, `Timing`, `Check`, `Quiet`.
