### Worked examples

```mathematica
In[1]:= TimeConstrained[Integrate[x^2 Exp[x], x], 10]  (* a body that finishes inside its budget returns its value unchanged *)
```

```mathematica
In[1]:= TimeConstrained[Do[Sin[1.0], {10^9}], 0.5]  (* an over-budget body is aborted, even inside an auto-compiled loop *)
```

```mathematica
In[1]:= TimeConstrained[While[True], 0.2, "timed out"]  (* the three-argument form returns failexpr on timeout *)
```

```mathematica
In[1]:= TimeConstrained[2 + 2, 1, Print["never printed"]; -1]  (* failexpr is held and evaluated only if the body is aborted *)
```

```mathematica
In[1]:= t = 0.1; TimeConstrained[While[True], 2 t, {"budget", 2 t}]  (* only the budget is evaluated up front, so it can be any expression *)
```

```mathematica
In[1]:= TimeConstrained[FactorInteger[2^64 + 1], Infinity]  (* an Infinity budget imposes no limit *)
```

```mathematica
In[1]:= TimeConstrained[Print["not run"], 0]  (* a non-positive budget aborts at once, without evaluating the body *)
```

```mathematica
In[1]:= {t1, r} = AbsoluteTiming[TimeConstrained[TimeConstrained[While[True], 30, "inner"], 0.3, "outer"]]; {r, t1 < 1}  (* nesting clamps: the 30 s inner budget is cut to the outer 0.3 s *)
```

```mathematica
In[1]:= TimeConstrained[{TimeConstrained[While[True], 0.1, "inner"], 7}, 5]  (* an inner timeout is local; the outer computation carries on *)
```

```mathematica
In[1]:= TimeConstrained[x^2, a]  (* a budget that is not a number is left unevaluated *)
```

### Notes

- `TimeConstrained[expr, t]` returns the value of `expr` if it finishes within `t` seconds, and `$Aborted` otherwise. `TimeConstrained[expr, t, failexpr]` returns `failexpr` instead, evaluated after the timeout with a fresh budget.
- Attributes `HoldAll`, `Protected`. Only `t` is evaluated up front. It may be an Integer, a Real, a big integer, a Rational or `Infinity`. Zero, negative and NaN budgets abort immediately. Budgets are capped near 2147483 s.
- On POSIX the budget is CPU time, measured with `ITIMER_PROF`/`SIGPROF`. This interrupts even a long-running C builtin or compiled loop. A cooperative wall-clock deadline is also polled at every evaluator step and by heavy kernels such as CAD, as a backstop on hosts with unreliable profiling timers (WSL 1). Native Windows has only the cooperative layer, so a single long C builtin there runs to completion before the deadline is noticed.
- Nesting clamps and never refunds (v0.253). An inner budget is `min(inner, outer remaining)`, and its CPU time is charged to the outer budget, so `TimeConstrained[expr, t]` bounds the whole call however many budgets `Integrate`, `Simplify`, `DSolve` or `NIntegrate` open inside it. When the clamped inner budget is what fires, the inner `failexpr` is the value you see.
- After a timeout the kernel restores the evaluator recursion depth, `Block` frames, message-suppression depth and every C counter registered with `tc_register_unwind_counter` (v0.451). An aborted call therefore cannot disable later behaviour in the session.
- The abort unwinds with `siglongjmp`, so expression nodes built by the abandoned computation are leaked by design. Results may differ between runs because the cut-off point depends on load and on cache state.
- See also: `AbsoluteTiming`, `Timing`, `Check`, `Quiet`.
