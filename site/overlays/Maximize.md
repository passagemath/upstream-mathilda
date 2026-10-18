### Worked examples

```mathematica
In[1]:= Maximize[-x^4 + 2 x^2 + 1, x]  (* a double-hump quartic: the two global maxima at x = -1 and x = 1 tie, and one is reported *)
```

```mathematica
In[1]:= Maximize[{x y, x + y == 10}, {x, y}]  (* an equality constraint: largest product for a fixed sum *)
```

```mathematica
In[1]:= Maximize[{x + 2 y, x^2 + y^2 <= 5}, {x, y}]  (* a linear objective on a disk: the maximum is on the boundary in the gradient direction *)
```

```mathematica
In[1]:= Maximize[{3 x + 2 y, x + y <= 4 && x + 3 y <= 6 && x >= 0 && y >= 0}, {x, y}]  (* a small linear program solved exactly: the optimum is at the vertex (4, 0) *)
```

```mathematica
In[1]:= Maximize[{x^2 y, x^2 + y^2 == 3}, {x, y}]  (* a Lagrange-multiplier problem with an irrational optimal point *)
```

```mathematica
In[1]:= Maximize[{Sqrt[x] + Sqrt[1 - x], 0 <= x <= 1}, x]  (* a radical objective, made polynomial with auxiliary variables *)
```

```mathematica
In[1]:= Maximize[{x y, x^2 + y^2 <= 8 && x >= 0}, {x, y}, Integers]  (* the Integers domain: the region is bounded by its continuous relaxation and then enumerated exactly *)
```

```mathematica
In[1]:= Maximize[x^3 - x, x]  (* unbounded above: returns Infinity and emits Maximize::natt *)
```

```mathematica
In[1]:= Maximize[{x, x^2 + 1 <= 0}, x]  (* an empty region: -Infinity and Maximize::infeas *)
```

```mathematica
In[1]:= Maximize[Sin[x], x]  (* transcendental objectives are not handled yet, so the call stays unevaluated *)
```

### Notes

- `Maximize` is `Minimize` applied to `-f` with the value negated. It has the same argument forms (`f`, `{f, cons}`, variable list, optional domain `Reals`/`Integers`) and returns `{f_max, {x -> x_max, ...}}` with exact values.
- The sentinels are swapped relative to Minimize: unbounded gives `{Infinity, {x -> Indeterminate, ...}}` with `Maximize::natt`; infeasible gives `{-Infinity, ...}` with `Maximize::infeas`. Both go through the message funnel, so `Quiet` and `Check` work.
- Inexact input (`Maximize[1.5 x - x^2, x]`) is passed to `NMaximize`. Compared with `NMaximize`, Maximize proves global optimality exactly (or declines), while `NMaximize` always returns a machine-precision point from a stochastic search.
- `TimeConstraint -> t` (CPU seconds, default 30) bounds the whole call. Attributes: `{Protected}`, not `HoldAll`.
- See also: `Minimize`, `NMaximize`, `FindMaximum`, `Reduce`.
