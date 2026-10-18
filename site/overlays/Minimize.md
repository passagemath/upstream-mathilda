### Worked examples

```mathematica
In[1]:= Minimize[x^4 - 3 x^2 + x, x]  (* a quartic with two wells: the global minimiser is the least real root of the cubic f' = 0, kept exact as a Root object *)
```

```mathematica
In[1]:= N[Minimize[x^4 - 3 x^2 + x, x]]  (* numerically -3.51391 at x = -1.30084, the deeper well; FindMinimum started at x = 2 stops in the shallower well at -1.07023 *)
```

```mathematica
In[1]:= Maximize[{x y z, x + y + z == 3 && x >= 0 && y >= 0 && z >= 0}, {x, y, z}]  (* mixed equality and inequality constraints: AM-GM on the simplex *)
```

```mathematica
In[1]:= Simplify[First[Minimize[{x^2 + y^2, (x - 2)^2 + (y - 3)^2 == 1}, {x, y}]]]  (* squared distance from the origin to a circle, (Sqrt[13] - 1)^2; the raw value is an unsimplified radical sum *)
```

```mathematica
In[1]:= Minimize[{(x^2 + 1)/(x + 2), x >= 0}, x]  (* rational objective: the denominator is proved positive on the region, giving 2 Sqrt[5] - 4 at x = Sqrt[5] - 2 *)
```

```mathematica
In[1]:= Minimize[Sum[Abs[x - i^2], {i, 1, 10}], x]  (* piecewise-linear objective minimised exactly over its sign cells; every x in [25, 36] is optimal *)
```

```mathematica
In[1]:= Minimize[{x + y, x^2 + y^2 == 25}, {x, y}, Integers]  (* the Integers domain: lattice points on a circle, found with Solve over the integers *)
```

```mathematica
In[1]:= Minimize[{x + y, x <= y^2}, {x, y}]  (* unbounded below: proved with quantifier elimination, returns -Infinity and emits Minimize::natt *)
```

```mathematica
In[1]:= Minimize[{x + y, x^2 < -1}, {x, y}]  (* an empty region: Infinity with Indeterminate rules and Minimize::infeas *)
```

```mathematica
In[1]:= Minimize[Exp[x] - x, x]  (* transcendental objectives are not handled yet, so the call stays unevaluated *)
```

### Notes

- `Minimize[f, x]`, `Minimize[f, {x, y, ...}]`, `Minimize[{f, cons}, vars]` and `Minimize[..., dom]` return `{f_opt, {x -> x_opt, ...}}` with **exact** values (Integer, Rational, radical, or `Root`). Each result is checked with exact oracles: `Solve` over the Reals, exact algebraic sign comparison, and `Reduce`/`Resolve` certificates. If any oracle cannot decide, the call is left unevaluated. Minimize does not guess.
- The minimum value is the objective evaluated at the minimiser. It can come back unreduced (for example a polynomial in a `Root`, or a sum of squared radicals), so wrap it in `Simplify` or `RootReduce` for a canonical form.
- **Inexact input** (`Minimize[2.5 x^2 - 3 x, x]`) is passed straight to `NMinimize`. `NMinimize` is a stochastic global search over machine numbers that always returns *some* point. Minimize either proves its answer or declines. For example, `NMinimize[{x + y, x <= y^2}, {x, y}]` returns a large finite negative number, while Minimize proves the problem is unbounded below (`-Infinity`).
- Special returns: unbounded gives `{-Infinity, {x -> Indeterminate, ...}}` with `Minimize::natt`; infeasible gives `{Infinity, {x -> Indeterminate, ...}}` with `Minimize::infeas`. Both messages go through the message funnel, so `Quiet` and `Check` see them. An infimum that is approached but not attained (for example `Minimize[{x, x > 0}, x]`) currently declines.
- Domains: `Reals` (the default) and `Integers`. The Integers domain can also be written as `Element[{vars}, Integers]` among the constraints, but it must cover every variable; mixed integer/continuous problems decline. `Complexes` with exact input stays unevaluated.
- Options: `TimeConstraint -> t` (seconds of CPU time, default 30; `Infinity` for no limit) bounds all internal `Solve`/`Reduce` probes together.
- Attributes: `{Protected}`. Minimize is not `HoldAll`, as in Mathematica, so a variable that has a value is substituted before Minimize sees it.
- Scaling: the CAD certificates are doubly exponential in the number of variables, and at most six inequality constraints are allowed. Problems in a few variables with moderate degree are the practical range.
- See also: `Maximize`, `NMinimize`, `FindMinimum`, `Reduce`, `Solve`, `FindInstance`.
