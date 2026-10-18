### Worked examples

```mathematica
In[1]:= DSolve[{y'[x] == x + y[x], y[0] == 1}, y, x]  (* an initial-value problem: C[1] is fitted and the answer is a pure Function *)
```

```mathematica
In[1]:= DSolve[y''[x] + 4 y[x] == Log[x], y[x], x]  (* variation of parameters; the particular needs SinIntegral and CosIntegral *)
```

```mathematica
In[1]:= DSolve[y[x] == x y'[x] + y'[x]^2, y[x], x, IncludeSingularSolutions -> True]  (* Clairaut: the line family plus the singular envelope -x^2/4 *)
```

```mathematica
In[1]:= DSolve[y'[x] == y[x]^Sqrt[2] + y[x], y[x], x]  (* a Bernoulli equation with an irrational exponent, answered as an implicit first integral *)
```

```mathematica
In[1]:= DSolve`InverseOrientation[y'[x] == y[x]/(x + y[x]^3), y[x], x]  (* linear in x(y) once the roles of x and y are swapped *)
```

```mathematica
In[1]:= DSolve[y''[x] - (x^2 + 3) y[x] == 0, y[x], x]  (* Kovacic: a Liouvillian solution x E^(x^2/2) and its partner *)
```

```mathematica
In[1]:= DSolve[x^2 y''[x] + (-x^2/4 + 2 x - 3/16) y[x] == 0, y[x], x]  (* Whittaker's equation, emitted in its confluent 1F1 form *)
```

```mathematica
In[1]:= DSolve[{y'[t] + y[t] == DiracDelta[t - 1], y[0] == 0}, y[t], t]  (* impulse forcing: the Green's function of the operator *)
```

```mathematica
In[1]:= DSolve[{x'[t] == y[t], y'[t] == -x[t], x[0] == 1, y[0] == 0}, {x[t], y[t]}, t]  (* a first-order linear system with initial data *)
```

```mathematica
In[1]:= DSolve[y''[x] + Sin[x] y[x] == 0, y[x], x]  (* no closed form: left unevaluated rather than truncated *)
In[2]:= DSolve[y''[x] + Sin[x] y[x] == 0, y[x], x, Method -> "PowerSeries"]  (* a series answer only on request *)
```

### Notes

- `DSolve[eqn, y[x], x]` returns rules for `y[x]`; `DSolve[eqn, y, x]` returns `y -> Function[...]`, which also replaces derivatives, so `eqn /. sol` checks the answer. Constants are `C[1], C[2], ...`; rename them with `GeneratedParameters`.
- Options: `GeneratedParameters -> C`, `Assumptions -> True`, `Method -> Automatic`, `IncludeSingularSolutions -> False`. Attributes: `Protected` only (not `HoldAll`).
- **An automatic answer is always a solution.** Every branch is checked by back-substitution. A truncated series, or a relation that still contains `y'` (a first integral), is never returned from `Method -> Automatic`. Those come only from a pinned method (`"PowerSeries"`/`"FrobeniusSeries"`, `"FirstOrderPowerSeries"`, `"ReducibleFirstIntegral"`). An equation the cascade cannot close stays unevaluated, as in Mathematica.
- Answers may be explicit rules, implicit equations `F(x, y[x]) == C[1]` (separable, exact, Bernoulli, Abel, Lie symmetry), or parametric (Lagrange, SolvableForY). A non-elementary quadrature can stay as `Inactive[Integrate]`.
- `Method -> "<name>"` pins one method with no fallback, e.g. `"Separable"`, `"Bernoulli"`, `"Kovacic"`, `"LieSymmetry"`, `"EulerCauchy"`, `"SpecialFunctionForm"`. The automatic-only late passes, such as `DSolve`InverseOrientation`, `DSolve`ExactMuU`, `DSolve`BernoulliExpanded` and `DSolve`HomogeneousG`, can also be called directly as qualified heads.
- Systems (`{y1, y2, ...}`) and PDEs (`u, {x, t}`) go through their own dispatcher branches: decoupling, triangular, constant- and variable-coefficient linear, 2-D autonomous; PDEs use characteristics, Charpit, d'Alembert and heat-kernel IVPs. `DSolve`EigenvalueProblem` handles Sturm–Liouville problems.
- Known gap: for Laguerre-type `x y'' + (1 - x) y' + n y == 0` with integer `n`, the `HypergeometricU[-n, 1, x]` partner is a polynomial multiple of the first solution, so the returned basis is not independent.
- See also: `NDSolve`, `Integrate`, `Series`, `PDEClassify`.
