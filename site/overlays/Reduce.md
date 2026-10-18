---
status: Partial
references:
  - "Collins & Hong, \"Partial Cylindrical Algebraic Decomposition for Quantifier Elimination\", J. Symbolic Computation 12 (1991)."
  - "McCallum, \"An Improved Projection Operation for Cylindrical Algebraic Decomposition\", in Quantifier Elimination and Cylindrical Algebraic Decomposition (1998)."
  - "Basu, Pollack & Roy, \"Algorithms in Real Algebraic Geometry\" (2nd ed., 2006), Ch. 5 & 11 (sign determination, CAD)."
  - "Schrijver, \"Theory of Linear and Integer Programming\" (1986), §12.2 (Fourier-Motzkin elimination)."
---
### Worked examples

```mathematica
In[1]:= Reduce[a x == b, x]  (* the degenerate branch a == 0 is kept, unlike Solve *)
```

```mathematica
In[1]:= Reduce[x^2 < 1]  (* variables inferred; the ordering inequality selects the Reals *)
```

```mathematica
In[1]:= Reduce[1/x < 1, x, Reals]  (* the pole at 0 is an excluded breakpoint of the sign diagram *)
```

```mathematica
In[1]:= Reduce[Sin[x] == 1/2 && 0 < x < 2 Pi, x]  (* periodic families filtered against a bounding region *)
```

```mathematica
In[1]:= Reduce[x^2 + y^2 < 1 && y > x^2, {x, y}, Reals]  (* two-variable CAD with algebraic breakpoints *)
```

```mathematica
In[1]:= Reduce[x y z > 0, {x, y, z}, Reals]  (* CAD in three variables, split by the sign of each coordinate *)
```

```mathematica
In[1]:= Reduce[x + y + z == 1 && x > 0 && y > 0 && z > 0 && x y + x z + y z < 9 x y z]  (* exact sum-of-squares certificate proves the region empty *)
```

```mathematica
In[1]:= Reduce[Exists[y, x^2 + y^2 < 1], {x}, Reals]  (* an embedded quantifier is eliminated *)
```

```mathematica
In[1]:= Reduce[2 x + 3 y == 1, {x, y}, Integers]  (* linear Diophantine family with an integer parameter *)
```

```mathematica
In[1]:= Reduce[x^2 + y^2 == 1, {x, y}]  (* honest decline: a positive-dimensional system over the complexes stays unevaluated *)
```

### Notes

- **Reduce and Solve answer different questions.** `Solve[eqns, vars]`
  returns generic solutions as replacement rules and drops the special
  cases. `Reduce[expr, vars]` returns a quantifier-free `And`/`Or` formula
  covering the whole solution set, every parametric and boundary branch
  included. Use `Reduce` when you need the exact conditions, especially
  with symbolic parameters or with inequalities.
- **Domain and variables.** The default domain is `Complexes`. A statement
  with an ordering inequality, or with a real-only function such as `Abs`,
  `Floor` or `Max`, defaults to `Reals`. `Reduce[expr]` and
  `Reduce[expr, dom]` infer the variables, taking every free symbol except
  named constants, in alphabetical order. A third argument (`Complexes`,
  `Reals`, `Integers` or `Rationals`) chooses the domain explicitly.
- **Engines over the reals.** One variable goes to an exact
  real-algebraic sign diagram, which handles poles, radicals, `Abs`,
  `Floor`, piecewise heads and bounded trig regions. Several variables go
  to Fourier–Motzkin for linear systems. Higher-degree regions go to an
  exact SOS/Positivstellensatz emptiness certificate, then an n-variable
  McCallum CAD, then a zero-dimensional solver. Rational atoms in several
  variables are cleared to polynomial form first.
- **Options.** `Cubics`, `Quartics`, `Modulus -> p`,
  `GeneratedParameters`, `Backsubstitution`, `WorkingPrecision` and
  `Method` are accepted. `TimeConstrained[Reduce[...], t]` interrupts the
  CAD cleanly, without leaking.
- **Soundness over completeness.** The engine is designed to return either
  the exact solution set or nothing. An undecidable sign, an unsupported
  construct, or an engine not yet wired (positive-dimensional nonlinear
  systems over the complexes) leaves the call unevaluated.
- `Protected`. Symbolic, with no NDArray or `Compile[]` path.
- See also: `Solve`, `Resolve`, `FindInstance`, `CylindricalDecomposition`,
  `Refine`.
