### Worked examples

```mathematica
In[1]:= Resolve[Exists[x, a x == 1], Reals]  (* an inverse for a exists exactly when a != 0 *)
```

```mathematica
In[1]:= Resolve[ForAll[x, x^2 + a x + 1 > 0], Reals]  (* discriminant band for a positive-definite quadratic *)
```

```mathematica
In[1]:= Resolve[Exists[x, x^4 + p x + q < 0], Reals]  (* parametric elimination of a quartic: q below the minimum value 3 (p/4)^(4/3) *)
```

```mathematica
In[1]:= Resolve[ForAll[x, Exists[y, y > x]]]  (* an alternating prefix, eliminated inner block first *)
```

```mathematica
In[1]:= Resolve[ForAll[x, Abs[x] < d, x^2 < 9], Reals]  (* a restricted quantifier whose Abs is case-split before the CAD *)
```

```mathematica
In[1]:= Resolve[ForAll[eps, eps > 0, Exists[del, del > 0, ForAll[x, 0 < Abs[x - 2] < del, Abs[(3 x - 1) - 5] < eps]]], Reals]  (* the epsilon-delta proof that 3 x - 1 tends to 5 at x = 2 *)
```

```mathematica
In[1]:= Resolve[ForAll[eps, eps > 0, Exists[del, del > 0, ForAll[x, 0 < Abs[x - 2] < del, Abs[(3 x - 1) - 6] < eps]]], Reals]  (* the wrong limit 6 is refuted *)
```

```mathematica
In[1]:= Resolve[ForAll[eps, eps > 0, Exists[nn, nn > 0, ForAll[x, x > nn, Abs[1/x] < eps]]], Reals]  (* a limit at infinity, decidable once 1/x is cleared to polynomial form *)
```

```mathematica
In[1]:= Resolve[Exists[{x, y, z}, x + y + z == 1 && x >= 0 && y >= 0 && z >= 0 && x^2 + y^2 + z^2 > 4], Reals]  (* a three-variable sentence decided False *)
```

```mathematica
In[1]:= Resolve[Exists[x, Sin[x] == a], Reals]  (* transcendental body: no polynomial CAD, so it is left unevaluated *)
```

### Notes

- `Resolve[expr]` and `Resolve[expr, Reals]` eliminate the quantifiers
  `Exists` and `ForAll` from `expr`. A fully quantified sentence gives
  `True` or `False`; otherwise the result is an equivalent quantifier-free
  condition on the free variables. The first argument must be quantified,
  so `Resolve[x > 0]` stays unevaluated. A quantified statement inside
  `Reduce` goes through the same engine.
- **Domain.** `Reals` is the default and the only domain supported. Any
  other domain declines.
- **What it handles.** It handles polynomial bodies, and also bodies with
  `Abs`, `Max`/`Min`, piecewise selectors, square roots, or rational atoms
  with variable denominators, all rewritten to polynomial branches first.
  The restricted form `ForAll[x, cond, body]` can appear at any level of
  an alternating prefix. Together these make ε–δ statements decidable for
  semialgebraic functions, including limits at infinity and infinite
  limits.
- **Soundness rule.** An undecidable sign, a non-`Reals` domain, a
  transcendental body, or a free leading coefficient that mixes a degree
  drop with an irrational boundary all leave the call unevaluated, rather
  than return a formula that might be wrong.
- `Protected`. Symbolic only.
- See also: `Reduce`, `Exists`, `ForAll`, `FindInstance`,
  `CylindricalDecomposition`.
