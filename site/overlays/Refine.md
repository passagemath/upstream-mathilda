### Worked examples

```mathematica
In[1]:= Refine[Sqrt[x^2], x > 0]  (* the square root of a square collapses once the sign is known *)
```

```mathematica
In[1]:= Refine[Sqrt[1 + b^2/a^2], a > 0]  (* radical reduction over a quotient: the denominator is exposed and |a| pulled out *)
```

```mathematica
In[1]:= Refine[Sqrt[(a + b)^2 (a - b)^2], a > b > 0]  (* each squared factor resolved by a CAD sign proof *)
```

```mathematica
In[1]:= Refine[Sign[x^3 - x], 0 < x < 1]  (* sign of a compound polynomial on an interval *)
```

```mathematica
In[1]:= Refine[{Floor[n], Mod[2 n, 2], (-1)^(2 n), Sin[n Pi]}, Element[n, Integers]]  (* integer-assumption structural rules shared with Simplify *)
```

```mathematica
In[1]:= Refine[Arg[x], x < 0]  (* the argument of a negative real *)
```

```mathematica
In[1]:= Refine[x^2 + y^2 >= 2 x y, Element[{x, y}, Reals]]  (* a predicate decided by Reduce/CAD entailment *)
```

```mathematica
In[1]:= Refine[x == y, x > y]  (* an equation refuted by the assumption *)
```

```mathematica
In[1]:= Assuming[a > 0 && b > 0, Refine[Sqrt[a^2 b^2]]]  (* with no positional assumption Refine reads $Assumptions *)
```

```mathematica
In[1]:= Refine[x^2 > 1, x > 0]  (* undecided: true for x > 1 only, so the predicate is returned unchanged *)
```

### Notes

- `Refine[expr, assum]` gives the form `expr` takes when its symbols are
  given values that satisfy `assum`. It shares the assumption engine of
  `Simplify` (`AssumeCtx`, `apply_assumption_rules`,
  `assume_structural_rewrite`), so every rewrite `Refine` makes is one
  `Simplify` can also make. `Refine` is the rewrite-and-decide pass without
  `Simplify`'s complexity search: it turns `Sqrt[x^2]` into `x` under
  `x > 0`, but it does not factor or expand. Conversely, a correct rewrite
  that `Simplify`'s scorer would reject (`ArcSin[Sin[x]] -> Pi - x` on
  `Pi/2 < x < Pi`) still comes through `Refine`.
- **Predicates.** `Element[...]`, equations, inequalities and logical
  combinations are decided to `True` or `False` when the assumptions prove
  or refute them. Domain membership goes to the assumption-aware prover,
  equations to the assumption-aware zero test plus equality substitution,
  and inequalities to a `Reduce`/CAD entailment over the reals: `P` is
  `True` exactly when `assum && !P` is unsatisfiable. A predicate that is
  neither proved nor refuted comes back unchanged.
- **Radicals.** Under the sign facts, a radical of a product or quotient of
  squares loses each squared factor (as its absolute value, then resolved).
  `Integrate` uses this pass to tidy a definite integral computed with
  `Assumptions`.
- **Defaults and budget.** `Refine[expr]` uses `$Assumptions`, including
  whatever an enclosing `Assuming` adds. `Assumptions -> X` replaces
  `$Assumptions`, and a positional assumption is conjoined with it.
  `TimeConstraint` (default 30 s) caps the time spent on a condition check.
  At most 24 `Reduce` entailments run per call, over at most 6 variables.
- `Protected`. Symbolic only, with no NDArray or `Compile[]` path.
- See also: `Simplify`, `Assuming`, `Reduce`, `PossibleZeroQ`.
