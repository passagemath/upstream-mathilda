---
status: Stable
references:
---
### Worked examples

```mathematica
In[1]:= Simplify[(x^2 - 1)/(x - 1)]  (* rational cancellation through the SHAPE_RATIONAL fast path *)
```

```mathematica
In[1]:= Simplify[Tan[x] + Tan[x + Pi/3] + Tan[x + 2 Pi/3] - 3 Tan[3 x]]  (* constant-phase angle expansion proves the triple-tangent identity *)
```

```mathematica
In[1]:= Simplify[Tan[Pi/7] Tan[2 Pi/7] Tan[3 Pi/7]]  (* trig at rational multiples of Pi is folded as an exact algebraic number *)
```

```mathematica
In[1]:= Simplify[4 ArcTan[1/5] - ArcTan[1/239] - Pi/4]  (* Machin's formula, certified by exp(i e) == 1 in qqbar *)
```

```mathematica
In[1]:= Simplify[{(-1)^(n^2), I^(4 n + 2), Floor[n/2] + Ceiling[n/2], Gamma[n + 2]/Gamma[n], GCD[n, n + 1]}, Element[n, Integers]]  (* integer-assumption rules: parity, roots of unity, floor division, Gamma ratios, coprimality *)
```

```mathematica
In[1]:= Simplify[Abs[x - 2] + Abs[x + 3], 0 < x < 1]  (* the Reduce-backed sign oracle resolves both absolute values *)
```

```mathematica
In[1]:= Simplify[ArcCos[1 - 2 x^2], 0 < x < 1]  (* multiple-angle collapse, emitted only after Reduce proves it on the region *)
```

```mathematica
In[1]:= Simplify[Sqrt[E^(2 I x)], -Pi/2 < x < Pi/2]  (* branch-gated power fold: Im of the exponent stays in the principal strip *)
```

```mathematica
In[1]:= Simplify[Abs[Cos[x]], 0 < x < Pi]  (* honest decline: Cos changes sign at Pi/2 *)
```

```mathematica
In[1]:= Simplify[Expand[(x + y)^20]/(x + y)^18, TimeConstraint -> 5]  (* a scalar TimeConstraint bounds the whole call *)
```

### Notes

- **Search, not rewriting.** `Simplify` seeds a candidate set with the input
  and the output of many exact rewriters. It then runs `Together`, `Cancel`,
  `Expand`, `Factor`, `Apart`, `TrigExpand`, `TrigReduce`, `TrigToExp` and
  related transforms on each candidate, and returns the lowest-scoring form
  (`SimplifyCount` by default, or any `ComplexityFunction`). A correct
  rewrite that does not lower the score can still lose: on
  `Pi/2 < x < Pi`, `Simplify[ArcSin[Sin[x]]]` stays as it is, while
  `Refine` returns `Pi - x`.
- **Assumptions.** The assumption argument is a positional second argument
  or `Assumptions -> …`; a positional one is conjoined with `$Assumptions`
  (and with `Assuming`). The rules that depend on assumptions are sound by
  construction. Each needs its assumption proved, by the structural provers
  or a `Reduce`/CAD entailment, and declines otherwise: `(-1)^(2 n)` does not
  reduce for real `n`, and `Abs[Cos[x]]` on `0 < x < Pi` is left alone.
  `(a^p)^q -> a^(p q)` needs both `a > 0` and `p` real. Integer-only folds such as `Gamma[n + 2]/Gamma[n] -> n (1 + n)` fire only
  under an integer assumption.
- **Not a decision procedure for predicates.** `Simplify[x > 0, x > 1]`
  stays `x > 0`. Use `Refine` (or `Reduce`) to decide a relation under
  assumptions.
- **TimeConstraint.** A scalar `t` caps the whole call. `{tLoc, tTot}` is
  also accepted, with `tTot` as the whole-call cap. The default is
  `Infinity`.
- **Threading.** `Simplify` threads over lists, equations and logical
  connectives itself. Its only attribute is `Protected`. Inexact input is rationalised, simplified, and
  numericalised again.
- **No numeric surfaces.** This head is symbolic, with no NDArray or
  `Compile[]` path.
- See also: `FullSimplify`, `Refine`, `PossibleZeroQ`, `Together`, `Factor`,
  `TrigReduce`, `Assuming`.
