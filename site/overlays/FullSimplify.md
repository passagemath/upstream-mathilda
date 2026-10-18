### Worked examples

```mathematica
In[1]:= FullSimplify[Gamma[x + 1]/Gamma[x]]  (* the gamma recurrence, loaded only because Gamma appears *)
```

```mathematica
In[1]:= FullSimplify[x Gamma[x]]  (* an equal-complexity repackaging kept by the final consolidation sweep *)
```

```mathematica
In[1]:= FullSimplify[Gamma[1/3] Gamma[2/3]]  (* reflection pair Gamma[z] Gamma[1 - z] == Pi/Sin[Pi z] *)
```

```mathematica
In[1]:= FullSimplify[Gamma[1 + I] Gamma[1 - I]]  (* conjugate pair on the line Re z == 1 *)
```

```mathematica
In[1]:= FullSimplify[Beta[a, b] Gamma[a + b]]  (* the Beta function is rewritten through Gamma, then cancelled *)
```

```mathematica
In[1]:= FullSimplify[a + Erf[x] + Erfc[x]]  (* complementary error functions inside a larger sum *)
```

```mathematica
In[1]:= FullSimplify[PolyLog[2, z] + PolyLog[2, -z]]  (* dilogarithm duplication formula *)
```

```mathematica
In[1]:= FullSimplify[Log[2 Sec[x] - 2 Tan[x]], -Pi/2 < x < Pi/2]  (* real-log Gudermannian family with a positive coefficient split off *)
```

```mathematica
In[1]:= FullSimplify[(-1)^(2 n + 1) Gamma[n + 1]/Gamma[n], Element[n, Integers]]  (* assumptions reach the Simplify core *)
```

```mathematica
In[1]:= FullSimplify[Gamma[x + 2]/Gamma[x]]  (* limitation: only a literal +1 shift is matched *)
```

### Notes

- **A wrapper around `Simplify`.** `FullSimplify[expr, assum, opts]` loads
  the identity libraries for the heads that appear in `expr` and passes
  them to `Simplify` as `TransformationFunctions -> {Automatic, …}`.
  Because `Automatic` is always kept, the result is never more complex than
  `Simplify[expr, assum]`. An input that contains none of the library heads
  (`Gamma`, `LogGamma`, `PolyGamma`, `Pochhammer`, `Beta`, `Factorial`,
  `Erf`, `Erfc`, `PolyLog`, `Surd`, `Log`) costs the same as `Simplify`.
- **Options.** `ComplexityFunction` and `TransformationFunctions` (user
  transforms are merged with the relevant ones) are passed through.
  `TimeConstraint -> {tLoc, tTot}` caps each transform at `tLoc` seconds and
  the whole call at `tTot`; a bare `t` means `{t, t}`. On a total timeout
  the plain `Simplify` result comes back.
- **The log rules need a positive argument.** The `Log[Sec ± Tan]`,
  `Log[Csc ± Cot]`, `Log[Coth ± Csch]` and half-angle `Log[Tan[x/2]]` rules
  are identities only where the log argument is positive. The examples
  above therefore state a region such as `-Pi/2 < x < Pi/2`. At present the
  rules also fire without that assumption, and where the argument is
  negative they then drop an `I Pi` term.
- **Coverage is growing.** The gamma recurrence matches a literal `+1`
  shift only, and each family has its own file under
  `src/internal/simp/transforms/`.
- `Protected`. There is no numeric surface.
- See also: `Simplify`, `Refine`, `FunctionExpand`, `Gamma`, `PolyLog`.
