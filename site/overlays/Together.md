---
status: Stable
references:
  - "Geddes, Czapor & Labahn, \"Algorithms for Computer Algebra\" (1992), on rational function arithmetic and common denominators."
  - "von zur Gathen & Gerhard, \"Modern Computer Algebra\", on polynomial GCDs used in cancellation."
---
### Worked examples

```mathematica
In[1]:= Together[1/x + 1/y]  (* distinct denominators are multiplied *)
```

```mathematica
In[1]:= Together[1/(x - 1) - 2/(x^2 - 1)]  (* a common factor hidden across the denominators is cancelled *)
```

```mathematica
In[1]:= Together[a/(a - b) + b/(b - a)]  (* a telescoping sum collapses to 1 *)
```

```mathematica
In[1]:= Together[E^x/(E^x + 1) + 1/(E^x + 1)]  (* the exponential kernel is a generator of the rational function field *)
```

```mathematica
In[1]:= Together[1/(x - Sqrt[2]) + Sqrt[2]/(x^2 - 2), Extension -> Sqrt[2]]  (* combination over Q(Sqrt[2]) *)
```

```mathematica
In[1]:= Together[Log[x^2] + 1/x]  (* logarithms stay opaque: Log[x^2] is not rewritten to 2 Log[x] *)
```

```mathematica
In[1]:= Together[E^Cos[x] Integrate[E^-Cos[x], x]]  (* an unevaluated integral is an atom, never re-evaluated with a symbol inside *)
```

```mathematica
In[1]:= Together[{1/x + 1, x == 1/y + 1}]  (* threads over lists and equations *)
```

```mathematica
In[1]:= Together[0.5/x + 1/(2 x)]  (* inexact input is rationalised, combined and numericalised again *)
```

### Notes

- `Together[expr]` puts a sum of terms over a single common denominator
  and cancels the GCD of the result. It is the inverse of `Apart`. Unlike a
  full expansion, it does not multiply out `Power[Plus, n]` factors that it
  does not need to.
- **Kernels.** Transcendental subexpressions (`Log[...]`, `E^f`, inverse
  trig) and algebraic constants (`I`, radicals) are treated as generators
  of the rational function field. Exponentials with commensurate exponents
  share one generator. Logarithms are left as they are, because rewriting
  `Log[x^2]` as `2 Log[x]` is false for `x < 0`.
- **Binding operators are atoms.** Unevaluated `Integrate`, `D`, `Sum`,
  `Limit`, `Function`, `Hold*` and `Inactive[op][...]` are never looked
  inside.
- **Extension.** `Extension -> α`, a list, or `Automatic` combines over the
  number field. Plain algebraic constants are often handled without it.
- `Listable`, `Protected`. Symbolic only.
- See also: `Cancel`, `Apart`, `Factor`, `Numerator`, `Denominator`,
  `Simplify`.
