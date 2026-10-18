---
status: Stable
references:
  - "von zur Gathen & Gerhard, \"Modern Computer Algebra\", on polynomial GCD computation."
  - "Geddes, Czapor & Labahn, \"Algorithms for Computer Algebra\" (1992), on rational function simplification."
---
### Worked examples

```mathematica
In[1]:= Cancel[(x^2 - 1)/(x^2 - 2 x + 1)]  (* the common factor x - 1 is divided out; the rest stays a reduced quotient *)
```

```mathematica
In[1]:= Cancel[(x^4 - 1)/(x^2 - 1)]  (* exact polynomial division by the GCD *)
```

```mathematica
In[1]:= Cancel[(x - y)^2/(x^2 - y^2) + a]  (* each term of a sum is cancelled; the sum is not combined *)
```

```mathematica
In[1]:= Cancel[(x^3 - 1)/(x - 1) == y]  (* threads over equations *)
```

```mathematica
In[1]:= Cancel[(x^2 - 2)/(x - Sqrt[2]), Extension -> Sqrt[2]]  (* cancellation over Q(Sqrt[2]) *)
```

```mathematica
In[1]:= Cancel[(x^2 + 1)/(x - I), Extension -> Automatic]  (* the extension Q(I) is detected automatically *)
```

```mathematica
In[1]:= Cancel[(E^(2 x) - 1)/(E^x - 1)]  (* both exponentials reduce to the single generator E^x *)
```

```mathematica
In[1]:= Cancel[(Sin[x]^2 - 1)/(Sin[x] - 1)]  (* a trig kernel is an opaque generator *)
```

```mathematica
In[1]:= Cancel[Integrate[f[x], x] (x^2 - 1)/(x - 1)]  (* an unevaluated integral is kept as an atom *)
```

```mathematica
In[1]:= Cancel[(x^2 - 1.)/(x - 1)]  (* inexact coefficients round-trip *)
```

### Notes

- `Cancel[expr]` divides out the common polynomial factors of a numerator
  and denominator. It does not split into partial fractions (`Apart`), and
  it does not combine a sum of fractions (`Together`): each term of a sum
  is cancelled on its own. The surviving factors are not expanded again.
- **Generators.** Exponential, trig, log and other kernels become
  independent generators, with `E^(a+b)` split so related exponentials
  agree. Algebraic constants carry their minimal polynomial, so
  `(x^2 - 2)/(x - Sqrt[2])` cancels to `Sqrt[2] + x` even without
  `Extension`. `Extension -> α`, a list, or `Automatic` makes the field
  explicit.
- **Value-preserving.** Since v0.434, unevaluated `Integrate`, `D`, `Sum`
  and similar heads are never substituted into, and logarithms are not
  expanded.
- `Listable`, `Protected`. Symbolic only.
- See also: `Together`, `Factor`, `PolynomialGCD`, `Apart`, `Simplify`.
