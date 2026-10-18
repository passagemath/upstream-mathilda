### Worked examples

```mathematica
In[1]:= FactorTerms[6 x^2 + 4 x]  (* the integer content is pulled out; the polynomial part is not factored *)
```

```mathematica
In[1]:= FactorTerms[x^2/2 + x/3]  (* rational coefficients: the content is a rational number *)
```

```mathematica
In[1]:= FactorTerms[3 x^2 y + 6 x y^2, x]  (* with a variable: every factor free of x is pulled out *)
```

```mathematica
In[1]:= FactorTerms[a x - a y + a, {x, y}]  (* a list of variables: content with respect to all of them *)
```

```mathematica
In[1]:= f = 2 x^2 y z + 2 x^2 y + 4 x^2 z + 4 x^2 + 4 y^2 z^2 + 4 z y^2 + 8 z^2 y + 2 z y - 6 y - 12 z - 12; FactorTerms[f, {x, y}]  (* contents peeled for the shrinking subsets of the variable list *)
```

```mathematica
In[1]:= FactorTerms[7 x + (14 y + 21)/z]  (* rational functions round-trip through Together *)
```

```mathematica
In[1]:= FactorTerms[2 x^2 - 4 == 6 y]  (* threads over equations *)
```

### Notes

- `FactorTerms[poly]` pulls out the overall numerical factor, the GCD of
  the coefficients, without touching the polynomial structure. Unlike
  `Factor`, it never splits `1 + 2 x + x^2` into `(1 + x)^2`.
- `FactorTerms[poly, x]` also pulls out every factor that does not depend
  on `x`, that is, the content of `poly` with respect to `x`.
  `FactorTerms[poly, {x1, …, xn}]` repeats this for the subsets
  `{x1, …, xn}`, …, `{x1}`.
- **Never a zero content.** Since v0.277, a nonzero input always has a unit
  content. A numerator that `Together` left with a Laurent term used to
  produce a content of `0`, which multiplied the whole result down to `0`.
- Threads over lists, equations, inequalities and logical connectives,
  though the attribute set is just `Protected`. Symbolic only.
- See also: `FactorTermsList`, `Factor`, `FactorSquareFree`, `Together`.
