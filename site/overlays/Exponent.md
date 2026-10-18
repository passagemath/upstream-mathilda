### Worked examples

```mathematica
In[1]:= Exponent[1 + x^2 + x^5, x]  (* the highest power of x *)
```

```mathematica
In[1]:= Exponent[(1 + x)^3, x]  (* the argument is expanded first, so it need not be given expanded *)
```

```mathematica
In[1]:= Exponent[1/x + x^2, x, List]  (* a custom h collects the whole exponent set, negative powers included *)
```

```mathematica
In[1]:= Exponent[1 + x + x^2, x, Min]  (* with h = Min, the lowest power *)
```

```mathematica
In[1]:= Exponent[x^(1/2) + x^(3/2), x]  (* fractional exponents *)
```

```mathematica
In[1]:= Exponent[x^n + x^2, x]  (* a symbolic exponent stays inside Max *)
```

```mathematica
In[1]:= Exponent[a x^2 + b x y^3, {x, y}]  (* a list of forms threads *)
```

```mathematica
In[1]:= Exponent[2 Sin[x]^3 + Sin[x], Sin[x]]  (* the form may be any kernel *)
```

```mathematica
In[1]:= Exponent[0, x]  (* the zero polynomial: Max of the empty set *)
```

### Notes

- `Exponent[expr, form]` gives the largest power of `form` in the expanded
  `expr`, applying `h` (default `Max`) to the set of exponents found.
  `Exponent[expr, form, h]` substitutes any other `h`, so `h = List`
  returns the sorted, de-duplicated exponent set. `form` may be a symbol, a
  kernel or a product (`Exponent[(x y)^2 + x^3 y, x y]` is `2`).
- The reading is **purely syntactic**, with no zero-coefficient
  recognition. The genuine zero polynomial has an empty exponent set, so
  `Exponent[0, x]` is `Max[]`, which is `-Infinity`.
- Since v0.278 a non-flat nested `Times` inside a term is searched
  recursively, so `Exponent[-a x^2 + I x, x]` is `2`.
- `Listable`, `Protected`. Symbolic, with no packed or `Compile[]` path.
- See also: `Coefficient`, `CoefficientList`, `Variables`, `Expand`.
