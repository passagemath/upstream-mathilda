### Worked examples

```mathematica
In[1]:= ProductLog[E]  (* since e = 1 * e^1, the principal value is 1 *)
```

```mathematica
In[1]:= ProductLog[-1/E]  (* the branch point where the k = 0 and k = -1 branches meet *)
```

```mathematica
In[1]:= ProductLog[-Pi/2]  (* an exact complex value *)
```

```mathematica
In[1]:= N[ProductLog[1], 30]  (* the omega constant *)
```

```mathematica
In[1]:= ProductLog[-1, -0.2]  (* the other real branch, k = -1 *)
```

```mathematica
In[1]:= ProductLog[-0.5]  (* below -1/e the principal value is complex *)
```

```mathematica
In[1]:= ProductLog[1. + 2. I]  (* complex argument, machine complex kernel *)
```

```mathematica
In[1]:= ProductLog[3, 1.5]  (* a higher complex branch *)
```

```mathematica
In[1]:= ProductLog[{0.5, 1., 10.}]  (* Listable *)
```

```mathematica
In[1]:= ProductLog[-1, -0.1`40]  (* arbitrary precision on the k = -1 branch *)
```

```mathematica
In[1]:= Series[ProductLog[x], {x, 0, 4}]  (* the Lagrange-inversion series *)
```

```mathematica
In[1]:= D[ProductLog[x], x]
```

### Notes

- `ProductLog[z]` is the Lambert W function, the principal branch `W_0(z)`: the
  solution `w` of `z = w e^w`. `ProductLog[k, z]` selects the `k`-th branch for
  an explicit integer `k`; `k = 0` is the principal branch and `k = -1` is the
  other real branch.
- A few exact values are returned directly (`0`, `E`, `-1/E`, `-Pi/2`). Every
  other exact argument stays symbolic.
- Machine-precision principal-branch values, real or complex, come from one
  raw-`double` kernel. It handles the branch point `-1/e` carefully (the shift
  `z + 1/e` is formed without cancellation) and checks every result. The same
  kernel serves packed arrays, `NDArray` arguments and `Compile[]` (scalar,
  complex and rank-1 array), so these agree bit for bit. Other branches, and
  arbitrary precision, use a complex-MPFR Halley iteration.
- The result is real on the real-valued domains (`W_0` for `x >= -1/e`, `W_{-1}`
  for `-1/e <= x < 0`) and `Complex` otherwise.
- The branch index must be an `Integer`. `N[ProductLog[-1, -1/10], 40]` turns
  the index into `-1.` and stays unevaluated. Write `ProductLog[-1, N[-1/10, 40]]`
  instead.
- Attributes `Listable`, `NumericFunction`, `Protected`.
- See also: `Log`, `Exp`, `Solve`.
