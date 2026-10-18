### Worked examples

```mathematica
In[1]:= Sign[-7]
```

```mathematica
In[1]:= Sign[{-2, 0, 5}]  (* Listable *)
```

```mathematica
In[1]:= Sign[-2.5]  (* the sign of an inexact real is still the exact integer -1 *)
```

```mathematica
In[1]:= Sign[-1/10^25]  (* a bignum rational *)
```

```mathematica
In[1]:= Sign[3 + 4 I]  (* a complex number gives the unit vector z/Abs[z] *)
```

```mathematica
In[1]:= Sign[(1 + I)^2]  (* (1 + I)^2 = 2 I points straight up *)
```

```mathematica
In[1]:= Sign[2 - 2 I]
```

```mathematica
In[1]:= Sign[N[1 + I, 25]]  (* MPFR components keep their precision *)
```

```mathematica
In[1]:= Sign[Interval[{-1, 2}]]  (* an interval maps monotonically *)
```

```mathematica
In[1]:= Sign[NDArray[{-1.5, 0., 2.}]]  (* the narrowing kernel writes an int64 buffer *)
```

```mathematica
In[1]:= Sign[Pi - 3]  (* not yet evaluated: exact real constants other than rationals stay as they are *)
```

### Notes

- For a real number `Sign` is the exact integer `-1`, `0` or `1`, even for an inexact input (`Sign[-2.5]` is `-1`, as in Mathematica). For a nonzero complex `z` it is `z/Abs[z]`, and `Sign[0]` is `0`.
- Limitation: an exact real constant that is not rational (`Sign[Pi]`, `Sign[Sqrt[2]]`, `Sign[Pi - 3]`) stays unevaluated, where Mathematica gives `1`. `Sign[N[expr]]` still returns an exact integer.
- `Listable`, `NumericFunction`, `Protected`. The packed/NDArray kernel narrows real and int64 data to int64 (`Sign` is on `pack.c`'s `INT64_OK` list), and `Compile[]` gives an Integer result for real arguments and a Complex one for complex arguments.
- See also: `Abs`, `UnitStep`, `Arg`, `Positive`.
