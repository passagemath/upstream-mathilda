### Worked examples

```mathematica
In[1]:= Arg[-1]  (* the principal value lies in (-Pi, Pi] *)
```

```mathematica
In[1]:= Arg[1 + I]
```

```mathematica
In[1]:= Table[Arg[(1 + I)^k], {k, 0, 8}]  (* powers of 1 + I go once round the circle and fold back at Pi *)
```

```mathematica
In[1]:= Arg[(-1)^(1/6)]  (* an exact root of unity gives an exact multiple of Pi *)
```

```mathematica
In[1]:= Arg[-(-1)^(5/6)]
```

```mathematica
In[1]:= Arg[-1/10^400]  (* read from the exact sign, so it does not underflow *)
```

```mathematica
In[1]:= Arg[2 + 3 I]  (* a generic Gaussian gives the two-argument ArcTan *)
```

```mathematica
In[1]:= N[Arg[2 + 3 I], 40]
```

```mathematica
In[1]:= Arg[N[-1 + I, 30]]  (* an MPFR argument goes through mpfr_atan2 at its own precision *)
```

```mathematica
In[1]:= Arg[NDArray[{1. + I, -1., -I}]]  (* the NDArray kernel *)
```

### Notes

- `Arg[z]` is the principal argument, in `(-Pi, Pi]`. It is `0` for a positive real, `Pi` for a negative real, and `0` at the origin.
- Exact Gaussians on an axis or a diagonal give exact multiples of `Pi`, and so do products of rational powers of negative numbers such as `(-1)^(1/6)`. Any other exact Gaussian gives `ArcTan[re, im]`, which equals the `ArcTan[im/re]` Mathematica prints when `re > 0`.
- Exact real constants that are not rational stay unevaluated: `Arg[Pi - 4]` should be `Pi`. Apply `N` first.
- `Listable`, `NumericFunction`, `Protected`. There is an NDArray/packed kernel (`atan2`, real result), and `Compile[]` lowers `Arg` at scalar and array shapes.
- See also: `Abs`, `Sign`, `ArcTan`, `AbsArg`.
