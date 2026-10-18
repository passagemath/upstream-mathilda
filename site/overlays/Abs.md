### Worked examples

```mathematica
In[1]:= Abs[-5]
```

```mathematica
In[1]:= Abs[3 + 4 I]  (* the modulus of a Gaussian integer, exactly *)
```

```mathematica
In[1]:= Abs[{-1, 2, -3/4}]  (* Listable *)
```

```mathematica
In[1]:= Abs[(1 + I)^10]  (* |1 + I|^10 = Sqrt[2]^10 *)
```

```mathematica
In[1]:= Abs[Sqrt[2] + Sqrt[3] I]  (* the modulus of an exact complex expression *)
```

```mathematica
In[1]:= Abs[-9223372036854775807 - 1]  (* |-2^63| does not fit an int64, so it becomes a bignum *)
```

```mathematica
In[1]:= Abs[Root[#^3 - 2 # + 1 &, 1]]  (* a negative real Root stays an exact algebraic number *)
```

```mathematica
In[1]:= Abs[Interval[{-2, 1}]]  (* interval arithmetic *)
```

```mathematica
In[1]:= Abs[N[Gamma[1/3 + 2 I], 40]]  (* an MPFR complex goes through mpfr_hypot at 40 digits *)
```

```mathematica
In[1]:= Abs[NDArray[{-1, 2, -3}, DataType -> "int64"]]  (* the int64 kernel keeps integers exact *)
```

```mathematica
In[1]:= Abs[Sqrt[2]]  (* not yet evaluated: exact real constants other than rationals stay as they are *)
```

### Notes

- `Abs[z]` is the absolute value of a real number and the modulus `Sqrt[Re[z]^2 + Im[z]^2]` of a complex one. Exact inputs give exact results.
- Limitation: an exact real constant that is not rational, such as `Sqrt[2]`, `-Pi` or `Pi - 4`, stays unevaluated (Mathematica gives `Sqrt[2]`, `Pi`, `4 - Pi`). `Abs[N[...]]` is numeric.
- Symbolic arguments stay unevaluated (`Abs[x]`). There is no `ComplexExpand`-style assumption that symbols are real.
- `Listable`, `NumericFunction`, `Protected`. There is a packed/NDArray kernel with an exact int64 arm (and a `hypot` arm for real and complex data), and `Compile[]` lowers it at scalar and array shapes (`Abs` of a Complex is Real).
- See also: `Sign`, `Arg`, `Norm`, `RealAbs`.
