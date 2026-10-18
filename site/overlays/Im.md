### Worked examples

```mathematica
In[1]:= Im[3 + 4 I]
```

```mathematica
In[1]:= Im[7]
```

```mathematica
In[1]:= Im[Sqrt[-4]]  (* Sqrt[-4] evaluates to 2 I *)
```

```mathematica
In[1]:= Im[Log[-1]]  (* the principal branch: Log[-1] = I Pi *)
```

```mathematica
In[1]:= Im[{1 + 2 I, 3 - 4 I, 5}]  (* Listable *)
```

```mathematica
In[1]:= Im[Abs[z]]  (* a real-valued head gives 0 *)
```

```mathematica
In[1]:= Im[Gamma[1 + I]]  (* no closed form, so it stays symbolic *)
```

```mathematica
In[1]:= N[Im[Gamma[1 + I]], 40]  (* but it evaluates numerically *)
```

```mathematica
In[1]:= Im[NDArray[{1, 2}, DataType -> "int64"]]  (* the narrowing kernel gives exact zeros for integer data *)
```

### Notes

- `Im[z]` is the imaginary part of a numeric `z`. It is `0` for real numbers, for real `Root` and `AlgebraicNumber` objects, and for calls to the real-valued heads `Re`, `Im`, `Abs`, `Arg`.
- Symbols are not assumed real (`Im[x]` stays unevaluated). An exact real constant with no complex part written out, such as `Im[Pi - 4]`, also stays unevaluated at present (Mathematica gives `0`).
- `Listable`, `NumericFunction`, `Protected`. The packed/NDArray kernel gives an int64 zero buffer for real or integer data and the imaginary parts for complex data. `Compile[]` lowers `Im` at scalar and array shapes.
- See also: `Re`, `ReIm`, `Conjugate`, `ComplexExpand`.
