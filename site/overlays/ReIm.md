### Worked examples

```mathematica
In[1]:= ReIm[3 + 4 I]
```

```mathematica
In[1]:= ReIm[(2 + I)^3]  (* (2 + I)^3 = 2 + 11 I *)
```

```mathematica
In[1]:= ReIm[(3 + 4 I)/(1 - 2 I)]  (* a complex quotient reduces to Gaussian integers *)
```

```mathematica
In[1]:= ReIm[N[E^(I Pi/4), 20]]  (* both parts are 1/Sqrt[2], to 20 digits *)
```

```mathematica
In[1]:= ReIm[{1 + 2 I, 3.5}]  (* Listable: one pair per element *)
```

```mathematica
In[1]:= ReIm[Abs[z]]  (* a real-valued head gives {f[z], 0} *)
```

```mathematica
In[1]:= ReIm[z]  (* symbolic z stays unevaluated *)
```

### Notes

- `ReIm[z]` is `{Re[z], Im[z]}`, and a real argument gives `{z, 0}`.
- Unlike `Re` and `Im`, `ReIm` has no packed/NDArray kernel and no `Compile[]` lowering. `ReIm[NDArray[...]]` stays unevaluated, and a compiled body that uses it runs in the interpreter (a recorded backlog item). Use `{Re[v], Im[v]}` on large arrays or inside `Compile[]`.
- `Listable`, `NumericFunction`, `Protected`.
- See also: `Re`, `Im`, `AbsArg`, `Conjugate`.
