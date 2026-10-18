### Worked examples

```mathematica
In[1]:= WhittakerM[2, 3, 4.5]  (* machine value *)
```

```mathematica
In[1]:= Exp[-2.25] 4.5^3.5 Hypergeometric1F1[1.5, 7, 4.5]  (* the same value from the defining 1F1 form *)
```

```mathematica
In[1]:= WhittakerM[0, 1/2, 2.] - 2 Sinh[1.]  (* checked against M_{0,1/2}(z) = 2 sinh(z/2) *)
```

```mathematica
In[1]:= N[WhittakerM[2, 1/3, 1/7], 30]  (* arbitrary precision through Hypergeometric1F1 *)
```

```mathematica
In[1]:= WhittakerM[1, 1/3, -2.]  (* negative z: the principal branch of z^(m+1/2) *)
```

```mathematica
In[1]:= WhittakerM[2, 3 - I, 4.5]  (* complex order *)
```

```mathematica
In[1]:= WhittakerM[1, 1, 0]  (* zero at the origin when Re m > -1/2 *)
```

```mathematica
In[1]:= D[WhittakerM[k, m, z], z]  (* DLMF 13.4.27 *)
```

```mathematica
In[1]:= Series[WhittakerM[2, 1/2, x], {x, 0, 6}]  (* origin series *)
```

```mathematica
In[1]:= Series[WhittakerM[2, 1/2, x], {x, Infinity, 5}]  (* at infinity only the terminating case is expanded: here the exact closed form *)
```

### Notes

- `Listable`, `NumericFunction`, `Protected`. The definition is `M_{k,m}(z) = e^(−z/2) z^(m+1/2) 1F1(m−k+1/2; 2m+1; z)`.
- Exact input stays symbolic, with one exception: `M_{k,m}(0) = 0` when Re m > −1/2.
- Machine input (k and m real; z real or complex) runs the shared double kernel `sf_machine_whittaker_m` or its complex twin, and so do buffers. A real z < 0 gives the principal complex value. If the kernel declines, or for complex k, m or precision above machine, the 1F1 form is evaluated at the input precision.
- `D` in z follows DLMF 13.4.27. Derivatives in k and m stay inert.
- `Series` at 0 gives the `z^(m+1/2)` Puiseux form, for symbolic m as well. At `Infinity` only the terminating (elementary) case is expanded. `FunctionExpand` does not expand `WhittakerM`, and an `Interval` argument stays unevaluated.
- See also: `WhittakerW`, `Hypergeometric1F1`.
