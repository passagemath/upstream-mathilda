### Worked examples

```mathematica
In[1]:= ChebyshevU[6, x]  (* the explicit degree-6 polynomial *)
```

```mathematica
In[1]:= {ChebyshevU[-1, x], ChebyshevU[-4, x]}  (* negative orders reflect: U_{-1} = 0, U_{-m} = -U_{m-2} *)
```

```mathematica
In[1]:= Simplify[TrigExpand[ChebyshevU[4, Cos[t]] Sin[t] - Sin[5 t]]]  (* the defining identity U_n(cos t) sin t = sin((n+1) t) *)
```

```mathematica
In[1]:= Expand[ChebyshevT[5, x] - (ChebyshevU[5, x] - ChebyshevU[3, x])/2]  (* the bridge T_n = (U_n - U_{n-2})/2 *)
```

```mathematica
In[1]:= {ChebyshevU[60, 0.7], Sin[61 ArcCos[0.7]]/Sin[ArcCos[0.7]]}  (* a high degree at machine precision via the stable recurrence *)
```

```mathematica
In[1]:= N[ChebyshevU[60, 7/10], 25]  (* the same value at 25 digits *)
```

```mathematica
In[1]:= {ChebyshevU[0.5, 0.3], Sin[1.5 ArcCos[0.3]]/Sin[ArcCos[0.3]]}  (* a non-integer order *)
```

```mathematica
In[1]:= ChebyshevU[3, 2. + I]  (* a machine complex argument uses the complex recurrence *)
```

```mathematica
In[1]:= Integrate[ChebyshevU[3, x]^2 Sqrt[1 - x^2], {x, -1, 1}]  (* the orthogonality norm with weight Sqrt[1 - x^2] *)
```

```mathematica
In[1]:= D[ChebyshevU[n, x], x]  (* the derivative for a symbolic order *)
```

```mathematica
In[1]:= cu = Compile[{{v, _Real, 1}}, ChebyshevU[9, v]]; cu[Range[-1., 1., 0.1]] === ChebyshevU[9, Range[-1., 1., 0.1]]  (* the compiled array result is bit-identical to the packed one *)
```

### Notes

- `ChebyshevU[n, x]` is `U_n(x)`, with `U_n(cos θ) sin θ = sin((n+1)θ)`; it is orthogonal on `[-1, 1]` with weight `Sqrt[1 - x^2]`. Negative orders reflect as in Mathematica.
- Special values hold for a symbolic order: `ChebyshevU[n, 1] = n + 1` and `ChebyshevU[n, 0] = Cos[n Pi/2]`.
- Numeric surfaces: the scalar, packed/`NDArray` (including complex buffers) and `Compile[]` paths share one `double` recurrence. A fractional order goes through `Hypergeometric2F1` or the `ArcCos` closed form.
- Attributes: `Listable`, `NumericFunction`, `Protected`.
- See also: `ChebyshevT`, `GegenbauerC`, `LegendreP`.
