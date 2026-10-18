### Worked examples

```mathematica
In[1]:= ChebyshevT[8, x]  (* the explicit degree-8 polynomial with integer coefficients *)
```

```mathematica
In[1]:= ChebyshevT[-3, x]  (* the order index is even: T_{-n} = T_n *)
```

```mathematica
In[1]:= Simplify[TrigExpand[ChebyshevT[5, Cos[t]] - Cos[5 t]]]  (* the defining identity T_n(cos t) = cos(n t) *)
```

```mathematica
In[1]:= Expand[ChebyshevT[3, ChebyshevT[4, x]] - ChebyshevT[12, x]]  (* the composition law T_m(T_n(x)) = T_{mn}(x) *)
```

```mathematica
In[1]:= {ChebyshevT[60, 0.7], Cos[60 ArcCos[0.7]]}  (* a high degree at machine precision via the stable recurrence *)
```

```mathematica
In[1]:= N[ChebyshevT[60, 7/10], 30]  (* the same value to 30 digits through the MPFR recurrence *)
```

```mathematica
In[1]:= {ChebyshevT[1/3, 0.5], Cos[1/3 ArcCos[0.5]]}  (* a non-integer order, through 2F1 / Cos[v ArcCos[x]] *)
```

```mathematica
In[1]:= ChebyshevT[2.5, 3.]  (* outside [-1, 1] a fractional order grows like Cosh[v ArcCosh[x]] *)
```

```mathematica
In[1]:= Integrate[ChebyshevT[4, x]^2/Sqrt[1 - x^2], {x, -1, 1}]  (* the orthogonality norm, pi/2 for n >= 1 *)
```

```mathematica
In[1]:= ChebyshevT[5, Range[0., 1., 0.25]]  (* a packed Range input stays on the buffer kernel *)
```

```mathematica
In[1]:= CompileDiagnostics[{{v, _Real, 1}}, ChebyshevT[7, v]]  (* the head lowers in Compile at array shape *)
```

```mathematica
In[1]:= ChebyshevT[1/2, x]  (* an exact fractional order with a symbolic argument stays unevaluated *)
```

### Notes

- `ChebyshevT[n, x]` is `T_n(x)`, defined by `T_n(cos θ) = cos(nθ)` and orthogonal on `[-1, 1]` with weight `1/Sqrt[1 - x^2]`. It agrees with Mathematica, including the branch cut on `(-∞, -1)` for non-integer orders.
- An integer order with an exact argument gives the exact polynomial; with an inexact argument it uses the stable recurrence, never the expanded polynomial. A fractional order is evaluated numerically only when some argument is inexact.
- Numeric surfaces: the machine `double` recurrence (real and complex) is shared by the scalar, packed/`NDArray` and `Compile[]` paths, so all of them agree bit for bit. An arbitrary-precision argument, a degree above 1000, or an overflow uses the MPFR recurrence.
- Attributes: `Listable`, `NumericFunction`, `Protected`.
- See also: `ChebyshevU`, `LegendreP`, `Hypergeometric2F1`.
