### Worked examples

```mathematica
In[1]:= LaguerreL[4, x]  (* the ordinary Laguerre polynomial with rational coefficients *)
```

```mathematica
In[1]:= LaguerreL[3, a, x]  (* the generalized polynomial with a symbolic parameter *)
```

```mathematica
In[1]:= LaguerreL[3, 1/2, 2]  (* an exact value with a rational parameter *)
```

```mathematica
In[1]:= LaguerreL[n, a, 0]  (* the value at the origin is a binomial coefficient *)
```

```mathematica
In[1]:= Simplify[x D[LaguerreL[4, a, x], {x, 2}] + (a + 1 - x) D[LaguerreL[4, a, x], x] + 4 LaguerreL[4, a, x]]  (* the Laguerre equation is satisfied identically *)
```

```mathematica
In[1]:= Integrate[LaguerreL[3, 2, x]^2 x^2 Exp[-x], {x, 0, Infinity}]  (* the norm Gamma[n + a + 1]/n! = 20 *)
```

```mathematica
In[1]:= Integrate[LaguerreL[2, x] LaguerreL[3, x] Exp[-x], {x, 0, Infinity}]  (* distinct degrees are orthogonal *)
```

```mathematica
In[1]:= {LaguerreL[30, 2, 10.5], N[LaguerreL[30, 2, 21/2], 25]}  (* a degree-30 value at machine and 25-digit precision *)
```

```mathematica
In[1]:= LaguerreL[0.5, 1.2]  (* a fractional order through Binomial[n + a, n] Hypergeometric1F1[-n, a + 1, x] *)
```

```mathematica
In[1]:= LaguerreL[3, 0.5, Range[0., 2., 0.25]]  (* the generalized form over a packed Range input *)
```

```mathematica
In[1]:= CompileDiagnostics[{{v, _Real, 1}}, LaguerreL[5, 0.5, v]]  (* the three-argument form does not yet lower at array shape *)
```

### Notes

- `LaguerreL[n, x]` equals `LaguerreL[n, 0, x]`. `LaguerreL[n, a, x]` is the generalized (associated) polynomial `L_n^(a)(x)`, orthogonal on `(0, ∞)` with weight `x^a e^(-x)`.
- A symbolic parameter `a` gives the expanded polynomial in `a` and `x` (up to degree 200). A numeric `a` gives exact rational coefficients (up to degree 2000).
- A fractional or complex order is evaluated as `Binomial[n + a, n] Hypergeometric1F1[-n, a + 1, x]` when some argument is inexact.
- Numeric surfaces: both forms stay on the buffer for packed and `NDArray` input. `Compile[]` lowers `LaguerreL[n, x]` at scalar and array shapes, but lowers `LaguerreL[n, a, x]` only at scalar shape.
- Attributes: `Listable`, `NumericFunction`, `Protected`.
- See also: `HermiteH`, `Hypergeometric1F1`, `Binomial`.
