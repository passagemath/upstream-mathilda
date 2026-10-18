### Worked examples

```mathematica
In[1]:= HermiteH[6, x]  (* the physicists' polynomial with integer coefficients *)
```

```mathematica
In[1]:= HermiteH[n, 0]  (* the value at the origin for a symbolic order *)
```

```mathematica
In[1]:= Simplify[D[HermiteH[5, x], {x, 2}] - 2 x D[HermiteH[5, x], x] + 10 HermiteH[5, x]]  (* the Hermite equation y'' - 2x y' + 2n y = 0 *)
```

```mathematica
In[1]:= Integrate[HermiteH[3, x]^2 Exp[-x^2], {x, -Infinity, Infinity}]  (* the norm 2^n n! Sqrt[Pi] *)
```

```mathematica
In[1]:= Integrate[HermiteH[3, x] HermiteH[4, x] Exp[-x^2], {x, -Infinity, Infinity}]  (* distinct degrees are orthogonal *)
```

```mathematica
In[1]:= HermiteH[40, 2.5]  (* a degree-40 value at machine precision *)
```

```mathematica
In[1]:= N[HermiteH[40, 5/2], 20]  (* the same value at 20 digits *)
```

```mathematica
In[1]:= HermiteH[1/2, 0.8]  (* a fractional order through the 1F1 form *)
```

```mathematica
In[1]:= {HermiteH[-1, 0.5], Sqrt[Pi] Exp[0.25] Erfc[0.5]/2}  (* a negative order: H_{-1}(x) = Sqrt[Pi] e^(x^2) Erfc[x]/2 *)
```

```mathematica
In[1]:= HermiteH[4, Range[0., 2., 0.25]]  (* a packed Range input runs on the buffer kernel *)
```

```mathematica
In[1]:= Series[HermiteH[n, x], {x, Infinity, 2}]  (* the algebraic asymptotic series at infinity *)
```

```mathematica
In[1]:= HermiteH[1/2, x]  (* a fractional order with a symbolic argument stays unevaluated *)
```

### Notes

- `HermiteH[n, x]` is the physicists' Hermite polynomial, orthogonal on `(-∞, ∞)` with weight `e^(-x^2)`. For the probabilists' `He_n(x)`, use `2^(-n/2) HermiteH[n, x/Sqrt[2]]`.
- A non-integer, complex or negative order with an inexact argument goes through a two-term `Hypergeometric1F1` form, evaluated with extra working precision that scales with `Re(x^2)` to cover the cancellation.
- Numeric surfaces: the scalar, packed/`NDArray` (real and complex) and `Compile[]` paths share one `double` recurrence for integer orders up to 1000.
- The asymptotic series at `Infinity` omits the subdominant `e^(x^2)` branch.
- Attributes: `Listable`, `NumericFunction`, `Protected`.
- See also: `LaguerreL`, `Hypergeometric1F1`, `Erfc`.
