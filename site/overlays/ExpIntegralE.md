### Worked examples

```mathematica
In[1]:= ExpIntegralE[0, z]  (* order 0 is elementary *)
```

```mathematica
In[1]:= ExpIntegralE[3, 0]  (* the value at the origin is 1/(n - 1) for n >= 2 *)
```

```mathematica
In[1]:= {ExpIntegralE[2, 3/2], N[ExpIntegralE[2, 3/2]]}  (* an exact argument stays symbolic until N *)
```

```mathematica
In[1]:= ExpIntegralE[1, 2.] == Gamma[0, 2.]  (* the order-1 case is the incomplete gamma Gamma[0, z] *)
```

```mathematica
In[1]:= ExpIntegralE[3, -1.]  (* a negative real argument lands on the upper side of the cut *)
```

```mathematica
In[1]:= ExpIntegralE[1, 2. + I]  (* a machine complex argument *)
```

```mathematica
In[1]:= D[ExpIntegralE[n, z], z]  (* the derivative lowers the order *)
```

```mathematica
In[1]:= Integrate[x^(s - 1) ExpIntegralE[2, x], {x, 0, Infinity}, Assumptions -> s > 0]  (* the Mellin transform Gamma[s]/(s + n - 1) *)
```

```mathematica
In[1]:= ce = Compile[{{v, _Real, 1}}, ExpIntegralE[3, v]]; ce[Range[0.5, 5., 0.5]] === ExpIntegralE[3, Range[0.5, 5., 0.5]]  (* the compiled array result matches the packed one bit for bit *)
```

```mathematica
In[1]:= {ExpIntegralE[5/2, 1.2], 1.2^(3/2) Gamma[-3/2, 1.2]}  (* a non-integer order stays unevaluated; E_n(z) = z^(n-1) Gamma[1 - n, z] gives the value *)
```

### Notes

- `ExpIntegralE[n, z] = ∫₁^∞ e^(-z t)/t^n dt`, so `ExpIntegralE[1, z] = Gamma[0, z]`. In general `E_n(z) = z^(n-1) Gamma[1-n, z]`.
- Only an integer order `n ≥ 0` is evaluated numerically, and only at machine precision. Mathematica also handles non-integer orders and arbitrary precision; here `ExpIntegralE[5/2, 1.2]` and `N[ExpIntegralE[2, 3/2], 30]` stay unevaluated, and `Gamma[1 - n, z]` gives those values.
- A machine-real underflow returns `0.`; for example `ExpIntegralE[2, 800.]` gives `0.0`, because the true value is near 10⁻³⁵¹.
- Numeric surfaces: one `double` kernel (real and complex) is shared by the scalar, packed/`NDArray` and `Compile[]` paths.
- Attributes: `Listable`, `NumericFunction`, `Protected`.
- See also: `ExpIntegralEi`, `Gamma`, `LogIntegral`.
