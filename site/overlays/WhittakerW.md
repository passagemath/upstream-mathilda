### Worked examples

```mathematica
In[1]:= WhittakerW[2, 3, 4.5]  (* machine value *)
```

```mathematica
In[1]:= Exp[-2.25] 4.5^3.5 HypergeometricU[1.5, 7, 4.5]  (* the same value from the defining U form *)
```

```mathematica
In[1]:= WhittakerW[0, 1/4, 2.] - Sqrt[2./Pi] BesselK[1/4, 1.]  (* checked against W_{0,m}(z) = √(z/π) K_m(z/2) *)
```

```mathematica
In[1]:= N[WhittakerW[3/5, 1/3, 1/7], 30]  (* arbitrary precision through HypergeometricU / arb *)
```

```mathematica
In[1]:= WhittakerW[1/5, 2 - I, 2.]  (* complex order *)
```

```mathematica
In[1]:= WhittakerW[0.5, 0.25, 1.] - WhittakerW[0.5, -0.25, 1.]  (* symmetry W_{k,m} = W_{k,-m} *)
```

```mathematica
In[1]:= {WhittakerW[0, 0, 0], WhittakerW[0, 1, 0], WhittakerW[k, 1/2, 0]}  (* values at the origin, decided by 1/2 - |Re m| *)
```

```mathematica
In[1]:= D[WhittakerW[k, m, z], z]  (* DLMF 13.14.31 *)
```

```mathematica
In[1]:= Series[WhittakerW[2, 1/2, x], {x, Infinity, 5}]  (* terminating case, via the Laguerre rewrite *)
```

```mathematica
In[1]:= WhittakerW[2, 3, {1.5, 2.5, 3.5}]  (* listable *)
```

### Notes

- `Listable`, `NumericFunction`, `Protected`. The definition is `W_{k,m}(z) = e^(−z/2) z^(m+1/2) U(m−k+1/2, 2m+1, z)`, which is even in m. The machine kernel folds m to |m|.
- Exact values at z = 0 depend on the leading exponent 1/2 − |Re m|: the result is `0`, `ComplexInfinity`, or `1/Gamma[1−k]` when m = ±1/2. Other exact input stays symbolic.
- Machine input runs `sf_machine_whittaker_w`, which is the `HypergeometricU` double kernel with the e^(−z/2) z^(m+1/2) prefactor. Buffers share it. Complex orders and precision above machine go through `HypergeometricU`, which uses arb.
- `D` in z follows DLMF 13.14.31. `Series` covers the terminating case only, at 0 and at `Infinity`, through `U(−n, b, z) = (−1)^n n! LaguerreL[n, b−1, z]`.
- See also: `WhittakerM`, `HypergeometricU`.
