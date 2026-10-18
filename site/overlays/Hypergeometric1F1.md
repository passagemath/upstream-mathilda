### Worked examples

```mathematica
In[1]:= Hypergeometric1F1[1, 2, z]  (* the elementary case (E^z - 1)/z *)
```

```mathematica
In[1]:= Hypergeometric1F1[a, a, z]  (* equal parameters cancel *)
```

```mathematica
In[1]:= Hypergeometric1F1[-2, 1, z]  (* terminates to a Laguerre polynomial *)
```

```mathematica
In[1]:= N[Hypergeometric1F1[1, 2, 3], 40]  (* 40 digits from the series *)
In[2]:= N[(E^3 - 1)/3, 40]  (* the closed form agrees *)
```

```mathematica
In[1]:= Hypergeometric1F1[1., 2., -40.]  (* an alternating series, kept accurate by Kummer's transformation *)
```

```mathematica
In[1]:= Hypergeometric1F1[1.5, 2.5, 100.]  (* large argument, the asymptotic expansion *)
```

```mathematica
In[1]:= Hypergeometric1F1[1., 2., 900.]  (* past the double range: promoted to MPFR *)
```

```mathematica
In[1]:= Hypergeometric1F1[1 + I, 2, 0.5]  (* complex parameter *)
```

```mathematica
In[1]:= N[Hypergeometric1F1[1/2, 3/2, -30], 30]  (* arbitrary precision for negative argument *)
```

### Notes

- `Hypergeometric1F1[a, b, z]` is Kummer's confluent hypergeometric function
  `M(a, b, z)`, equal to `HypergeometricPFQ[{a}, {b}, z]`; it is entire in `z`.
  An unevaluated call comes back as `HypergeometricPFQ[{a}, {b}, z]`.
- Machine-precision calls, real or complex, go to one raw-`double` kernel: closed
  forms, the terminating polynomial, the DLMF 13.7.2 expansion for `|z| > 20`,
  Kummer's transformation for `Re z < 0`, else the series, re-summed in
  double-double when it cancels. The same kernel serves packed arrays, `NDArray`
  arguments and `Compile[]`, so all of these agree bit for bit.
- A result outside the double range is recomputed in MPFR instead of
  overflowing.
- Attributes `Listable`, `NumericFunction`, `Protected`. `Compile[]` lowers real
  and complex scalars; rank-1 array arguments are not lowered yet.
- See also: `HypergeometricPFQ`, `HypergeometricU`, `WhittakerM`, `LaguerreL`,
  `Erf`.
