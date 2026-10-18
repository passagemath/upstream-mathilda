### Worked examples

```mathematica
In[1]:= Hypergeometric2F1[1, 1, 2, z]  (* the logarithm *)
```

```mathematica
In[1]:= Hypergeometric2F1[-3, 1, 1, z]  (* terminates to (1 - z)^3 *)
```

```mathematica
In[1]:= Hypergeometric2F1[1, 1, 3/2, z]  (* the central-binomial ArcSin family *)
```

```mathematica
In[1]:= N[Hypergeometric2F1[1/2, 1/2, 3/2, 1/4]/2, 40]  (* since z 2F1(1/2,1/2;3/2;z^2) = ArcSin[z], here at z = 1/2 *)
In[2]:= N[ArcSin[1/2], 40]  (* the same value from ArcSin *)
```

```mathematica
In[1]:= Hypergeometric2F1[0.3, 0.7, 1.9, -50.]  (* far outside the unit disc, by a connection formula *)
```

```mathematica
In[1]:= Hypergeometric2F1[0.5, 0.5, 1.5, 3.]  (* on the cut z > 1: the limit from below, ArcSin[Sqrt[3]]/Sqrt[3] *)
```

```mathematica
In[1]:= Hypergeometric2F1[0.3, 0.7, 1.9, 2. + 1. I]  (* complex argument beyond |z| = 1 *)
```

```mathematica
In[1]:= N[Hypergeometric2F1[1/3, 2/3, 5/6, 99/100], 50]  (* arbitrary precision near the unit circle, via arb *)
```

```mathematica
In[1]:= Hypergeometric2F1[1, 1, 2, {0.5, -0.5}]  (* Listable *)
```

```mathematica
In[1]:= Hypergeometric2F1[1/2, 1/3, 2, 1]  (* exact Gauss sum at z = 1 is not applied: stays unevaluated *)
```

### Notes

- `Hypergeometric2F1[a, b, c, z]` is the Gauss hypergeometric function, equal to
  `HypergeometricPFQ[{a, b}, {c}, z]`; an unevaluated call comes back in the
  `HypergeometricPFQ` form.
- At machine precision it is evaluated on the whole complex plane (real
  parameters): direct and Pfaff series, the `1-z` and `1/(1-z)` connection
  formulas, and continuation of the hypergeometric ODE where all of these are
  slow. On the cut `z > 1` a real argument gives the value from below, the same
  convention as Mathematica. With complex parameters only `|z| <= 0.9` is
  evaluated.
- Arbitrary precision is limited to `|z| < 1`, where arb's `acb_hypgeom_2f1`
  is used. An exact or high-precision `z` outside the disc stays unevaluated, and
  so does the exact value at `z = 1`. A machine `z = 1.` uses Gauss's sum.
- Attributes `Listable`, `NumericFunction`, `Protected`. Packed arrays and
  `NDArray` arguments map the machine kernel over the buffer (the coefficient
  table is built once per buffer); a buffer with an element on the cut comes back
  complex. `Compile[]` lowers real and complex scalars; rank-1 array arguments are
  not lowered yet.
- See also: `HypergeometricPFQ`, `Hypergeometric1F1`, `LegendreP`, `EllipticK`.
