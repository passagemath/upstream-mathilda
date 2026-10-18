### Worked examples

```mathematica
In[1]:= LerchPhi[1, 2, 1]  (* at z = 1 the Lerch transcendent is a Hurwitz zeta *)
```

```mathematica
In[1]:= LerchPhi[z, -1, 1]  (* a negative integer s gives a rational function *)
```

```mathematica
In[1]:= LerchPhi[-1, 2, 1/2]  (* the alternating half-integer sum is 4 times Catalan's constant *)
```

```mathematica
In[1]:= LerchPhi[1/2, 2, 0]  (* a = 0 drops the singular term and reduces to PolyLog[2, 1/2] *)
```

```mathematica
In[1]:= LerchPhi[2, 3, 1]  (* a positive integer a shifts onto PolyLog *)
```

```mathematica
In[1]:= LerchPhi[{0.1, 0.5, 0.9}, 2.0, 1.5]  (* machine reals on the double-double series kernel *)
```

```mathematica
In[1]:= LerchPhi[0.5 I, 2.0, 1.5]  (* complex z on the machine kernel *)
```

```mathematica
In[1]:= N[LerchPhi[1/3, 3/2, 2/5], 40]  (* arbitrary precision via arb *)
```

```mathematica
In[1]:= LerchPhi[2.0 + 1.0 I, 1.5, 0.5]  (* |z| > 1 off the cut: the Erdelyi continuation *)
```

```mathematica
In[1]:= LerchPhi[2.0, 2.0, 0.5]  (* on the cut with integer s there is no continuation, so it stays unevaluated *)
```

### Notes

- `LerchPhi[z, s, a] = Σ_{k≥0} zᵏ/(k+a)ˢ` generalizes `Zeta`, `HurwitzZeta`
  and `PolyLog`: `LerchPhi[1, s, a] = Zeta[s, a]` and
  `z LerchPhi[z, s, 1] = PolyLog[s, z]`. For `Re a < 0` it uses the symmetric
  power `((k+a)²)^(−s/2)` and drops any `k + a = 0` term, as Mathematica does.
- Many shapes reduce exactly: `z = 0`, `s = 0`, `z = ±1`, integer `a`, and
  negative integer `s`.
- Options: `DoublyInfinite -> True` gives the two-sided sum.
  `IncludeSingularTerm -> True` keeps the `k + a = 0` term, which makes the
  result `ComplexInfinity` at a non-positive integer `a`.
- Numeric surfaces: for `|z| < 1`, machine-precision real and complex arguments
  run on one `double` kernel, shared by scalar, packed/`NDArray` and scalar
  `Compile[]`. Rank-1 `Compile[]` does not lower. Higher precision with
  `|z| ≤ 1` uses arb's `acb_dirichlet_lerch_phi`. For `|z| > 1`, the Erdélyi
  continuation applies only within its domain.
- Known gap: values on the negative real axis with `z < -1` (e.g.
  `LerchPhi[-3.0, 1.5, 0.5]`) are currently inaccurate, so do not rely on them.
- Attributes: `Listable`, `NumericFunction`, `Protected`.
- See also: `Zeta`, `HurwitzZeta`, `PolyLog`.
