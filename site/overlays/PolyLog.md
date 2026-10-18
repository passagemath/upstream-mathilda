### Worked examples

```mathematica
In[1]:= PolyLog[2, 1]  (* at z = 1 the polylogarithm is a zeta value *)
```

```mathematica
In[1]:= PolyLog[3, 1/2]  (* Euler's closed form mixing Log[2], Pi^2 and Zeta[3] *)
```

```mathematica
In[1]:= PolyLog[2, I]  (* the dilogarithm at I: Catalan's constant appears *)
```

```mathematica
In[1]:= PolyLog[-2, z]  (* negative integer orders are rational functions *)
```

```mathematica
In[1]:= N[PolyLog[3, 1/2 + I/2], 30]  (* a complex argument at 30 digits *)
```

```mathematica
In[1]:= N[PolyLog[1/3, 9/10], 40]  (* a fractional order at high precision via arb *)
```

```mathematica
In[1]:= PolyLog[2, 3.0]  (* on the cut the value is taken from below: Im = -Pi Log[3] *)
```

```mathematica
In[1]:= PolyLog[2, -1000.0]  (* real for z < -1, through the 1/z inversion *)
```

```mathematica
In[1]:= PolyLog[2, {0.25, 0.5, 0.75}]  (* a machine list runs through the packed dilogarithm kernel *)
```

```mathematica
In[1]:= PolyLog[2.5, 1000.0]  (* a non-integer order far from the origin is outside every regime, so it stays unevaluated *)
```

### Notes

- `PolyLog[n, z] = Li_n(z) = Σ_{k≥1} zᵏ/kⁿ`, with a branch cut along
  `[1, ∞)`. Values on the cut are taken from below, as in Mathematica.
- Exact closed forms: `PolyLog[n, 0] = 0`, `PolyLog[1, z] = -Log[1 - z]`, and
  `PolyLog[0, z] = z/(1 - z)`. Negative integer orders give rational functions.
  `PolyLog[n, ±1]` reduces to `Zeta[n]`. `Li_2(1/2)`, `Li_3(1/2)` and `Li_2(±I)`
  have closed forms.
- Numeric surfaces: a real order at machine precision runs on one `double`
  kernel (the Bernoulli-series dilogarithm, plus power series, Jonquière and
  inversion for other orders). The scalar, packed/`NDArray` and `Compile[]`
  paths (scalar and array, real and complex `z`) all share it. A complex order
  and arbitrary precision use complex MPFR; real `-1 < z < 1` at high precision
  uses arb's `acb_polylog`.
- `PolyLog[n, p, z]` (Nielsen) is accepted but left unevaluated.
- Attributes: `Listable`, `NumericFunction`, `Protected`.
- See also: `LerchPhi`, `Zeta`, `Log`.
