### Worked examples

```mathematica
In[1]:= StruveH[1/2, x]  (* no closed form is substituted for exact symbolic input *)
```

```mathematica
In[1]:= StruveH[0, 1.5]  (* machine value from the double series kernel *)
```

```mathematica
In[1]:= NIntegrate[2/Pi Sin[1.5 Cos[t]], {t, 0, Pi/2}]  (* the Poisson-type integral (2/π)∫ sin(z cos θ) dθ gives the same number *)
```

```mathematica
In[1]:= StruveH[1/2, 2.] - Sqrt[2/(Pi 2.)] (1 - Cos[2.])  (* agreement with the elementary half-integer form, to rounding *)
```

```mathematica
In[1]:= N[StruveH[1, 13/10], 40]  (* arbitrary precision through the regularized 1F2 in arb *)
```

```mathematica
In[1]:= N[StruveH[-3/2, 2], 20] - N[-BesselJ[3/2, 2], 20]  (* the arb path at ν = -3/2, checked against H_{-3/2} = -J_{3/2} *)
```

```mathematica
In[1]:= StruveH[0, 30.]  (* large argument: the Y_ν + K_ν asymptotic regime *)
```

```mathematica
In[1]:= NIntegrate[2/Pi Sin[30 Cos[t]], {t, 0, Pi/2}]  (* independent check of the previous line *)
```

```mathematica
In[1]:= StruveH[1, {0.5, 1., 2.}]  (* listable; a packed list runs the same kernel *)
```

```mathematica
In[1]:= N[StruveH[0, 1 + I], 20]  (* complex z needs N: the machine kernel is real-only *)
```

### Notes

- `Listable`, `NumericFunction`, `Protected`. Exact input is never rewritten: `StruveH[1/2, x]` stays as it is.
- The only exact value is `StruveH[nu, 0] = 0` for a real ν > −1.
- Machine-precision real input runs `sf_machine_struve_h`: the convergent 1F2 series, or for large z the `Y_ν + K_ν` asymptotic, whichever predicts the smaller error. The same kernel backs packed lists, `NDArray` and `Compile[]` (scalar and rank-1), so all of them agree bit for bit. Between about z = 17 and z = 24 (small ν) the double answer is good to about 1e-9 absolute; `N[..., p]` gives full precision there.
- With precision above machine (`N[StruveH[ν, z], p]`, or an MPFR argument) the value comes from FLINT/arb as `(z/2)^(ν+1) · 1F2~(1; 3/2, ν+3/2; −z²/4)`, using the regularized 1F2, so it holds for any complex ν and z.
- A machine complex argument (`StruveH[0, 1. + I]`) stays unevaluated. Use `N` on exact input to get a complex value.
- No derivative rule exists yet: `D[StruveH[n, x], x]` returns an inert `Derivative[0, 1][StruveH][n, x]`, so `Series` at a regular point yields only `Derivative` coefficients.
- See also: `BesselJ`, `BesselY`, `HypergeometricPFQ`.
