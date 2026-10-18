### Worked examples

```mathematica
In[1]:= {FresnelC[0], FresnelC[Infinity], FresnelC[I Infinity]}  (* odd and entire, with limits 1/2 and I/2 *)
In[2]:= FresnelC[1.8]  (* machine real: auxiliary f, g tables with an exactly reduced phase *)
In[3]:= FresnelC[123456.789]  (* a huge argument is still accurate because the phase Pi x^2/2 is reduced exactly *)
In[4]:= FresnelC[2.5 + I]  (* machine complex through the exactly rotated erf kernel *)
In[5]:= FresnelC[1.*^10 + 1.*^-10 I]  (* the scaled axis series handles a huge argument just off the real axis *)
In[6]:= N[FresnelC[2], 50]  (* arbitrary precision through the guarded MPFR series *)
In[7]:= FresnelC[{0.5, 1.0, 1.5}]  (* threads element-wise over a packed real list *)
In[8]:= D[FresnelC[x^2], x]  (* the integrand, with the chain rule *)
In[9]:= FresnelC[-x]  (* the odd-symmetry fold *)
```

### Notes

- `FresnelC[z] = ∫₀ᶻ Cos[Pi t²/2] dt` (Mathematica's π/2 convention), entire and odd; `FresnelC[±Infinity] = ±1/2`, `FresnelC[±I Infinity] = ±I/2`.
- Machine real and complex arguments use the `sf_machine_error.c` kernels — Maclaurin series, auxiliary-function tables with exact phase reduction, the axis Taylor series, and the complex erf core fed exactly rotated arguments — shared bit for bit with packed lists, `NDArray` and `Compile[]`.
- Arbitrary precision uses the convergent series with `(π/2)|z|²/ln2` guard bits, or (real arguments only) the DLMF 7.12 asymptotic expansion.
- Attributes: `Listable`, `NumericFunction`, `Protected`. See also `FresnelS`, `Erf`, `Erfi`.
