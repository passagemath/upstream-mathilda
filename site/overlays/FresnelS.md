### Worked examples

```mathematica
In[1]:= {FresnelS[0], FresnelS[Infinity], FresnelS[-I Infinity]}  (* note the sign flip on the imaginary directions *)
In[2]:= FresnelS[1.8]  (* machine real on the shared Fresnel kernel *)
In[3]:= FresnelS[123456.789]  (* the exactly reduced phase keeps a huge argument accurate *)
In[4]:= FresnelS[2.5 + I]  (* machine complex *)
In[5]:= FresnelS[1.*^10 + 1.*^-10 I]  (* the tiny imaginary component is accurate relative to itself *)
In[6]:= N[FresnelS[2], 30]  (* arbitrary precision *)
In[7]:= FresnelS[{0.5, 1., 1.5}] == Map[FresnelS, {0.5, 1., 1.5}]  (* the packed list and the scalar calls agree *)
In[8]:= Series[FresnelS[x], {x, 0, 7}]  (* the Maclaurin series begins with Pi x^3/6 *)
In[9]:= D[FresnelS[x], x]  (* the integrand *)
```

### Notes

- `FresnelS[z] = ∫₀ᶻ Sin[Pi t²/2] dt`, entire and odd; `FresnelS[±Infinity] = ±1/2` but `FresnelS[±I Infinity] = ∓I/2`.
- `FresnelS` and `FresnelC` share every kernel: the pair (C, S) is computed together. Machine arguments run on the `sf_machine_error.c` kernels, shared bit for bit with packed lists, `NDArray` and `Compile[]`; arbitrary precision uses the guarded MPFR series or the real-axis asymptotic expansion.
- Attributes: `Listable`, `NumericFunction`, `Protected`. See also `FresnelC`, `Erf`, `Erfi`.
