### Worked examples

```mathematica
In[1]:= Beta[3, 5]  (* a positive-integer argument collapses to an exact rational *)
In[2]:= {Beta[1/2, 1/2], Beta[1/3, 1/3]}  (* half-integers fold to Pi; other rationals keep a Gamma quotient *)
In[3]:= {Beta[-2, 1], Beta[0, b]}  (* a cancelling pair of gamma poles has a finite limit; a surviving pole is ComplexInfinity *)
In[4]:= Beta[2.3, 3.2]  (* machine real, computed directly in double *)
In[5]:= Beta[2.5 + I, 1 - I]  (* machine complex: the three Gammas under one exponential *)
In[6]:= N[Beta[22/10, 33/10], 50]  (* arbitrary precision through the Gamma ratio *)
In[7]:= Beta[1/2, 3, 4]  (* the incomplete beta terminates when b is a positive integer *)
In[8]:= Beta[0.3, 0.5, 2.5]  (* a numeric incomplete beta through Hypergeometric2F1 *)
In[9]:= cb = Compile[{{v, _Real, 1}}, Beta[v, 2.5]]; cb[{1.5, 2.5}] == Beta[{1.5, 2.5}, 2.5]  (* rank-1 Compile[] runs the same kernel as the packed map *)
In[10]:= Beta[a, b]  (* symbolic arguments stay unevaluated, as in Mathematica *)
```

### Notes

- `Beta[a, b] = Gamma[a] Gamma[b]/Gamma[a+b]`; `Beta[z, a, b]` is the incomplete beta `∫₀ᶻ t^(a-1) (1-t)^(b-1) dt` (the integration limit comes first), and `Beta[z0, z1, a, b]` the generalized form `Beta[z1, a, b] - Beta[z0, a, b]`.
- Machine-precision `Beta[a, b]` (real or complex) is computed in double by `sfg_beta`/`sfg_beta_c`, correcting the rounding of `a + b` through `PolyGamma`; poles and `|a|, |b|, |a+b| > 170` fall back to the exact/MPFR path.
- Packed and `NDArray` arguments run element-wise — including both arguments arrays — and `Compile[]` lowers scalar and array-with-scalar calls.
- Attributes: `Listable`, `NumericFunction`, `Protected`. See also `Gamma`, `Pochhammer`, `Hypergeometric2F1`.
