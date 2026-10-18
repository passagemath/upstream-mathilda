### Worked examples

```mathematica
In[1]:= CosIntegral[{-Infinity, Infinity, -I Infinity, I Infinity}]  (* the limits along the four axis directions *)
In[2]:= CosIntegral[2.8]  (* machine real *)
In[3]:= CosIntegral[3.5]  (* inside [2, 4] the real kernel uses the factored form (x - x0) Q(x) through the first zero x0 = 3.3842... *)
In[4]:= CosIntegral[-2.]  (* a negative real lies on the cut: the value from above, Ci(2) + I Pi *)
In[5]:= CosIntegral[3. I]  (* on the imaginary axis, CoshIntegral[3] + I Pi/2 *)
In[6]:= CosIntegral[2.5 + I]  (* machine complex *)
In[7]:= N[CosIntegral[2], 50]  (* arbitrary precision *)
In[8]:= Series[CosIntegral[x], {x, 0, 6}]  (* the logarithmic Taylor series at the origin *)
In[9]:= Normal[Series[CosIntegral[x], {x, Infinity, 3}]]  (* the trig-prefactored asymptotic expansion *)
In[10]:= CosIntegral[-x]  (* no symmetry fold *)
```

### Notes

- `CosIntegral[z] = -∫_z^∞ Cos[t]/t dt` has a logarithmic singularity at 0 and a branch cut on `(-∞, 0]`; a negative real gives the from-above value `Ci(|x|) + I Pi`.
- Machine arguments use the `sf_machine_expint.c` kernels (Cephes `sici` structure on the real line with a factored form through the first zero; series or E1 continued fraction in the complex plane), shared bit for bit with packed lists, `NDArray` and `Compile[]`; arbitrary precision uses the guarded series or the Stokes-corrected asymptotic expansion.
- `D[CosIntegral[z], z] = Cos[z]/z`. Attributes: `Listable`, `NumericFunction`, `Protected`. See also `SinIntegral`, `CoshIntegral`, `ExpIntegralEi`.
