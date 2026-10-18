### Worked examples

```mathematica
In[1]:= SinIntegral[{-Infinity, Infinity, -I Infinity, I Infinity}]  (* limits along the axes *)
In[2]:= SinIntegral[2.8]  (* machine real *)
In[3]:= SinIntegral[1.*^6]  (* approaches Pi/2 for large x *)
In[4]:= SinIntegral[3. I]  (* on the imaginary axis, I SinhIntegral[3] *)
In[5]:= SinIntegral[2.5 + I]  (* machine complex *)
In[6]:= N[SinIntegral[2], 50]  (* arbitrary precision *)
In[7]:= SinIntegral[-x]  (* the odd-symmetry fold *)
In[8]:= Series[SinIntegral[x], {x, 0, 7}]  (* the Taylor series at the origin *)
In[9]:= Normal[Series[SinIntegral[x], {x, Infinity, 3}]]  (* the asymptotic expansion *)
In[10]:= D[SinIntegral[x], x]  (* the derivative is Sinc *)
```

### Notes

- `SinIntegral[z] = ∫₀ᶻ Sin[t]/t dt` is entire and odd, with `SinIntegral[±Infinity] = ±Pi/2`.
- Machine arguments use the `sf_machine_expint.c` kernels (Cephes `sici` structure on the real line; series, backward E1 continued fraction, or asymptotic series in the complex plane), shared bit for bit with packed lists, `NDArray` and `Compile[]`; arbitrary precision uses the guarded Maclaurin series or the asymptotic expansion.
- Attributes: `Listable`, `NumericFunction`, `Protected`. See also `CosIntegral`, `SinhIntegral`, `Sinc`, `ExpIntegralEi`.
