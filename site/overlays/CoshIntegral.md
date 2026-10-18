### Worked examples

```mathematica
In[1]:= {CoshIntegral[0], CoshIntegral[Infinity], CoshIntegral[I Infinity]}  (* the logarithmic singularity and the limits *)
In[2]:= CoshIntegral[2.8]  (* machine real *)
In[3]:= CoshIntegral[0.5238]  (* beside the zero x0 = 0.52382..., still relatively accurate *)
In[4]:= CoshIntegral[-2.]  (* a negative real is on the cut: the value from above, Chi(2) + I Pi *)
In[5]:= CoshIntegral[3. I]  (* on the imaginary axis, CosIntegral[3] + I Pi/2 *)
In[6]:= CoshIntegral[3.5 + I]  (* machine complex *)
In[7]:= N[CoshIntegral[2], 50]  (* arbitrary precision *)
In[8]:= CoshIntegral[{1.0, 2.0, 3.0}]  (* threads element-wise over a packed real list *)
In[9]:= Series[CoshIntegral[x], {x, 0, 4}]  (* the logarithmic Taylor series *)
In[10]:= D[CoshIntegral[x], x]  (* the derivative is Cosh[x]/x *)
```

### Notes

- `CoshIntegral[z] = EulerGamma + Log[z] + ∫₀ᶻ (Cosh[t] - 1)/t dt` has a logarithmic singularity at 0 and a branch cut on `(-∞, 0]`; it is the imaginary-axis sibling of `CosIntegral` (`Chi(z) = Ci(I z) - I Pi/2`).
- Machine arguments use the `sf_machine_expint.c` kernels (Cephes `shichi` on the real line, `(Ei - E1)/2` beyond 4; `Ci(iz)` with an exactly placed constant in the complex plane), shared bit for bit with packed lists, `NDArray` and `Compile[]`; overflow and arbitrary precision go to MPFR.
- Attributes: `Listable`, `NumericFunction`, `Protected`. See also `SinhIntegral`, `CosIntegral`, `ExpIntegralEi`.
