### Worked examples

```mathematica
In[1]:= {SinhIntegral[0], SinhIntegral[Infinity], SinhIntegral[I Infinity]}  (* exact special values *)
In[2]:= SinhIntegral[-2 x]  (* odd in its argument *)
In[3]:= SinhIntegral[2.8]  (* machine real *)
In[4]:= SinhIntegral[3. I]  (* on the imaginary axis, I SinIntegral[3] *)
In[5]:= SinhIntegral[2.5 + I]  (* machine complex: an exact rotation of the Si kernel *)
In[6]:= SinhIntegral[10.^6]  (* past the double range the result is an extended-exponent number, not Infinity *)
In[7]:= N[SinhIntegral[2], 50]  (* arbitrary precision *)
In[8]:= Series[SinhIntegral[x], {x, 0, 7}]  (* all-positive Taylor coefficients *)
In[9]:= Integrate[Sinh[a x]/x, x]  (* the integral that defines it *)
```

### Notes

- `SinhIntegral[z] = ∫₀ᶻ Sinh[t]/t dt` is entire and odd, the imaginary-axis sibling of `SinIntegral` (`Shi(z) = -I Si(I z)`); `SinhIntegral[±I Infinity] = ±I Pi/2`.
- Machine arguments use the `sf_machine_expint.c` kernels (Cephes `shichi` on the real line, `(Ei + E1)/2` beyond 4), shared bit for bit with packed lists, `NDArray` and `Compile[]`; overflow past `x ≈ 720` and arbitrary precision go to MPFR.
- Attributes: `Listable`, `NumericFunction`, `Protected`. See also `CoshIntegral`, `SinIntegral`, `ExpIntegralEi`.
