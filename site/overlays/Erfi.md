### Worked examples

```mathematica
In[1]:= {Erfi[0], Erfi[I Infinity], Erfi[-x]}  (* the imaginary-axis limit is finite, and the function is odd *)
In[2]:= Erfi[2.5]  (* machine real: e^(x^2) times a tabulated Dawson factor *)
In[3]:= Erfi[26.]  (* still a machine double near the top of the range *)
In[4]:= Erfi[30.]  (* past the double range the MPFR series returns an extended-exponent number *)
In[5]:= Erfi[1.5 - I]  (* machine complex through the complex erf kernel *)
In[6]:= N[Erfi[1/2], 50]  (* the all-positive Maclaurin series at arbitrary precision *)
In[7]:= N[Erfi[1], 30]  (* another arbitrary-precision value *)
In[8]:= Series[Erfi[x], {x, 0, 7}]  (* like Erf's series without the alternating signs *)
In[9]:= D[Erfi[z], z]  (* the derivative has the positive exponent *)
```

### Notes

- `Erfi[z] = -I Erf[I z] = (2/Sqrt[Pi]) ∫₀ᶻ e^(t²) dt` is entire and odd; `Erfi[±I Infinity] = ±I`.
- There is no libm or MPFR `erfi`: machine reals use a Dawson-factor table times `e^(x²)` (exact square), machine complex the complex erf kernel, both shared bit for bit with packed lists, `NDArray` and `Compile[]`; arbitrary precision sums the all-positive Maclaurin series (real) or the guarded erf series (complex).
- Attributes: `Listable`, `NumericFunction`, `Protected`. See also `Erf`, `Erfc`, `FresnelC`.
