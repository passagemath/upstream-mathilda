### Worked examples

```mathematica
In[1]:= {InverseErf[0], InverseErf[1], InverseErf[-1]}  (* exact values at the centre and the endpoints *)
In[2]:= InverseErf[0.6]  (* machine real on the piecewise-polynomial kernel *)
In[3]:= InverseErf[1/{2., 3., 4., 5.}]  (* listable, element by element *)
In[4]:= {InverseErf[0.999999], Erf[InverseErf[0.999999]]}  (* deep in the tail, and the round trip back through Erf *)
In[5]:= N[InverseErf[33/100], 50]  (* arbitrary precision by Newton's method with precision doubling *)
In[6]:= N[Sqrt[2] InverseErf[2 (95/100) - 1], 30]  (* the 95% quantile of the standard normal distribution *)
In[7]:= InverseErf[0.4, 0.2]  (* the generalized form: the z with Erf[0.4, z] = 0.2 *)
In[8]:= Series[InverseErf[x], {x, 0, 7}]  (* the Taylor series through the derivative rule *)
In[9]:= {InverseErf[2], InverseErf[0.5 + 0.1 I]}  (* outside [-1, 1], and complex arguments, stay unevaluated *)
```

### Notes

- `InverseErf[s]` solves `s = Erf[z]`; as in Mathematica it gives explicit values only for real `s` in `[-1, 1]`, and `InverseErf[z0, s]` inverts the generalized `Erf[z0, z]`.
- Machine reals use `sf_inverse_erf_r` (polynomials in `s` and in Giles' variable `w = -log(1 - s²)`, one Halley step in the far tail), shared bit for bit with packed lists, `NDArray` and `Compile[]`; arbitrary precision uses Newton's method on `mpfr_erf`.
- `D[InverseErf[z], z] = (Sqrt[Pi]/2) E^(InverseErf[z]^2)`. Exact rationals such as `InverseErf[1/2]` stay unevaluated.
- Attributes: `Listable`, `NumericFunction`, `Protected`. See also `InverseErfc`, `Erf`, `Quantile`.
