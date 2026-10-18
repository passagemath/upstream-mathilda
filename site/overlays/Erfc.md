### Worked examples

```mathematica
In[1]:= {Erfc[0], Erfc[Infinity], Erfc[-Infinity]}  (* exact special values *)
In[2]:= Erfc[10.]  (* computed as e^(-x^2) erfcx(x), with no cancellation against 1 *)
In[3]:= 1 - Erf[10.]  (* the naive complement loses every digit *)
In[4]:= Erfc[27.]  (* the tail is rounded once into the subnormal range *)
In[5]:= N[Erfc[3/2], 50]  (* arbitrary-precision real through mpfr_erfc *)
In[6]:= Erfc[1.5 - I]  (* machine complex on the shared kernel *)
In[7]:= N[Erfc[1 + I], 25]  (* arbitrary-precision complex *)
In[8]:= Series[Erfc[x], {x, 0, 5}]  (* the Maclaurin series of 1 - erf *)
In[9]:= Erfc[-x]  (* no symmetry fold: erfc(-x) = 2 - erfc(x) is not applied *)
```

### Notes

- `Erfc[z] = 1 - Erf[z]` is entire; unlike `Erf` it is not odd, and `Erfc[-x]` is left as written.
- Machine arguments use the `sf_machine_error.c` kernels shared bit for bit with packed lists, `NDArray` and `Compile[]`: `e^(-x²) erfcx(x)` on the real line (accurate far into the tail), and `e^(-z²) w(iz)` or Salzer's series in the complex plane.
- Arbitrary precision uses `mpfr_erfc` (real) and the guarded complex erf series with the complement formed at working precision.
- Attributes: `Listable`, `NumericFunction`, `Protected`. See also `Erf`, `Erfi`, `InverseErfc`.
