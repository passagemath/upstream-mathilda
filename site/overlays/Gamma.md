### Worked examples

```mathematica
In[1]:= Gamma[5]  (* the shifted factorial, 4! *)
In[2]:= Gamma[7/2]  (* half-integers close to rational multiples of Sqrt[Pi] *)
In[3]:= {Gamma[-1/2], Gamma[0]}  (* the functional equation continues past the poles at the non-positive integers *)
In[4]:= N[Gamma[1/3], 40]  (* arbitrary-precision real through mpfr_gamma *)
In[5]:= Gamma[2.3 + I]  (* machine complex on the shared double kernel *)
In[6]:= cg = Compile[{{z, _Complex}}, Gamma[z]]; cg[2.3 + I] === Gamma[2.3 + I]  (* the compiled call and the scalar builtin agree bit for bit *)
In[7]:= N[Gamma[3 + 4 I], 20]  (* arbitrary-precision complex through Spouge's approximation *)
In[8]:= Gamma[3, x]  (* a positive-integer first argument gives the finite closed form *)
In[9]:= {Gamma[1.5, 7.5], Gamma[1, 1.1, 2.2]}  (* upper and generalized incomplete gamma, numerically *)
In[10]:= Gamma[3/2, I]  (* an exact complex incomplete gamma has no closed form and stays unevaluated *)
```

### Notes

- `Gamma[n] = (n-1)!` for positive integers; half-integers reduce exactly to rational multiples of `Sqrt[Pi]`; non-positive integers are poles (`ComplexInfinity`).
- Machine reals use libm `tgamma`; past `Gamma[171.6]` the value is recomputed in MPFR and returned as a machine-precision number with an extended exponent instead of overflowing.
- Machine complex arguments use the `sf_machine_gamma.c` kernel (Taylor tiles near the real axis, shifted Stirling elsewhere). Scalar, packed list, `NDArray` and `Compile[]` share it bit for bit; arbitrary-precision complex uses Spouge's approximation at the requested precision.
- `Gamma[a, z]` is the upper incomplete gamma and `Gamma[a, z0, z1] = Gamma[a, z0] - Gamma[a, z1]`; inexact arguments go through `mpfr_gamma_inc` (real) or a series/continued-fraction split (complex).
- Attributes: `Listable`, `NumericFunction`, `Protected`. See also `LogGamma`, `PolyGamma`, `Beta`, `Pochhammer`, `Factorial`.
