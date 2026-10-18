### Worked examples

```mathematica
In[1]:= {Erf[0], Erf[Infinity], Erf[-z]}  (* exact special values and the odd-symmetry fold *)
In[2]:= Erf[0.95]  (* machine real on the shared erfcx-table kernel *)
In[3]:= Erf[1.*^-300 + 2. I]  (* next to the imaginary axis the tiny real part is still accurate to every digit *)
In[4]:= N[Erf[1], 40]  (* arbitrary-precision real through mpfr_erf *)
In[5]:= N[Erf[1 + I], 20]  (* arbitrary-precision complex through the guarded Maclaurin series *)
In[6]:= Erf[1.5, 2]  (* the generalized form erf(z1) - erf(z0) *)
In[7]:= Erf[{0.5, 1.5, 2.5}] == Table[Erf[x], {x, {0.5, 1.5, 2.5}}]  (* the packed list runs the same kernel as the scalar call *)
In[8]:= Series[Erf[x], {x, 0, 7}]  (* the odd Maclaurin series *)
In[9]:= D[Erf[x^2], x]  (* the chain rule on the Gaussian derivative *)
In[10]:= Erf[2]  (* an exact non-special argument stays unevaluated *)
```

### Notes

- `Erf[z] = (2/Sqrt[Pi]) ∫₀ᶻ e^(-t²) dt` is entire and odd; `Erf[z0, z1] = erf(z1) - erf(z0)` reduces only when both values are concrete.
- Machine real and complex arguments use the `sf_machine_error.c` kernels (erfcx tables; Maclaurin, axis Taylor series, a Faddeeva-function table and Salzer's series in the complex plane), shared bit for bit by the scalar call, packed lists, `NDArray` and `Compile[]`.
- Arbitrary precision uses `mpfr_erf` (real) and a guarded DLMF 7.6.2 series (complex).
- Known gap: a machine complex argument with `|Im z| > 26` whose value leaves the double range answers `Overflow[]` (e.g. `Erf[3. + 30. I]`), where `Erfc` and `Erfi` return extended-exponent numbers.
- Attributes: `Listable`, `NumericFunction`, `Protected`. See also `Erfc`, `Erfi`, `InverseErf`, `FresnelC`, `FresnelS`.
