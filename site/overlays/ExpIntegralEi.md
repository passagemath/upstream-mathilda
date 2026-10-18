### Worked examples

```mathematica
In[1]:= {ExpIntegralEi[0], ExpIntegralEi[Infinity], ExpIntegralEi[I Infinity]}  (* exact special values *)
In[2]:= ExpIntegralEi[2.]  (* machine real on the Cody-Thacher-style kernel *)
In[3]:= ExpIntegralEi[0.3725]  (* close to the zero x0 = 0.37250741..., still relatively accurate *)
In[4]:= ExpIntegralEi[-1.]  (* on the branch cut the principal value is real *)
In[5]:= {ExpIntegralEi[-1. + 10^-10 I], ExpIntegralEi[-1. - 10^-10 I]}  (* just above and below the cut the imaginary part jumps by 2 Pi *)
In[6]:= ExpIntegralEi[2. + I]  (* machine complex *)
In[7]:= ExpIntegralEi[800.]  (* past the double range the MPFR path returns an extended-exponent number *)
In[8]:= N[ExpIntegralEi[2], 50]  (* arbitrary precision through mpfr_eint *)
In[9]:= N[ExpIntegralEi[I], 30]  (* arbitrary-precision complex *)
In[10]:= D[ExpIntegralEi[x^2], x]  (* the derivative E^z/z with the chain rule *)
```

### Notes

- `ExpIntegralEi[z] = -PV ∫_{-z}^∞ e^(-t)/t dt`, with a branch cut on the negative real axis where the principal (real) value is returned. Exact non-zero arguments stay unevaluated.
- Machine real and complex arguments use the `sf_machine_expint.c` kernels (Cody–Thacher pieces on the real line; Horner series, a backward-evaluated E1 continued fraction and the asymptotic series in the complex plane), shared bit for bit with packed lists, `NDArray` and `Compile[]`.
- Arbitrary precision: `mpfr_eint` for `x > 0`, the guarded convergent series otherwise, switching to the asymptotic expansion for large `|z|`.
- Attributes: `Listable`, `NumericFunction`, `Protected`. See also `ExpIntegralE`, `LogIntegral`, `CosIntegral`, `SinIntegral`, `CoshIntegral`, `SinhIntegral`.
