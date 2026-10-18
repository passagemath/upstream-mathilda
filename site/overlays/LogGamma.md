### Worked examples

```mathematica
In[1]:= LogGamma[5]  (* exact at positive integers: Log[4!] *)
In[2]:= LogGamma[-3/2]  (* negative half-integers carry the branch term -Ceiling[-z] Pi I *)
In[3]:= LogGamma[-2.5]  (* a negative machine real gives the complex value from above the cut *)
In[4]:= LogGamma[1000.]  (* finite where Gamma[1000.] would leave the double range *)
In[5]:= N[LogGamma[100], 40]  (* arbitrary precision *)
In[6]:= {LogGamma[2. + 30. I], Log[Gamma[2. + 30. I]]}  (* the continuous branch: the imaginary part keeps growing where Log[Gamma] wraps *)
In[7]:= LogGamma[-4.5 + 3. I] - Log[Gamma[-4.5 + 3. I]]  (* the two differ by a multiple of 2 Pi I *)
In[8]:= N[LogGamma[1 + I], 30]  (* arbitrary-precision complex through Stirling's series *)
In[9]:= D[LogGamma[z], z]  (* the derivative is the digamma function *)
```

### Notes

- `LogGamma[z]` is the analytic continuation of `log Γ(z)` with a single branch cut on the negative real axis — not `Log[Gamma[z]]`, whose imaginary part is confined to `(-Pi, Pi]`.
- Exact integers give `Log[(n-1)!]`; half-integers give the log of the exact `Sqrt[Pi]` form (plus the branch term below zero); non-positive integers give `Infinity`.
- Machine complex arguments use the `sf_machine_gamma.c` kernel (Hare's winding-counted continuous branch), shared bit for bit with packed lists, `NDArray` and `Compile[]`; machine reals use libm `lgamma`.
- A machine pole (`LogGamma[0.]`) currently answers `ComplexInfinity`, while the exact pole `LogGamma[0]` gives `Infinity`.
- Attributes: `Listable`, `NumericFunction`, `Protected`. See also `Gamma`, `PolyGamma` (`PolyGamma[-1, z]` is `LogGamma[z]`).
