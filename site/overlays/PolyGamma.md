### Worked examples

```mathematica
In[1]:= PolyGamma[5]  (* the digamma at a positive integer: a harmonic number minus EulerGamma *)
In[2]:= {PolyGamma[1, 1], PolyGamma[3, 1]}  (* odd orders at integers close through zeta at even arguments *)
In[3]:= PolyGamma[0, 1/2]  (* the Gauss digamma closed form at a half-integer *)
In[4]:= PolyGamma[2.5]  (* a machine argument is computed directly by the shared kernel *)
In[5]:= PolyGamma[1.461632144968362]  (* relatively accurate right at the positive zero of the digamma *)
In[6]:= PolyGamma[3, 2.5 + I]  (* machine complex, order 3 *)
In[7]:= N[PolyGamma[2, 1 + I], 30]  (* arbitrary-precision complex through the Bernoulli asymptotic *)
In[8]:= PolyGamma[-1, z]  (* order -1 is the log-gamma function *)
In[9]:= {PolyGamma[2, 1], PolyGamma[0, 1/4]}  (* even orders at integers, and non-half-integer rationals, stay symbolic *)
```

### Notes

- `PolyGamma[z]` is the digamma ψ(z) (`PolyGamma[0, z]`); `PolyGamma[n, z]` is ψ⁽ⁿ⁾(z). Non-positive integers are poles (`ComplexInfinity`); `PolyGamma[-1, z]` returns `LogGamma[z]`.
- Exact closed forms: ψ(m) = `H[m-1] - EulerGamma`, odd-order ψ⁽ⁿ⁾(m) as rational plus rational·`Pi^(n+1)`, and the Gauss digamma theorem at positive half-integers. `ζ(odd)` values (even order at integers) are left symbolic.
- Machine real or complex `z` with integer order `0 ≤ n ≤ 40` runs on the `sf_machine_gamma.c` kernels, shared bit for bit with packed lists, `NDArray` and `Compile[]` (both `PolyGamma[x]` and `PolyGamma[n, x]` compile); otherwise MPFR (`mpfr_digamma`, or the recurrence-shifted Bernoulli asymptotic).
- Complex or fractional orders stay unevaluated. Attributes: `Listable`, `NumericFunction`, `Protected`. See also `Gamma`, `LogGamma`, `HarmonicNumber`, `Zeta`.
