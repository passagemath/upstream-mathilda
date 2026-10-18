### Worked examples

```mathematica
In[1]:= EllipticK[0]  (* at m = 0 the integrand is 1, so K is the quarter period Pi/2 *)
Out[1]= 1/2 Pi
```

```mathematica
In[1]:= EllipticK[1]  (* the integrand loses its last denominator at t = Pi/2 *)
Out[1]= ComplexInfinity
```

```mathematica
In[1]:= EllipticK[-1]  (* a lemniscatic singular value: closed form in Gamma, not a decimal *)
Out[1]= (1/4 Gamma[1/4]^2)/Sqrt[2 Pi]
```

```mathematica
In[1]:= EllipticK[1/2]  (* the second singular value; EllipticE has no closed form at either *)
Out[1]= (8 Pi^(3/2))/Gamma[-1/4]^2
```

```mathematica
In[1]:= N[EllipticK[1/2], 30]  (* the same number to 30 digits, through Arb *)
Out[1]= 1.854074677301371918433850347195
```

```mathematica
In[1]:= Series[EllipticK[m], {m, 0, 3}]  (* a dedicated kernel: the generic Taylor path hits 0/0 here *)
Out[1]= 1/2 Pi + 1/8 Pi m + 9/128 Pi m^2 + 25/512 Pi m^3 + O[m]^4
```

```mathematica
In[1]:= N[EllipticK[9/5], 20]  (* past the branch point at m = 1 the value is complex *)
Out[1]= 1.41933775128651291929 - 1.348846512193268578*I
```

```mathematica
In[1]:= EllipticK[2.]  (* a machine m past the cut gets a machine Complex: the limit from below, from the complex AGM kernel *)
```

```mathematica
In[1]:= N[(EllipticK[1/2] - I EllipticK[1/2])/Sqrt[2], 20]  (* cross-check by the reciprocal-parameter identity K(m) = (K(1/m) - I K(1 - 1/m))/Sqrt[m] at m = 2 *)
```

```mathematica
In[1]:= EllipticK[0.3 + 0.4 I]  (* complex parameter at machine precision: AGM in raw double re/im *)
```

```mathematica
In[1]:= NIntegrate[1/Sqrt[1 - (0.3 + 0.4 I) Sin[t]^2], {t, 0, Pi/2}]  (* the defining integral, by quadrature, agrees *)
```

```mathematica
In[1]:= EllipticK[NDArray[{0.5 + 0.5 I, 2. + 0. I, -3. + I}]]  (* a complex NDArray runs the same kernel element-wise *)
```

```mathematica
In[1]:= CompileDiagnostics[{{z, _Complex}}, EllipticK[z]]  (* and it lowers in Compile[] at a complex scalar *)
```

```mathematica
In[1]:= N[EllipticK[1 - 10^-17], 20]  (* here dK/dm ~ 2^56, and all 20 digits are still right *)
Out[1]= 20.9582676515692789829
```

```mathematica
In[1]:= D[EllipticK[m], m]  (* the parameter derivative, in terms of E and K themselves *)
Out[1]= (1/2 (EllipticE[m] - EllipticK[m] (1 - m)))/(m (1 - m))
```

```mathematica
In[1]:= EllipticK[Interval[{1/4, 1/2}]]  (* certified: K is increasing in m below 1 *)
Out[1]= Interval[{EllipticK[1/4], (8 Pi^(3/2))/Gamma[-1/4]^2}]
```

```mathematica
In[1]:= Precision[EllipticK[0.5]]  (* Real in, Real out -- so the result packs *)
Out[1]= MachinePrecision
```

### Notes

`EllipticK[m]` takes the **parameter** `m = k²`, not the modulus `k` — the Wolfram
convention, and the one place a reader of these signatures goes silently wrong.

Three argument classes are answered by three different routes. An exact `0`, `1`, `∞`,
`1/2` or `−1` reduces symbolically; any other exact argument stays symbolic, since a
decimal would be a loss of information rather than an evaluation. A machine argument
is answered in `double` (`src/special_functions/sf_machine_ellk.c`): a real `m < 1` by a
Cephes-form `P(t) − log(t) Q(t)` fit, about 1 ulp, and a complex `m` or a real `m > 1`
by the arithmetic-geometric mean, which returns a machine `Complex` on the principal
branch (cut on `[1, ∞)`, value on the cut taken from below, as in Mathematica). Any
request above machine precision goes to Arb's `acb_elliptic_k`, rigorous to the digits
asked for.

- Attributes: `Listable`, `NumericFunction`, `Protected`. The same machine kernels serve
  the scalar call, packed `List`s, real and complex `NDArray`s and `Compile[]` (`_Real`,
  `_Complex`, rank-1), bit for bit.
- `EllipticK[1.]` and `EllipticK[SetPrecision[1, 30]]` are `ComplexInfinity`, like the
  exact `1`: the pole is tested before any numeric path.
- See also: `EllipticE`, `EllipticF`, `EllipticPi`, `Hypergeometric2F1`
  (`K(m) = (π/2) ₂F₁(½, ½; 1; m)`).

The conditioning example is worth reading twice: `dK/dm ≈ 1/(2(1−m))`, which at
`m = 1 − 10⁻¹⁷` is about `2⁵⁶`, so an arbitrary-precision request has to be *planned* for
the loss rather than merely issued at the precision asked for. How much a function loses
is set by its conditioning, not by how its argument is spelled.
