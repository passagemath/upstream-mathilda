### Worked examples

```mathematica
In[1]:= EllipticE[0]  (* complete, m = 0: the quarter period again *)
Out[1]= 1/2 Pi
```

```mathematica
In[1]:= EllipticE[1]  (* the integrand is Abs[Cos[t]], whose quarter-period integral is 1 *)
Out[1]= 1
```

```mathematica
In[1]:= EllipticE[Pi/2, 1/2]  (* one argument is complete, two incomplete -- and a full quarter period closes the gap *)
Out[1]= EllipticE[1/2]
```

```mathematica
In[1]:= N[EllipticE[1/2], 30]  (* unlike K, E has no closed form at 1/2 or -1 *)
Out[1]= 1.350643881047675502520174735339
```

```mathematica
In[1]:= EllipticE[1/2, 1]  (* E(phi|1) = Integrate[Abs[Cos[t]]], which is Sin[phi] on |phi| <= Pi/2 *)
Out[1]= Sin[1/2]
```

```mathematica
In[1]:= EllipticE[2, 1]  (* past Pi/2 that identity is false, so the call stays symbolic rather than answering Sin[2] *)
Out[1]= EllipticE[2, 1]

In[1]:= N[EllipticE[2, 1], 20]
Out[1]= 1.09070257317431830461

In[1]:= N[2 - Sin[2], 20]
Out[1]= 1.09070257317431830461
```

```mathematica
In[1]:= D[EllipticE[m], m]  (* the parameter derivative of the complete integral *)
Out[1]= (1/2 (EllipticE[m] - EllipticK[m]))/m
```

```mathematica
In[1]:= D[EllipticE[phi, m], phi]  (* the amplitude derivative is the integrand *)
Out[1]= Sqrt[1 - m Sin[phi]^2]
```

```mathematica
In[1]:= Series[EllipticE[m], {m, 0, 3}]  (* same dedicated kernel as K, with the 1/(1-2k) factor *)
Out[1]= 1/2 Pi + -1/8 Pi m + -3/128 Pi m^2 + -5/512 Pi m^3 + O[m]^4
```

```mathematica
In[1]:= EllipticE[Interval[{1/4, 1/2}]]  (* certified decreasing in m, so the endpoints come back swapped *)
Out[1]= Interval[{EllipticE[1/2], EllipticE[1/4]}]
```

```mathematica
In[1]:= EllipticE[2.]  (* past m = 1 a machine argument gets a machine Complex, the limit from below the cut *)
```

```mathematica
In[1]:= N[EllipticE[2], 20]  (* the arbitrary-precision Arb path agrees *)
```

```mathematica
In[1]:= EllipticE[0.3 + 0.4 I]  (* complex parameter at machine precision, by the AGM *)
```

```mathematica
In[1]:= NIntegrate[Sqrt[1 - (0.3 + 0.4 I) Sin[t]^2], {t, 0, Pi/2}]  (* the defining integral by quadrature *)
```

```mathematica
In[1]:= f = Compile[{{z, _Complex}}, EllipticE[z]]; f[0.3 + 0.4 I] == EllipticE[0.3 + 0.4 I]  (* compiled and interpreted calls share one kernel *)
```

```mathematica
In[1]:= EllipticE[{0.1, 0.2, 0.3}]  (* Listable; the unary kernel runs element-wise on the buffer *)
Out[1]= {1.53076, 1.48904, 1.44536}
```

### Notes

`EllipticE` is arity-overloaded exactly as in the Wolfram Language: **one** argument is the
complete integral, **two** the incomplete one. A wrong count emits `EllipticE::argt` and
leaves the call unevaluated.

The `m = 1` pair of examples is the interesting one. `E(φ|1) = ∫₀^φ |cos t| dt` equals
`Sin[φ]` only for `|φ| ≤ π/2`; applied unconditionally it made `EllipticE[2, 1]` answer
`0.909297` where the value is `2 − Sin[2] = 1.090703`, and jump 0.18 against its own
neighbour at `m = 1 − 10⁻¹⁸`. The reduction is now gated on the principal strip, so the
exact call stays symbolic — correct but not closed-form — and `N[]` routes it to Arb, where
it agrees with `2 − Sin[2]` to every digit shown.

Both arities have `double` kernels, but different ones. The **complete** `E[m]` uses
`src/special_functions/sf_machine_ellk.c`: a Cephes-form `P(t) − t log(t) Q(t)` fit for
real `m ≤ 1` (about 1 ulp), and the arithmetic-geometric mean for complex `m` or real
`m > 1`, which returns a machine `Complex` on the principal branch (cut `[1, ∞)`, value on
the cut from below). The same kernels serve packed `List`s, real and complex `NDArray`s and
`Compile[]`, bit for bit. The **incomplete** `E[φ, m]` uses Carlson's `R_F` and `R_D`,
computed from one shared duplication loop, on the real principal domain only; a complex
amplitude or parameter is answered by Arb even at machine precision.

- Attributes: `Listable`, `NumericFunction`, `Protected`.
- See also: `EllipticK`, `EllipticF`, `EllipticPi`, `Hypergeometric2F1`
  (`E(m) = (π/2) ₂F₁(−½, ½; 1; m)`).
