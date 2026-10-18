### Worked examples

```mathematica
In[1]:= AiryAi[0]  (* closed-form value at the origin *)
In[2]:= AiryAi[2]  (* other exact arguments stay symbolic *)
In[3]:= AiryAi[{-2., 0., 2.}]  (* oscillatory on the left, decaying on the right *)
```

```mathematica
In[1]:= D[AiryAi[z], {z, 2}]  (* the Airy equation y'' = z y *)
In[2]:= Series[AiryAi[x], {x, Infinity, 2}]  (* asymptotic expansion with the essential factor kept symbolic *)
```

```mathematica
In[1]:= AiryAi[2.5 + I]  (* complex argument, machine kernel *)
In[2]:= AiryAi[-50.]  (* far into the oscillatory region *)
In[3]:= AiryAi[120.]  (* below the double range; promoted instead of flushed to zero *)
In[4]:= AiryAi[150. I]  (* overflow on the imaginary axis returns an arbitrary-precision value *)
```

```mathematica
In[1]:= N[AiryAi[2], 50]  (* fifty digits through arb *)
In[2]:= AiryAi[4.8`40]  (* precision follows the input *)
In[3]:= N[AiryAi[-3] AiryBiPrime[-3] - AiryAiPrime[-3] AiryBi[-3], 40]  (* the Wronskian Ai Bi' - Ai' Bi = 1/Pi *)
```

```mathematica
In[1]:= AiryAi[NDArray[{1., 2., 3.}]]  (* NDArray input runs element-wise *)
In[2]:= cf = Compile[{{z, _Complex}}, AiryAi[z]]; cf[2.5 + I]  (* complex Compile matches the scalar value *)
```

### Notes

- `AiryAi[z]` is the solution of `y'' == z y` that decays as `z -> +Infinity`. It is entire, `AiryAi[±Infinity] = 0`, and the origin value is `1/(3^(2/3) Gamma[2/3])`.
- Machine real and complex arguments use the raw-double kernel in `sf_machine_bessel.c`. On `[-10, 10]` it evaluates piecewise Chebyshev tables; beyond that the Poincaré expansions apply, with `zeta` in double-double. For complex `z` it uses `K_{1/3}`. The same helper serves packed lists, `NDArray` and `Compile[]` (scalar real, scalar complex and rank-1). Arbitrary precision goes through arb's `acb_hypgeom_airy`, and an overflow or underflow is promoted to an MPFR value.
- `D[AiryAi[z], z] = AiryAiPrime[z]`, which shares the same evaluators.
- Attributes: `Listable`, `NumericFunction`, `Protected`.
- See also: `AiryBi`, `AiryAiPrime`, `BesselK`.
