### Worked examples

```mathematica
In[1]:= AiryBi[0]  (* closed-form value at the origin *)
In[2]:= AiryBi[Infinity]  (* the dominant solution grows without bound *)
In[3]:= AiryBi[{-2., 0., 2.}]  (* oscillatory on the left, growing on the right *)
```

```mathematica
In[1]:= D[AiryBi[z], z]  (* derivative head *)
In[2]:= Series[AiryBi[x], {x, 0, 4}]  (* Maclaurin series; every third coefficient vanishes *)
```

```mathematica
In[1]:= AiryBi[2.5 + I]  (* complex argument, machine kernel *)
In[2]:= AiryBi[-50.]  (* far into the oscillatory region *)
In[3]:= AiryBi[1000.]  (* overflow past the double range returns an arbitrary-precision value *)
```

```mathematica
In[1]:= N[AiryBi[2], 50]  (* fifty digits through arb *)
In[2]:= AiryBi[1.8`30]  (* precision follows the input *)
In[3]:= N[AiryBi[-4 + 3 I], 30]  (* complex arbitrary precision *)
```

```mathematica
In[1]:= cf = Compile[{{v, _Real, 1}}, AiryBi[v]]; cf[{-2., 0., 2.}]  (* compiled, bit-identical to the scalar values *)
```

### Notes

- `AiryBi[z]` is the solution of `y'' == z y` that grows as `z -> +Infinity`. It is entire, `AiryBi[-Infinity] = 0`, and the origin value is `1/(3^(1/6) Gamma[2/3])`.
- Machine real and complex arguments use the raw-double kernel in `sf_machine_bessel.c`. On `[-10, 10]` it evaluates piecewise Chebyshev tables; beyond that the Poincaré expansions apply. For complex `z` it uses `I_{±1/3}`. That helper also serves packed lists, `NDArray` and `Compile[]`. Arbitrary precision goes through arb's `acb_hypgeom_airy`, and the MPFR core is the fallback.
- `D[AiryBi[z], z] = AiryBiPrime[z]`. `AiryBiPrime[-Infinity]` is left unevaluated.
- Attributes: `Listable`, `NumericFunction`, `Protected`.
- See also: `AiryAi`, `AiryBiPrime`, `BesselI`.
