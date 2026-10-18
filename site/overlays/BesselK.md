### Worked examples

```mathematica
In[1]:= BesselK[0, 0]  (* a directed divergence at the origin *)
In[2]:= BesselK[I, 0]  (* pure-imaginary order oscillates boundedly, so there is no limit *)
In[3]:= BesselK[{1, 2, 3}, 1.0]  (* listable over the order *)
```

```mathematica
In[1]:= BesselK[3/2, x]  (* half-integer orders give decaying elementary forms *)
In[2]:= BesselK[-2, x]  (* even in the order *)
In[3]:= Series[BesselK[0, x], {x, 0, 3}]  (* logarithmic expansion at the origin *)
In[4]:= Limit[BesselK[1, x], x -> Infinity]  (* exponential decay *)
```

```mathematica
In[1]:= BesselK[1/3, 2.0 + 1.0 I]  (* fractional order, complex argument, machine kernel *)
In[2]:= BesselK[0, -1.0]  (* across the cut: K_0(1) - i Pi I_0(1) *)
In[3]:= BesselK[1 + I, 3.0 - 2 I]  (* complex order takes the MPFR core *)
In[4]:= N[BesselK[0, 4], 50]  (* fifty digits *)
In[5]:= BesselK[0, 800.]  (* an underflowing double is promoted rather than flushed to zero *)
```

```mathematica
In[1]:= BesselK[0, NDArray[{0.5, 1., 2.}]]  (* an NDArray runs through the element-wise kernel *)
In[2]:= BesselK[1, -z]  (* no parity fold across the branch cut *)
```

### Notes

- `BesselK[n, z]` is the modified Bessel function of the second kind. It decays like `e^{-z}`, is even in the order, diverges at the origin, and has a cut along the negative real axis.
- Machine-number calls with a real order use the raw-double kernel in `sf_machine_bessel.c`. Integer orders take table-fitted `K_0`/`K_1` and upward recurrence. Fractional orders take Temme's series, tabulated Taylor nodes or Steed's CF2, plus Hankel. Packed lists, `NDArray` and `Compile[]` (scalar real, scalar complex and rank-1) call the same helper, bit for bit. Complex order, arbitrary precision and underflow go to the MPFR core: the logarithmic series, the connection formula or the asymptotic series.
- `D[BesselK[n, z], z] = -(BesselK[n-1, z] + BesselK[n+1, z])/2`.
- Attributes: `Listable`, `NumericFunction`, `Protected`.
- See also: `BesselI`, `BesselJ`, `BesselY`, `AiryAi`.
