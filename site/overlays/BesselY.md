### Worked examples

```mathematica
In[1]:= BesselY[0, 0]  (* logarithmic singularity at the origin *)
In[2]:= BesselY[-1/2, 0]  (* negative half-odd orders vanish there: Y_{-1/2} = J_{1/2} *)
In[3]:= BesselY[0, {1.0, 2.0, 3.0}]  (* listable, machine kernel *)
```

```mathematica
In[1]:= BesselY[3/2, x]  (* half-integer orders close in elementary functions *)
In[2]:= BesselY[-3, x]  (* negative integer order reflects: Y_{-n} = (-1)^n Y_n *)
In[3]:= Series[BesselY[0, x], {x, 0, 3}]  (* the logarithmic series DLMF 10.8.1 *)
```

```mathematica
In[1]:= BesselY[0, -1.0]  (* across the branch cut the value is complex *)
In[2]:= BesselY[2.5, 4.0 + 1.0 I]  (* fractional order, complex argument *)
In[3]:= BesselY[0.5 I, 3 - I]  (* complex order goes to the MPFR core *)
In[4]:= N[BesselY[0, 1], 50]  (* fifty digits *)
```

```mathematica
In[1]:= N[BesselJ[1, 5] BesselY[0, 5] - BesselJ[0, 5] BesselY[1, 5] - 2/(5 Pi), 30]  (* the Wronskian J_1 Y_0 - J_0 Y_1 = 2/(Pi z) holds to working precision *)
```

```mathematica
In[1]:= BesselY[1, -z]  (* no parity fold: the negative axis is a branch cut *)
In[2]:= BesselY[11/2, 1]  (* exact numeric arguments stay unevaluated *)
```

### Notes

- `BesselY[n, z]` is the Bessel function of the second kind. It is singular at the origin (`Y_0(0) = -Infinity`, `ComplexInfinity` for most other orders) and has a cut along the negative real axis.
- Machine-number calls with a real order use the shared raw-double kernel (`sf_machine_bessel.c`): libm `yn` for integer order with `x > 0`, and Amos's construction from `K_nu(-iz)` and `I_nu(-iz)` otherwise. The same helper serves packed lists, `NDArray` and `Compile[]` (scalar real, scalar complex and rank-1), bit for bit. Complex order, arbitrary precision and `|nu| > 1000` use MPFR: `mpfr_yn`, the logarithmic series, the connection formula or Hankel's expansion.
- Half-integer orders with a symbolic argument give spherical-Bessel closed forms.
- `D[BesselY[n, z], z] = (BesselY[n-1, z] - BesselY[n+1, z])/2`.
- Attributes: `Listable`, `NumericFunction`, `Protected`.
- See also: `BesselJ`, `BesselK`, `BesselI`.
