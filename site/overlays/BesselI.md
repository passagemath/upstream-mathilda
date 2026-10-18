### Worked examples

```mathematica
In[1]:= BesselI[0, 0]  (* regular at the origin *)
In[2]:= BesselI[-1/2, 0]  (* a negative fractional order diverges there *)
In[3]:= BesselI[{0, 1, 2}, 1.]  (* listable over the order *)
```

```mathematica
In[1]:= BesselI[-1/2, x]  (* half-integer orders give hyperbolic closed forms; I is not even at half-integer order *)
In[2]:= BesselI[-3, x]  (* even in an integer order *)
In[3]:= BesselI[1, -z]  (* integer-order parity *)
In[4]:= Series[BesselI[0, x], {x, 0, 8}]  (* all-positive Frobenius series *)
```

```mathematica
In[1]:= BesselI[2.5, -3.0]  (* on the cut the principal value e^{i nu pi} I_nu(3) is purely imaginary *)
In[2]:= BesselI[3 + I, 1.5 - I]  (* complex order and argument through the MPFR core *)
In[3]:= N[BesselI[0, 1], 50]  (* fifty digits *)
In[4]:= BesselI[0, 800.]  (* beyond the double range the value comes back at arbitrary precision *)
```

```mathematica
In[1]:= N[BesselI[0, 2] BesselK[1, 2] + BesselI[1, 2] BesselK[0, 2], 30]  (* the Wronskian I_0 K_1 + I_1 K_0 = 1/z, here 1/2 *)
In[2]:= Integrate[BesselI[1, x], x]  (* from the integral tables *)
In[3]:= Limit[BesselI[0, x], x -> Infinity]  (* exponential growth *)
```

### Notes

- `BesselI[n, z]` is the modified Bessel function of the first kind. It grows like `e^z`, is even in an integer order, and has a cut along the negative real axis for non-integer order.
- Machine-number calls with a real order use the raw-double kernel in `sf_machine_bessel.c`: Cephes-style `I_0`/`I_1` on the real line, and otherwise the ascending series or Miller's backward recurrence. Packed lists, `NDArray` and `Compile[]` (scalar real, scalar complex and rank-1) call the same helper. Overflow, complex order and arbitrary precision fall through to the MPFR power or asymptotic series.
- `D[BesselI[n, z], z] = (BesselI[n-1, z] + BesselI[n+1, z])/2`.
- Attributes: `Listable`, `NumericFunction`, `Protected`.
- See also: `BesselK`, `BesselJ`, `BesselY`.
