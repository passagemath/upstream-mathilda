### Worked examples

```mathematica
In[1]:= BesselJ[0, 0]  (* regular at the origin: J_0(0) = 1 *)
In[2]:= BesselJ[{0, 1, 2}, 2.5]  (* listable over the order, machine kernel *)
```

```mathematica
In[1]:= BesselJ[5/2, x]  (* half-integer orders close in elementary functions *)
In[2]:= BesselJ[1, -2 x]  (* integer-order parity J_n(-z) = (-1)^n J_n(z) *)
In[3]:= Series[BesselJ[0, x], {x, 0, 6}]  (* the alternating Frobenius series *)
```

```mathematica
In[1]:= BesselJ[1/2, -2.]  (* on the branch cut a real argument gives a complex value *)
In[2]:= BesselJ[1.5, 2.0 + 3.0 I]  (* complex argument at machine precision *)
In[3]:= BesselJ[7/3 + I, 4.5 - I]  (* complex order takes the MPFR core *)
In[4]:= N[BesselJ[0, 4], 50]  (* fifty digits from the arbitrary-precision path *)
```

```mathematica
In[1]:= BesselJ[1, NDArray[{0.5, 1., 1.5}]]  (* an NDArray runs through the element-wise kernel *)
In[2]:= cf = Compile[{{v, _Real, 1}}, BesselJ[1, v]]; cf[{0.5, 1., 1.5}]  (* compiled, bit-identical to the scalar values *)
```

```mathematica
In[1]:= BesselJ[11/2, 1]  (* an exact numeric argument stays unevaluated, as in Mathematica *)
In[2]:= N[BesselJ[11/2, 1]]  (* numericalise explicitly instead *)
```

### Notes

- `BesselJ[n, z]` is the Bessel function of the first kind, regular at the origin. `J_0(0) = 1`, `J_n(0) = 0` for integer `n != 0`, and a non-integer order with `Re nu < 0` gives `ComplexInfinity` at 0.
- Machine-number calls with a **real** order go to the raw-double kernel in `sf_machine_bessel.c`: libm `jn` for real integer orders, and Amos's I/K scheme otherwise. A packed list, a visible `NDArray` and `Compile[]` (scalar real, scalar complex and rank-1) call the same helper, so all five surfaces agree bit for bit. A complex order, an argument with more than machine precision, or `|nu| > 1000` uses the MPFR core (`mpfr_jn`, power series or Hankel asymptotic).
- Half-integer orders with a symbolic argument rewrite to `Sin`/`Cos` forms, and negative integer orders reflect by `J_{-n} = (-1)^n J_n`. These are C rules now, no longer `bessel.m` DownValues.
- Non-integer order has a branch cut along the negative real axis.
- `D[BesselJ[n, z], z] = (BesselJ[n-1, z] - BesselJ[n+1, z])/2`.
- Attributes: `Listable`, `NumericFunction`, `Protected`.
- See also: `BesselY`, `BesselI`, `BesselK`, `BesselJZero`, `AiryAi`.
