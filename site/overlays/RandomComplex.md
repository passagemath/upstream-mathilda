### Worked examples

```mathematica
In[1]:= SeedRandom[1]; RandomComplex[]  (* uniform in the unit square with corners 0 and 1 + I *)
```

```mathematica
In[1]:= SeedRandom[1]; RandomComplex[2 + 2 I]  (* a bare corner z gives the rectangle from 0 to z; same draw, scaled by 2 *)
```

```mathematica
In[1]:= SeedRandom[1]; RandomComplex[{-1 - I, 1 + I}, {2, 2}]  (* a 2x2 array of points in the square centred at the origin *)
```

```mathematica
In[1]:= SeedRandom[5]; RandomComplex[{-Pi - I, Pi + I}]  (* symbolic corners are reduced to numbers first *)
```

```mathematica
In[1]:= SeedRandom[7]; pts = RandomComplex[{-1 - 2 I, 3 + I}, 10^5]; {PackedArrayQ[pts], MinMax[Re[pts]], MinMax[Im[pts]]}  (* a packed complex64 buffer whose parts stay inside the rectangle *)
```

```mathematica
In[1]:= SeedRandom[7]; Mean[RandomComplex[{-1 - 2 I, 3 + I}, 10^5]]  (* the sample mean is close to the centre 1 - I/2 *)
```

```mathematica
In[1]:= SeedRandom[11]; Count[RandomComplex[1 + I, 10^4], w_ /; Abs[w] < 1] / 10.^4  (* about Pi/4 = 0.785 of the unit square lies inside the unit circle *)
```

```mathematica
In[1]:= SeedRandom[3]; a = RandomComplex[{0, 1 + I}, {2, 2}]; SeedRandom[3]; b = Table[RandomComplex[{0, 1 + I}], 2, 2]; a == b  (* the packed fill draws in the same order as repeated scalar calls *)
```

```mathematica
In[1]:= SeedRandom[4]; Precision /@ Re[RandomComplex[{0, 1 + I}, 3, WorkingPrecision -> 20]]  (* above machine precision the parts are MPFR numbers *)
```

```mathematica
In[1]:= RandomComplex[x]  (* a non-numeric range stays unevaluated *)
```

### Notes

- The range is an **axis-aligned rectangle**, not a disk. The real and imaginary parts are two independent uniform draws, real part first. A bare `z` means the rectangle with corners `0` and `z`; `{z1, z2}` means the rectangle spanned by `z1` and `z2`. Corners can be any numeric expression, for example `-Pi - I`.
- Machine draws come from the same xoshiro256++ stream that `RandomReal` uses (53 bits per part), so `SeedRandom` makes them reproducible. `WorkingPrecision -> d` above machine precision draws each part with `mpfr_urandomb` from the GMP generator, which is seeded together with the machine stream.
- `RandomComplex[range, n]` and `RandomComplex[range, {n1, ...}]` with a truly complex range fill a packed `complex64` buffer directly, so `Sin`, `Exp`, `HermiteH` and the reductions run on the buffer. A range with zero imaginary extent (for example `{0, 5}`) folds to Reals and returns an ordinary unpacked List, even at large sizes.
- `RandomComplex` has no `Compile[]` lowering: `CompileDiagnostics` reports it as the unlowered subexpression. Draw the array outside `Compile` and pass it in.
- The precision of a `WorkingPrecision -> d` draw is reported from the MPFR bit count, so it is slightly above `d` (`20.169` for `d = 20`).
- Attributes: `Protected`. See also: `RandomReal`, `RandomInteger`, `RandomVariate`, `SeedRandom`.
