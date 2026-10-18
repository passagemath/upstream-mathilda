### Worked examples

```mathematica
In[1]:= Zeta[6]  (* even positive integers: a rational multiple of a power of Pi *)
```

```mathematica
In[1]:= Table[Zeta[-n], {n, 1, 6}]  (* negative integers give -BernoulliB[n+1]/(n+1); the even ones are the trivial zeros *)
```

```mathematica
In[1]:= N[Zeta[5], 50]  (* odd values stay symbolic; N reaches them through arb *)
```

```mathematica
In[1]:= Zeta[0.5 + 14.134725141734693 I]  (* the first nontrivial zero, to machine precision *)
```

```mathematica
In[1]:= N[Zeta[2 + 3 I], 30]  (* a complex argument at 30 digits *)
```

```mathematica
In[1]:= Zeta[{0.5, 2.0, 3.0, -1.5}]  (* a machine-real list runs through the packed kernel *)
```

```mathematica
In[1]:= Series[Zeta[s], {s, 1, 2}]  (* the Laurent expansion at the pole defines the Stieltjes constants *)
```

```mathematica
In[1]:= Zeta[4, 5]  (* Hurwitz form at an integer shift: Zeta[4] minus a finite power sum *)
```

```mathematica
In[1]:= {Zeta[-2.5, -0.5], Zeta[3, -1.5]}  (* the two-argument form uses ((k+a)^2)^(-s/2), so it is real for negative a *)
```

### Notes

- `Zeta[s]` is the Riemann zeta function, `Zeta[s, a]` the generalized zeta
  Σ_{k≥0} ((k+a)²)^(−s/2). That is Mathematica's convention: for `Re a > 0` it
  equals `HurwitzZeta[s, a]`, and for negative real `a` it stays real where the
  principal-branch `HurwitzZeta` has poles or complex values.
- Exact closed forms: `Zeta[0] = -1/2`, `Zeta[1] = ComplexInfinity`, even positive
  integers give rational·π²ⁿ, negative integers give rationals, and
  `Zeta[Infinity] = 1`. Odd positive integers stay symbolic.
- Numeric surfaces: machine-precision real and complex `s` (with `Re s ≥ 1/2` for
  complex `s`) use one `double` kernel, shared bit-for-bit by the scalar,
  packed/`NDArray` and `Compile[]` paths. Arbitrary precision uses arb
  (`acb_dirichlet_zeta`) for `Zeta[s]` and a complex-MPFR Euler–Maclaurin kernel
  for `Zeta[s, a]`.
- `Compile[]` lowers `Zeta[s]` (scalar and array, real and complex) but not the
  two-argument form. Inside compiled code, use `HurwitzZeta[s, a]` instead.
- Attributes: `Listable`, `NumericFunction`, `Protected`.
- See also: `HurwitzZeta`, `LerchPhi`, `PolyLog`, `StieltjesGamma`, `BernoulliB`.
