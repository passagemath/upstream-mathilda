### Worked examples

```mathematica
In[1]:= HurwitzZeta[2, 1]  (* reduces to Zeta[2] when the second argument is 1 *)
```

```mathematica
In[1]:= HurwitzZeta[s, 1/2]  (* the half-integer shift: (2^s - 1) Zeta[s] *)
```

```mathematica
In[1]:= HurwitzZeta[4, 5]  (* positive integer second argument: Zeta[4] minus a finite power sum *)
```

```mathematica
In[1]:= HurwitzZeta[-3, -2]  (* at a non-positive integer a, -BernoulliB[4, -2]/4 *)
```

```mathematica
In[1]:= HurwitzZeta[2, -1]  (* principal powers keep the pole at a = -1 *)
```

```mathematica
In[1]:= HurwitzZeta[3.0, {0.5, 1.0, 2.5}]  (* machine reals, threaded over a list *)
```

```mathematica
In[1]:= HurwitzZeta[2.0 + 1.0 I, 0.5 + 0.5 I]  (* complex s and a on the machine kernel *)
```

```mathematica
In[1]:= N[HurwitzZeta[1/3, 2/7], 40]  (* arbitrary precision via arb *)
```

```mathematica
In[1]:= {HurwitzZeta[3, -1.5], Zeta[3, -1.5]}  (* principal versus symmetric powers for negative a *)
```

### Notes

- `HurwitzZeta[s, a] = Σ_{k≥0} (k+a)^(−s)`, continued to all `s ≠ 1`. Unlike
  two-argument `Zeta`, it uses *principal-branch* powers. It therefore keeps
  poles at `a = 0, -1, -2, …` and differs from `Zeta[s, a]` for negative real `a`.
- Exact reductions: `s == 1` gives `ComplexInfinity`, `a == 1` gives `Zeta[s]`,
  `a == 1/2` gives `(2^s - 1) Zeta[s]`, and a positive integer `a` gives
  `Zeta[s]` minus a finite power sum. At a non-positive integer `a`, the result
  is `ComplexInfinity` for positive integer `s` and the Bernoulli-polynomial
  value for non-positive integer `s`.
- Numeric surfaces: machine-precision real and complex arguments with
  `Re s ≥ 1/2` and `Re a > 0` use one Euler–Maclaurin `double` kernel, shared by
  the scalar, packed/`NDArray` and `Compile[]` paths (scalar and array, real and
  complex). Higher precision goes to arb's `acb_dirichlet_hurwitz`, and other
  inputs go to a complex-MPFR Euler–Maclaurin kernel.
- Known gap: for a negative non-integer real `a` with non-integer real `s`
  (e.g. `HurwitzZeta[-2.5, -0.5]`), the principal value is complex, but only its
  real part is returned.
- Attributes: `Listable`, `NumericFunction`, `Protected`.
- See also: `Zeta`, `LerchPhi`, `PolyGamma`, `BernoulliB`.
