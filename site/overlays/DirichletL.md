### Worked examples

```mathematica
In[1]:= DirichletL[4, 2, 1]  (* the Leibniz series 1 - 1/3 + 1/5 - ... = Pi/4 *)
```

```mathematica
In[1]:= DirichletL[3, 2, 1]  (* the real character mod 3 *)
```

```mathematica
In[1]:= DirichletL[8, 2, 1]  (* an even character: the logarithmic class-number form *)
```

```mathematica
In[1]:= DirichletL[12, 3, 1]  (* an imprimitive character: Pi/4 times the Euler factor 1 + 1/3 *)
```

```mathematica
In[1]:= {DirichletL[4, 2, 0], DirichletL[4, 2, -2], DirichletL[4, 2, -3]}  (* non-positive integers give generalized Bernoulli numbers; odd characters vanish at negative odd s *)
```

```mathematica
In[1]:= DirichletL[5, 2, -4]  (* a complex character gives a value in Q(I) *)
```

```mathematica
In[1]:= {N[DirichletL[4, 2, 2], 30], N[Catalan, 30]}  (* the value at s = 2 for the character mod 4 is Catalan's constant *)
```

```mathematica
In[1]:= DirichletL[4, 2, 0.5 + 6.0209489 I]  (* near the first zero on the critical line *)
```

```mathematica
In[1]:= N[DirichletL[5, 2, 5/2 + I], 25]  (* a complex character at complex s, 25 digits through arb *)
```

```mathematica
In[1]:= DirichletL[3, 2, Range[0.5, 2., 0.5]]  (* a packed Range of s runs the machine kernel over the buffer *)
```

```mathematica
In[1]:= DirichletL[5, 1, Interval[{1.23, 1.24}]]  (* a rigorous enclosure for a real character *)
```

```mathematica
In[1]:= DirichletL[4, 2, 2]  (* an exact s >= 2 stays unevaluated, as in Mathematica *)
```

### Notes

- `DirichletL[k, j, s]` is `L(χ, s) = Σ χ(n) n^(-s)`, continued analytically to the whole plane, for `χ = DirichletCharacter[k, j, ·]`. It uses the same labelling as `DirichletCharacter`.
- Exact values are given at `s ≤ 0`, where they are generalized Bernoulli numbers, and at `s = 1` for every non-principal character. This includes imprimitive characters, which Mathematica leaves unevaluated at `s = 1`. Exact `s ≥ 2` stays unevaluated, as in Mathematica.
- Numerics: a machine `s` uses an Euler–Maclaurin / functional-equation kernel (max relative error about 6e-14 against arb in testing). Arbitrary precision uses arb's `acb_dirichlet_l`, and an `Interval` `s` with a real character gives a rigorous enclosure. `N` keeps `k` and `j` exact.
- Surfaces: packed and `NDArray` `s` with scalar labels, and `Compile[]` at scalar and array shapes, all share the machine kernel bit for bit.
- Attributes: `Listable`, `NumericFunction`, `Protected`. `DirichletL[1, 1, s]` is `Zeta[s]` for numeric `s`.
- See also: `DirichletCharacter`, `Zeta`, `HurwitzZeta`, `Catalan`.
