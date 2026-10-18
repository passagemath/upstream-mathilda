### Worked examples

```mathematica
In[1]:= PossibleZeroQ[(x - 1) (x + 1) - (x^2 - 1)]  (* a rational identity, decided exactly by Together and Cancel *)
```

```mathematica
In[1]:= PossibleZeroQ[Sqrt[2] + Sqrt[3] - Sqrt[5 + 2 Sqrt[6]]]  (* a nested-radical denesting, confirmed on the MPFR ladder *)
```

```mathematica
In[1]:= PossibleZeroQ[Tan[Pi/2 - x] - Cot[x]]  (* constant-phase angle expansion in the exact trig/exp test *)
```

```mathematica
In[1]:= PossibleZeroQ[Cosh[14]^30 - (1 - Tanh[14]^2)^(-15)]  (* huge operands: the cancellation is judged by its trend up to 1000 bits *)
```

```mathematica
In[1]:= PossibleZeroQ[Cosh[14]^30 - (1 - Tanh[14]^2)^(-15) + 1]  (* the same scale with a genuine +1 is rejected *)
```

```mathematica
In[1]:= PossibleZeroQ[PolyLog[2, x^2] + 1]  (* numeric-function declines are re-drawn, not taken as undecidable *)
```

```mathematica
In[1]:= PossibleZeroQ[-E^(-a x)]  (* structural certificate: an exponential never vanishes *)
```

```mathematica
In[1]:= PossibleZeroQ[Sqrt[x^2] - x, Assumptions -> x > 0]  (* sampling is confined to the assumed region *)
```

```mathematica
In[1]:= PossibleZeroQ[{Max[x, y] - x, Max[x, y] - y}, Assumptions -> x > y]  (* a coupling inequality is enforced by rejection sampling *)
```

```mathematica
In[1]:= PossibleZeroQ[Exp[2 Pi I k] - 1, Assumptions -> Element[k, Integers]]  (* zero on the integer lattice only *)
```

### Notes

- **Heuristic, as in Mathematica.** `PossibleZeroQ` runs exact symbolic
  stages first (rational normalisation, the exact trig/exp kernel test,
  structural non-zero certificates). Then come numeric evaluation on an MPFR
  precision ladder and Schwartz–Zippel sampling at random points. `False`
  means a definite non-zero was found. `True` means the expression was
  proved zero, or behaved like zero at every point tested; when unsure it
  answers `True` and prints no message.
- **Assumptions.** `PossibleZeroQ[e, Assumptions -> a]` and an enclosing
  `Assuming` restrict the sample points to the assumed region. That covers
  sign and interval bounds, `Element[x, Integers]`,
  `Element[x, Complexes]`, and inequalities coupling several variables. An
  option replaces `$Assumptions`, and a positional assumption is conjoined
  with it.
- **Real sampling by default.** Without assumptions, free symbols are
  sampled as reals. Functions that are not analytic, such as `Conjugate`,
  `Re` and `Im`, then look like the identity:
  `PossibleZeroQ[Conjugate[x] - x]` is `True`. Pass
  `Assumptions -> Element[x, Complexes]` to sample off the real line, which
  gives `False`.
- **Threading.** A list in the first argument threads with the options
  broadcast. The head is `Protected` and not `Listable`, so a list
  assumption is read as one conjunction.
- See also: `Simplify`, `Refine`, `Equal`, `ComplexExpand`.
