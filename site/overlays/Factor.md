---
status: Stable
references:
  - "B. M. Trager, \"Algebraic factoring and rational function integration\", SYMSAC 1976 — the norm / sqfr_norm / alg_factor approach used for the Extension path."
  - "Geddes, Czapor & Labahn, \"Algorithms for Computer Algebra\" (1992), Ch. 8 (polynomial factorization)."
---
### Worked examples

```mathematica
In[1]:= Factor[6 x^2 + 7 x + 2]  (* factorisation over the integers with a non-monic leading coefficient *)
```

```mathematica
In[1]:= Factor[x^10 - 1]  (* cyclotomic factors of x^10 - 1 *)
```

```mathematica
In[1]:= Length[Factor[x^99 - y^99]]  (* the multivariate FLINT path: one factor for each divisor of 99 *)
```

```mathematica
In[1]:= Factor[(x^3 + 2 x^2)/(x^2 - 4 y^2) - (x + 2)/(x^2 - 4 y^2)]  (* a rational function: combined first, then numerator and denominator factored *)
```

```mathematica
In[1]:= Factor[x^4 + 1, Extension -> Sqrt[2]]  (* an irreducible quartic splits over Q(Sqrt[2]) *)
```

```mathematica
In[1]:= Factor[x^4 - 5 x^2 + 6, Extension -> {Sqrt[2], Sqrt[3]}]  (* a tower of extensions, reduced to one primitive element *)
```

```mathematica
In[1]:= Factor[x^2 - 2 Sqrt[2] x + 2, Extension -> Automatic]  (* the extension is detected from the coefficients *)
```

```mathematica
In[1]:= Factor[x^4 + 1, Modulus -> 3]  (* factorisation in F_3[x] *)
```

```mathematica
In[1]:= Factor[E^(2 x) - 1]  (* a polynomial in the kernel E^x *)
```

```mathematica
In[1]:= Factor[Integrate[f[x], x] x^2 - Integrate[f[x], x]]  (* an unevaluated integral is carried as an opaque factor *)
```

### Notes

- Over the integers, `Factor` returns irreducible factors with the integer
  content pulled out front. `x^2 + 1` stays as it is because it has no
  rational factors. Factors print in canonical term order, constant term
  first (`-1 + x`, not `x - 1`).
- **Engines.** With FLINT, univariate integer polynomials go to
  `fmpz_poly_factor` and multivariate rational polynomials to
  `fmpq_mpoly_factor`. Otherwise the native Berlekamp–Zassenhaus and Hensel
  code runs, with heuristics for three or more variables. `Extension -> α`
  (or a list, or `Automatic`) factors over the number field with Trager's
  algorithm. `Modulus -> p` factors over `F_p`. A composite modulus is left
  unevaluated.
- **Rational functions and kernels.** A rational expression is combined
  with `Together` first. Exponentials, trig functions and fractional powers
  are factored as polynomials in their kernels. Unevaluated `Integrate`,
  `D` and `Sum` are opaque atoms.
- **Threading.** `Listable`, and threads over equations, inequalities and
  logical connectives as well.
- `Protected`. Symbolic only.
- See also: `FactorList`, `FactorSquareFree`, `FactorTerms`, `Expand`,
  `Together`, `FactorInteger`.
