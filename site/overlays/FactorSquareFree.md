### Worked examples

```mathematica
In[1]:= FactorSquareFree[x^5 - x^4 - x + 1]  (* the square-free factor 1 + x + x^2 + x^3 is reducible but is not split *)
```

```mathematica
In[1]:= FactorSquareFree[x^8 + 4 x^6 + 6 x^4 + 4 x^2 + 1]  (* a perfect fourth power found without full factorisation *)
```

```mathematica
In[1]:= FactorSquareFree[x^3 y^2 + 2 x^2 y^3 + x y^4]  (* multivariate: square-free parts x and y (x + y), the latter squared *)
```

```mathematica
In[1]:= FactorSquareFree[x^2 - 2 Sqrt[2] x + 2, Extension -> Sqrt[2]]  (* a square over Q(Sqrt[2]) *)
```

```mathematica
In[1]:= FactorSquareFreeList[(x^2 + 1)^3 (x - 1)^2]  (* the companion list of {factor, multiplicity} pairs *)
```

```mathematica
In[1]:= FactorSquareFree[x^4 - 2 x^2 + 1.0]  (* inexact coefficients are rationalised and restored *)
```

```mathematica
In[1]:= FactorSquareFree[(x^2 - 1)/(x^2 + 2 x + 1)]  (* limitation: a rational function is returned unchanged *)
```

### Notes

- `FactorSquareFree[poly]` writes `poly` as a product of pairwise coprime,
  square-free factors raised to their multiplicities. It uses Yun's
  algorithm, built on GCDs of the polynomial with its derivative. It does
  not split a square-free factor that is reducible, which makes it cheaper
  than `Factor` when only the multiplicities matter.
- The grouping of factors is a valid square-free decomposition but not
  always the one Mathematica prints: `x (x y + y^2)^2` is equal to
  `x y^2 (x + y)^2`.
- **Extension.** `Extension -> α` (or a list, or `Automatic`) decomposes
  over the number field by way of `Factor`'s Trager engine.
- **Soundness.** If the computed factors do not reproduce the input
  exactly, the expanded polynomial is returned unfactored rather than an
  unequal product.
- **Not supported.** Rational functions and `Modulus` are not handled; such
  calls stay as they are.
- `Listable`, `Protected`. Symbolic only.
- See also: `FactorSquareFreeList`, `Factor`, `FactorTerms`,
  `PolynomialGCD`, `D`.
