---
status: Stable
references:
  - "Knuth, \"The Art of Computer Programming, Vol. 2: Seminumerical Algorithms\", on arbitrary-precision integer addition."
  - "Geddes, Czapor & Labahn, \"Algorithms for Computer Algebra\" (1992), on normal forms for sums."
---
### Worked examples

```mathematica
In[1]:= 2^100 + 3^50  (* integer overflow promotes to a bignum *)
```

```mathematica
In[1]:= 1/2 + 1/3 + 1/6  (* exact rationals reduce to an integer *)
```

```mathematica
In[1]:= a + a + b + 2 a  (* like terms collect *)
```

```mathematica
In[1]:= a + b - (a + b)  (* a negated sum is distributed so it cancels *)
```

```mathematica
In[1]:= 1. + Pi  (* an inexact summand numericalises the exact ones *)
```

```mathematica
In[1]:= Plus @@ Table[x^Mod[k, 3], {k, 10^5}]  (* hash-based collection: 10^5 terms in linear time *)
```

```mathematica
In[1]:= Series[Exp[x], {x, 0, 4}] + Series[Sin[x], {x, 0, 4}]  (* power series add through the series algebra *)
```

```mathematica
In[1]:= Interval[{1, 2}] + Interval[{-1, 3}]  (* interval arithmetic *)
```

```mathematica
In[1]:= NDArray[{1., 2., 3.}] + 10  (* element-wise on the raw buffer, with scalar broadcasting *)
```

```mathematica
In[1]:= Infinity - Infinity  (* indeterminate, with the Infinity::indet message *)
```

### Notes

- `Plus` is `Flat`, `Listable`, `NumericFunction`, `OneIdentity`, `Orderless` and `Protected`. `Plus[]` is `0` and `Plus[x]` is `x`.
- Like terms are collected structurally, by an exact match of the non-numeric part, through a hash table. The result is canonically sorted, and only when it is out of order. A collected term is a flat `Times`, so `2 a b` is `Times[2, a, b]`.
- A sum is not distributed over a product (`c (a + b)` stays as it is), and no algebraic identities are applied: `Sin[x]^2 + Cos[x]^2` needs `Simplify`.
- Integer sums overflow into GMP bignums. An inexact summand makes the whole numeric part inexact (`1. + Pi`), at the lowest precision present.
- Infinities follow Mathematica: a real infinity absorbs finite terms; `Infinity - Infinity` and `ComplexInfinity + ComplexInfinity` give `Indeterminate` with `Infinity::indet`. A directed complex infinity keeps its direction.
- Numeric surfaces. Packed lists add without unpacking, NDArrays add element-wise with scalar broadcasting (mismatched shapes give `NDArray::shape` and the sum stays unevaluated), and `Compile[]` lowers `Plus` for integer, real, complex and array arguments.
- See also: `Total`, `Times`, `Subtract`, `Accumulate`.
