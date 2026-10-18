# FactorSquareFree

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`FactorSquareFree[poly]`**

writes poly as a product of pairwise-coprime square-free factors, collecting repeated factors into powers.

**`FactorSquareFree[poly, Extension -> alpha]`**

performs the decomposition

<details>
<summary>Notes</summary>

Computed via the Yun / Musser square-free decomposition using polynomial GCDs of poly with its derivative; cheaper than full Factor and sufficient when only multiplicities are needed. over Q(alpha) (alpha = Sqrt\[c\], c^(1/n), I, or a list of generators), grouping the irreducible-over-Q(alpha) factors by multiplicity; the default Extension -\> None factors over Q.

</details>

## Examples (12)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= FactorSquareFree[x^5 - x^3 - x^2 + 1]
Out[1]= (-1 + x)^2 (1 + 2 x + 2 x^2 + x^3)

In[2]:= FactorSquareFree[x^4 - 9x^3 + 29x^2 - 39x + 18]
Out[2]= (-3 + x)^2 (2 - 3 x + x^2)

In[3]:= FactorSquareFree[x^5 - x^3 y^2 - x^2 y^3 + y^5]
Out[3]= (x - y)^2 (x^3 + 2 x^2 y + 2 x y^2 + y^3)

In[4]:= FactorSquareFree[{(x^2 - 1)(x - 1), (x^4 - 1)(x^2 - 1)}]
Out[4]= {(1 + x) (-1 + x)^2, (1 + x^2) (-1 + x^2)^2}
```

### Worked examples (1)

```mathematica
In[5]:= FactorSquareFree[x^2 - 2 Sqrt[2] x + 2, Extension -> Sqrt[2]]
Out[5]= (-Sqrt[2] + x)^2
```

### Applications (7)

The square-free factor 1 + x + x^2 + x^3 is reducible but is not split

```mathematica
In[6]:= FactorSquareFree[x^5 - x^4 - x + 1]
Out[6]= (-1 + x)^2 (1 + x + x^2 + x^3)
```

A perfect fourth power found without full factorisation

```mathematica
In[7]:= FactorSquareFree[x^8 + 4 x^6 + 6 x^4 + 4 x^2 + 1]
Out[7]= (1 + x^2)^4
```

Multivariate: square-free parts x and y (x + y), the latter squared

```mathematica
In[8]:= FactorSquareFree[x^3 y^2 + 2 x^2 y^3 + x y^4]
Out[8]= x (x y + y^2)^2
```

A square over Q(Sqrt[2])

```mathematica
In[9]:= FactorSquareFree[x^2 - 2 Sqrt[2] x + 2, Extension -> Sqrt[2]]
Out[9]= (-Sqrt[2] + x)^2
```

The companion list of {factor, multiplicity} pairs

```mathematica
In[10]:= FactorSquareFreeList[(x^2 + 1)^3 (x - 1)^2]
Out[10]= {{1, 1}, {-1 + x, 2}, {1 + x^2, 3}}
```

Inexact coefficients are rationalised and restored

```mathematica
In[11]:= FactorSquareFree[x^4 - 2 x^2 + 1.0]
Out[11]= (-1.0 + x^2)^2
```

Limitation: a rational function is returned unchanged

```mathematica
In[12]:= FactorSquareFree[(x^2 - 1)/(x^2 + 2 x + 1)]
Out[12]= (-1 + x^2)/(1 + 2 x + x^2)
```

## Implementation notes

**Algorithm.** `builtin_factorsquarefree` (in `src/poly/facpoly_squarefree.inc`, compiled into `facpoly.c`) computes the square-free decomposition `p = ∏ p_i^i` with each `p_i` square-free and pairwise coprime, using **Yun's algorithm**. A trailing `Extension -> α` (an explicit generator, a list, or `Automatic`; anything other than `None`) goes to `factor_square_free_over_extension`. That routine calls `Factor[poly, opts…]`, the Trager engine over `Q(α)`, and rebuilds the square-free decomposition by grouping irreducible factors of equal multiplicity, so it needs no separate square-free primitive over `Q(α)`. If `Factor` declines, it falls back to the path over Q. With FLINT, a free-symbol radicand `Sqrt[k]` goes to `flint_parametric_sqrt_factor_squarefree`. Any other trailing argument, `Modulus` included, leaves the call unevaluated. After a `Rationalize`/`Numericalize` round-trip for inexact inputs, the dispatcher `factor_square_free_dispatcher` expands and routes to `factor_square_free_poly`. `FactorSquareFreeList` splits the result into `{factor, multiplicity}` pairs.

`factor_square_free_poly` works recursively on the variable list. It first extracts the polynomial content (`poly_content`) and recurses on it; on the primitive part `pp` it runs Yun's recurrence: set `A = pp`, `B = gcd(A, A')`, `C = A/B`, `D = A'/B − C'`, then iterate `P_i = gcd(C, D)` (the degree-`i` square-free factor), `C ← C/P_i`, `D ← D/P_i − (C/P_i)'`, accumulating `P_i^i`. Derivatives are computed by the local `poly_deriv` and GCDs by `poly_gcd_internal`; exact divisions by `exact_poly_div`. An F4 "cheap squarefree" pre-check (`sqfree_cheap_check`, a univariate-substitution probe) short-circuits the expensive `gcd(pp, pp')` when `pp` is provably square-free. A final `exact_poly_div` recovers any leading or missing scalar factor, so the product round-trips to the input. Since v0.458 that division is a soundness check as well. If the computed factors do not divide `P` exactly, as can happen when the "variables" are transcendental atoms (`Log[…]`, `Arg[…]`) over a coefficient ring the classical PRS does not model, the routine returns `Expand[P]` unfactored. Before that change it substituted `1` for the cofactor and returned a product unequal to the input.

**Data structures.** `Expr*` polynomials manipulated through the generic evaluator (`eval_and_free`) and the poly helpers (`get_degree_poly`, `poly_gcd_internal`, `exact_poly_div`); factors collected in a growable `Expr**` buffer with their multiplicities encoded as `Power[P_i, i]`.

**Complexity / limits.** Dominated by the GCD chain; the cheap pre-check avoids one full multivariate GCD in the common square-free case. Square-free decomposition is a prerequisite stage for full `Factor`. Only polynomials are decomposed: a rational function comes back unchanged, and `Modulus` is not supported. `FactorSquareFree` is `Listable`, `Protected`, and symbolic only.

- `Listable`, `Protected`.
- Automatically threads over lists, as well as equations, inequalities and logic functions.
- Works on both univariate and multivariate polynomials.
- Bigint coefficients (above 2^63) are coefficients, not variables: `FactorSquareFree[Expand[3^45 (1 + t)^2]]` gives `2954312706550833698643 (1 + t)^2` (v0.458; it used to come back unfactored).
- Self-checking: if the computed content and square-free factors do not reproduce the input exactly, which can happen when the "variables" are transcendental atoms over a coefficient ring the classical PRS does not model, the input is returned unfactored, never as an unequal product (v0.458).
- `Extension` runs the decomposition over an algebraic number field Q(alpha).  It reduces to `Factor[poly, Extension -> alpha]` (Trager engine) and groups the irreducible-over-Q(alpha) factors by multiplicity, so it also accepts inputs with algebraic-number coefficients, e.g. `FactorSquareFree[x^2 - 2 Sqrt[2] x + 2, Extension -> Sqrt[2]]` → `(-Sqrt[2] + x)^2`.  Because square-free structure is field-independent for rational inputs, an extension leaves a rational-coefficient answer unchanged; the default `Extension -> None` uses the fast over-Q Yun path.
- Multivariate inputs use a cheap squarefree pre-check (F4 Stage 1): after content extraction in the main variable, `sqfree_cheap_check` substitutes integer values from `{1, -1, 2, -2, 3, -3, 4}` for the other variables and tests `gcd(image, image')` over `Z[x]`.  If any image is squarefree at an alpha that preserves the leading-x degree, the pre-check proves `pp` is squarefree in x and the expensive multivariate `gcd(pp, pp')` is skipped.  Soundness comes from content extraction guaranteeing any repeated factor of `pp` involves the main variable nontrivially.  Measured 6.6× speedup on 4-variable squarefree inputs (6.27 s → 0.95 s); non-squarefree inputs fall through to the original Yun loop with negligible overhead.
- The cheap pre-check's univariate `gcd(image, image')` runs through `zupoly_gcd` (subresultant PRS, GMP `mpz_t` coefficients).  The previous implementation used `poly_gcd_internal` (Knuth-style primitive PRS at the Expr level) which suffers exponential coefficient growth on the intermediate pseudo-remainders; on a degree-31 univariate image (e.g. `Factor[Expand[x^2 (z^13 - x^12)(z^4 + 3 x^9 - y^13)(17 - 5 y - z^14)]]`) it ran for >120 s.  Routing the same gcd through subresultant PRS keeps coefficient sizes polynomially bounded and runs in sub-millisecond time on the same input, bringing the full Factor call to under 1 s.

**Attributes:** `Listable`, `Protected`.

## References

**See also:** [I](../../mathematical-constants/I/)

- D. Y. Y. Yun, "On square-free decomposition algorithms", SYMSAC 1976.
- K. O. Geddes, S. R. Czapor, G. Labahn, *Algorithms for Computer Algebra* (Kluwer, 1992).
- Source: [`src/poly/facpoly_squarefree.inc`](https://github.com/stblake/mathilda/blob/main/src/poly/facpoly_squarefree.inc)
- Specification: [`docs/spec/builtins/structural-manipulation.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/structural-manipulation.md)
- Tests: [`tests/test_facpoly.c`](https://github.com/stblake/mathilda/blob/main/tests/test_facpoly.c)
- Tests: [`tests/test_factorsquarefreelist.c`](https://github.com/stblake/mathilda/blob/main/tests/test_factorsquarefreelist.c)
- Tests: [`tests/test_flint_bridge.c`](https://github.com/stblake/mathilda/blob/main/tests/test_flint_bridge.c)
- Tests: [`tests/test_inexact_dispatch.c`](https://github.com/stblake/mathilda/blob/main/tests/test_inexact_dispatch.c)

## Notes & additional examples

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
