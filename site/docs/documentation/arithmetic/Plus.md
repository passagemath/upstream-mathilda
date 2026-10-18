# Plus

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

x + y + ... or Plus\[x, y, ...\] represents a sum of terms. Plus is Flat, Orderless, OneIdentity, Listable, and NumericFunction: nested Plus is auto-flattened, terms are sorted into canonical order, like terms are combined, and integer arguments are summed exactly via GMP (with int64 fast path and BigInt overflow promotion).

## Examples (14)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= 1 + 2 + x + 2*x
Out[1]= 3 + 3 x

In[2]:= a + b - (a + b)
Out[2]= 0

In[3]:= Interval[{1, 2}] + Interval[{-1, 3}]
Out[3]= Interval[{0, 5}]

In[4]:= NDArray[{1., 2., 3.}] + 10
Out[4]= NDArray[{11.0, 12.0, 13.0}]
```

### Applications (10)

Integer overflow promotes to a bignum

```mathematica
In[5]:= 2^100 + 3^50
Out[5]= 1267651318126217093349291975625
```

Exact rationals reduce to an integer

```mathematica
In[6]:= 1/2 + 1/3 + 1/6
Out[6]= 1
```

Like terms collect

```mathematica
In[7]:= a + a + b + 2 a
Out[7]= 4 a + b
```

A negated sum is distributed so it cancels

```mathematica
In[8]:= a + b - (a + b)
Out[8]= 0
```

An inexact summand numericalises the exact ones

```mathematica
In[9]:= 1. + Pi
Out[9]= 4.14159
```

Hash-based collection: 10^5 terms in linear time

```mathematica
In[10]:= Plus @@ Table[x^Mod[k, 3], {k, 10^5}]
Out[10]= 33333 + 33334 x + 33333 x^2
```

Power series add through the series algebra

```mathematica
In[11]:= Series[Exp[x], {x, 0, 4}] + Series[Sin[x], {x, 0, 4}]
Out[11]= 1 + 2 x + 1/2 x^2 + 1/24 x^4 + O[x]^5
```

Interval arithmetic

```mathematica
In[12]:= Interval[{1, 2}] + Interval[{-1, 3}]
Out[12]= Interval[{0, 5}]
```

Element-wise on the raw buffer, with scalar broadcasting

```mathematica
In[13]:= NDArray[{1., 2., 3.}] + 10
Out[13]= NDArray[{11.0, 12.0, 13.0}]
```

Indeterminate, with the Infinity::indet message

```mathematica
In[14]:= Infinity - Infinity
Out[14]= Indeterminate
```

## Implementation notes

**Algorithm.** `Plus` is `Flat | Listable | NumericFunction | OneIdentity | Orderless | Protected`, so the evaluator has already flattened nested sums and threaded over `List` arguments before `builtin_plus` runs. `Plus[]` is `0` and `Plus[x]` is `x`. The builtin then works in a fixed order. (1) **Machine-integer fast path.** When every argument is an `int64`, they are summed with overflow-checked `ci_add_i64`. On overflow it falls through, and the generic path promotes the sum to a GMP bigint. (2) `algnum_plus_combine` adds same-generator `AlgebraicNumber` terms in the number field, folding loose rationals in. (3) A single **fused guard pass** checks whether any argument is an NDArray, a `-(...)` sum, a `SeriesData`, an inexact number or an infinity. For an ordinary symbolic or rational sum, which is the common case, it jumps straight to term grouping. Otherwise the special passes run in order. NDArray operands add element-wise over raw buffers (`ndarray_elementwise`, which broadcasts numeric scalars); mismatched shapes give `NDArray::shape` and the sum stays unevaluated. `Times[-1, Plus[...]]` is distributed, so `a + b - (a + b)` cancels. `SeriesData` summands fold through `series_combine_plus`. Inexact contagion (`numeric_contagion_args`) numericalises the exact numeric parts, so `1. + Pi` is `4.14159`. Finally, Infinity and Indeterminate follow Mathematica's table: `Infinity - Infinity` and any `ComplexInfinity` clash give `Infinity::indet` and `Indeterminate`, and a finite part is absorbed by a real infinity. (4) **Grouping.** `get_coeff_base` splits each term into a numeric coefficient and a base. Bare numbers go to a running numeric total via `add_numbers`, which also handles Rational, bigint, Complex, MPFR and `Interval` (`interval_add`). Bases are collected through an open-addressing hash table keyed by `expr_hash`, and like coefficients are summed. (5) **Reassembly.** Zero groups are dropped, and each `c·base` is emitted as a *flat* `Times`, with the coefficient spliced in front of the base's factors (never `Times[2, Times[a, b]]`). The result is sorted with `expr_orderless_sort` only when an O(n) check finds it out of order.

**Data structures.** There is a `TermGroup {base, coeff}` array and a power-of-two open-addressing table (`slot_group`/`slot_hash`, at least 2n slots). For sums of at most 8 terms both live on the stack, as does the result vector, so a small `Plus` in a tight loop allocates nothing beyond its result. Overflow is propagated as `Overflow[]`. The packed and NDArray surfaces work through `pack.c`: `Plus` is on the `AWARE` and broadcast lists, so packed integer or real lists add without unpacking (`Range[10^6] + Range[10^6]` stays packed). `Compile[]` lowers `Plus` at scalar and array shapes, with overflow-checked integer arithmetic and real, complex and array arms (`compile_infer.c`, `compile_emit_arith.c`).

**Complexity / limits.** Like-term collection is O(n) expected, through the hash table plus one canonical sort of the collapsed terms only, so `Plus @@` of 10^5 monomials does not pay an O(n²) scan or an O(n log n) sort of the raw input. Two caveats. A product times a sum, `c (a + b)`, is not distributed, as in Mathematica. And like-term collection is purely structural, so `Sin[x]^2 + Cos[x]^2` is left alone: that is `Simplify`'s job.

- `Flat`, `Orderless`, `Listable`.
- Combines numeric constants and collects like terms (e.g., `x + 2x` becomes `3x`).
- Like-term collection is hash-based (O(N) expected) and canonically sorts only
  the collapsed result, so a large flat sum (e.g. `Total` of 10^5 monomials) does
  not pay the generic `Orderless` sort over all N raw terms. `Times` collects
  factors the same way.
- A collected term is a flat `Times`: `2 a b + x` holds `Times[2, a, b]`, never
  `Times[2, Times[a, b]]`. The coefficient is spliced in front of the
  monomial's factors.
- Returns `0` if no arguments are provided.
- Returns `Overflow[]` if integer addition overflows or if any argument is `Overflow[]`.
- A finite term plus a real `±Infinity` absorbs the finite part (`Infinity + 3 → Infinity`),
  but a finite term plus a **directed (complex) infinity** keeps its direction rather
  than collapsing to a real `Infinity`: `1/2 - I Infinity` stays
  `Plus[1/2, Times[Complex[0, -1], Infinity]]` (so a Mellin–Barnes / Bromwich
  vertical line `s = c ± I Infinity` remains recognisable to `Integrate`).

**Attributes:** `Flat`, `Listable`, `NumericFunction`, `OneIdentity`, `Orderless`, `Protected`.

## References

**See also:** [Flat](../../expression-information/Flat/), [Orderless](../../expression-information/Orderless/), [Total](../../arithmetic/Total/), [Times](../../arithmetic/Times/), [Integrate](../../calculus/Integrate/)

- Knuth, "The Art of Computer Programming, Vol. 2: Seminumerical Algorithms", on arbitrary-precision integer addition.
- Geddes, Czapor & Labahn, "Algorithms for Computer Algebra" (1992), on normal forms for sums.
- Source: [`src/plus.c`](https://github.com/stblake/mathilda/blob/main/src/plus.c)
- Specification: [`docs/spec/builtins/arithmetic.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/arithmetic.md)
- Tests: [`tests/test_assoc_ops.c`](https://github.com/stblake/mathilda/blob/main/tests/test_assoc_ops.c)
- Tests: [`tests/test_association.c`](https://github.com/stblake/mathilda/blob/main/tests/test_association.c)
- Tests: [`tests/test_bigint.c`](https://github.com/stblake/mathilda/blob/main/tests/test_bigint.c)
- Tests: [`tests/test_characteristicpolynomial.c`](https://github.com/stblake/mathilda/blob/main/tests/test_characteristicpolynomial.c)

## Notes & additional examples

### Notes

- `Plus` is `Flat`, `Listable`, `NumericFunction`, `OneIdentity`, `Orderless` and `Protected`. `Plus[]` is `0` and `Plus[x]` is `x`.
- Like terms are collected structurally, by an exact match of the non-numeric part, through a hash table. The result is canonically sorted, and only when it is out of order. A collected term is a flat `Times`, so `2 a b` is `Times[2, a, b]`.
- A sum is not distributed over a product (`c (a + b)` stays as it is), and no algebraic identities are applied: `Sin[x]^2 + Cos[x]^2` needs `Simplify`.
- Integer sums overflow into GMP bignums. An inexact summand makes the whole numeric part inexact (`1. + Pi`), at the lowest precision present.
- Infinities follow Mathematica: a real infinity absorbs finite terms; `Infinity - Infinity` and `ComplexInfinity + ComplexInfinity` give `Indeterminate` with `Infinity::indet`. A directed complex infinity keeps its direction.
- Numeric surfaces. Packed lists add without unpacking, NDArrays add element-wise with scalar broadcasting (mismatched shapes give `NDArray::shape` and the sum stays unevaluated), and `Compile[]` lowers `Plus` for integer, real, complex and array arguments.
- See also: `Total`, `Times`, `Subtract`, `Accumulate`.
