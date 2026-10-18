# Abs

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Abs[z]`**

gives the absolute value (modulus) of numeric z, Sqrt\[Re\[z\]^2 + Im\[z\]^2\] for complex z.

## Examples (18)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (7)

```mathematica
In[1]:= Re[1/(2 + 3 I)]
Out[1]= 2/13

In[2]:= Im[Log[-1]]
Out[2]= Pi

In[3]:= ReIm[(2 + I)^3]
Out[3]= {2, 11}

In[4]:= Abs[-9223372036854775807 - 1]
Out[4]= 9223372036854775808

In[5]:= Sign[3 + 4 I]
Out[5]= 3/5 + 4/5*I

In[6]:= Conjugate[(2 + I)/(1 - 3 I)]
Out[6]= -1/10 - 7/10*I

In[7]:= Arg[(-1)^(1/6)]
Out[7]= 1/6 Pi
```

### Applications (11)

```mathematica
In[8]:= Abs[-5]
Out[8]= 5
```

The modulus of a Gaussian integer, exactly

```mathematica
In[9]:= Abs[3 + 4 I]
Out[9]= 5
```

Listable

```mathematica
In[10]:= Abs[{-1, 2, -3/4}]
Out[10]= {1, 2, 3/4}
```

|1 + I|^10 = Sqrt[2]^10

```mathematica
In[11]:= Abs[(1 + I)^10]
Out[11]= 32
```

The modulus of an exact complex expression

```mathematica
In[12]:= Abs[Sqrt[2] + Sqrt[3] I]
Out[12]= Sqrt[5]
```

|-2^63| does not fit an int64, so it becomes a bignum

```mathematica
In[13]:= Abs[-9223372036854775807 - 1]
Out[13]= 9223372036854775808
```

A negative real Root stays an exact algebraic number

```mathematica
In[14]:= Abs[Root[#^3 - 2 # + 1 &, 1]]
Out[14]= -Root[#1^3 - 2 #1 + 1 &, 1]
```

Interval arithmetic

```mathematica
In[15]:= Abs[Interval[{-2, 1}]]
Out[15]= Interval[{0, 2}]
```

An MPFR complex goes through mpfr_hypot at 40 digits

```mathematica
In[16]:= Abs[N[Gamma[1/3 + 2 I], 40]]
Out[16]= 0.096659594257326641410227978598687975893249
```

The int64 kernel keeps integers exact

```mathematica
In[17]:= Abs[NDArray[{-1, 2, -3}, DataType -> "int64"]]
Out[17]= NDArray[{1, 2, 3}]
```

Not yet evaluated: exact real constants other than rationals stay as they are

```mathematica
In[18]:= Abs[Sqrt[2]]
Out[18]= Abs[Sqrt[2]]
```

## Performance

Against other systems, from the benchmark suite (same input, results cross-checked for agreement):

| case | Mathilda | Wolfram | Python |
|---|---:|---:|---:|
| NI 50-digit Gaussian | 6.87 s | 2.11 s | 1.1 s |
| NI 2-D ridge | 2.54 s | 43.5 s | 0.685 s |
| NI oscillatory k=40 | 0.515 s | 4.05 s | 0.002 s |
| NI oscillatory k=200 | 0.477 s | 21.5 s | 0.002 s |
| NI oscillatory k=1000 | 0.389 s | 140 s | 0.002 s |
| NI oscillatory k=1001 nonzero | 0.382 s | 1.03 s | 0.626 s |

## Implementation notes

**Algorithm.** `builtin_abs` tries the following cases in order. (1) An `Interval` goes to `interval_apply_function` (`iv_abs`), so `Abs[Interval[{-2, 1}]]` is `Interval[{0, 2}]`. (2) A real `AlgebraicNumber` or `Root` object is compared with 0 through FLINT `qqbar` (`flint_qqbar_compare`) and returned as itself or its negation, so it stays an exact algebraic object. (3) A `Complex[re, im]` literal, or an expression that `complex_decompose` splits into concretely numeric real and imaginary parts (a recursive Plus/Times walk that propagates Complex factors), becomes the modulus. When either part is MPFR this goes through `mpfr_hypot` at the combined precision, which is also stable when the two magnitudes are very different. Otherwise it builds `(re^2 + im^2)^(1/2)` exactly, so `Abs[3 + 4 I]` is `5` and `Abs[Sqrt[2] + Sqrt[3] I]` is `Sqrt[5]`. (4) Real leaves: `mpz_abs` for a bigint, a sign flip for an int64 (`INT64_MIN` is promoted to the bignum 2^63 instead of wrapping), `fabs` for a machine Real, `mpfr_abs` at the input precision, and `|n|/d` for a Rational, including one with bignum parts. Anything else returns `NULL`, so the call stays unevaluated.

**Data structures.** The routine works directly on `Expr` leaves. MPFR scratch values are sized by `numeric_combined_bits`. The packed and NDArray surfaces use `NDKU_AbsInt` (`ndkernels.c`): an exact int64 arm (it abandons the whole buffer on `INT64_MIN`, so the List path answers with a bignum), and a `hypot` arm for real and complex buffers that writes a real buffer. `Abs` is on `pack.c`'s `AWARE` list. In `Compile[]`, `Abs` of an Integer stays an Integer, of a Real a Real, and of a Complex a Real, at scalar and array shapes.

**Complexity / limits.** Each case costs O(1) numeric work, except the `complex_decompose` walk, which is linear in the size of the expression. Symbolic arguments stay unevaluated (`Abs[x]`), as in Mathematica. A gap: a *real* exact constant that is neither a rational nor a `Complex` literal, such as `Abs[Sqrt[2]]`, `Abs[-Pi]` or `Abs[Pi - 4]`, also stays unevaluated at present, whereas Mathematica returns `Sqrt[2]`, `Pi` and `4 - Pi`. Apply `N` first for a numeric value.

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## References

**See also:** [Re](../../arithmetic/Re/), [Im](../../arithmetic/Im/), [ReIm](../../arithmetic/ReIm/), [Sign](../../arithmetic/Sign/), [Conjugate](../../arithmetic/Conjugate/), [Arg](../../arithmetic/Arg/), [Real](../../other-advanced/Real/), [Rational](../../arithmetic/Rational/)

- Source: [`src/complex.c`](https://github.com/stblake/mathilda/blob/main/src/complex.c)
- Specification: [`docs/spec/builtins/arithmetic.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/arithmetic.md)
- Tests: [`tests/test_accuracygoal.c`](https://github.com/stblake/mathilda/blob/main/tests/test_accuracygoal.c)
- Tests: [`tests/test_airybi.c`](https://github.com/stblake/mathilda/blob/main/tests/test_airybi.c)
- Tests: [`tests/test_algebraicnumber.c`](https://github.com/stblake/mathilda/blob/main/tests/test_algebraicnumber.c)
- Tests: [`tests/test_assoc_ops.c`](https://github.com/stblake/mathilda/blob/main/tests/test_assoc_ops.c)

## Notes & additional examples

### Notes

- `Abs[z]` is the absolute value of a real number and the modulus `Sqrt[Re[z]^2 + Im[z]^2]` of a complex one. Exact inputs give exact results.
- Limitation: an exact real constant that is not rational, such as `Sqrt[2]`, `-Pi` or `Pi - 4`, stays unevaluated (Mathematica gives `Sqrt[2]`, `Pi`, `4 - Pi`). `Abs[N[...]]` is numeric.
- Symbolic arguments stay unevaluated (`Abs[x]`). There is no `ComplexExpand`-style assumption that symbols are real.
- `Listable`, `NumericFunction`, `Protected`. There is a packed/NDArray kernel with an exact int64 arm (and a `hypot` arm for real and complex data), and `Compile[]` lowers it at scalar and array shapes (`Abs` of a Complex is Real).
- See also: `Sign`, `Arg`, `Norm`, `RealAbs`.
