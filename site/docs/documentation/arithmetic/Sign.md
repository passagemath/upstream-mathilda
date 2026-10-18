# Sign

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Sign[x]`**

gives -1, 0, or 1 for real numeric x according to its sign, and z/Abs\[z\] for a nonzero numeric complex z.

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
In[8]:= Sign[-7]
Out[8]= -1
```

Listable

```mathematica
In[9]:= Sign[{-2, 0, 5}]
Out[9]= {-1, 0, 1}
```

The sign of an inexact real is still the exact integer -1

```mathematica
In[10]:= Sign[-2.5]
Out[10]= -1
```

A bignum rational

```mathematica
In[11]:= Sign[-1/10^25]
Out[11]= -1
```

A complex number gives the unit vector z/Abs[z]

```mathematica
In[12]:= Sign[3 + 4 I]
Out[12]= 3/5 + 4/5*I
```

(1 + I)^2 = 2 I points straight up

```mathematica
In[13]:= Sign[(1 + I)^2]
Out[13]= I
```

```mathematica
In[14]:= Sign[2 - 2 I]
Out[14]= (1/2 - 1/2*I) Sqrt[2]
```

MPFR components keep their precision

```mathematica
In[15]:= Sign[N[1 + I, 25]]
Out[15]= 0.70710678118654752440084433 + 0.70710678118654752440084433*I
```

An interval maps monotonically

```mathematica
In[16]:= Sign[Interval[{-1, 2}]]
Out[16]= Interval[{-1, 1}]
```

The narrowing kernel writes an int64 buffer

```mathematica
In[17]:= Sign[NDArray[{-1.5, 0., 2.}]]
Out[17]= NDArray[{-1, 0, 1}]
```

Not yet evaluated: exact real constants other than rationals stay as they are

```mathematica
In[18]:= Sign[Pi - 3]
Out[18]= Sign[-3 + Pi]
```

## Implementation notes

**Algorithm.** `builtin_sign` maps an `Interval` monotonically through `interval_thread_monotone`, so `Sign[Interval[{-1, 2}]]` is `Interval[{-1, 1}]`. For a real number it returns the **exact** integer -1, 0 or 1 whatever the input type: a direct comparison for an int64 or a machine Real (so `Sign[-2.5]` is `-1`, not `-1.`), `mpz_sgn` for a bigint, `mpfr_sgn` for MPFR, and `expr_numeric_sign` for a Rational, including one with bignum parts such as `1/10^25`. For a `Complex[re, im]` with numeric parts it returns the unit vector `z/Abs[z]`, and `0` for `0 + 0 I`. Exact parts build `z · Abs[z]^-1` and let the evaluator fold it, so `Sign[3 + 4 I]` is `3/5 + 4/5 I` and `Sign[2 - 2 I]` is `(1/2 - 1/2 I) Sqrt[2]`. MPFR parts are divided by `mpfr_hypot` at the combined precision. Anything else returns `NULL`.

**Data structures.** The routine allocates no temporaries beyond the MPFR scratch values. On the packed and NDArray surfaces `NDKU_Sign` is a *narrowing* kernel: a real or int64 buffer gives an int64 buffer of -1, 0 and 1 (`-0.` gives 0; a non-finite element abandons the buffer), and a complex buffer gives `z/|z|`. `Sign` is on `pack.c`'s `AWARE` and `INT64_OK` lists. In `Compile[]` the result type is Integer for an Integer or Real argument and Complex for a Complex one, at scalar and array shapes.

**Complexity / limits.** Each case is O(1). `Sign[x]` stays unevaluated for symbolic `x`. A gap: an exact real constant that is not a rational, such as `Sign[Pi]`, `Sign[Sqrt[2]]` or `Sign[Pi - 3]`, also stays unevaluated at present, whereas Mathematica returns `1`. Use `Sign[N[expr]]`, which is still exact.

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## References

**See also:** [Re](../../arithmetic/Re/), [Im](../../arithmetic/Im/), [ReIm](../../arithmetic/ReIm/), [Abs](../../arithmetic/Abs/), [Conjugate](../../arithmetic/Conjugate/), [Arg](../../arithmetic/Arg/), [Real](../../other-advanced/Real/), [Rational](../../arithmetic/Rational/)

- Source: [`src/complex.c`](https://github.com/stblake/mathilda/blob/main/src/complex.c)
- Specification: [`docs/spec/builtins/arithmetic.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/arithmetic.md)
- Tests: [`tests/test_autocompile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_autocompile.c)
- Tests: [`tests/test_bignum_rational_numeric.c`](https://github.com/stblake/mathilda/blob/main/tests/test_bignum_rational_numeric.c)
- Tests: [`tests/test_compile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile.c)
- Tests: [`tests/test_compile_coverage.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile_coverage.c)

## Notes & additional examples

### Notes

- For a real number `Sign` is the exact integer `-1`, `0` or `1`, even for an inexact input (`Sign[-2.5]` is `-1`, as in Mathematica). For a nonzero complex `z` it is `z/Abs[z]`, and `Sign[0]` is `0`.
- Limitation: an exact real constant that is not rational (`Sign[Pi]`, `Sign[Sqrt[2]]`, `Sign[Pi - 3]`) stays unevaluated, where Mathematica gives `1`. `Sign[N[expr]]` still returns an exact integer.
- `Listable`, `NumericFunction`, `Protected`. The packed/NDArray kernel narrows real and int64 data to int64 (`Sign` is on `pack.c`'s `INT64_OK` list), and `Compile[]` gives an Integer result for real arguments and a Complex one for complex arguments.
- See also: `Abs`, `UnitStep`, `Arg`, `Positive`.
