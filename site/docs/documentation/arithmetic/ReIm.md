# ReIm

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`ReIm[z]`**

gives {Re\[z\], Im\[z\]}, the real and imaginary parts of numeric z as a list; real-valued arguments give {z, 0}.

## Examples (14)

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

### Applications (7)

```mathematica
In[8]:= ReIm[3 + 4 I]
Out[8]= {3, 4}
```

(2 + I)^3 = 2 + 11 I

```mathematica
In[9]:= ReIm[(2 + I)^3]
Out[9]= {2, 11}
```

A complex quotient reduces to Gaussian integers

```mathematica
In[10]:= ReIm[(3 + 4 I)/(1 - 2 I)]
Out[10]= {-1, 2}
```

Both parts are 1/Sqrt[2], to 20 digits

```mathematica
In[11]:= ReIm[N[E^(I Pi/4), 20]]
Out[11]= {0.707106781186547524402, 0.707106781186547524402}
```

Listable: one pair per element

```mathematica
In[12]:= ReIm[{1 + 2 I, 3.5}]
Out[12]= {{1, 2}, {3.5, 0}}
```

A real-valued head gives {f[z], 0}

```mathematica
In[13]:= ReIm[Abs[z]]
Out[13]= {Abs[z], 0}
```

Symbolic z stays unevaluated

```mathematica
In[14]:= ReIm[z]
Out[14]= ReIm[z]
```

## Implementation notes

**Algorithm.** `builtin_reim` returns `{Re[z], Im[z]}` as a two-element `List` and follows the same order of cases as `builtin_re` and `builtin_im`. A call to a head that is real by construction (`Re`, `Im`, `Abs`, `Arg`) gives `{f[z], 0}`. A `Complex[re, im]` literal gives `{re, im}`. A real number leaf (an Integer, bigint, machine Real, MPFR, or Rational with int64 or bignum parts) gives `{x, 0}`. Otherwise the shared `complex_decompose` Plus/Times walk supplies both parts when they are concretely numeric. Anything else returns `NULL`. Unlike `Re` and `Im`, `ReIm` has no special case for real `AlgebraicNumber` or `Root` objects.

**Data structures.** The result is a fresh `List` of the two parts. Because `ReIm` is `Listable`, a list argument gives a list of pairs.

**Complexity / limits.** The cost is the same as one `Re` plus one `Im`, and the decomposition is done only once. `ReIm` has no packed or NDArray kernel and no `Compile[]` lowering. `ReIm[NDArray[...]]` stays unevaluated, a packed complex list is unpacked before the call, and a compiled body containing `ReIm` falls back to the interpreter. The head is recorded in the fast-path sweep's `OFF_BUFFER` backlog (`tools/nd_fastpath_sweep.py`). Use `{Re[v], Im[v]}` inside `Compile[]` or on large arrays. Symbolic `z` stays unevaluated, as do exact real constants such as `ReIm[Pi]`.

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## References

**See also:** [Re](../../arithmetic/Re/), [Im](../../arithmetic/Im/), [Abs](../../arithmetic/Abs/), [Sign](../../arithmetic/Sign/), [Conjugate](../../arithmetic/Conjugate/), [Arg](../../arithmetic/Arg/), [Real](../../other-advanced/Real/), [Rational](../../arithmetic/Rational/)

- Source: [`src/complex.c`](https://github.com/stblake/mathilda/blob/main/src/complex.c)
- Specification: [`docs/spec/builtins/arithmetic.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/arithmetic.md)
- Tests: [`tests/test_bignum_rational_numeric.c`](https://github.com/stblake/mathilda/blob/main/tests/test_bignum_rational_numeric.c)
- Tests: [`tests/test_complexexpand.c`](https://github.com/stblake/mathilda/blob/main/tests/test_complexexpand.c)
- Tests: [`tests/test_core.c`](https://github.com/stblake/mathilda/blob/main/tests/test_core.c)

## Notes & additional examples

### Notes

- `ReIm[z]` is `{Re[z], Im[z]}`, and a real argument gives `{z, 0}`.
- Unlike `Re` and `Im`, `ReIm` has no packed/NDArray kernel and no `Compile[]` lowering. `ReIm[NDArray[...]]` stays unevaluated, and a compiled body that uses it runs in the interpreter (a recorded backlog item). Use `{Re[v], Im[v]}` on large arrays or inside `Compile[]`.
- `Listable`, `NumericFunction`, `Protected`.
- See also: `Re`, `Im`, `AbsArg`, `Conjugate`.
