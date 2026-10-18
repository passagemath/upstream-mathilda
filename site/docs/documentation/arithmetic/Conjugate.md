# Conjugate

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Conjugate[z]`**

gives the complex conjugate Re\[z\] - I Im\[z\] of numeric z; real and real-valued (Re/Im/Abs/Arg) arguments are returned unchanged.

## Examples (16)

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

### Applications (9)

```mathematica
In[8]:= Conjugate[3 + 4 I]
Out[8]= 3 - 4*I
```

Listable; real numbers are fixed

```mathematica
In[9]:= Conjugate[{1 + I, 2 - 3 I, 5}]
Out[9]= {1 - I, 2 + 3*I, 5}
```

The conjugate of an exact quotient

```mathematica
In[10]:= Conjugate[(2 + I)/(1 - 3 I)]
Out[10]= -1/10 - 7/10*I
```

The squared modulus

```mathematica
In[11]:= z Conjugate[z] /. z -> 3 + 4 I
Out[11]= 25
```

An involution, even for symbolic z

```mathematica
In[12]:= Conjugate[Conjugate[z]]
Out[12]= z
```

An exact real constant is its own conjugate

```mathematica
In[13]:= Conjugate[Sqrt[2] Log[3]]
Out[13]= Log[3] Sqrt[2]
```

```mathematica
In[14]:= Conjugate[Sqrt[-2]]
Out[14]= -I Sqrt[2]
```

The NDArray kernel

```mathematica
In[15]:= Conjugate[NDArray[{1. + 2. I, 3. - I}]]
Out[15]= NDArray[{1.0 - 2.0*I, 3.0 + 1.0*I}]
```

Symbols are not assumed real

```mathematica
In[16]:= Conjugate[z]
Out[16]= Conjugate[z]
```

## Implementation notes

**Algorithm.** `builtin_conjugate` first checks its arity. Any argument count other than one emits `Conjugate::argx` (through `mth_message`) and leaves the call unevaluated. Then, in order: `Conjugate[Conjugate[z]]` gives `z` (it is an involution). A call to a head that is real by construction (`Re`, `Im`, `Abs`, `Arg`) is a fixed point. A `Complex[re, im]` literal gives `Complex[re, -im]`. An Integer, machine Real or Rational is returned unchanged. So is any expression that `is_numeric_real` accepts, meaning it reduces to a real number under machine `N`; this covers `Sqrt[2] Log[3]` and `Pi - 4`. Finally, an expression that `complex_decompose` splits into numeric real and imaginary parts is rebuilt as `re - I im` and evaluated, which is how `Conjugate[(2 + I)/(1 - 3 I)]` becomes `-1/10 - 7/10 I`. Anything else returns `NULL`.

**Data structures.** The routine allocates no state beyond the decomposition's fresh `re`/`im` trees. On the packed and NDArray surfaces `Conjugate` is a closed kernel (`UK_CLOSED(Conjugate, x, conj(z))`): the identity on real buffers and `conj` on complex ones. In `Compile[]` its result type is its argument's type, at scalar and array shapes.

**Complexity / limits.** The leaf cases are O(1), and the reality test costs one machine `N`. Symbols are not assumed real, so `Conjugate[z]` stays unevaluated, as in Mathematica. Note that the reality test for non-literal expressions is numerical, at machine precision.

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## References

**See also:** [Re](../../arithmetic/Re/), [Im](../../arithmetic/Im/), [ReIm](../../arithmetic/ReIm/), [Abs](../../arithmetic/Abs/), [Sign](../../arithmetic/Sign/), [Arg](../../arithmetic/Arg/), [Real](../../other-advanced/Real/), [Rational](../../arithmetic/Rational/)

- Source: [`src/complex.c`](https://github.com/stblake/mathilda/blob/main/src/complex.c)
- Specification: [`docs/spec/builtins/arithmetic.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/arithmetic.md)
- Tests: [`tests/test_bignum_rational_numeric.c`](https://github.com/stblake/mathilda/blob/main/tests/test_bignum_rational_numeric.c)
- Tests: [`tests/test_blas.c`](https://github.com/stblake/mathilda/blob/main/tests/test_blas.c)
- Tests: [`tests/test_compile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile.c)
- Tests: [`tests/test_complexexpand.c`](https://github.com/stblake/mathilda/blob/main/tests/test_complexexpand.c)

## Notes & additional examples

### Notes

- `Conjugate[z]` is `Re[z] - I Im[z]`. Real numbers, exact real constants (tested by evaluating them numerically at machine precision) and calls to the real-valued heads `Re`, `Im`, `Abs`, `Arg` are fixed points, and `Conjugate[Conjugate[z]]` is `z`.
- Symbols are not assumed real, so `Conjugate[z]` stays unevaluated. Calling it with any number of arguments other than one emits `Conjugate::argx`.
- `Listable`, `NumericFunction`, `Protected`. There is a packed/NDArray kernel (`conj`), and `Compile[]` lowers `Conjugate` at scalar and array shapes.
- See also: `Re`, `Im`, `Abs`, `ConjugateTranspose`.
