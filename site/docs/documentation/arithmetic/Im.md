# Im

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Im[z]`**

gives the imaginary part of numeric z, and 0 for real or real-valued (Re/Im/Abs/Arg) arguments.

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
In[8]:= Im[3 + 4 I]
Out[8]= 4

In[9]:= Im[7]
Out[9]= 0
```

Sqrt[-4] evaluates to 2 I

```mathematica
In[10]:= Im[Sqrt[-4]]
Out[10]= 2
```

The principal branch: Log[-1] = I Pi

```mathematica
In[11]:= Im[Log[-1]]
Out[11]= Pi
```

Listable

```mathematica
In[12]:= Im[{1 + 2 I, 3 - 4 I, 5}]
Out[12]= {2, -4, 0}
```

A real-valued head gives 0

```mathematica
In[13]:= Im[Abs[z]]
Out[13]= 0
```

No closed form, so it stays symbolic

```mathematica
In[14]:= Im[Gamma[1 + I]]
Out[14]= Im[Gamma[1 + I]]
```

But it evaluates numerically

```mathematica
In[15]:= N[Im[Gamma[1 + I]], 40]
Out[15]= -0.15494982830181068512495513048388660519589
```

The narrowing kernel gives exact zeros for integer data

```mathematica
In[16]:= Im[NDArray[{1, 2}, DataType -> "int64"]]
Out[16]= NDArray[{0, 0}]
```

## Implementation notes

**Algorithm.** `builtin_im` mirrors `builtin_re`. It returns `0` for a call to a head that is real by construction (`Re`, `Im`, `Abs`, `Arg`), for a real `AlgebraicNumber` or `Root` (FLINT `qqbar` reality test), and for any real number leaf (an Integer, bigint, machine Real, MPFR, or Rational with int64 or bignum parts). A `Complex[re, im]` literal gives `im`. Otherwise `complex_decompose`, a recursive Plus/Times walk that multiplies Complex factors through, supplies the imaginary part when both decomposed parts are concretely numeric (`is_numeric_real`). This is how `Im[Sqrt[-4]]` is `2` and `Im[Log[-1]]` is `Pi`: those arguments have already evaluated to `2 I` and `I Pi`. Anything else returns `NULL`, including special-function values with no closed form such as `Im[Gamma[1 + I]]`, whose `N` still evaluates.

**Data structures.** The routine works on fresh `re`/`im` trees from the decomposition. On the packed and NDArray surfaces the kernel `NDKU_ImInt` narrows: a real or int64 buffer gives an int64 buffer of exact zeros (matching the List path's exact `0`), and a complex buffer gives the real buffer of imaginary parts. In `Compile[]` `Im` lowers at scalar and array shapes, and it is `0` over integers.

**Complexity / limits.** The leaf cases are O(1), and the decomposition is linear in the size of the expression. Symbols are not assumed real (`Im[x]` stays unevaluated). A gap: an exact real constant with no Complex literal, such as `Im[Pi - 4]`, stays unevaluated at present (Mathematica gives `0`).

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## References

**See also:** [Re](../../arithmetic/Re/), [ReIm](../../arithmetic/ReIm/), [Abs](../../arithmetic/Abs/), [Sign](../../arithmetic/Sign/), [Conjugate](../../arithmetic/Conjugate/), [Arg](../../arithmetic/Arg/), [Real](../../other-advanced/Real/), [Rational](../../arithmetic/Rational/)

- Source: [`src/complex.c`](https://github.com/stblake/mathilda/blob/main/src/complex.c)
- Specification: [`docs/spec/builtins/arithmetic.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/arithmetic.md)
- Tests: [`tests/test_accuracygoal.c`](https://github.com/stblake/mathilda/blob/main/tests/test_accuracygoal.c)
- Tests: [`tests/test_airyai.c`](https://github.com/stblake/mathilda/blob/main/tests/test_airyai.c)
- Tests: [`tests/test_airybi.c`](https://github.com/stblake/mathilda/blob/main/tests/test_airybi.c)
- Tests: [`tests/test_algebraicnumber.c`](https://github.com/stblake/mathilda/blob/main/tests/test_algebraicnumber.c)

## Notes & additional examples

### Notes

- `Im[z]` is the imaginary part of a numeric `z`. It is `0` for real numbers, for real `Root` and `AlgebraicNumber` objects, and for calls to the real-valued heads `Re`, `Im`, `Abs`, `Arg`.
- Symbols are not assumed real (`Im[x]` stays unevaluated). An exact real constant with no complex part written out, such as `Im[Pi - 4]`, also stays unevaluated at present (Mathematica gives `0`).
- `Listable`, `NumericFunction`, `Protected`. The packed/NDArray kernel gives an int64 zero buffer for real or integer data and the imaginary parts for complex data. `Compile[]` lowers `Im` at scalar and array shapes.
- See also: `Re`, `ReIm`, `Conjugate`, `ComplexExpand`.
