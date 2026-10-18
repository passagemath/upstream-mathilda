# Re

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Re[z]`**

gives the real part of numeric z; Re\[Re\[z\]\], Re\[Im\[z\]\], Re\[Abs\[z\]\], Re\[Arg\[z\]\] fold since those heads are real-valued.

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
In[8]:= Re[3 + 4 I]
Out[8]= 3
```

Listable

```mathematica
In[9]:= Re[{1 + I, 2 - 3 I, 5}]
Out[9]= {1, 2, 5}
```

(1 + I)^10 = 32 I is purely imaginary

```mathematica
In[10]:= Re[(1 + I)^10]
Out[10]= 0
```

The exact real part of a quotient

```mathematica
In[11]:= Re[1/(2 + 3 I)]
Out[11]= 2/13
```

A bignum is its own real part

```mathematica
In[12]:= Re[10^30]
Out[12]= 1000000000000000000000000000000
```

Re, Im, Abs and Arg are real-valued by construction

```mathematica
In[13]:= Re[Re[z]]
Out[13]= Re[z]
```

A real Root object is returned unchanged

```mathematica
In[14]:= Re[Root[#^3 - 2 &, 1]]
Out[14]= Root[#1^3 - 2 &, 1]
```

The NDArray projection kernel

```mathematica
In[15]:= Re[NDArray[{1. + 2. I, 3. - I}]]
Out[15]= NDArray[{1.0, 3.0}]
```

Symbols are not assumed real: use ComplexExpand for that

```mathematica
In[16]:= Re[x + I y]
Out[16]= Re[x + I y]
```

## Implementation notes

**Algorithm.** `builtin_re` tries the following cases in order. A call to a head that is real by construction (`Re`, `Im`, `Abs`, `Arg`) is returned unchanged, so `Re[Re[z]]` is `Re[z]`. A real `AlgebraicNumber` or `Root` (tested with FLINT `qqbar`, `flint_qqbar_is_real`) is returned as itself; a non-real one stays wrapped, as in Mathematica. A `Complex[re, im]` literal gives `re`. A real number leaf (`is_real_number_leaf`: an Integer, bigint, machine Real, MPFR, or Rational with int64 or bignum parts) gives itself. Otherwise `complex_decompose` walks the Plus/Times structure, multiplying Complex factors through. When both parts it finds are concretely numeric (`is_numeric_real`, judged by machine `N`), the real part is returned; this is what reduces `Re[1/(2 + 3 I)]` to `2/13`. Anything else returns `NULL` and stays unevaluated.

**Data structures.** The decomposition builds and evaluates fresh `re`/`im` trees. No other state is kept. On the packed and NDArray surfaces `Re` is a projection kernel (`UK_PROJ(Re, ar)`): a complex buffer gives a real buffer, and real data passes through. In `Compile[]` it lowers at scalar and array shapes, and over integers it is the identity.

**Complexity / limits.** The leaf cases are O(1). The decomposition is linear in the size of the expression, plus one machine-`N` reality check per part. Symbols are not assumed real, so `Re[x + I y]` stays unevaluated (use `ComplexExpand` for that convention). A gap: an exact real constant with no Complex literal inside, such as `Re[Pi]`, is not yet recognised and stays unevaluated (Mathematica gives `Pi`).

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## References

**See also:** [Im](../../arithmetic/Im/), [ReIm](../../arithmetic/ReIm/), [Abs](../../arithmetic/Abs/), [Sign](../../arithmetic/Sign/), [Conjugate](../../arithmetic/Conjugate/), [Arg](../../arithmetic/Arg/), [Real](../../other-advanced/Real/), [Rational](../../arithmetic/Rational/)

- Source: [`src/complex.c`](https://github.com/stblake/mathilda/blob/main/src/complex.c)
- Specification: [`docs/spec/builtins/arithmetic.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/arithmetic.md)
- Tests: [`tests/test_accuracygoal.c`](https://github.com/stblake/mathilda/blob/main/tests/test_accuracygoal.c)
- Tests: [`tests/test_airyai.c`](https://github.com/stblake/mathilda/blob/main/tests/test_airyai.c)
- Tests: [`tests/test_airybi.c`](https://github.com/stblake/mathilda/blob/main/tests/test_airybi.c)
- Tests: [`tests/test_algebraicnumber.c`](https://github.com/stblake/mathilda/blob/main/tests/test_algebraicnumber.c)

## Notes & additional examples

### Notes

- `Re[z]` is the real part of a numeric `z`. A real number comes back unchanged.
- Symbols are not assumed real, so `Re[x + I y]` stays unevaluated; `ComplexExpand[Re[x + I y]]` gives `x`. A non-real `Root` or `AlgebraicNumber` also stays wrapped, as in Mathematica.
- Limitation: an exact real constant with no complex part written out, such as `Re[Pi]`, is not yet recognised and stays unevaluated (Mathematica gives `Pi`).
- `Listable`, `NumericFunction`, `Protected`. There is a packed/NDArray projection kernel, and `Compile[]` lowers `Re` at scalar and array shapes.
- See also: `Im`, `ReIm`, `Conjugate`, `ComplexExpand`.
