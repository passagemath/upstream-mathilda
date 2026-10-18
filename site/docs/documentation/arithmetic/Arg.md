# Arg

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Arg[z]`**

gives the argument (phase angle in (-Pi, Pi\]) of numeric z; 0 for nonnegative reals, Pi for negative reals.

## Examples (17)

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

### Applications (10)

The principal value lies in (-Pi, Pi]

```mathematica
In[8]:= Arg[-1]
Out[8]= Pi
```

```mathematica
In[9]:= Arg[1 + I]
Out[9]= 1/4 Pi
```

Powers of 1 + I go once round the circle and fold back at Pi

```mathematica
In[10]:= Table[Arg[(1 + I)^k], {k, 0, 8}]
Out[10]= {0, 1/4 Pi, 1/2 Pi, 3/4 Pi, Pi, -3/4 Pi, -1/2 Pi, -1/4 Pi, 0}
```

An exact root of unity gives an exact multiple of Pi

```mathematica
In[11]:= Arg[(-1)^(1/6)]
Out[11]= 1/6 Pi
```

```mathematica
In[12]:= Arg[-(-1)^(5/6)]
Out[12]= -1/6 Pi
```

Read from the exact sign, so it does not underflow

```mathematica
In[13]:= Arg[-1/10^400]
Out[13]= Pi
```

A generic Gaussian gives the two-argument ArcTan

```mathematica
In[14]:= Arg[2 + 3 I]
Out[14]= ArcTan[2, 3]
```

```mathematica
In[15]:= N[Arg[2 + 3 I], 40]
Out[15]= 0.98279372324732906798571061101466601449686
```

An MPFR argument goes through mpfr_atan2 at its own precision

```mathematica
In[16]:= Arg[N[-1 + I, 30]]
Out[16]= 2.356194490192344928846982537459
```

The NDArray kernel

```mathematica
In[17]:= Arg[NDArray[{1. + I, -1., -I}]]
Out[17]= NDArray[{0.785398, 3.14159, -1.5708}]
```

## Implementation notes

**Algorithm.** `builtin_arg` returns the principal argument in `(-Pi, Pi]`, trying the following cases in order. (1) A pure MPFR real gives the exact `0` or `Pi` by sign. (2) An exact real leaf (an Integer, bigint or Rational, bignum parts included) gives `0` or `Pi` from its exact sign, so even `Arg[-1/10^400]` is `Pi`. (3) **Root-of-unity constants.** `arg_turns_of` reads a product of real-rational powers and pure-imaginary Gaussians, such as `(-1)^(1/6)`, `-(-1)^(5/6)` or `I (-1)^(1/4)`. It adds their half-turns mod 2 (`[c<0] + Σ eᵢ[bᵢ<0]`), with denominators capped at 2^40, and `arg_from_turns` reduces the total into `(-Pi, Pi]` as an exact rational multiple of `Pi`. (4) A `Complex[re, im]`: with an MPFR part it uses `mpfr_atan2` at the combined precision. With exact parts it returns exact multiples of `Pi` for the axes and diagonals (`0`, `Pi`, `±Pi/2`, `±Pi/4`, `±3 Pi/4`), and otherwise the symbolic two-argument `ArcTan[re, im]`. With a machine part it uses libm `atan2(im, re)`. (5) A machine Real gives `atan2(0, x)`. Anything else, including symbols, returns `NULL`.

**Data structures.** The routine allocates no structures. A rational turn count, an int64 pair, carries the root-of-unity case. On the packed and NDArray surfaces `Arg` is a projection kernel (`UK_PROJ(Arg, atan2(ai, ar))`) that writes a real buffer. In `Compile[]` it lowers at scalar and array shapes: an integer argument stays integral, and the `Pi` case is handed back to the interpreter.

**Complexity / limits.** Each case is O(1), and the root-of-unity reader is linear in the number of factors. Generic exact Gaussians stay as `ArcTan[re, im]`, which Mathematica would print as `ArcTan[3/2]` for `Arg[2 + 3 I]`; `N` of either gives the same number. A gap: an exact real constant that is not rational, such as `Arg[Pi - 4]`, stays unevaluated (Mathematica gives `Pi`). A product such as `-Sqrt[2]` is covered by the power-product reader.

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## References

**See also:** [Re](../../arithmetic/Re/), [Im](../../arithmetic/Im/), [ReIm](../../arithmetic/ReIm/), [Abs](../../arithmetic/Abs/), [Sign](../../arithmetic/Sign/), [Conjugate](../../arithmetic/Conjugate/), [Real](../../other-advanced/Real/), [Rational](../../arithmetic/Rational/)

- Source: [`src/complex.c`](https://github.com/stblake/mathilda/blob/main/src/complex.c)
- Specification: [`docs/spec/builtins/arithmetic.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/arithmetic.md)
- Tests: [`tests/test_bignum_rational_numeric.c`](https://github.com/stblake/mathilda/blob/main/tests/test_bignum_rational_numeric.c)
- Tests: [`tests/test_compile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile.c)
- Tests: [`tests/test_compile_coverage.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile_coverage.c)
- Tests: [`tests/test_compiledfunction.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compiledfunction.c)

## Notes & additional examples

### Notes

- `Arg[z]` is the principal argument, in `(-Pi, Pi]`. It is `0` for a positive real, `Pi` for a negative real, and `0` at the origin.
- Exact Gaussians on an axis or a diagonal give exact multiples of `Pi`, and so do products of rational powers of negative numbers such as `(-1)^(1/6)`. Any other exact Gaussian gives `ArcTan[re, im]`, which equals the `ArcTan[im/re]` Mathematica prints when `re > 0`.
- Exact real constants that are not rational stay unevaluated: `Arg[Pi - 4]` should be `Pi`. Apply `N` first.
- `Listable`, `NumericFunction`, `Protected`. There is an NDArray/packed kernel (`atan2`, real result), and `Compile[]` lowers `Arg` at scalar and array shapes.
- See also: `Abs`, `Sign`, `ArcTan`, `AbsArg`.
