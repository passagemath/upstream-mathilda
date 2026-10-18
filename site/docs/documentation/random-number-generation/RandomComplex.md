# RandomComplex

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`RandomComplex[]`**

gives a pseudorandom complex number with real and imaginary parts in the range 0 to 1.

**`RandomComplex[{zmin, zmax}]`**

gives a pseudorandom complex number in the rectangle with corners given by the complex numbers zmin and zmax.

**`RandomComplex[zmax]`**

gives a pseudorandom complex number in the rectangle whose corners are the origin and zmax.

**`RandomComplex[range, n]`**

gives a list of n pseudorandom complex numbers.

**`RandomComplex[range, {n1, n2, ...}]`**

gives an n1 x n2 x ... array of pseudorandom complex numbers.

**`RandomComplex[spec, WorkingPrecision -> n]`**

yields complex numbers whose real and imaginary parts have n digits of precision. Leading or trailing digits of the generated parts can be 0. n may be MachinePrecision (the default) or a positive number of decimal digits.

## Examples (17)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (6)

```mathematica
In[1]:= SeedRandom[42]; Head[RandomComplex[]]
Out[1]= Complex

In[2]:= SeedRandom[42]; RandomComplex[2 + 3 I]
Out[2]= 1.62861 + 0.956463*I

In[3]:= SeedRandom[42]; Length[RandomComplex[{0, 1 + I}, 5]]
Out[3]= 5

In[4]:= SeedRandom[42]; Dimensions[RandomComplex[{0, 1 + I}, {3, 4}]]
Out[4]= {3, 4}

In[5]:= RandomComplex[x]
Out[5]= RandomComplex[x]

In[6]:= SeedRandom[7]; pts = RandomComplex[{-1 - 2 I, 3 + I}, 10^5]; {PackedArrayQ[pts], MinMax[Re[pts]], MinMax[Im[pts]]}
Out[6]= {True, {-0.999983, 2.99996}, {-1.99997, 0.999965}}
```

### Options (1)

```mathematica
In[7]:= SeedRandom[42]; z = RandomComplex[1 + I, WorkingPrecision -> 30]; {Precision[Re[z]], Precision[Im[z]]}
Out[7]= {30.103, 30.103}
```

### Applications (10)

Uniform in the unit square with corners 0 and 1 + I

```mathematica
In[8]:= SeedRandom[1]; RandomComplex[]
Out[8]= 0.811612 + 0.747105*I
```

A bare corner z gives the rectangle from 0 to z; same draw, scaled by 2

```mathematica
In[9]:= SeedRandom[1]; RandomComplex[2 + 2 I]
Out[9]= 1.62322 + 1.49421*I
```

A 2x2 array of points in the square centred at the origin

```mathematica
In[10]:= SeedRandom[1]; RandomComplex[{-1 - I, 1 + I}, {2, 2}]
Out[10]= {{0.623224 + 0.494209*I, -0.799698 + 0.492434*I}, {-0.630643 + 0.180958*I, 0.973748 + 0.0468337*I}}
```

Symbolic corners are reduced to numbers first

```mathematica
In[11]:= SeedRandom[5]; RandomComplex[{-Pi - I, Pi + I}]
Out[11]= -1.30676 + 0.222879*I
```

A packed complex64 buffer whose parts stay inside the rectangle

```mathematica
In[12]:= SeedRandom[7]; pts = RandomComplex[{-1 - 2 I, 3 + I}, 10^5]; {PackedArrayQ[pts], MinMax[Re[pts]], MinMax[Im[pts]]}
Out[12]= {True, {-0.999983, 2.99996}, {-1.99997, 0.999965}}
```

The sample mean is close to the centre 1 - I/2

```mathematica
In[13]:= SeedRandom[7]; Mean[RandomComplex[{-1 - 2 I, 3 + I}, 10^5]]
Out[13]= 0.999867 - 0.50121*I
```

About Pi/4 = 0.785 of the unit square lies inside the unit circle

```mathematica
In[14]:= SeedRandom[11]; Count[RandomComplex[1 + I, 10^4], w_ /; Abs[w] < 1] / 10.^4
Out[14]= 0.7883
```

The packed fill draws in the same order as repeated scalar calls

```mathematica
In[15]:= SeedRandom[3]; a = RandomComplex[{0, 1 + I}, {2, 2}]; SeedRandom[3]; b = Table[RandomComplex[{0, 1 + I}], 2, 2]; a == b
Out[15]= True
```

Above machine precision the parts are MPFR numbers

```mathematica
In[16]:= SeedRandom[4]; Precision /@ Re[RandomComplex[{0, 1 + I}, 3, WorkingPrecision -> 20]]
Out[16]= {20.169, 20.169, 20.169}
```

A non-numeric range stays unevaluated

```mathematica
In[17]:= RandomComplex[x]
Out[17]= RandomComplex[x]
```

## Implementation notes

**Algorithm.** `builtin_randomcomplex` first strips a trailing `WorkingPrecision -> d` option (`parse_working_precision_opt`). If `d` is above machine precision, the call goes to `randomcomplex_mpfr` (in `USE_MPFR` builds); otherwise it goes to `randomcomplex_machine`. `parse_complex_range` turns the range into the corners of an axis-aligned rectangle. A bare `z` means `[0, Re z] × [0, Im z]`; `{z1, z2}` means `[Re z1, Re z2] × [Im z1, Im z2]`. `expr_to_complex_parts` reads each corner from a `Complex` atom or a real literal. For a symbolic corner such as `-Pi - I` or `Exp[I Pi/4]` it falls back to `numericalize` at machine precision. A non-numeric range leaves the call unevaluated. A single point is `random_complex_range`: two `random_uniform_01()` draws (real part first, then imaginary), mapped affinely onto the rectangle. `random_uniform_01` takes the top 53 bits of one xoshiro256++ output (`xs_next`, seeded by splitmix64), which is the stream `RandomReal` uses. Since v0.402, `RandomComplex[range, n]` and `RandomComplex[range, {n1, …}]` with a truly complex range (imaginary extent not identically zero) allocate a packed `NDT_COMPLEX64` buffer with `ndbuild_open` and fill it row-major. The draws are in the same re-then-im order, so a seeded stream gives the same values as repeated scalar calls. A range with zero imaginary extent keeps the boxed path (`random_complex_array`), where `make_complex` folds each `im = 0` point to a Real. The MPFR path parses the corners with `random_get_approx_mpfr` at the target bit count. It draws each part with `mpfr_urandomb` from the GMP Mersenne-Twister state, which `SeedRandom` seeds together with the xoshiro stream. Its arrays are boxed `List`s (`random_complex_array_mpfr`).

**Data structures.** The machine scalar is a `Complex[Real, Real]`. Machine arrays are a packed `complex64` NDArray buffer (interleaved re/im doubles), which `Sin`, `Exp`, `HermiteH`, reductions and the other buffer kernels consume directly. The degenerate real-valued case and every MPFR result are nested boxed `List`s of `Complex` (or `Real`) leaves. The generator state is four `uint64_t` words for xoshiro256++ plus one `gmp_randstate_t` for the MPFR draws.

**Complexity / limits.** The cost is `O(N)` for `N` points: two xoshiro steps, about 2.4 ns each, plus one buffer store per point on the packed path. The boxed and MPFR paths add one allocation per element. Dimension specs must be non-negative machine integers, else the call is unevaluated. The distribution is uniform on the rectangle only; there is no disk or other region form. `RandomComplex` has no `Compile[]` lowering: `CompileDiagnostics` names it as the unlowered head. The packed producer is the array surface; there is no element-wise ND kernel because the head takes no numeric array argument. Attributes: `Protected`.

- `Protected`.
- `RandomComplex[{zmin, zmax}]` chooses complex numbers uniformly in the rectangle with corners at `zmin` and `zmax`.
- RandomComplex gives a different sequence of pseudorandom complex numbers whenever you run Mathilda. You can start with a particular seed using SeedRandom.
- Uses 53 bits of randomness per component for full double-precision mantissa coverage.
- Accepts integer, real, rational, and complex range arguments, as well as symbolic-but-numeric corners that `N[]` can reduce to a number, e.g. `RandomComplex[{-Pi - I, Pi + I}]` (a `Plus` that only collapses to a `Complex` once `Pi` is numeric). When the range has no imaginary component, the result simplifies to a real.
- `WorkingPrecision -> n` accepts `MachinePrecision` (the default) or a positive number of decimal digits. Digit counts above MachinePrecision route generation through MPFR, so the real and imaginary parts are MPFR atoms at the requested precision.
- Machine-precision `RandomComplex[range, n]` / `RandomComplex[range, {n1, ...}]` with a genuinely complex range fill a packed `complex64` buffer directly (v0.402), in the same re-then-im draw order as scalar calls, so a seeded stream is unchanged. A zero-imaginary range returns an unpacked List of Reals; MPFR results are boxed.

**Attributes:** `Protected`.

## References

**See also:** [Plus](../../arithmetic/Plus/), [Complex](../../arithmetic/Complex/), [Pi](../../mathematical-constants/Pi/)

- Source: [`src/random.c`](https://github.com/stblake/mathilda/blob/main/src/random.c)
- Specification: [`docs/spec/builtins/random-number-generation.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/random-number-generation.md)
- Tests: [`tests/test_erf.c`](https://github.com/stblake/mathilda/blob/main/tests/test_erf.c)
- Tests: [`tests/test_erfc.c`](https://github.com/stblake/mathilda/blob/main/tests/test_erfc.c)
- Tests: [`tests/test_erfi.c`](https://github.com/stblake/mathilda/blob/main/tests/test_erfi.c)
- Tests: [`tests/test_fresnelc.c`](https://github.com/stblake/mathilda/blob/main/tests/test_fresnelc.c)

## Notes & additional examples

### Notes

- The range is an **axis-aligned rectangle**, not a disk. The real and imaginary parts are two independent uniform draws, real part first. A bare `z` means the rectangle with corners `0` and `z`; `{z1, z2}` means the rectangle spanned by `z1` and `z2`. Corners can be any numeric expression, for example `-Pi - I`.
- Machine draws come from the same xoshiro256++ stream that `RandomReal` uses (53 bits per part), so `SeedRandom` makes them reproducible. `WorkingPrecision -> d` above machine precision draws each part with `mpfr_urandomb` from the GMP generator, which is seeded together with the machine stream.
- `RandomComplex[range, n]` and `RandomComplex[range, {n1, ...}]` with a truly complex range fill a packed `complex64` buffer directly, so `Sin`, `Exp`, `HermiteH` and the reductions run on the buffer. A range with zero imaginary extent (for example `{0, 5}`) folds to Reals and returns an ordinary unpacked List, even at large sizes.
- `RandomComplex` has no `Compile[]` lowering: `CompileDiagnostics` reports it as the unlowered subexpression. Draw the array outside `Compile` and pass it in.
- The precision of a `WorkingPrecision -> d` draw is reported from the MPFR bit count, so it is slightly above `d` (`20.169` for `d = 20`).
- Attributes: `Protected`. See also: `RandomReal`, `RandomInteger`, `RandomVariate`, `SeedRandom`.
