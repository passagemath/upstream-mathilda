# RealDigits

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`RealDigits[x]`**

gives a list {digits, exp} of the digits in the approximate real number x together with the exponent such that the first digit is the coefficient of 10^(exp - 1).

**`RealDigits[x, b]`**

gives base-b digits.

**`RealDigits[x, b, len]`**

gives len digits.

**`RealDigits[x, b, len, n]`**

gives len digits starting from the coefficient of b^n.

<details>
<summary>Notes</summary>

For rationals with non-terminating expansions the digit list ends in a nested list of the recurring block.  For inexact (machine or MPFR) reals, digits beyond the available precision are returned as Indeterminate.  The sign of x is discarded.

</details>

## Examples (17)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (7)

```mathematica
In[1]:= RealDigits[123.55555]
Out[1]= {{1, 2, 3, 5, 5, 5, 5, 5, 0, 0, 0, 0, 0, 0, 0, 0}, 3}

In[2]:= RealDigits[Pi, 10, 25]
Out[2]= {{3, 1, 4, 1, 5, 9, 2, 6, 5, 3, 5, 8, 9, 7, 9, 3, 2, 3, 8, 4, 6, 2, 6, 4, 3}, 1}

In[3]:= RealDigits[19/7]
Out[3]= {{2, {7, 1, 4, 2, 8, 5}}, 1}

In[4]:= RealDigits[5.635, 10, 20]
Out[4]= {{5, 6, 3, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, Indeterminate, Indeterminate, Indeterminate, Indeterminate}, 1}

In[5]:= RealDigits[Pi, 10, 20, -5]
Out[5]= {{9, 2, 6, 5, 3, 5, 8, 9, 7, 9, 3, 2, 3, 8, 4, 6, 2, 6, 4, 3}, -4}

In[6]:= RealDigits[1.234, 2, 15]
Out[6]= {{1, 0, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 0, 1}, 1}

In[7]:= {RealDigits[1000.], RealDigits[0.999]}
Out[7]= {{{1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, 4}, {{9, 9, 9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, 0}}
```

### Applications (10)

A machine real: sixteen digits and the exponent 3

```mathematica
In[8]:= RealDigits[123.456]
Out[8]= {{1, 2, 3, 4, 5, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, 3}
```

Exact rationals end in a nested list holding the recurring block

```mathematica
In[9]:= RealDigits[22/7]
Out[9]= {{3, {1, 4, 2, 8, 5, 7}}, 1}
```

In base 3, 1/12 = 0.0020202..., a pure cycle starting at 3^-3

```mathematica
In[10]:= RealDigits[1/12, 3]
Out[10]= {{{2, 0}}, -2}
```

Symbolic constants are numericalised once a length fixes the precision

```mathematica
In[11]:= RealDigits[Pi, 10, 25]
Out[11]= {{3, 1, 4, 1, 5, 9, 2, 6, 5, 3, 5, 8, 9, 7, 9, 3, 2, 3, 8, 4, 6, 2, 6, 4, 3}, 1}
```

Hexadecimal digits of E: 2.B7E151628...

```mathematica
In[12]:= RealDigits[E, 16, 10]
Out[12]= {{2, 11, 7, 14, 1, 5, 1, 6, 2, 8}, 1}
```

Five digits starting at the coefficient of 10^-3

```mathematica
In[13]:= RealDigits[1/7, 10, 5, -3]
Out[13]= {{2, 8, 5, 7, 1}, -2}
```

Past the 50 available digits the tail is Indeterminate

```mathematica
In[14]:= RealDigits[N[E, 50], 10, 60]
Out[14]= {{2, 7, 1, 8, 2, 8, 1, 8, 2, 8, 4, 5, 9, 0, 4, 5, 2, 3, 5, 3, 6, 0, 2, 8, 7, 4, 7, 1, 3, 5, 2, 6, 6, 2, 4, 9, 7, 7, 5, 7, 2, 4, 7, 0, 9, 3, 7, 0, 0, 0, Indeterminate, Indeterminate, Indeterminate, Indeterminate, Indeterminate, Indeterminate, Indeterminate, Indeterminate, Indeterminate, Indeterminate}, 1}
```

At and just below a power of the base the exponent stays consistent with the digits

```mathematica
In[15]:= {RealDigits[1000.], RealDigits[0.999]}
Out[15]= {{{1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, 4}, {{9, 9, 9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, 0}}
```

Any integer base works: here the digits are groups of three decimals

```mathematica
In[16]:= RealDigits[2^100, 1000]
Out[16]= {{1, 267, 650, 600, 228, 229, 401, 496, 703, 205, 376}, 11}
```

Without a length there is no precision to pick, so the call stays unevaluated

```mathematica
In[17]:= RealDigits[Pi]
Out[17]= RealDigits[Pi]
```

## Algorithm

real.c

RealDigits builtin -- positional-notation digit expansion.

```text
  RealDigits[x]               default base 10, length set by Precision[x].
  RealDigits[x, b]            base b, length set by Precision[x] / Log10[b].
  RealDigits[x, b, len]       exactly `len` digits, MSD-first.
  RealDigits[x, b, len, n]    `len` digits, first one = coefficient of b^n.

Result form is `{ digits-list, exp }`.  The first element of digits-list is
the coefficient of b^(exp-1).  Sign of x is discarded.  Exact rationals
```

with non-terminating base-b expansions return a list ending in a nested

```text
list of the recurring block.  Inexact (machine or MPFR) reals get
```

Indeterminate for any requested digit beyond the available precision.

```text
  x can be: Integer, BigInt, Rational[n,d], Real (machine), or
            EXPR_MPFR (arbitrary precision; USE_MPFR builds only).
```

The general algorithm scales |x| by base^(-low) where `low` is the lowest digit position we need, floors to an integer N, and reads off the base-b

```text
digits of N (padding with leading zeros as needed).  This single GMP /
MPFR shift handles every numeric type uniformly.  For the special case of
```

an exact rational with no explicit `len`, a remainder-tracked long division detects terminating vs recurring expansions and emits the nested-list form.

```text
Implementation only supports integer bases b >= 2.  Non-integer bases
```

(e.g. GoldenRatio) emit a `::ibase` diagnostic and leave the call unevaluated -- adding them requires MPFR floor-iteration and has been deferred.

## Implementation notes

**Algorithm.** `builtin_realdigits` accepts 1–4 arguments (any other count gives `RealDigits::argb`). It first classifies `x` with `rd_classify` into Integer/BigInt, Rational, machine Real or MPFR. Complex input, or any other numeric non-real, gives `RealDigits::realx`. It then validates the base: an integer `b >= 2` that fits in `unsigned long`, else `::ibase`. It validates the length (a non-negative integer up to `REAL_MAX_LEN = 10^6`, or `Automatic`, else `::intnn`) and the start position (a machine integer, else `::int`). A symbolic `x` such as `Pi`, `E` or `Sqrt[2] + Pi` is numericalised by `rd_try_numericalize` only when a length is given. It asks `numericalize` for MPFR at `len·log2 b + 32` guard bits, padded by `|n|` when a start position is given. With no length the call is left unevaluated. Three paths follow:

1. **Exact input, no length or start position.** `rd_rational_exact_form` runs remainder-tracked long division. It takes the integer-part digits by repeated `mpz_tdiv_qr`, then generates fractional digits and records each remainder. A repeated remainder marks the start of the period, which is emitted as a nested sublist. A zero remainder means the expansion terminates.
2. **Exact input with a length or start position.** `rd_scaled_floor_exact` computes `N = floor(|x|·b^(-low))` in exact `mpq` arithmetic. `rd_digits_lsd` reads off the base-`b` digits of `N`, and the window is padded with leading zeros.
3. **Inexact input.** `rd_inexact_canonical_digits` calls `mpfr_get_str` once, at the available digit count `rd_available_digits` (`prec_bits / log2 b`, rounded), with round-to-nearest. Those digits are placed at positions `exp_p-1 … exp_p-avail`. Positions above that range are `0` and positions below it are `Indeterminate`. The reported exponent comes from `rd_inexact_natural_exp`, which reads it from the same `mpfr_get_str` call at the same digit count. Because the exponent and the digit layout come from one source, they cannot disagree. Before v0.358 the exponent came from a `floor(log x / log b)` quotient and a 2-digit rounding, and both were wrong at a power of the base: `RealDigits[1000.]` lost its leading `1`, and `RealDigits[0.999]` reported exponent `1`. Without MPFR, a double-precision `floor(v·b^-pos) mod b` fallback is used.

Zero is handled separately. Exact `0` gives `{{0}, 0}`. Inexact zero gives `{{0}, -Floor[Accuracy]}`, which is `-323` for a machine `0.`. Both are padded to `len` when a length is given.

**Data structures.** Base, integer part and scaled floor are GMP `mpz_t`; rationals are `mpq_t`. Inexact values are copied into an `mpfr_t` at their own precision (53 bits for a machine real). Cycle detection keeps a growable array of `(remainder, position)` pairs. The result is a plain boxed `List`: the digit list (with a trailing nested `List` for a period), paired with an integer exponent.

**Complexity / limits.** Path 3 costs one `mpfr_get_str` per call, which is quasi-linear in the precision. Path 2 costs one big-integer scaling and a digit loop of `O(len)` divisions. Cycle detection does a linear scan of the remainders seen so far, so a period of length `p` costs `O(p²)` comparisons. It is capped at `REAL_MAX_CYCLE_DIGITS = 10^5` digits; past the cap the expansion is returned flat. Only integer bases are implemented: Real/Rational bases such as `GoldenRatio` emit `::ibase`. `RealDigits` is `Listable` and `Protected`. It returns a nested structure rather than an element-wise array, so there is no ND kernel or `Compile[]` lowering. `tools/nd_fastpath_sweep.py` records it on the `OFF_BUFFER` list, so a packed machine-real list threads element by element.

- `Protected`, `Listable`. Threads over lists in any argument position.
- Works for `Integer`, `BigInt`, `Rational`, machine `Real`, and (under
  `USE_MPFR`) arbitrary-precision `MPFR` numbers. The sign of `x` is
  discarded.
- For integers and rationals with terminating base-`b` expansions, the
  digit list is flat. For rationals with non-terminating expansions, the
  list ends in a nested list giving the recurring block:
  `RealDigits[19/7]` returns `{{2, {7, 1, 4, 2, 8, 5}}, 1}`.
- For inexact (machine or MPFR) reals, the default `len` is
  `Round[Precision[x] / Log10[b]]`. Requesting more digits than the
  precision allows produces `Indeterminate` at the LSB end. The digits
  themselves use the canonical round-to-nearest representation supplied
  by MPFR -- so `RealDigits[123.55555]` returns the literal decimal
  digits, not the binary IEEE tail.
- Symbolic numeric constants such as `Pi`, `E`, and `GoldenRatio` are
  numericalized to MPFR at the matching precision when an explicit `len`
  is given. `RealDigits[Pi]` (no `len`) is left unevaluated.
- `RealDigits[0]` returns `{{0}, 0}`. `RealDigits[0.]` returns
  `{{0}, -Floor[Accuracy[0.]]}` — `{{0}, -323}` for machine precision
  (`Accuracy[0.] ≈ 323.607`), and `{{0}, -p}` for an MPFR zero
  `0``p` of precision `p` digits.
- Bases must be integers `>= 2`. Non-integer (Real / Rational) bases
  trigger `RealDigits::ibase` and leave the call unevaluated.
  Non-integer-base expansions (e.g. `GoldenRatio`) are not yet supported.
- `FromDigits` can be used as the inverse of `RealDigits` for the
  integer / terminating-rational case.

**Attributes:** `Listable`, `Protected`.

## References

**See also:** [Rational](../../arithmetic/Rational/), [Real](../../other-advanced/Real/), [Pi](../../mathematical-constants/Pi/), [E](../../mathematical-constants/E/), [GoldenRatio](../../mathematical-constants/GoldenRatio/), [FromDigits](../../arithmetic/FromDigits/)

- Source: [`src/real.c`](https://github.com/stblake/mathilda/blob/main/src/real.c)
- Specification: [`docs/spec/builtins/arithmetic.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/arithmetic.md)
- Tests: [`tests/test_real_digits.c`](https://github.com/stblake/mathilda/blob/main/tests/test_real_digits.c)

## Notes & additional examples

### Notes

- `RealDigits[x]` returns `{digits, e}`; the first digit is the coefficient of `b^(e-1)`. The sign of `x` is discarded. `Listable`, `Protected`.
- Exact `Integer`/`Rational` input with no length gives the full expansion; a non-terminating one ends in a nested list of the period (`{{3, {1, 4, 2, 8, 5, 7}}, 1}`).
- For inexact input the default length is the available precision in base `b` (16 decimal digits for a machine real, more for an MPFR number). Positions below that precision come back as `Indeterminate`. Digits are the round-to-nearest digits of the stored value, so `RealDigits[123.55555]` gives the decimal digits you typed, not the binary tail.
- `RealDigits[0.]` is `{{0}, -323}`, i.e. `-Floor[Accuracy[0.]]`. `RealDigits[0]` is `{{0}, 0}`.
- The 4-argument form `RealDigits[x, b, len, n]` starts at the coefficient of `b^n`; `len = Automatic` keeps only the digits within the available precision.
- Only integer bases `b >= 2` are supported. A Real or Rational base gives `RealDigits::ibase` and the call stays unevaluated, unlike Mathematica, which accepts non-integer bases. Complex input gives `RealDigits::realx`. All messages go through the message funnel, so `Quiet` and `Check` work on them.
- `RealDigits` returns a nested `{digits, e}` structure, not an element-wise array, so it has no packed-array or `Compile[]` path. A machine-real list threads element by element through `Listable`.
- To rebuild a number from its digits, use `FromDigits[d] b^(e - Length[d])` on a flat digit list `d`. The Mathematica form `FromDigits[{d, e}]` is not supported yet.
- See also: `IntegerDigits`, `FromDigits`, `MantissaExponent`, `DigitCount`.
