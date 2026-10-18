### Worked examples

```mathematica
In[1]:= RealDigits[123.456]  (* a machine real: sixteen digits and the exponent 3 *)
```

```mathematica
In[1]:= RealDigits[22/7]  (* exact rationals end in a nested list holding the recurring block *)
```

```mathematica
In[1]:= RealDigits[1/12, 3]  (* in base 3, 1/12 = 0.0020202..., a pure cycle starting at 3^-3 *)
```

```mathematica
In[1]:= RealDigits[Pi, 10, 25]  (* symbolic constants are numericalised once a length fixes the precision *)
```

```mathematica
In[1]:= RealDigits[E, 16, 10]  (* hexadecimal digits of E: 2.B7E151628... *)
```

```mathematica
In[1]:= RealDigits[1/7, 10, 5, -3]  (* five digits starting at the coefficient of 10^-3 *)
```

```mathematica
In[1]:= RealDigits[N[E, 50], 10, 60]  (* past the 50 available digits the tail is Indeterminate *)
```

```mathematica
In[1]:= {RealDigits[1000.], RealDigits[0.999]}  (* at and just below a power of the base the exponent stays consistent with the digits *)
```

```mathematica
In[1]:= RealDigits[2^100, 1000]  (* any integer base works: here the digits are groups of three decimals *)
```

```mathematica
In[1]:= RealDigits[Pi]  (* without a length there is no precision to pick, so the call stays unevaluated *)
```

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
