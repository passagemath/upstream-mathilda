### Worked examples

```mathematica
In[1]:= HypergeometricPFQ[{1, 1}, {2}, z]  (* the 2F1 shape with a logarithmic closed form *)
```

```mathematica
In[1]:= HypergeometricPFQ[{a, b, c}, {a, d, e}, z]  (* common upper and lower parameters cancel *)
```

```mathematica
In[1]:= HypergeometricPFQ[{-3, 1}, {1}, z]  (* non-positive integer upper parameter: a polynomial *)
```

```mathematica
In[1]:= N[HypergeometricPFQ[{1, 2, 3}, {4, 5}, 1/2], 40]  (* 3F2 summed at 40 digits *)
```

```mathematica
In[1]:= HypergeometricPFQ[{1, 2}, {3, 4}, 0.5 + 0.5 I]  (* complex argument, machine complex kernel *)
```

```mathematica
In[1]:= HypergeometricPFQ[{1, 1}, {3, 3, 3}, {0.5, 1., 2.}]  (* threads over a list of arguments *)
```

```mathematica
In[1]:= D[HypergeometricPFQ[{a}, {b}, x], x]  (* parameters shift by one *)
```

```mathematica
In[1]:= HypergeometricPFQ[{1, 1, 1}, {2, 2}, 1.5]  (* 3F2 outside its unit disc: no continuation, stays unevaluated *)
```

### Notes

- `HypergeometricPFQ[{a1, ...}, {b1, ...}, z]` is the generalized hypergeometric
  function `pFq`. The series converges for all `z` when `p <= q` and for `|z| < 1`
  when `p == q + 1`; only there is the general series summed. Mathematica also
  continues `p == q + 1` outside the disc; Mathilda does that only for the `2F1`
  shape.
- The `0F1`, `1F1` and `2F1` shapes go to dedicated machine kernels, so
  `HypergeometricPFQ[{a}, {b}, z]` and `Hypergeometric1F1[a, b, z]` return the
  same bits.
- Machine-precision input (exact parameters with a machine `z` included) is
  evaluated by one raw-`double` kernel, real or complex, that measures its own
  cancellation and re-sums in double-double when needed. Higher precision uses
  the MPFR series, or arb for `2F1`.
- Attributes `NumericFunction`, `Protected` — not `Listable` (the parameter
  arguments are lists), but a `List` third argument is threaded. Packed arrays
  and `NDArray` arguments map the machine kernel over the buffer. `Compile[]`
  lowers scalar calls (real and complex); the rank-1 array shape is not lowered
  yet.
- See also: `Hypergeometric0F1`, `Hypergeometric1F1`, `Hypergeometric2F1`,
  `HypergeometricU`, `MeijerG`.
