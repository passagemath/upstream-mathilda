# HypergeometricPFQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`HypergeometricPFQ[{a1, ...}, {b1, ...}, z]`**

is the generalized hypergeometric function pFq(a;b;z), the series

<details>
<summary>Notes</summary>

Sum (prod\_i Pochhammer\[a\_i, k\] / prod\_j Pochhammer\[b\_j, k\]) z^k / k!. Common upper/lower parameters cancel; a non-positive integer upper parameter terminates the series to a polynomial. Evaluates to machine, arbitrary-precision (MPFR), and complex numbers in the convergent regime (p\<=q for all z; p==q+1 for |z|\<1), and for 2F1 in the whole plane.

</details>

## Examples (15)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (7)

```mathematica
In[1]:= HypergeometricPFQ[{a1, a2, a3}, {b1, b2, b3}, 0]
Out[1]= 1

In[2]:= HypergeometricPFQ[{}, {}, z]
Out[2]= E^z

In[3]:= HypergeometricPFQ[{a, b, c}, {a, d, e}, z]
Out[3]= HypergeometricPFQ[{b, c}, {d, e}, z]

In[4]:= HypergeometricPFQ[{1, 1}, {3, 3, 3}, 2.]
Out[4]= 1.07893

In[5]:= HypergeometricPFQ[{1, 2, 3, 4}, {5, 6, 7}, {0.1, 0.3, 0.5}]
Out[5]= {1.01164, 1.03627, 1.06296}

In[6]:= N[HypergeometricPFQ[{1, 1, 1}, {3/2, 3/2, 3/2}, 10], 50]
Out[6]= 530.19188827362590438855961685444087792733053398358

In[7]:= D[HypergeometricPFQ[{a1, a2}, {b1, b2, b3}, x], x]
Out[7]= (a1 a2 HypergeometricPFQ[{1 + a1, 1 + a2}, {1 + b1, 1 + b2, 1 + b3}, x])/(b1 b2 b3)
```

### Applications (8)

The 2F1 shape with a logarithmic closed form

```mathematica
In[8]:= HypergeometricPFQ[{1, 1}, {2}, z]
Out[8]= -Log[1 - z]/z
```

Common upper and lower parameters cancel

```mathematica
In[9]:= HypergeometricPFQ[{a, b, c}, {a, d, e}, z]
Out[9]= HypergeometricPFQ[{b, c}, {d, e}, z]
```

Non-positive integer upper parameter: a polynomial

```mathematica
In[10]:= HypergeometricPFQ[{-3, 1}, {1}, z]
Out[10]= 1 - 3 z + 3 z^2 - z^3
```

3F2 summed at 40 digits

```mathematica
In[11]:= N[HypergeometricPFQ[{1, 2, 3}, {4, 5}, 1/2], 40]
Out[11]= 1.1898747542564229318256831180919799547257
```

Complex argument, machine complex kernel

```mathematica
In[12]:= HypergeometricPFQ[{1, 2}, {3, 4}, 0.5 + 0.5 I]
Out[12]= 1.0824 + 0.0966608*I
```

Threads over a list of arguments

```mathematica
In[13]:= HypergeometricPFQ[{1, 1}, {3, 3, 3}, {0.5, 1., 2.}]
Out[13]= {1.01881, 1.03822, 1.07893}
```

Parameters shift by one

```mathematica
In[14]:= D[HypergeometricPFQ[{a}, {b}, x], x]
Out[14]= (a HypergeometricPFQ[{1 + a}, {1 + b}, x])/b
```

3F2 outside its unit disc: no continuation, stays unevaluated

```mathematica
In[15]:= HypergeometricPFQ[{1, 1, 1}, {2, 2}, 1.5]
Out[15]= HypergeometricPFQ[{1, 1, 1}, {2, 2}, 1.5]
```

## Options & behaviour

> **Packed arrays.** With a real or complex machine buffer as `z` and every
> parameter a machine Real or exact Integer / Rational, the machine kernel is
> mapped over the buffer — the same kernel the scalar call uses, so the results
> agree bit-for-bit. For real parameters over a float64 buffer a vector driver
> hoists the parameter work out of the loop (for `2F1`, the series coefficients
> are tabulated once per buffer; the scalar path forms the identical
> coefficients on the fly). A `2F1` element on the cut widens the result to
> complex; an int64 buffer, an MPFR parameter, or a declined element takes the
> element-wise route. For `2F1` series that are one-signed and never grow
> (decided from the parameters), the buffer sums four elements in lock-step with
> the scalar loop's exact arithmetic. 10^6 elements: `1F1(1;2;x)` 11 ms (was
> 1.3 s; SciPy 27 ms), `2F1(1,2;3;x)` 0.14 s (was 2.6 s; SciPy 0.16 s),
> `2F2(1,2;3,4;x)` 0.19 s (was 1.3 s).

## Implementation notes

**Algorithm.** `builtin_hypergeometric_pfq` handles `pFq(a;b;z)`, the series
`Sum_k (prod_i (a_i)_k / prod_j (b_j)_k) z^k / k!`. Dispatch order: validate
(arity 3, first two arguments `List`); a packed / `NDArray` `z` goes to
`hg_pfq_ndarray` (else `ndarray_delist_and_reeval`); a plain `List` `z` is
threaded by hand; then the **machine gate** `hg_machine_pfq` — every parameter and
`z` a machine-scalar leaf (Integer, Rational, Real, or `Complex` of those) with at
least one machine `Real` — answers from the raw-`double` kernel *before* any exact
rewrite, so exact parameters with a machine `z` take the same path and bits as the
buffer and `Compile[]` surfaces. Only after that come the exact steps: **generic
cancellation** of a parameter common to both lists (`try_cancel`), `z == 0 -> 1`,
**termination** to the explicit polynomial for a non-positive integer upper
parameter (`try_terminate`), **reduction** to closed forms (`try_reduce`: `0F0 =
Exp`, `1F0 = (1-z)^-a`, `0F1(1/2)`/`0F1(3/2)` as `Cosh`/`Sinh`, `1F1(1;2;z)`, the
`2F1` logarithm and `ArcSin` families, and a table-backed very-well-poised `4F3`
Ramanujan `1/Pi` value), and finally `try_numeric` for non-machine inexact input.
There the governing precision is the *minimum* inexact precision of `z` and every
parameter; above 53 bits a `2F1` with `|z| < 1` goes to arb
(`flint_num_hyp2f1`, `acb_hypgeom_2f1`), anything else to the MPFR `(re, im)`
series `hgpfq_mpfr_sum`. A machine partial sum that overflows promotes the whole
call to MPFR (`numeric_promote_machine_call`).

**Data structures.** The machine kernel lives in
`src/special_functions/sf_machine_hyper.c` (`sf_machine_pfq` / `sf_machine_pfq_c`,
plain C99 `double` with complex values as explicit re/im pairs, no
`<complex.h>`). `sf_machine_pfq` sends `{}/{b}`, `{a}/{b}` and `{a,b}/{c}` to the
`0F1` / `1F1` / `2F1` kernels, so every spelling of one function reaches the same
code; other shapes use a **loss-measuring** series: summed once in double while
tracking its peak term, and re-summed in double-double (`hdd`, Dekker /
Hida–Li–Bailey) when `peak/|sum|` exceeds 2 or the peak lies past the 6th term,
declining beyond `2^45` of cancellation (`HY_LOSS_MAX`). Parameters are collected
in `HgParams` (at most 16 upper and 16 lower). The buffer path runs the vector
driver `sf_machine_pfq_vec` for real parameters over a float64 buffer (the `2F1`
coefficients are tabulated once per buffer and four elements are summed in
lock-step, bit-identical to the scalar loop), and otherwise a general loop that
widens the output dtype to complex when a parameter is complex or a `2F1` element
falls on the cut. ND: n-ary kernel `NDKN_HypergeometricPFQ` (`REG_N`), which
receives `p` as a leading operand; all four hypergeometric heads are on `pack.c`'s
`AWARE` list. Attributes: `NumericFunction`, `Protected` — not `Listable`, since
the first two arguments are themselves lists.

**Complexity / limits.** Series term caps: 20000 (double), 6000 (double-double),
`HGPFQ_MACHINE_MAX_TERMS = 200000` / `HGPFQ_MPFR_MAX_TERMS = 1000000` on the
generic path. The general series is summed only where it converges — `p <= q`
for any `z`, `p == q+1` for `|z| < 1`; other shapes stay unevaluated, except that
the `0F1`/`1F1`/`2F1` shapes have their own continuations (`2F1` is continued to
the whole plane at machine precision; see `Hypergeometric2F1`). `Compile[]`
lowers the head at real and complex scalar shapes (the two literal lists are
destructured in `compile_emit_mathfn.c` into one `OP_KERNN` / `OP_KERNN_C`
operand block, capped at `KERNN_MAX_OPERANDS = 24`); the rank-1 array shape
reports `Compiled -> False`, while the packed / `NDArray` buffer path still runs
the kernel at the REPL.

- Attributes `NumericFunction`, `Protected`.
- `HypergeometricPFQ[a, b, 0]` is `1`; threads over a `List` third argument.
- Parameters common to the upper and lower lists cancel; a non-positive
  integer upper parameter terminates the series to a polynomial (valid for
  symbolic `z`).
- Reduces to elementary functions for simple parameters: `0F0 -> E^z`,
  `1F0(a) -> (1-z)^(-a)`, `0F1(1/2) -> Cosh[2 Sqrt[z]]`,
  `0F1(3/2) -> Sinh[2 Sqrt[z]]/(2 Sqrt[z])`, `1F1(1;2) -> (E^z-1)/z`,
  `2F1(1,1;2) -> -Log[1-z]/z`. The central-binomial / arcsin family also closes:
  `2F1(1,1;1/2;z) -> 1/(1-z) + Sqrt[z] ArcSin[Sqrt[z]]/(1-z)^(3/2)`,
  `2F1(1,1;3/2;z) -> ArcSin[Sqrt[z]]/(Sqrt[z] Sqrt[1-z])`,
  `2F1(2,1;3/2;z) -> 1/(2(1-z)) + ArcSin[Sqrt[z]]/(2 Sqrt[z] (1-z)^(3/2))`.
  A table-backed very-well-poised reduction closes the classic
  central-binomial-cubed Ramanujan `1/Pi` series:
  `4F3({1/2,1/2,1/2,5/4}, {1/4,1,1}, -1) -> 2/Pi` (general `1/Pi` summation is an
  open problem, so this is a recognizer for the known class, not universal).
- Numeric evaluation at machine, arbitrary (MPFR), and complex precision, with
  output precision tracking the input. The general series is summed where it
  converges — `p <= q` (entire) and `p == q+1` with `|z| < 1`; for `p > q+1`,
  and `p == q+1` with `|z| >= 1` other than `2F1`, the call stays unevaluated.
  `2F1` (`{a, b}, {c}`) is analytically continued to the whole plane (see
  `Hypergeometric2F1`).
- **Machine precision: one raw-`double` kernel** (`sf_machine_hyper.c`), real and
  complex (raw re/im), answering every machine-precision call — exact
  (Integer / Rational) parameters included — *before* the exact rewrites, so the
  scalar result is bit-identical to the packed / NDArray buffer and `Compile[]`
  results. `{}/{b}`, `{a}/{b}` and `{a, b}/{c}` dispatch to the
  `Hypergeometric0F1` / `1F1` / `2F1` kernels, so every spelling of the same
  function returns the same bits.
- **Cancellation and rounding are measured, not ignored.** Each series is
  summed once in double with compensated (TwoSum) summation while tracking its
  largest term and where it occurs. If the largest term exceeds the sum by more
  than 2× (cancellation — e.g. the alternating series at negative `z`, whose
  terms reach `e^|z|` times the result) or lies deeper than the 6th term (each
  term then carries the rounding of the ratio products before it), the same
  series is re-summed in double-double arithmetic (~106 bits); beyond 2^45 of
  cancellation the kernel declines and the MPFR series answers.
  `1F1(1;2;-40)` agrees with `(E^z-1)/z` to one ulp; `2F1(1,2;3;0.89)` (380
  terms) is within 1.4 ulp where the plain sum was 13 ulp off.
- Compiles (`Compile`, `Plot`, `NIntegrate`, `NDSolve`, …) at real and complex
  scalar shapes; a complex argument lowers to the kernel's complex arm
  (`OP_KERNN_C`). The rank-1 array shape of the n-ary heads is still open
  (`COMPILE_MISSING.md`).
- `D[HypergeometricPFQ[{a},{b},x], x]
   = (prod a_i / prod b_j) HypergeometricPFQ[{a_i+1},{b_j+1},x]`.

**Attributes:** `NumericFunction`, `Protected`.

## References

**See also:** [List](../../other-advanced/List/), [Hypergeometric2F1](../../special-functions/Hypergeometric2F1/), [Hypergeometric0F1](../../special-functions/Hypergeometric0F1/), [Compile](../../control-flow/Compile/), [Plot](../../graphics/Plot/), [NIntegrate](../../numerical-calculus/NIntegrate/), [NDSolve](../../numerical-calculus/NDSolve/)

- DLMF §16 — generalized hypergeometric functions pFq.
- T. J. Dekker, A floating-point technique for extending the available precision, Numer. Math. 18 (1971) 224-242 — the double-double re-sum.
- Y. Hida, X. S. Li and D. H. Bailey, Algorithms for quad-double precision floating point arithmetic, ARITH-15 (2001).
- Source: [`src/special_functions/hypergeopfq.c`](https://github.com/stblake/mathilda/blob/main/src/special_functions/hypergeopfq.c)
- Specification: [`docs/spec/builtins/special-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/special-functions.md)
- Tests: [`tests/test_compile_coverage.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile_coverage.c)
- Tests: [`tests/test_hypergeopfq.c`](https://github.com/stblake/mathilda/blob/main/tests/test_hypergeopfq.c)
- Tests: [`tests/test_integrate_ramanujan.c`](https://github.com/stblake/mathilda/blob/main/tests/test_integrate_ramanujan.c)
- Tests: [`tests/test_meijerg.c`](https://github.com/stblake/mathilda/blob/main/tests/test_meijerg.c)

## Notes & additional examples

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
