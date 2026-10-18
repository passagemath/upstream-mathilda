# WhittakerM

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`WhittakerM[k, m, z]`**

gives the Whittaker function M\_{k,m}(z).

<details>
<summary>Notes</summary>

Related to the Kummer confluent hypergeometric function by M\_{k,m}(z) = E^(-z/2) z^(m+1/2) Hypergeometric1F1\[m-k+1/2, 2m+1, z\]. With an inexact argument it evaluates numerically at machine or arbitrary (MPFR) precision, real or complex, through that confluent form; a real machine argument uses a shared double kernel that also drives the NDArray/packed fast path and Compile\[\]. Special value: M\_{k,m}(0) = 0 for Re(m) \> -1/2. D\[WhittakerM\[k,m,z\],z\] = (1/2-k/z) WhittakerM\[k,m,z\] + ((1/2+k+m)/z) WhittakerM\[k+1,m,z\]. Series at the origin (the z^(1/2+m) Puiseux form) is supported; Series at Infinity only for the terminating (elementary) case. Listable.

</details>

## Examples (13)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= WhittakerM[2, 3, 4.5]
Out[1]= 68.967

In[2]:= D[WhittakerM[k, m, z], z]
Out[2]= ((1/2 + k + m) WhittakerM[1 + k, m, z])/z + (1/2 - k/z) WhittakerM[k, m, z]

In[3]:= Series[WhittakerM[2, 1/2, x], {x, 0, 6}]
Out[3]= x - x^2 + 3/8 x^3 - 1/12 x^4 + 5/384 x^5 - 1/640 x^6 + O[x]^7
```

### Applications (10)

Machine value

```mathematica
In[4]:= WhittakerM[2, 3, 4.5]
Out[4]= 68.967
```

The same value from the defining 1F1 form

```mathematica
In[5]:= Exp[-2.25] 4.5^3.5 Hypergeometric1F1[1.5, 7, 4.5]
Out[5]= 68.967
```

Checked against M_{0,1/2}(z) = 2 sinh(z/2)

```mathematica
In[6]:= WhittakerM[0, 1/2, 2.] - 2 Sinh[1.]
Out[6]= 4.44089e-16
```

Arbitrary precision through Hypergeometric1F1

```mathematica
In[7]:= N[WhittakerM[2, 1/3, 1/7], 30]
Out[7]= 0.165649771906678290310364828072
```

Negative z: the principal branch of z^(m+1/2)

```mathematica
In[8]:= WhittakerM[1, 1/3, -2.]
Out[8]= -4.8379 + 2.79316*I
```

Complex order

```mathematica
In[9]:= WhittakerM[2, 3 - I, 4.5]
Out[9]= -18.2886 - 74.6279*I
```

Zero at the origin when Re m > -1/2

```mathematica
In[10]:= WhittakerM[1, 1, 0]
Out[10]= 0
```

DLMF 13.4.27

```mathematica
In[11]:= D[WhittakerM[k, m, z], z]
Out[11]= ((1/2 + k + m) WhittakerM[1 + k, m, z])/z + (1/2 - k/z) WhittakerM[k, m, z]
```

Origin series

```mathematica
In[12]:= Series[WhittakerM[2, 1/2, x], {x, 0, 6}]
Out[12]= x - x^2 + 3/8 x^3 - 1/12 x^4 + 5/384 x^5 - 1/640 x^6 + O[x]^7
```

At infinity only the terminating case is expanded: here the exact closed form

```mathematica
In[13]:= Series[WhittakerM[2, 1/2, x], {x, Infinity, 5}]
Out[13]= x E^(-1/2 x) (1 - 1/2 x)
```

## Algorithm

Mathilda -- WhittakerM, the Whittaker M function M_{k,m}(z).

```text
  WhittakerM[k, m, z] = e^(-z/2) z^(m+1/2) Hypergeometric1F1[m-k+1/2, 2m+1, z].
```

Evaluation is layered so each argument shape takes the cheapest numeric route (mirroring the LaguerreL / HermiteH confluent-hypergeometric design):

```text
  z == 0 (exact), Re(m) > -1/2   ->  exact 0 (M_{k,m} ~ z^(m+1/2) near 0).
  machine-precision call         ->  the shared double helpers
  (k, m exact or machine real;         sf_machine_whittaker_m (real z > 0) /
   z machine real or complex)          sf_machine_whittaker_m_c (complex z, or
                                       real z < 0 on the principal branch),
                                     bit-for-bit with the NDArray element-wise
                                     buffer loop and the Compile[] VM.
                                     Declines (2m+1 pole, cancellation,
                                      overflow) -> template.
  some inexact argument          ->  the confluent-hypergeometric form built as
                                       e^(-z/2) z^(m+1/2) Hypergeometric1F1[...],
                                     evaluated so it inherits 1F1's machine +
                                     MPFR, real + complex numerics. Kept only
                                     if it collapses to a number.
  everything else (symbolic)     ->  stays symbolic (return NULL).
```

WhittakerM is a single growing PRODUCT (prefactor * 1F1), not a difference, so it adds no cancellation over 1F1 itself and needs no guard digits (unlike HermiteH): bad 1F1 cancellation is declined inside sf_machine_pfq / re-summed in MPFR by the 1F1 evaluator.

Attributes: Listable, NumericFunction, Protected.

Deferred / best-effort, as the sibling special functions: tight Interval enclosures (the generic NumericFunction + D[] interval threading in eval.c declines here because D[WhittakerM,z] mentions WhittakerM itself -- the self-referential-derivative guard), the general asymptotic Series at infinity (a two-scale Stokes transseries; only the terminating/elementary case is emitted, in calculus/series.c), symbolic FunctionExpand, and general-n SeriesCoefficient. The D[] rule and the origin Series live in src/calculus/deriv.c and src/calculus/series.c respectively.

## Implementation notes

**Algorithm.** `builtin_whittaker_m` runs three stages in order.

1. `whit_ndarray` threads the machine kernels over a buffer z, holding k and m constant.
2. `whit_machine` handles a scalar machine call. It applies when k and m are exact or machine reals and z is a machine real or complex number, and it calls `sf_machine_whittaker_m`. That kernel computes `e^(−z/2) z^(m+1/2) · sf_machine_hyper1f1(m−k+1/2, 2m+1, z)` and declines inside the 1F1 when cancellation is too large. For a real z < 0, or a complex z, the call goes to the complex twin `sf_machine_whittaker_m_c`, which takes the principal branch.
3. `whittaker_m_eval` returns the exact `0` at z = 0 when Re m > −1/2. For any other inexact input, `whittaker_m_numeric_hyp` builds `Exp[−z/2] z^(m+1/2) Hypergeometric1F1[m−k+1/2, 2m+1, z]` and evaluates it. The result inherits 1F1's MPFR and complex numerics, and it is kept only if it collapses to a number.

`D` in z comes from DLMF 13.4.27, implemented in `deriv.c`. `series.c` emits the Puiseux origin series `z^(1/2+m)(...)`, and at `Infinity` it emits only the terminating closed form, which arises when `m−k+1/2` is a non-positive integer.

**Data structures.** The machine path keeps no state; it takes doubles in and returns doubles or complex doubles. The function pointer typedefs `WhitReal` and `WhitCplx` let WhittakerM and WhittakerW share the scalar and buffer drivers. A real buffer gets a fast pass that stays real. A point on the cut promotes the whole output buffer to complex. `REG_N(WhittakerM)` registers the n-ary kernel for the Compile VM.

**Complexity / limits.** Machine cost is one 1F1 evaluation plus `exp` and `pow`. The documented timings are 0.2 µs for a real argument and 0.4 µs for a complex one. Precision above machine costs whatever `Hypergeometric1F1` costs. `Compile[]` lowers scalar arguments (real and complex). The rank-1 array shape is not lowered: `CompileDiagnostics[{{v,_Real,1}}, WhittakerM[1,1/2,v]]` reports `Compiled -> False`. The general Stokes asymptotic at infinity is not implemented, and neither is `FunctionExpand`. An `Interval` argument stays unevaluated, because the derivative refers back to WhittakerM.

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## References

**See also:** [SeriesCoefficient](../../power-series/SeriesCoefficient/), [FunctionExpand](../../special-functions/FunctionExpand/)

- DLMF §13.14 — Whittaker functions, M_{k,m}(z) = e^(−z/2) z^(m+1/2) M(m−k+1/2, 1+2m, z) (§13.14.2).
- DLMF §13.4.27 (the derivative recurrence used for D).
- Source: [`src/special_functions/whittaker.c`](https://github.com/stblake/mathilda/blob/main/src/special_functions/whittaker.c)
- Specification: [`docs/spec/builtins/special-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/special-functions.md)
- Tests: [`tests/test_dsolve_m19_stress.c`](https://github.com/stblake/mathilda/blob/main/tests/test_dsolve_m19_stress.c)
- Tests: [`tests/test_dsolve_m55_stress.c`](https://github.com/stblake/mathilda/blob/main/tests/test_dsolve_m55_stress.c)
- Tests: [`tests/test_sf_machine_hyper.c`](https://github.com/stblake/mathilda/blob/main/tests/test_sf_machine_hyper.c)
- Tests: [`tests/test_whittaker.c`](https://github.com/stblake/mathilda/blob/main/tests/test_whittaker.c)

## Notes & additional examples

### Notes

- `Listable`, `NumericFunction`, `Protected`. The definition is `M_{k,m}(z) = e^(−z/2) z^(m+1/2) 1F1(m−k+1/2; 2m+1; z)`.
- Exact input stays symbolic, with one exception: `M_{k,m}(0) = 0` when Re m > −1/2.
- Machine input (k and m real; z real or complex) runs the shared double kernel `sf_machine_whittaker_m` or its complex twin, and so do buffers. A real z < 0 gives the principal complex value. If the kernel declines, or for complex k, m or precision above machine, the 1F1 form is evaluated at the input precision.
- `D` in z follows DLMF 13.4.27. Derivatives in k and m stay inert.
- `Series` at 0 gives the `z^(m+1/2)` Puiseux form, for symbolic m as well. At `Infinity` only the terminating (elementary) case is expanded. `FunctionExpand` does not expand `WhittakerM`, and an `Interval` argument stays unevaluated.
- See also: `WhittakerW`, `Hypergeometric1F1`.
