# WhittakerW

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`WhittakerW[k, m, z]`**

gives the Whittaker function W\_{k,m}(z).

**`D[WhittakerW[k,m,z],z] = (1/2-k/z) WhittakerW[k,m,z] - WhittakerW[1+k,m,z]/z.`**

<details>
<summary>Notes</summary>

Related to the Tricomi confluent hypergeometric function by W\_{k,m}(z) = E^(-z/2) z^(m+1/2) HypergeometricU\[m-k+1/2, 2m+1, z\]. With an inexact argument it evaluates numerically at machine or arbitrary (MPFR) precision, real or complex, through that confluent form; a real machine argument uses a shared double kernel that also drives the NDArray/packed fast path and Compile\[\]. Values at z = 0: 0 for |Re m| \< 1/2, 1/Gamma\[1-k\] for m = +-1/2, ComplexInfinity for |Re m| \> 1/2. Series (origin and Infinity) for the terminating (elementary) case via LaguerreL. Listable.

</details>

## Examples (13)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= WhittakerW[2, 3, 4.5]
Out[1]= 11.3063

In[2]:= D[WhittakerW[k, m, z], z]
Out[2]= -WhittakerW[1 + k, m, z]/z + (1/2 - k/z) WhittakerW[k, m, z]

In[3]:= WhittakerW[k, 1/2, 0]
Out[3]= 1/Gamma[1 - k]
```

### Applications (10)

Machine value

```mathematica
In[4]:= WhittakerW[2, 3, 4.5]
Out[4]= 11.3063
```

The same value from the defining U form

```mathematica
In[5]:= Exp[-2.25] 4.5^3.5 HypergeometricU[1.5, 7, 4.5]
Out[5]= 11.3063
```

Checked against W_{0,m}(z) = √(z/π) K_m(z/2)

```mathematica
In[6]:= WhittakerW[0, 1/4, 2.] - Sqrt[2./Pi] BesselK[1/4, 1.]
Out[6]= -5.55112e-17
```

Arbitrary precision through HypergeometricU / arb

```mathematica
In[7]:= N[WhittakerW[3/5, 1/3, 1/7], 30]
Out[7]= 0.3970969447006262740475282170027
```

Complex order

```mathematica
In[8]:= WhittakerW[1/5, 2 - I, 2.]
Out[8]= 0.314759 - 1.27345*I
```

Symmetry W_{k,m} = W_{k,-m}

```mathematica
In[9]:= WhittakerW[0.5, 0.25, 1.] - WhittakerW[0.5, -0.25, 1.]
Out[9]= 0.0
```

Values at the origin, decided by 1/2 - |Re m|

```mathematica
In[10]:= {WhittakerW[0, 0, 0], WhittakerW[0, 1, 0], WhittakerW[k, 1/2, 0]}
Out[10]= {0, ComplexInfinity, 1/Gamma[1 - k]}
```

DLMF 13.14.31

```mathematica
In[11]:= D[WhittakerW[k, m, z], z]
Out[11]= -WhittakerW[1 + k, m, z]/z + (1/2 - k/z) WhittakerW[k, m, z]
```

Terminating case, via the Laguerre rewrite

```mathematica
In[12]:= Series[WhittakerW[2, 1/2, x], {x, Infinity, 5}]
Out[12]= -x (2 - x) E^(-1/2 x)
```

Listable

```mathematica
In[13]:= WhittakerW[2, 3, {1.5, 2.5, 3.5}]
Out[13]= {86.3829, 32.8767, 17.992}
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

**Algorithm.** `builtin_whittaker_w` uses the same three stages as WhittakerM, with the U kernels in place of the 1F1 kernels.

1. `whit_ndarray` threads the kernels over a buffer z.
2. `whit_machine` handles a scalar machine call with `sf_machine_whittaker_w`. That kernel folds m to |m|, which is valid because W is even in m, and so b = 2|m|+1 ≥ 1. It then returns `e^(−z/2) z^(|m|+1/2) · sf_machine_hyperu(|m|−k+1/2, 2|m|+1, z)`, with `sf_machine_whittaker_w_c` for z off the positive axis.
3. `whittaker_w_eval` settles z = 0 from the leading exponent `1/2 − |Re m|`: the value is `0` below 1/2, `ComplexInfinity` above it, and `1/Gamma[1−k]` for a real m = ±1/2. Any other inexact input is evaluated through the template `Exp[−z/2] z^(m+1/2) HypergeometricU[m−k+1/2, 2m+1, z]`, so it inherits the HypergeometricU arb path.

In the common case of integer or half-integer m, the U kernel's log series advances its three digammas by recurrence instead of calling digamma at every term. `D` comes from DLMF 13.14.31. The origin and infinity series in `series.c` cover only the terminating case `U(−n, b, z) = (−1)^n n! LaguerreL[n, b−1, z]`.

**Data structures.** The machine scalar and buffer drivers are shared with WhittakerM through the `WhitReal` and `WhitCplx` function pointers. Buffers stay real when every point is real and are promoted to complex otherwise. `REG_N(WhittakerW)` registers the kernel for the Compile VM.

**Complexity / limits.** Machine cost is one `sf_machine_hyperu` call plus the prefactor. The docs report the packed path running about 1.8× faster than SciPy's `exp·z^(m+½)·hyperu`. `Compile[]` lowers the scalar shape but not rank-1 arrays: `CompileDiagnostics[{{v,_Real,1}}, WhittakerW[1,1/2,v]]` reports `Compiled -> False`. The general two-branch origin series and the general asymptotic are not produced. The HypergeometricU terminating-case defect carries over: `WhittakerW[2, 1/2, 1.]` returns `0.0`, while the correct value, −e^(−1/2) = −0.606531, is obtained with `N[WhittakerW[2, 1/2, 1]]`.

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## References

**See also:** [HypergeometricU](../../special-functions/HypergeometricU/), [WhittakerM](../../special-functions/WhittakerM/)

- DLMF §13.14 — Whittaker functions, W_{k,m}(z) = e^(−z/2) z^(m+1/2) U(m−k+1/2, 1+2m, z) (§13.14.3).
- DLMF §13.14.31 (the derivative used for D); §13.2.9 (the integer-b log series behind the machine kernel).
- Source: [`src/special_functions/whittaker.c`](https://github.com/stblake/mathilda/blob/main/src/special_functions/whittaker.c)
- Specification: [`docs/spec/builtins/special-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/special-functions.md)
- Tests: [`tests/test_sf_machine_hyper.c`](https://github.com/stblake/mathilda/blob/main/tests/test_sf_machine_hyper.c)
- Tests: [`tests/test_whittaker.c`](https://github.com/stblake/mathilda/blob/main/tests/test_whittaker.c)

## Notes & additional examples

### Notes

- `Listable`, `NumericFunction`, `Protected`. The definition is `W_{k,m}(z) = e^(−z/2) z^(m+1/2) U(m−k+1/2, 2m+1, z)`, which is even in m. The machine kernel folds m to |m|.
- Exact values at z = 0 depend on the leading exponent 1/2 − |Re m|: the result is `0`, `ComplexInfinity`, or `1/Gamma[1−k]` when m = ±1/2. Other exact input stays symbolic.
- Machine input runs `sf_machine_whittaker_w`, which is the `HypergeometricU` double kernel with the e^(−z/2) z^(m+1/2) prefactor. Buffers share it. Complex orders and precision above machine go through `HypergeometricU`, which uses arb.
- `D` in z follows DLMF 13.14.31. `Series` covers the terminating case only, at 0 and at `Infinity`, through `U(−n, b, z) = (−1)^n n! LaguerreL[n, b−1, z]`.
- See also: `WhittakerM`, `HypergeometricU`.
