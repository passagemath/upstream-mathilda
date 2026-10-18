# HypergeometricU

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`HypergeometricU[a, b, z]`**

gives the Tricomi confluent hypergeometric function U(a, b, z).

**`D[HypergeometricU[a,b,z],z] = -a HypergeometricU[a+1,b+1,z]. Listable.`**

<details>
<summary>Notes</summary>

The second solution of z w'' + (b-z) w' - a w = 0, with a branch cut on z from -Infinity to 0. With an inexact argument it evaluates numerically at machine or arbitrary (MPFR) precision, real or complex: a real machine argument uses a shared double kernel (also driving the NDArray/packed fast path and Compile\[\]); MPFR/complex/z\<=0 use FLINT's rigorous acb\_hypgeom\_u. Exact reductions U(0,b,z) = 1 and U(a,a+1,z) = z^(-a).

</details>

## Examples (14)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= HypergeometricU[1, 2, z]
Out[1]= 1/z

In[2]:= HypergeometricU[1.5, 7, 4.5]
Out[2]= 0.554931

In[3]:= HypergeometricU[1/2, 1, 2.]
Out[3]= 0.645694

In[4]:= D[HypergeometricU[a, b, z], z]
Out[4]= -a HypergeometricU[1 + a, 1 + b, z]
```

### Applications (10)

The reduction U(a, a+1, z) = z^-a

```mathematica
In[5]:= HypergeometricU[1, 2, z]
Out[5]= 1/z
```

Machine value

```mathematica
In[6]:= HypergeometricU[1.5, 7, 4.5]
Out[6]= 0.554931
```

The Laplace-type integral U = (1/Γ(a))∫ e^(-zt) t^(a-1) (1+t)^(b-a-1) dt gives the same

```mathematica
In[7]:= NIntegrate[Exp[-4.5 t] t^0.5 (1 + t)^4.5, {t, 0, Infinity}]/Gamma[1.5]
Out[7]= 0.554931
```

Arbitrary precision through FLINT's acb_hypgeom_u

```mathematica
In[8]:= N[HypergeometricU[3/2, 7, 9/2], 30]
Out[8]= 0.5549311339066449498555946197268
```

Integer b: the log series, checked against U(1/2,1,z) = e^(z/2) K_0(z/2)/√π

```mathematica
In[9]:= HypergeometricU[1/2, 1, 2.] - Exp[1.] BesselK[0, 1.]/Sqrt[Pi]
Out[9]= 1.11022e-16
```

Real z < 0 lies on the cut: the principal complex value

```mathematica
In[10]:= HypergeometricU[0.7, 1.3, -2.]
Out[10]= -0.354426 - 0.642952*I
```

Complex parameters go to arb, demoted to machine precision

```mathematica
In[11]:= HypergeometricU[3 - I, 2 I, 3.0]
Out[11]= -0.00526539 + 0.000793166*I
```

DLMF 13.3.22

```mathematica
In[12]:= D[HypergeometricU[a, b, z], z]
Out[12]= -a HypergeometricU[1 + a, 1 + b, z]
```

Listable over a packed buffer

```mathematica
In[13]:= HypergeometricU[1, 2, {2., 4., 5.}]
Out[13]= {0.5, 0.25, 0.2}
```

The Laguerre-polynomial case is not expanded for symbolic arguments

```mathematica
In[14]:= HypergeometricU[-2, b, z]
Out[14]= HypergeometricU[-2, b, z]
```

## Algorithm

Mathilda -- HypergeometricU, the Tricomi confluent hypergeometric U(a, b, z).

```text
  HypergeometricU[a, b, z]  is the second solution of the confluent
  hypergeometric equation z w'' + (b - z) w' - a w = 0, with the branch cut on
  z from -infinity to 0 (the z^(1-b) term in its connection to Hypergeometric1F1).
```

Evaluation is layered so each argument shape takes the cheapest accurate route (mirroring Hypergeometric1F1 / WhittakerM):

```text
  U(0, b, z) = 1,  U(a, a+1, z) = z^(-a)   -> exact (symbolic or numeric).
  machine-precision call (a, b exact or  -> the shared double helpers
  machine real; z machine real/complex)     sf_machine_hyperu (real z > 0) /
                                            sf_machine_hyperu_c (complex z, or
                                            real z < 0 on the principal branch),
                                            bit-for-bit with the NDArray buffer
                                            loop and the Compile[] VM; checked
                                            BEFORE the exact reductions below.
                                            Declines -> FLINT.
  some inexact argument                  -> FLINT's rigorous acb_hypgeom_u via
                                            flint_num_hypergeometric_u (MPFR,
                                            complex, z <= 0). Kept iff it
                                            collapses to a number.
  everything else (symbolic)             -> stays symbolic (return NULL).
```

Attributes: Listable, NumericFunction, Protected.

Deferred, as the confluent siblings: the general two-branch Series at the origin (z^(1/2 +- ...) / log structure -- U is not analytic at 0) and FunctionExpand. The D[] rule (DLMF 13.3.22) lives in src/calculus/deriv.c.

## Implementation notes

**Algorithm.** `builtin_hypergeometric_u` checks argument shapes in the following order.

1. A buffer z (`NDArray` or packed list) goes to `hypergeometricu_ndarray`. That function never declines a buffer: it either returns a packed result or hands off to `ndarray_delist_and_reeval`.
2. A machine-scalar call goes to `hu_machine`. The call qualifies when a and b are real, z is real or complex, and at least one argument is a machine Real. This gate comes before the exact reductions.
3. `hypergeometricu_eval` applies the exact rules `U(0, b, z) = 1` and `U(a, a+1, z) = z^(−a)`. The second rule tests whether `b − a − 1` evaluates to the literal 0. Any other inexact input goes to FLINT's `acb_hypgeom_u` through `flint_num_hypergeometric_u`. When every inexact argument was machine precision, the result is demoted back to machine doubles.

The real kernel `sf_machine_hyperu` requires z > 0 and tries four regimes in turn:

- `hyperu_asymptotic`, the series z^(−a) 2F0(a, a−b+1; ; −1/z) with TwoSum compensation. It returns exactly when the series terminates. A truncated sum is accepted only if its least term is below half an ulp of the sum.
- `hyperu_integer_b`, the DLMF 13.2.9 log series for integer b ≥ 1. Its three digammas are each seeded once and then advanced with ψ(x+1) = ψ(x) + 1/x. For b ≤ 0 the kernel first applies the fold `U = z^(1−b) U(a−b+1, 2−b, z)`.
- `hyperu_connection`, the Kummer formula with two 1F1 terms, used for non-integer b.
- `sf_machine_hyperu_quad`, an exp-sinh double-exponential quadrature of the integral representation (step 0.125), used for a > 0. For a ≤ 0 it is applied to the Kummer image instead.

The integer-b and connection regimes decline when cancellation amplifies rounding by more than a factor of 4. A real z < 0 and any complex z go to `sf_machine_hyperu_c`, which returns the value on the principal branch.

**Data structures.** Both kernels work on scalar doubles only. `hu_point` dispatches one point between the real and complex kernels. In the buffer loop, a real `NDT_FLOAT64` buffer gets a fast pass. When any point lands on the branch cut, the loop switches to an interleaved re/im scratch buffer and promotes the output to complex64 (or complex32). `REG_N(HypergeometricU)` in `src/ndkernels.c` exposes the real kernel to the Compile VM.

**Complexity / limits.** The asymptotic regime runs up to about z terms (cap 400). The log series is linear in its term count. The quadrature uses 50–90 nodes and gives up when a > 160. `Compile[]` lowers the scalar shape only: `CompileDiagnostics[{{v,_Real,1}}, HypergeometricU[1,2,v]]` reports `Compiled -> False`. The terminating Laguerre case (a a non-positive integer) is not expanded symbolically. There is no Series at the origin and no `FunctionExpand` rule. One defect: when a is a non-positive integer, b is an integer and z is small (for example `HypergeometricU[-1, 2, 1.]`), the asymptotic regime stops at its growing first term. `hyperu_integer_b` then returns 0 because both of its 1/Γ prefactors vanish. That result is wrong.

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## References

**See also:** [Hypergeometric1F1](../../special-functions/Hypergeometric1F1/), [FunctionExpand](../../special-functions/FunctionExpand/)

- DLMF §13.2 — the Kummer and Tricomi functions; §13.2.9 — U at integer b (the logarithmic case); §13.2.42 — the connection between U and 1F1.
- DLMF §13.4.4 — the integral representation U = (1/Γ(a)) ∫₀^∞ e^(−zt) t^(a−1) (1+t)^(b−a−1) dt; §13.7 — the asymptotic expansion; §13.3.22 — the z-derivative.
- Source: [`src/special_functions/hypergeometricu.c`](https://github.com/stblake/mathilda/blob/main/src/special_functions/hypergeometricu.c)
- Specification: [`docs/spec/builtins/special-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/special-functions.md)
- Tests: [`tests/test_hypergeometricu.c`](https://github.com/stblake/mathilda/blob/main/tests/test_hypergeometricu.c)
- Tests: [`tests/test_sf_machine_hyper.c`](https://github.com/stblake/mathilda/blob/main/tests/test_sf_machine_hyper.c)
- Tests: [`tests/test_whittaker.c`](https://github.com/stblake/mathilda/blob/main/tests/test_whittaker.c)

## Notes & additional examples

### Notes

- `Listable`, `NumericFunction`, `Protected`. The function is the Tricomi solution of `z w'' + (b − z) w' − a w = 0`. Its branch cut runs along z ∈ (−∞, 0].
- Exact reductions: `U(0, b, z) = 1` and `U(a, a+1, z) = z^(−a)`. All other exact input stays symbolic, including the Laguerre-polynomial case of a non-positive integer a.
- Machine input (real a and b; z real or complex) runs a raw-double kernel. It tries, in order: the asymptotic 2F0 series, accepted only when its least term is below half an ulp; the DLMF 13.2.9 log series for integer b, whose digammas advance by recurrence; the Kummer two-1F1 connection for non-integer b; and an exp-sinh quadrature of the integral form. Packed and `NDArray` buffers share the kernel. A declined point goes to FLINT's `acb_hypgeom_u`, and the result is rounded back to machine precision.
- Complex parameters and precision above machine always use `acb_hypgeom_u`.
- `D` in z follows DLMF 13.3.22. Derivatives in a and b stay inert. There is no origin series and no `FunctionExpand` rule.
- See also: `Hypergeometric1F1`, `WhittakerW`, `MeijerG`.
