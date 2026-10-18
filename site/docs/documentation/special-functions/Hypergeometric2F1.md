# Hypergeometric2F1

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Hypergeometric2F1[a, b, c, z]`**

is the Gauss hypergeometric 2F1, equal to HypergeometricPFQ\[{a, b}, {c}, z\]; continued analytically to the whole

<details>
<summary>Notes</summary>

plane, with the cut z \> 1 taking the limit from below.

</details>

## Examples (15)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= Hypergeometric2F1[1, 1, 2, z]
Out[1]= -Log[1 - z]/z

In[2]:= Hypergeometric2F1[1, 1, 3/2, z]
Out[2]= ArcSin[Sqrt[z]]/(Sqrt[z] Sqrt[1 - z])

In[3]:= Hypergeometric2F1[1, 2, 3, 2.5]
Out[3]= -0.929749 - 1.00531*I

In[4]:= N[Hypergeometric2F1[1/3, 1/2, 1/5, 99/100], 50]
Out[4]= 24.8565393133858972090022500136712230259988212774851
```

### Applications (11)

The logarithm

```mathematica
In[5]:= Hypergeometric2F1[1, 1, 2, z]
Out[5]= -Log[1 - z]/z
```

Terminates to (1 - z)^3

```mathematica
In[6]:= Hypergeometric2F1[-3, 1, 1, z]
Out[6]= 1 - 3 z + 3 z^2 - z^3
```

The central-binomial ArcSin family

```mathematica
In[7]:= Hypergeometric2F1[1, 1, 3/2, z]
Out[7]= ArcSin[Sqrt[z]]/(Sqrt[z] Sqrt[1 - z])
```

Since z 2F1(1/2,1/2;3/2;z^2) = ArcSin[z], here at z = 1/2

```mathematica
In[8]:= N[Hypergeometric2F1[1/2, 1/2, 3/2, 1/4]/2, 40]
Out[8]= 0.52359877559829887307710723054658381403285
```

The same value from ArcSin

```mathematica
In[9]:= N[ArcSin[1/2], 40]
Out[9]= 0.52359877559829887307710723054658381403285
```

Far outside the unit disc, by a connection formula

```mathematica
In[10]:= Hypergeometric2F1[0.3, 0.7, 1.9, -50.]
Out[10]= 0.487749
```

On the cut z > 1: the limit from below, ArcSin[Sqrt[3]]/Sqrt[3]

```mathematica
In[11]:= Hypergeometric2F1[0.5, 0.5, 1.5, 3.]
Out[11]= 0.9069 - 0.661768*I
```

Complex argument beyond |z| = 1

```mathematica
In[12]:= Hypergeometric2F1[0.3, 0.7, 1.9, 2. + 1. I]
Out[12]= 1.03477 + 0.307616*I
```

Arbitrary precision near the unit circle, via arb

```mathematica
In[13]:= N[Hypergeometric2F1[1/3, 2/3, 5/6, 99/100], 50]
Out[13]= 2.95897980754556447425507657958878991881629488154278
```

Listable

```mathematica
In[14]:= Hypergeometric2F1[1, 1, 2, {0.5, -0.5}]
Out[14]= {1.38629, 0.81093}
```

Exact Gauss sum at z = 1 is not applied: stays unevaluated

```mathematica
In[15]:= Hypergeometric2F1[1/2, 1/3, 2, 1]
Out[15]= HypergeometricPFQ[{1/2, 1/3}, {2}, 1]
```

## Implementation notes

**Algorithm.** `builtin_hypergeometric_2f1` (`hypergeopfq.c`) first tries the
machine gate `hg_machine_head` — every argument a machine-scalar leaf
(`sf_machine_scalar_args`), at least one a machine `Real` — and otherwise rewrites
to `HypergeometricPFQ[{a, b}, {c}, z]`, inheriting its exact steps (cancellation,
`z == 0`, termination for non-positive integer `a`/`b`, the closed forms
`-Log[1-z]/z`, `(1-z)^-a` and the `ArcSin[Sqrt[z]]` family). The real machine
kernel `sf_machine_hyper2f1` builds a `Hy2Plan` from the parameters (polynomial,
`(1-z)^-b` / `(1-z)^-a`, or generic) and evaluates: the direct series on
`(0, 0.9]`, Pfaff's `(1-z)^-a F(a, c-b; c; z/(z-1))` on `[-9, 0)`, Gauss's sum
`G(c)G(c-a-b)/(G(c-a)G(c-b))` at `z = 1`, and the shared complex core
`hy2f1_core` elsewhere. That core picks the cheapest of direct / Pfaff series
(`|w| <= 0.75`), the `1-z` and `1/(1-z)` connection formulas (`|w| <= 0.6`,
accepted only while the combination's condition number is at most 2, with the
degenerate logarithmic cases when `c-a-b` is an integer), and — around
`z = e^(±i π/3)` or for ill-conditioned parameters — continuation of the
hypergeometric ODE by Taylor steps from `z0 = 1/2 ± 0.45 i` (`hy2f1_continue`).
A machine-real `z > 1` returns the complex value from below the cut (`zi = -0`).
Complex parameters use the series only, for `|z| <= 0.9`. At arbitrary precision
`|z| < 1` goes to arb's `acb_hypgeom_2f1` (`flint_num_hyp2f1` in
`src/flint_num_bridge.c`), with the MPFR series as fallback.

**Data structures.** Plain C99 `double`, complex values as re/im pairs, with
double-double (`hdd`) re-sums for a series that loses more than one bit or peaks
late; `Hy2Plan` / `Hy2Tab` hold the coefficient recurrence
`c_{k+1} = c_k (p0+k)(p1+k)/((q0+k)(k+1))`, which a buffer tabulates once. When
the terms are one-signed and never grow ("lean", decided from the parameters)
`sf_machine_pfq_vec` sums four buffer elements in lock-step, bit-identical to the
scalar loop. ND: n-ary kernel `NDKN_Hypergeometric2F1` (arity 4, `REG_N`); on
`pack.c`'s `AWARE` list with the family; a buffer containing a cut element comes
back complex. Attributes: `Listable`, `NumericFunction`, `Protected`.

**Complexity / limits.** Scalar about 0.3 µs real and 1.7 µs complex; 10^6 packed
elements about 0.14 s. Term caps 20000 (double) and 6000 (double-double); a
case where nothing converges declines to the exact / MPFR route. The arbitrary
precision path still requires `|z| < 1`, so an exact or MPFR `z` outside the unit
disc stays unevaluated, as does the exact Gauss value at `z = 1`.
`Compile[]` lowers at real and complex scalar shapes (through
`HypergeometricPFQ`'s `OP_KERNN` / `OP_KERNN_C`); the rank-1 array shape reports
`Compiled -> False`, while the packed / `NDArray` path still serves the REPL.

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## References

- DLMF §15 — the Gauss hypergeometric function 2F1.
- M. Abramowitz and I. A. Stegun, Handbook of Mathematical Functions (1964), 15.3.3-15.3.12 — the linear transformations and the degenerate logarithmic cases.
- R. C. Forrey, Computing the hypergeometric function, J. Comput. Phys. 137 (1997) 79-100.
- N. Michel and M. V. Stoitsov, Fast computation of the Gauss hypergeometric function with all its parameters complex, Comput. Phys. Commun. 178 (2008) 535-551.
- W. Bühring, An analytic continuation of the hypergeometric series, SIAM J. Math. Anal. 18 (1987) 884-889.
- Source: [`src/special_functions/hypergeopfq.c`](https://github.com/stblake/mathilda/blob/main/src/special_functions/hypergeopfq.c)
- Specification: [`docs/spec/builtins/special-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/special-functions.md)
- Tests: [`tests/test_compile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile.c)
- Tests: [`tests/test_compile_coverage.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile_coverage.c)
- Tests: [`tests/test_compiledfunction.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compiledfunction.c)
- Tests: [`tests/test_dsolve.c`](https://github.com/stblake/mathilda/blob/main/tests/test_dsolve.c)

## Notes & additional examples

### Notes

- `Hypergeometric2F1[a, b, c, z]` is the Gauss hypergeometric function, equal to
  `HypergeometricPFQ[{a, b}, {c}, z]`; an unevaluated call comes back in the
  `HypergeometricPFQ` form.
- At machine precision it is evaluated on the whole complex plane (real
  parameters): direct and Pfaff series, the `1-z` and `1/(1-z)` connection
  formulas, and continuation of the hypergeometric ODE where all of these are
  slow. On the cut `z > 1` a real argument gives the value from below, the same
  convention as Mathematica. With complex parameters only `|z| <= 0.9` is
  evaluated.
- Arbitrary precision is limited to `|z| < 1`, where arb's `acb_hypgeom_2f1`
  is used. An exact or high-precision `z` outside the disc stays unevaluated, and
  so does the exact value at `z = 1`. A machine `z = 1.` uses Gauss's sum.
- Attributes `Listable`, `NumericFunction`, `Protected`. Packed arrays and
  `NDArray` arguments map the machine kernel over the buffer (the coefficient
  table is built once per buffer); a buffer with an element on the cut comes back
  complex. `Compile[]` lowers real and complex scalars; rank-1 array arguments are
  not lowered yet.
- See also: `HypergeometricPFQ`, `Hypergeometric1F1`, `LegendreP`, `EllipticK`.
