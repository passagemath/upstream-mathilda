# Hypergeometric1F1

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Hypergeometric1F1[a, b, z]`**

is Kummer's confluent hypergeometric 1F1, equal to HypergeometricPFQ\[{a}, {b}, z\].

## Examples (14)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= Hypergeometric1F1[a, b, 0]
Out[1]= 1

In[2]:= Hypergeometric1F1[-2, 1, z]
Out[2]= 1 - 2 z + 1/2 z^2

In[3]:= Hypergeometric1F1[1., 2., -40.]
Out[3]= 0.025

In[4]:= Hypergeometric1F1[1 + I, 2, 0.5]
Out[4]= 1.27058 + 0.322122*I
```

### Applications (10)

The elementary case (E^z - 1)/z

```mathematica
In[5]:= Hypergeometric1F1[1, 2, z]
Out[5]= (-1 + E^z)/z
```

Equal parameters cancel

```mathematica
In[6]:= Hypergeometric1F1[a, a, z]
Out[6]= E^z
```

Terminates to a Laguerre polynomial

```mathematica
In[7]:= Hypergeometric1F1[-2, 1, z]
Out[7]= 1 - 2 z + 1/2 z^2
```

40 digits from the series

```mathematica
In[8]:= N[Hypergeometric1F1[1, 2, 3], 40]
Out[8]= 6.3618456410625559136428432181939059656629
```

The closed form agrees

```mathematica
In[9]:= N[(E^3 - 1)/3, 40]
Out[9]= 6.3618456410625559136428432181939059656629
```

An alternating series, kept accurate by Kummer's transformation

```mathematica
In[10]:= Hypergeometric1F1[1., 2., -40.]
Out[10]= 0.025
```

Large argument, the asymptotic expansion

```mathematica
In[11]:= Hypergeometric1F1[1.5, 2.5, 100.]
Out[11]= 4.01191e+41
```

Past the double range: promoted to MPFR

```mathematica
In[12]:= Hypergeometric1F1[1., 2., 900.]
Out[12]= 8.143126913674914e+387
```

Complex parameter

```mathematica
In[13]:= Hypergeometric1F1[1 + I, 2, 0.5]
Out[13]= 1.27058 + 0.322122*I
```

Arbitrary precision for negative argument

```mathematica
In[14]:= N[Hypergeometric1F1[1/2, 3/2, -30], 30]
Out[14]= 0.1618021593796400696905131930601
```

## Implementation notes

**Algorithm.** `builtin_hypergeometric_1f1` first tries the machine gate
`hg_machine_head` (every argument a machine-scalar leaf, at least one a machine
`Real`), which calls the raw-`double` kernel `sf_machine_hyper1f1` /
`sf_machine_hyper1f1_c` in `src/special_functions/sf_machine_hyper.c`; otherwise
it rewrites to `HypergeometricPFQ[{a}, {b}, z]` and inherits that engine's exact
steps (cancellation, so `a == b` gives `Exp[z]`; `z == 0 -> 1`; termination to a
polynomial for non-positive integer `a`; `1F1(1;2;z) = (E^z - 1)/z`) and its
arbitrary-precision MPFR series. The machine kernel tries, in order: the closed
forms for `z = 0` or `a = 0`, `a = b` (`exp`) and `(a, b) = (1, 2)` (`expm1(z)/z`,
with a complex form free of cancellation at small `|z|`); the terminating
polynomial for non-positive integer `a` (also past a lower-parameter pole it never
reaches); for `|z| > 20`, DLMF 13.7.2 with both exponential sectors kept (the
subdominant one with the Stokes-line `cos(π a)` factor), accepted only when its
least term is below half an ulp and the two sectors do not cancel; for
`Re z < 0`, Kummer's transformation `e^z 1F1(b-a; b; -z)` with `b - a` carried in
double-double, so the terms are one-signed; otherwise the plain series. A machine
overflow (e.g. `Hypergeometric1F1[1., 2., 900.]`) is promoted to MPFR by the
`PFQ` builtin.

**Data structures.** Plain C99 `double` with complex values as re/im pairs (no
`<complex.h>`); every series measures its own loss (peak term against the sum)
and re-sums in double-double (`hdd`) when it loses more than a bit. ND: n-ary
kernel `NDKN_Hypergeometric1F1` (arity 3, `REG_N`); a packed or `NDArray` `z`
also runs through `PFQ`'s buffer driver `sf_machine_pfq_vec`, the same kernel as
the scalar call, so the surfaces agree bit for bit. On `pack.c`'s `AWARE` list
with the family. Attributes: `Listable`, `NumericFunction`, `Protected`.

**Complexity / limits.** Entire in `z`. Scalar about 0.15 µs real and 0.35 µs
complex; 10^6 packed elements of `1F1(1;2;x)` about 11 ms. The asymptotic branch
needs `|a|, |b|, |b-a| <= 160`; a declined machine point falls to the exact /
MPFR route. A symbolic call is returned in its `HypergeometricPFQ` form.
`Compile[]` lowers at real and complex scalar shapes (via `HypergeometricPFQ`'s
`OP_KERNN` / `OP_KERNN_C`); the rank-1 array shape reports `Compiled -> False`,
while the buffer path still serves the REPL.

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## References

- DLMF §13 — confluent hypergeometric functions (Kummer's 1F1); 13.7.2 for the large-|z| expansion.
- M. Abramowitz and I. A. Stegun, Handbook of Mathematical Functions (1964), 13.1.27 — Kummer's transformation.
- T. J. Dekker, A floating-point technique for extending the available precision, Numer. Math. 18 (1971) 224-242.
- Source: [`src/special_functions/hypergeopfq.c`](https://github.com/stblake/mathilda/blob/main/src/special_functions/hypergeopfq.c)
- Specification: [`docs/spec/builtins/special-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/special-functions.md)
- Tests: [`tests/test_compile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile.c)
- Tests: [`tests/test_dsolve.c`](https://github.com/stblake/mathilda/blob/main/tests/test_dsolve.c)
- Tests: [`tests/test_dsolve_m19_stress.c`](https://github.com/stblake/mathilda/blob/main/tests/test_dsolve_m19_stress.c)
- Tests: [`tests/test_hypergeopfq.c`](https://github.com/stblake/mathilda/blob/main/tests/test_hypergeopfq.c)

## Notes & additional examples

### Notes

- `Hypergeometric1F1[a, b, z]` is Kummer's confluent hypergeometric function
  `M(a, b, z)`, equal to `HypergeometricPFQ[{a}, {b}, z]`; it is entire in `z`.
  An unevaluated call comes back as `HypergeometricPFQ[{a}, {b}, z]`.
- Machine-precision calls, real or complex, go to one raw-`double` kernel: closed
  forms, the terminating polynomial, the DLMF 13.7.2 expansion for `|z| > 20`,
  Kummer's transformation for `Re z < 0`, else the series, re-summed in
  double-double when it cancels. The same kernel serves packed arrays, `NDArray`
  arguments and `Compile[]`, so all of these agree bit for bit.
- A result outside the double range is recomputed in MPFR instead of
  overflowing.
- Attributes `Listable`, `NumericFunction`, `Protected`. `Compile[]` lowers real
  and complex scalars; rank-1 array arguments are not lowered yet.
- See also: `HypergeometricPFQ`, `HypergeometricU`, `WhittakerM`, `LaguerreL`,
  `Erf`.
