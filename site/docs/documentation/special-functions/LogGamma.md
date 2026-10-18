# LogGamma

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`LogGamma[z]`**

gives the log-gamma function log(Gamma(z)), analytic except for a branch

<details>
<summary>Notes</summary>

cut on the negative reals. Exact at integer and half-integer z (with the negative-axis branch term), divergent (Infinity) at non-positive integers, and evaluated numerically for real and complex z at machine or arbitrary (MPFR) precision. D\[LogGamma\[z\], z\] is PolyGamma\[0, z\]. Listable.

</details>

## Examples (12)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= LogGamma[-3/2]
Out[1]= -(2*I) Pi + Log[4/3 Sqrt[Pi]]

In[2]:= LogGamma[-2.5]
Out[2]= -0.0562437 - 9.42478*I

In[3]:= LogGamma[2. + 30. I]
Out[3]= -41.1026 + 74.356*I
```

### Applications (9)

Exact at positive integers: Log[4!]

```mathematica
In[4]:= LogGamma[5]
Out[4]= Log[24]
```

Negative half-integers carry the branch term -Ceiling[-z] Pi I

```mathematica
In[5]:= LogGamma[-3/2]
Out[5]= -(2*I) Pi + Log[4/3 Sqrt[Pi]]
```

A negative machine real gives the complex value from above the cut

```mathematica
In[6]:= LogGamma[-2.5]
Out[6]= -0.0562437 - 9.42478*I
```

Finite where Gamma[1000.] would leave the double range

```mathematica
In[7]:= LogGamma[1000.]
Out[7]= 5905.22
```

Arbitrary precision

```mathematica
In[8]:= N[LogGamma[100], 40]
Out[8]= 359.13420536957539877604401046028690961264
```

The continuous branch: the imaginary part keeps growing where Log[Gamma] wraps

```mathematica
In[9]:= {LogGamma[2. + 30. I], Log[Gamma[2. + 30. I]]}
Out[9]= {-41.1026 + 74.356*I, -41.1026 - 1.04221*I}
```

The two differ by a multiple of 2 Pi I

```mathematica
In[10]:= LogGamma[-4.5 + 3. I] - Log[Gamma[-4.5 + 3. I]]
Out[10]= 0.0 - 12.5664*I
```

Arbitrary-precision complex through Stirling's series

```mathematica
In[11]:= N[LogGamma[1 + I], 30]
Out[11]= -0.6509231993018563388852168315042 - 0.3016403204675331978875316577968*I
```

The derivative is the digamma function

```mathematica
In[12]:= D[LogGamma[z], z]
Out[12]= PolyGamma[0, z]
```

## Algorithm

Mathilda -- LogGamma[z], the log-gamma function log(Gamma(z)).

LogGamma is analytic throughout the complex plane except for a branch cut on the negative real axis. It is the analytic continuation of log(Gamma(z)) and is *not* Log[Gamma[z]] (which has a more involved branch structure). Because it never overflows where Gamma does, it is the natural primitive for factorial ratios and asymptotics.

Evaluation is layered, mirroring Gamma (gamma.c), so each kind of argument takes the cheapest exact or fastest numeric route:

```text
  exact integer n >= 1          ->  Log[(n-1)!]   (exact, via Gamma's machinery)
  exact half-integer            ->  Log of the exact Sqrt[Pi] form, with the
                                    branch term -Ceiling[-z] Pi I for z < 0
  non-positive integer          ->  Infinity      (pole)
  symbolic infinities           ->  Infinity / Indeterminate / ComplexInfinity
  machine real z > 0            ->  libm   lgamma
  machine real z < 0 (non-int)  ->  lgamma + branch term (complex result)
  arbitrary real                ->  MPFR   mpfr_lgamma  (+ branch term)
  machine complex               ->  sfg_loggamma_c (sf_machine_gamma.c), shared
                                    with the NDArray kernel and Compile[]
  arbitrary complex             ->  Spouge log-gamma   (MPFR, runtime coeffs)
  everything else               ->  stays symbolic (return NULL)
```

The imaginary branch term on the negative real axis is Im = -Pi Ceiling[-z], taken from above (Mathematica's convention); this is exact for the symbolic half-integer path and is reused by the numeric real paths.

Attributes: Listable, NumericFunction, Protected.

## Implementation notes

**Algorithm.** `builtin_loggamma` evaluates the analytic continuation of `log Γ(z)`. This is the continuous branch, with one cut on the negative real axis, and it is not `Log[Gamma[z]]`. `loggamma_one_arg` checks the cases in this order:

1. **Machine real.** For `x > 0` it uses libm `lgamma`; an overflow leaves the call unevaluated. At a non-positive integer it returns `ComplexInfinity`. For `x < 0` it returns `lgamma(x) - i π Ceiling[-x]`, the value approached from above.
2. **Complex with an MPFR part** → `lg_mpfr_complex`: Stirling's series in the file-local `lcx` toolkit. The argument is shifted up to `Re ~ K/3` with `K ~ wp/2` exact `mpq` Bernoulli terms, and the reflection uses a continued `log sin(π z)`.
3. **Machine complex** → `loggamma_machine_complex`, which calls `sfg_loggamma_c` in `sf_machine_gamma.c`. Taylor tiles about 30 centres in the strip `1 <= Re < 5/2` cover `-5 <= Re z < 10`, `|Im z| < 19/4`. On the real-axis row the tiled series is for `lnΓ/((z-1)(z-2))`, so both zeros are carried exactly. Elsewhere the kernel uses Stirling's series after shifting to `Re >= 7`. The shift's `N` logarithms are summed as one log of their product, with the `2πi` windings counted as the product crosses the negative axis (Hare's continuous branch). The kernel declines `Re z < -30`; those inputs fall back to a `g = 7` Lanczos log-gamma with a winding-correct reflection.
4. **Exact** (`lg_exact`). A positive integer gives `Log[(n-1)!]`. A half-integer gives the `Log` of the exact `Sqrt[Pi]` form, plus `-Ceiling[-z] π I` when `z < 0`. A non-positive integer gives `Infinity`.
5. **Symbolic infinities.** `Infinity → Infinity`, `-Infinity → Indeterminate`, and `ComplexInfinity` or a directed complex infinity → `ComplexInfinity`.
6. **Arbitrary-precision real** → `mpfr_lgamma` (which gives `log|Γ|` and its sign), plus the same branch term.

An `Interval` argument routes through `interval_apply_function`.

**Data structures.** `Expr`; `lcx` (a pair of `mpfr_t`) with an exact `mpq` Bernoulli generator for the arbitrary-precision complex path; `(re, im)` doubles in the machine kernel, with the tile constants stored as double-doubles. The ND kernel is `NDKU_LogGamma = { loggamma_machine_complex, ndk_LogGamma_r }`, registered with `REG_U`, so it is packed-aware. The real arm uses `lgamma` and declines `x <= 0`, where the value is complex; a declined buffer drops to the `List` path. The complex arm is the same continued log-gamma the scalar call uses. `Compile[]` lowers `LogGamma` through that kernel for real and complex scalars and for rank-1 arrays, and a compiled complex call agrees bit for bit with the interpreter.

**Complexity / limits.** The machine paths are `O(1)`. The complex kernel has a median error of 0.5 ulp and a maximum of 3.3 ulp on `0.5 <= Re z <= 4`, `|Im z| <= 3`; the Lanczos series it replaced was off by up to 689 ulp near the zeros at 1 and 2. The arbitrary-precision Stirling path costs more as the precision rises. The imaginary part grows without bound as `|Im z|` grows, which is what separates this function from `Log[Gamma[z]]`. One known inconsistency at the poles: the exact path returns `Infinity`, but a machine pole (`LogGamma[0.]`) returns `ComplexInfinity`. Attributes: `Listable`, `NumericFunction`, `Protected`.

- **Exact closed forms.** Integers reduce as `LogGamma[n] = Log[(n-1)!]`
  (`LogGamma[5] = Log[24]`); positive half-integers give `Log` of the exact
  `Sqrt[Pi]` form; negative half-integers carry the branch term
  `-Ceiling[-z] Pi I`, e.g. `LogGamma[-3/2] = -2 I Pi + Log[(4 Sqrt[Pi])/3]`.
- **Poles.** Non-positive integers diverge: `LogGamma[0] = LogGamma[-1] = … =
  Infinity`.
- **Symbolic infinities.** `LogGamma[Infinity] = Infinity`,
  `LogGamma[-Infinity] = Indeterminate`, `LogGamma[I Infinity] =
  LogGamma[ComplexInfinity] = ComplexInfinity`.
- **Branch continuity across the reflection.** For `Re[z] < 1/2` the
  arbitrary-precision path (and the machine fallback, `Re z < -30`) use
  `LogGamma[z] = Log[Pi] - Log[Sin[Pi z]] - LogGamma[1-z]`, where the
  log of the sine must be the *continued* one. Factoring
  `Sin[Pi z] = E^(-I Pi z) (E^(2 I Pi z) - 1)/(2 I)` puts the whole winding in
  the exact `-I Pi z` term and leaves a factor near `-1` whose principal log is
  safe. At `Im[z] = 0` this is the limit from above, matching the
  `Im = -Ceiling[-z] Pi` convention used on the real axis. Before 2026-07-27 a
  principal `Log[Pi/Sin[Pi z]]` was used instead: correct only in the strip
  `-1 < Re[z] < 0`, and short by a multiple of `2 Pi I` outside it, so
  `LogGamma[-4.5 + 3. I]` equalled `Log[Gamma[-4.5 + 3. I]]` rather than the
  continuation.
- **Numerics.** Machine real via `lgamma`; arbitrary-precision real via MPFR
  `lgamma`; machine complex via `sfg_loggamma_c` (`sf_machine_gamma.c`, shared
  with the packed / `NDArray` maps and `Compile[]`): Taylor tiles about 30
  centres in the strip `1 ≤ Re < 5/2` reached by the exact recurrence
  (`-5 ≤ Re z < 10`, `|Im z| < 19/4`; on the real-axis row the series is for
  lnΓ/((z−1)(z−2)), so both zeros are carried exactly), shifted Stirling
  elsewhere. The recurrence's N logs are one log of the product, with the
  2πi windings counted as the running product crosses the negative real
  axis — the same continuous branch (Hare 1997). Median 0.5 / max 3.3 ulp on
  `0.5 ≤ Re z ≤ 4, |Im z| ≤ 3` (the Lanczos series it replaces: median 7 /
  max 689 ulp, worst near the zeros at 1 and 2). `Re z < -30` falls back to the
  Lanczos/reflection path described above; arbitrary-precision
  complex via a Stirling series with argument reduction (the continuous
  branch — its imaginary part grows past π where `Log[Gamma]` would wrap).
  Negative real arguments return the complex value `log|Γ(z)| - Pi Ceiling[-z] I`.
- **Derivative.** `D[LogGamma[z], z] = PolyGamma[0, z]`; higher derivatives
  raise the `PolyGamma` order. Produced by `PolyGamma[-1, z]`.

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## References

**See also:** [Gamma](../../special-functions/Gamma/), [Log](../../elementary-functions/Log/), [NDArray](../../linear-algebra/NDArray/), [PolyGamma](../../special-functions/PolyGamma/)

- DLMF §5.11 — Stirling's series for ln Gamma.
- D. E. G. Hare, *Computing the principal branch of log-Gamma*, J. Algorithms 25 (1997) 221–236.
- C. Lanczos, *A precision approximation of the gamma function*, SIAM J. Numer. Anal. Ser. B 1 (1964) 86–96.
- T. J. Dekker, *A floating-point technique for extending the available precision*, Numer. Math. 18 (1971) 224–242.
- Source: [`src/special_functions/loggamma.c`](https://github.com/stblake/mathilda/blob/main/src/special_functions/loggamma.c)
- Specification: [`docs/spec/builtins/special-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/special-functions.md)
- Tests: [`tests/test_autocompile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_autocompile.c)
- Tests: [`tests/test_compile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile.c)
- Tests: [`tests/test_fullsimplify.c`](https://github.com/stblake/mathilda/blob/main/tests/test_fullsimplify.c)
- Tests: [`tests/test_gamma_machine.c`](https://github.com/stblake/mathilda/blob/main/tests/test_gamma_machine.c)

## Notes & additional examples

### Notes

- `LogGamma[z]` is the analytic continuation of `log Γ(z)` with a single branch cut on the negative real axis — not `Log[Gamma[z]]`, whose imaginary part is confined to `(-Pi, Pi]`.
- Exact integers give `Log[(n-1)!]`; half-integers give the log of the exact `Sqrt[Pi]` form (plus the branch term below zero); non-positive integers give `Infinity`.
- Machine complex arguments use the `sf_machine_gamma.c` kernel (Hare's winding-counted continuous branch), shared bit for bit with packed lists, `NDArray` and `Compile[]`; machine reals use libm `lgamma`.
- A machine pole (`LogGamma[0.]`) currently answers `ComplexInfinity`, while the exact pole `LogGamma[0]` gives `Infinity`.
- Attributes: `Listable`, `NumericFunction`, `Protected`. See also `Gamma`, `PolyGamma` (`PolyGamma[-1, z]` is `LogGamma[z]`).
