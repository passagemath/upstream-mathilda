# PolyGamma

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`PolyGamma[z]`**

gives the digamma function psi(z) (rewritten as PolyGamma\[0, z\]).

**`PolyGamma[n, z]`**

gives the n-th derivative of the digamma function, psi^(n)(z).

**`LogGamma[z]. Listable.`**

<details>
<summary>Notes</summary>

Positive-integer arguments reduce to exact values: psi(m) to a rational minus EulerGamma, and psi^(n)(m) for odd n to a rational plus a rational multiple of Pi^(n+1); even orders stay symbolic. Non-positive integer arguments give ComplexInfinity. Inexact real and complex arguments evaluate numerically at machine or arbitrary (MPFR) precision. PolyGamma\[-1, z\] gives

</details>

## Examples (16)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= PolyGamma[5]
Out[1]= 25/12 - EulerGamma

In[2]:= PolyGamma[0, 1/2]
Out[2]= -EulerGamma - Log[4]

In[3]:= PolyGamma[3, 2.5 + I]
Out[3]= 0.0415431 - 0.163057*I
```

### Worked examples (4)

```mathematica
In[4]:= PolyGamma[n, -k]
Out[4]= PolyGamma[n, -k]

In[5]:= PolyGamma[0, Infinity]
Out[5]= Infinity

In[6]:= PolyGamma[n, Infinity]
Out[6]= PolyGamma[n, Infinity]

In[7]:= PolyGamma[n, Indeterminate]
Out[7]= PolyGamma[n, Indeterminate]
```

### Applications (9)

The digamma at a positive integer: a harmonic number minus EulerGamma

```mathematica
In[8]:= PolyGamma[5]
Out[8]= 25/12 - EulerGamma
```

Odd orders at integers close through zeta at even arguments

```mathematica
In[9]:= {PolyGamma[1, 1], PolyGamma[3, 1]}
Out[9]= {1/6 Pi^2, 1/15 Pi^4}
```

The Gauss digamma closed form at a half-integer

```mathematica
In[10]:= PolyGamma[0, 1/2]
Out[10]= -EulerGamma - Log[4]
```

A machine argument is computed directly by the shared kernel

```mathematica
In[11]:= PolyGamma[2.5]
Out[11]= 0.703157
```

Relatively accurate right at the positive zero of the digamma

```mathematica
In[12]:= PolyGamma[1.461632144968362]
Out[12]= -3.07279e-16
```

Machine complex, order 3

```mathematica
In[13]:= PolyGamma[3, 2.5 + I]
Out[13]= 0.0415431 - 0.163057*I
```

Arbitrary-precision complex through the Bernoulli asymptotic

```mathematica
In[14]:= N[PolyGamma[2, 1 + I], 30]
Out[14]= 0.3685529315879351717366345429807 + 0.7666528503450662124026953776316*I
```

Order -1 is the log-gamma function

```mathematica
In[15]:= PolyGamma[-1, z]
Out[15]= LogGamma[z]
```

Even orders at integers, and non-half-integer rationals, stay symbolic

```mathematica
In[16]:= {PolyGamma[2, 1], PolyGamma[0, 1/4]}
Out[16]= {PolyGamma[2, 1], PolyGamma[0, 1/4]}
```

## Algorithm

Mathilda -- the PolyGamma function family.

```text
  PolyGamma[z]      digamma psi(z) = Gamma'(z)/Gamma(z)
                      (always rewrites to the two-argument form PolyGamma[0, z])
  PolyGamma[n, z]   n-th polygamma psi^(n)(z) = d^n/dz^n psi(z)
```

Evaluation is layered so each kind of argument takes the cheapest exact or fastest numeric route available:

```text
  order n = -1                -> LogGamma[z]   (inert; psi^(-1) = log-gamma)
  z a non-positive integer    -> ComplexInfinity (pole of every psi^(n))
  z a positive integer        -> exact closed form:
                                   n = 0          : H_{z-1} - EulerGamma
                                   n >= 1 odd     : rational + rational*Pi^(n+1)
                                                    (via zeta(n+1), even -> Bernoulli)
                                   n >= 1 even    : stays symbolic (zeta(odd), no
                                                    closed form)
  machine real / complex z,   -> sfg_polygamma / sfg_polygamma_c (sf_machine_gamma.c),
  integer 0 <= n <= 40           shared with the NDArray kernel and Compile[];
                                 the one-argument form computes directly
  z inexact real (Real/MPFR)  -> numeric: mpfr_digamma for n = 0, otherwise the
                                   recurrence-shift + Bernoulli asymptotic series
  z inexact complex           -> the same asymptotic, in complex arithmetic
  everything else             -> stays symbolic (return NULL)
```

MPFR provides a digamma (n = 0) but no higher polygamma and no Hurwitz zeta, so the n >= 1 numeric kernel and the complex digamma are implemented here from the classical asymptotic expansion (Abramowitz & Stegun 6.4.11):

```text
  psi^(n)(w) ~ (-1)^(n-1) [ (n-1)!/w^n + n!/(2 w^(n+1))
                            + Sum_{k>=1} B_{2k} (2k+n-1)!/(2k)! w^(-(2k+n)) ]   (n>=1)
  psi(w)     ~ ln w - 1/(2w) - Sum_{k>=1} B_{2k}/(2k) w^(-2k)                   (n=0)
```

valid for large |w|; the recurrence psi^(n)(z) = psi^(n)(z+1) - (-1)^n n! z^(-(n+1)) shifts a small or negative argument up into that regime.

Attributes: Listable, NumericFunction, Protected.

## Performance

Against other systems, from the benchmark suite (same input, results cross-checked for agreement):

| case | Mathilda | Wolfram | Python |
|---|---:|---:|---:|
| BesselJ[0, .] over 10^6 | 3.46e+03 s | 1.74e+03 s | 53 s |
| Zeta over 10^6 | 8.17 s | 4.31e+03 s | 5.7 s |
| AiryAi over 10^6 | 3.14 s | 107 s | 58.9 s |
| PolyGamma[0, .] over 10^6 | 1.55 s | 151 s | 8.19 s |
| Gamma over 10^6 | 1.16 s | 1.35 s | 7.34 s |
| Erf over 10^6 | 0.968 s | 1.23 s | 7.43 s |

## Implementation notes

**Algorithm.** `builtin_polygamma` evaluates `PolyGamma[n, z] = ψ⁽ⁿ⁾(z)`. It first tries `pg_machine`, the machine fast path. That path applies when `sf_machine_scalar_args` accepts the call and the order is an integer `0 <= n <= 40` (`SFG_POLYGAMMA_MAX_N`); the one-argument form counts as `n = 0`. It dispatches to `sf_machine_gamma.c`:

- **Real digamma** (`sfg_digamma`). On `[1, 2]`, a degree-22 Chebyshev interpolant of `ψ(x)/(x - x₀)` is multiplied back by `(x - x₀)`, with the zero `x₀ = 1.4616…` held as a double-double, so the answer stays relatively accurate right through the zero (the Cody–Strecok–Thacher technique). From `(0, 10)` the kernel reaches `[1, 2]` by integer shifts, which are exact in double. For `x >= 10` it uses the asymptotic series. For `x < 0` it uses reflection with an exactly reduced `π cot(π x)`.
- **Real ψ⁽ⁿ⁾, n >= 1, x > 0** (`sfg_polygamma`). Upward recurrence to `x >= 10 + n/2`, then a Hurwitz-zeta asymptotic tail. Every term has one sign, each gets a first-order correction for the rounding of `x + j`, and the sum is compensated.
- **Complex** (`sfg_polygamma_c`). For ψ, the kernel uses the Taylor tiles shared with `Gamma`/`LogGamma` (`-5 <= Re z < 10`, `|Im z| < 19/4`), a Taylor series about `x₀` within `0.35`, the asymptotic series after shifting to `Re z >= 10`, and reflection when `Re z < -5`. For `n >= 1`, it uses the recurrence plus the tail on `Re z >= 0`.

A machine `PolyGamma[x]` is answered directly and is identical to `PolyGamma[0, x]`. If the kernel declines (a pole, or `n >= 1` with `Re z < 0`), the layered path runs:

1. `PolyGamma[z]` rewrites to `PolyGamma[0, z]`.
2. `n = -1` gives the inert `LogGamma[z]`; `n <= -2` stays unevaluated.
3. Special values: `Infinity` gives `Infinity` (n = 0) or `0` (n >= 1). `ComplexInfinity` and `Indeterminate` give `Indeterminate`. A non-positive integer `z` is a pole, `ComplexInfinity`.
4. Positive-integer `z`: `n = 0` gives `H_{z-1} - EulerGamma`. Odd `n` gives a rational plus a rational multiple of `π^(n+1)` (through ζ(n+1) and Bernoulli numbers). Even `n` stays symbolic.
5. `n = 0` at a positive half-integer: the Gauss digamma closed form, e.g. `-EulerGamma - Log[4]` at `1/2`. Other rationals stay unevaluated.
6. Numeric: an inexact real uses `mpfr_digamma` for `n = 0`, otherwise a recurrence shift plus the A&S 6.4.11 Bernoulli asymptotic in MPFR. An inexact complex `z` uses the same asymptotic in the `pcx` toolkit.

An `Interval` argument with an integer order goes to `interval_polygamma`.

**Data structures.** `Expr`; a process-lifetime exact `mpq` Bernoulli cache; `pcx` (a pair of `mpfr_t`) for the complex asymptotic; `(re, im)` doubles and generated coefficient tables (`tools/sf_accuracy/gen_gamma_tables.py`) in the machine kernel. There are two ND kernels, both packed-aware. The unary `NDKU_PolyGamma = { ndk_PolyGamma_u, ndk_PolyGamma_r }` runs the real digamma or `sfg_polygamma_c(0, ·)`. The binary `NDKB_PolyGamma` runs `sfg_polygamma_c(n, ·)` and declines a non-integer, negative, or `> 40` order. `Compile[]` lowers `PolyGamma[x]` and `PolyGamma[n, x]` through these kernels, for real and complex scalars and for rank-1 arrays. The scalar, packed, `NDArray` and `Compile[]` answers agree bit for bit.

**Complexity / limits.** The machine kernels are `O(1)`, about 0.14 µs real and 0.35 µs complex. Real ψ has a median error of 0.5 ulp and a maximum of 3.5 ulp; complex has a median of 0.7 ulp, with maxima beside ψ's real zero. The MPFR asymptotic needs the argument shifted up to `Re ~ 0.13 wp + n`, so its cost grows with the precision. Limits: exact argument cap `100000`, exact order cap `256`, numeric order cap `1024`. Complex or fractional orders stay unevaluated. Attributes: `Listable`, `NumericFunction`, `Protected`.

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## References

**See also:** [EulerGamma](../../mathematical-constants/EulerGamma/), [NDArray](../../linear-algebra/NDArray/), [Gamma](../../special-functions/Gamma/), [Series](../../power-series/Series/), [D](../../calculus/D/)

- DLMF §5.15 — the polygamma functions.
- M. Abramowitz & I. A. Stegun, Handbook of Mathematical Functions, §6.4.11 — the asymptotic expansion.
- W. J. Cody, A. J. Strecok and H. C. Thacher, *Chebyshev approximations for the psi function*, Math. Comp. 27 (1973) 123–127.
- DLMF §25.11 — the Hurwitz zeta tail used for psi^(n), n >= 1.
- Source: [`src/special_functions/polygamma.c`](https://github.com/stblake/mathilda/blob/main/src/special_functions/polygamma.c)
- Specification: [`docs/spec/builtins/special-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/special-functions.md)
- Tests: [`tests/test_beta.c`](https://github.com/stblake/mathilda/blob/main/tests/test_beta.c)
- Tests: [`tests/test_compile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile.c)
- Tests: [`tests/test_flint_bridge.c`](https://github.com/stblake/mathilda/blob/main/tests/test_flint_bridge.c)
- Tests: [`tests/test_fullsimplify.c`](https://github.com/stblake/mathilda/blob/main/tests/test_fullsimplify.c)

## Notes & additional examples

### Notes

- `PolyGamma[z]` is the digamma ψ(z) (`PolyGamma[0, z]`); `PolyGamma[n, z]` is ψ⁽ⁿ⁾(z). Non-positive integers are poles (`ComplexInfinity`); `PolyGamma[-1, z]` returns `LogGamma[z]`.
- Exact closed forms: ψ(m) = `H[m-1] - EulerGamma`, odd-order ψ⁽ⁿ⁾(m) as rational plus rational·`Pi^(n+1)`, and the Gauss digamma theorem at positive half-integers. `ζ(odd)` values (even order at integers) are left symbolic.
- Machine real or complex `z` with integer order `0 ≤ n ≤ 40` runs on the `sf_machine_gamma.c` kernels, shared bit for bit with packed lists, `NDArray` and `Compile[]` (both `PolyGamma[x]` and `PolyGamma[n, x]` compile); otherwise MPFR (`mpfr_digamma`, or the recurrence-shifted Bernoulli asymptotic).
- Complex or fractional orders stay unevaluated. Attributes: `Listable`, `NumericFunction`, `Protected`. See also `Gamma`, `LogGamma`, `HarmonicNumber`, `Zeta`.
