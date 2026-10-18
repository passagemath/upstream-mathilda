# Gamma

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Gamma[z]`**

is the Euler gamma function Gamma(z).

**`Gamma[a, z]`**

is the upper incomplete gamma function Gamma(a, z).

**`Gamma[a, z0, z1]`**

is the generalized incomplete gamma Gamma(a, z0) - Gamma(a, z1).

<details>
<summary>Notes</summary>

Integer and half-integer arguments reduce to exact values ((z-1)!, and rational multiples of Sqrt\[Pi\]); non-positive integers give ComplexInfinity. Machine and arbitrary-precision (MPFR) real inputs evaluate numerically, as do machine-precision complex inputs. Listable.

</details>

## Examples (17)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= Gamma[7/2]
Out[1]= 15/8 Sqrt[Pi]

In[2]:= Gamma[2.3 + I]
Out[2]= 0.719141 + 0.540614*I

In[3]:= Gamma[3, x]
Out[3]= E^(-x) (2 + 2 x + x^2)
```

### Worked examples (4)

```mathematica
In[4]:= Gamma[-n]
Out[4]= Gamma[-n]

In[5]:= Gamma[Infinity]
Out[5]= Infinity

In[6]:= Gamma[-Infinity]
Out[6]= Indeterminate

In[7]:= Gamma[ComplexInfinity]
Out[7]= ComplexInfinity
```

### Applications (10)

The shifted factorial, 4!

```mathematica
In[8]:= Gamma[5]
Out[8]= 24
```

Half-integers close to rational multiples of Sqrt[Pi]

```mathematica
In[9]:= Gamma[7/2]
Out[9]= 15/8 Sqrt[Pi]
```

The functional equation continues past the poles at the non-positive integers

```mathematica
In[10]:= {Gamma[-1/2], Gamma[0]}
Out[10]= {-2 Sqrt[Pi], ComplexInfinity}
```

Arbitrary-precision real through mpfr_gamma

```mathematica
In[11]:= N[Gamma[1/3], 40]
Out[11]= 2.6789385347077476336556929409746776441286
```

Machine complex on the shared double kernel

```mathematica
In[12]:= Gamma[2.3 + I]
Out[12]= 0.719141 + 0.540614*I
```

The compiled call and the scalar builtin agree bit for bit

```mathematica
In[13]:= cg = Compile[{{z, _Complex}}, Gamma[z]]; cg[2.3 + I] === Gamma[2.3 + I]
Out[13]= True
```

Arbitrary-precision complex through Spouge's approximation

```mathematica
In[14]:= N[Gamma[3 + 4 I], 20]
Out[14]= 0.00522553847136921419473 - 0.172547079294300187719*I
```

A positive-integer first argument gives the finite closed form

```mathematica
In[15]:= Gamma[3, x]
Out[15]= E^(-x) (2 + 2 x + x^2)
```

Upper and generalized incomplete gamma, numerically

```mathematica
In[16]:= {Gamma[1.5, 7.5], Gamma[1, 1.1, 2.2]}
Out[16]= {0.00160996, 0.222068}
```

An exact complex incomplete gamma has no closed form and stays unevaluated

```mathematica
In[17]:= Gamma[3/2, I]
Out[17]= Gamma[3/2, I]
```

## Algorithm

Mathilda -- the Gamma function family.

```text
  Gamma[z]          Euler gamma function   Gamma(z) = Int_0^Inf t^(z-1) e^-t dt
  Gamma[a, z]       upper incomplete gamma Gamma(a,z) = Int_z^Inf t^(a-1) e^-t dt
  Gamma[a, z0, z1]  generalized incomplete = Gamma(a,z0) - Gamma(a,z1)
```

Evaluation is layered so each kind of argument takes the cheapest exact or fastest numeric route available:

```text
  exact integer / half-integer  ->  (z-1)! via the Factorial machinery
                                     (exact, BigInt, or rational*Sqrt[Pi])
  machine real        ->  libm   tgamma
  machine complex     ->  sfg_gamma_c (sf_machine_gamma.c), shared with the
                          NDArray kernel and Compile[]
  arbitrary real      ->  MPFR   mpfr_gamma / mpfr_gamma_inc
  everything else     ->  stays symbolic (return NULL)
```

Arbitrary-precision *complex* gamma is deliberately left symbolic: a fixed-coefficient Lanczos series only carries ~15 correct digits, so emitting it as an MPFR value would advertise a precision it does not have. Reporting the input unevaluated is the honest behaviour.

Attributes: Listable, NumericFunction, Protected.

## Implementation notes

**Algorithm.** `builtin_gamma` serves `Gamma[z]`, `Gamma[a, z]` (upper incomplete) and `Gamma[a, z0, z1] = Gamma[a, z0] - Gamma[a, z1]`; an `Interval` argument goes to `interval_apply_function`. `gamma_one_arg` tests the cheap numeric cases first, then the exact ones:

1. **Machine real** → libm `tgamma`. A pole (`NaN`) gives `ComplexInfinity`. An overflow (`Gamma[172.]` and beyond) is recomputed by `mpfr_gamma` at 53 bits, so the answer is a machine-precision number with an extended exponent and not `Infinity`.
2. **Complex with an MPFR part** → `gamma_mpfr_complex`: Spouge's approximation in the file-local complex-MPFR toolkit `gcx`. Its coefficients are computed at run time for the working precision (output precision + 64 guard bits), so the requested digits are correct.
3. **Machine complex** → `gamma_machine_complex`, which calls `sfg_gamma_c` in `sf_machine_gamma.c`. For `-5 <= Re z < 10`, `|Im z| < 19/4`, the exact recurrence moves `z` into the strip `1 <= Re < 5/2`, and lnΓ comes from a Taylor series about the nearest of 30 precomputed centres. Elsewhere the kernel shifts to `Re w >= 7` and sums Stirling's series. The `N` recurrence logarithms are taken as one log of the running product, with its windings counted (Hare's continuous branch). Γ is then `exp` of a double-double exponent, built with Dekker products, times the product correction. A real-valued result (`Im = 0`) comes back as a `Real`. Inputs the kernel declines (`Re z < -30`, overflow) fall back to the old `g = 7`, 9-term Lanczos series with reflection.
4. **Exact integer / half-integer** → `(z-1)!` through the `Factorial` machinery. This gives an exact integer, a `BigInt`, or a rational multiple of `Sqrt[Pi]`, and works for negative half-integers too. Non-positive integers are poles (`ComplexInfinity`).
5. **Symbolic infinities**: `Infinity → Infinity`, `-Infinity → Indeterminate`, `ComplexInfinity → ComplexInfinity`.
6. **Arbitrary-precision real** → `mpfr_gamma` at the input precision.

The incomplete form does the exact rewrites first: `Gamma[a, 0] = Gamma[a]` and `Gamma[a, Infinity] = 0`. A positive-integer `a <= 1000` with exact or symbolic `z` expands to `e^{-z} Sum_k ((n-1)!/k!) z^k`, with the coefficients built in GMP (`gamma_incomplete_int`). For a numeric real it calls `mpfr_gamma_inc` at the minimum input precision, keeping a 53-bit MPFR value when a double would overflow or underflow. For numeric complex it uses `gcx_inc_gamma`: the lower series `gamma(a,z) = z^a e^{-z} Sum z^n/(a)_{n+1}` when `Re z < Re a + 1`, otherwise a Lentz continued fraction, then multiplies by `e^{-z} z^a`.

**Data structures.** `Expr` trees; `gcx` (a pair of `mpfr_t`; alias-safe, explicit precision) for the complex and incomplete MPFR work; GMP for the factorial and closed-form coefficients. The machine kernel works on `(re, im)` doubles, without `<complex.h>`. Its tables (tile Taylor coefficients, Stirling `B_2k/(2k(2k-1))`) come from `tools/sf_accuracy/gen_gamma_tables.py` (mpmath), and `tools/sf_accuracy/gamma.py` gates their accuracy. The ND kernel is `NDKU_Gamma = { gamma_machine_complex, ndk_Gamma_r }`, registered with `REG_U`, so packed arrays and visible `NDArray`s are mapped directly: real elements use `tgamma` (a pole or overflow declines to the `List` path), complex elements use the same `gamma_machine_complex` as the scalar call. `Compile[]` lowers `Gamma` through that kernel for real scalars, complex scalars and rank-1 arrays. An integer argument uses the integer-closed opcode `OP_GAMMA_I`, so `Gamma[5]` stays the exact `24`. The scalar, packed, `NDArray` and `Compile[]` answers agree bit for bit.

**Complexity / limits.** The machine paths are `O(1)` per element. The complex kernel has a median error of 0.8 ulp and a maximum of about 4 ulp in the box `0.5 <= Re z <= 4`, `|Im z| <= 3` (the Lanczos series it replaced reached 50 ulp). The MPFR paths scale their term counts and guard bits with the precision. Exact complex arguments (`Gamma[3/2, I]`), exact non-integer `a` in `Gamma[a, z]`, and symbolic `z` stay unevaluated. The incomplete integer expansion is capped at `a <= 1000`. Attributes: `Listable`, `NumericFunction`, `Protected`.

- Exact reductions for `Gamma[z]`:
  - Positive integers: `Gamma[n] = (n-1)!` (exact, with GMP BigInt for large `n`).
  - Non-positive integers are poles: `Gamma[0]`, `Gamma[-n]` → `ComplexInfinity`.
  - Half-integers reduce to rational multiples of `Sqrt[Pi]`, e.g.
    `Gamma[1/2] = Sqrt[Pi]`, `Gamma[5/2] = 3/4 Sqrt[Pi]`,
    `Gamma[-1/2] = -2 Sqrt[Pi]` (via the Factorial functional equation).
  - `Gamma[Infinity]` → `Infinity`, `Gamma[-Infinity]` → `Indeterminate`,
    `Gamma[ComplexInfinity]` → `ComplexInfinity`.
- Exact reductions for the incomplete form:
  - `Gamma[a, 0] = Gamma[a]`, `Gamma[a, Infinity] = 0`.
  - **Positive integer first argument** reduces to its finite closed form
    `Gamma[n, z] = (n-1)! e^-z Σ_{k=0}^{n-1} z^k/k!` for symbolic or exact `z`,
    e.g. `Gamma[1, z] = E^-z`, `Gamma[2, x] = (1 + x) E^-x`,
    `Gamma[3, x] = (2 + 2 x + x^2) E^-x`, and exact `Gamma[2, 3] = 4/E^3`.
- Numeric evaluation:
  - Machine-precision real → libm `tgamma`. Machine-precision complex →
    the Gamma-family double kernel `sfg_gamma_c`
    (`src/special_functions/sf_machine_gamma.c`), shared bit-for-bit by the
    scalar call, packed / `NDArray` maps and `Compile[]`, e.g.
    `Gamma[2.3 + I] = 0.719141 + 0.540614 I`. For `-5 ≤ Re z < 10`,
    `|Im z| < 19/4` the exact recurrence moves z into the strip `1 ≤ Re < 5/2`
    and a Taylor series of lnΓ about the nearest of 30 precomputed centres is
    summed (the constant term carried as a double-double); elsewhere shifted
    Stirling (Re w ≥ 7) with the exponent formed in double-double (Dekker
    products, an extra-precise log|w|), then ONE complex exponential. Accuracy:
    median 0.8 / max 3.5 ulp on `0.5 ≤ Re z ≤ 4, |Im z| ≤ 3`; median 3.6 ulp
    over `-10..30 × ±20i`, where the conditioning |z ψ(z)| itself reaches 10².
    (The former Lanczos g = 7 series was median 9 / max 53 ulp on the same box.)
    `Re z < -30` and overflow fall back to the Lanczos/reflection path.
  - Arbitrary precision (MPFR) real → `mpfr_gamma`, output precision tracking
    the input, e.g. `N[Gamma[22/10], 50]` and `Gamma[2.2`200]`.
  - **Arbitrary-precision complex** `Gamma[z]` → Spouge's approximation,
    whose coefficients are computed at runtime to the requested precision
    (reflection for `Re(z) < 1/2`), e.g.
    `N[Gamma[I], 50] = -0.15494982830181068512… − 0.49801566811835604271… I`.
  - Incomplete real (machine or MPFR) → `mpfr_gamma_inc`, e.g.
    `Gamma[1.5, 7.5] = 0.00160996`, `Gamma[1, 1.1, 2.2] = 0.222068`.
  - **Incomplete complex** `Gamma[a, z]` (machine or arbitrary precision) →
    a lower-incomplete series (`Re(z) < Re(a)+1`) or a Lentz continued
    fraction otherwise, e.g. `Gamma[2.0, 1 + I] = (2 + I) e^-(1+I)` and
    `N[Gamma[3/2, 2 + I], 30] = 0.160487401929263240… − 0.176588715957602346… I`.
- Derivatives: `D[Gamma[a, z], z] = -z^(a-1) E^-z` (chain rule on both
  arguments; the `a`-derivative is the generic `Derivative[1,0][Gamma][a,z]`).
  `D[Gamma[z], z] = Gamma[z] PolyGamma[0, z]`, so higher derivatives compose
  through `PolyGamma`, e.g.
  `D[Gamma[z], {z, 2}] = Gamma[z] PolyGamma[1, z] + Gamma[z] PolyGamma[0, z]^2`.
- All other arguments (e.g. `Gamma[1/3]`, `Gamma[x]`, `Gamma[a, z]`,
  exact non-integer `Gamma[3/2, z]`, exact complex `Gamma[3/2, I]`) stay
  unevaluated.

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## References

**See also:** [NDArray](../../linear-algebra/NDArray/), [PolyGamma](../../special-functions/PolyGamma/)

- DLMF §5.2, §5.5 — the Euler gamma function and the reflection formula.
- DLMF §5.11 — Stirling's series for ln Gamma.
- DLMF §8.2, §8.11 — the (upper) incomplete gamma function Gamma(a,z).
- J. L. Spouge, *Computation of the gamma, digamma, and trigamma functions*, SIAM J. Numer. Anal. 31 (1994) 931–944.
- C. Lanczos, *A precision approximation of the gamma function*, SIAM J. Numer. Anal. Ser. B 1 (1964) 86–96.
- D. E. G. Hare, *Computing the principal branch of log-Gamma*, J. Algorithms 25 (1997) 221–236.
- T. J. Dekker, *A floating-point technique for extending the available precision*, Numer. Math. 18 (1971) 224–242.
- Source: [`src/special_functions/gamma.c`](https://github.com/stblake/mathilda/blob/main/src/special_functions/gamma.c)
- Specification: [`docs/spec/builtins/special-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/special-functions.md)
- Tests: [`tests/test_airyai.c`](https://github.com/stblake/mathilda/blob/main/tests/test_airyai.c)
- Tests: [`tests/test_airybi.c`](https://github.com/stblake/mathilda/blob/main/tests/test_airybi.c)
- Tests: [`tests/test_beta.c`](https://github.com/stblake/mathilda/blob/main/tests/test_beta.c)
- Tests: [`tests/test_compile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile.c)

## Notes & additional examples

### Notes

- `Gamma[n] = (n-1)!` for positive integers; half-integers reduce exactly to rational multiples of `Sqrt[Pi]`; non-positive integers are poles (`ComplexInfinity`).
- Machine reals use libm `tgamma`; past `Gamma[171.6]` the value is recomputed in MPFR and returned as a machine-precision number with an extended exponent instead of overflowing.
- Machine complex arguments use the `sf_machine_gamma.c` kernel (Taylor tiles near the real axis, shifted Stirling elsewhere). Scalar, packed list, `NDArray` and `Compile[]` share it bit for bit; arbitrary-precision complex uses Spouge's approximation at the requested precision.
- `Gamma[a, z]` is the upper incomplete gamma and `Gamma[a, z0, z1] = Gamma[a, z0] - Gamma[a, z1]`; inexact arguments go through `mpfr_gamma_inc` (real) or a series/continued-fraction split (complex).
- Attributes: `Listable`, `NumericFunction`, `Protected`. See also `LogGamma`, `PolyGamma`, `Beta`, `Pochhammer`, `Factorial`.
