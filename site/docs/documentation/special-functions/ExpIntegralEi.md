# ExpIntegralEi

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`ExpIntegralEi[z]`**

gives the exponential integral Ei(z), the principal value of -Integral\_{-z}^Infinity e^-t/t dt, with a branch cut on (-Infinity, 0).

**`ExpIntegralEi[0] = -Infinity, ExpIntegralEi[Infinity] = Infinity,`**

**`ExpIntegralEi[-Infinity] = 0, ExpIntegralEi[+-I Infinity] = +-I Pi. Real and`**

**`D[ExpIntegralEi[z], z] = E^z/z. Listable.`**

<details>
<summary>Notes</summary>

complex inputs evaluate numerically at machine or arbitrary (MPFR) precision;

</details>

## Examples (13)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= ExpIntegralEi[2.]
Out[1]= 4.95423

In[2]:= ExpIntegralEi[-1.]
Out[2]= -0.219384

In[3]:= ExpIntegralEi[2. + I]
Out[3]= 4.06998 + 3.40094*I
```

### Applications (10)

Exact special values

```mathematica
In[4]:= {ExpIntegralEi[0], ExpIntegralEi[Infinity], ExpIntegralEi[I Infinity]}
Out[4]= {-Infinity, Infinity, I Pi}
```

Machine real on the Cody-Thacher-style kernel

```mathematica
In[5]:= ExpIntegralEi[2.]
Out[5]= 4.95423
```

Close to the zero x0 = 0.37250741..., still relatively accurate

```mathematica
In[6]:= ExpIntegralEi[0.3725]
Out[6]= -2.88742e-05
```

On the branch cut the principal value is real

```mathematica
In[7]:= ExpIntegralEi[-1.]
Out[7]= -0.219384
```

Just above and below the cut the imaginary part jumps by 2 Pi

```mathematica
In[8]:= {ExpIntegralEi[-1. + 10^-10 I], ExpIntegralEi[-1. - 10^-10 I]}
Out[8]= {-0.219384 + 3.14159*I, -0.219384 - 3.14159*I}
```

Machine complex

```mathematica
In[9]:= ExpIntegralEi[2. + I]
Out[9]= 4.06998 + 3.40094*I
```

Past the double range the MPFR path returns an extended-exponent number

```mathematica
In[10]:= ExpIntegralEi[800.]
Out[10]= 3.412238865448377e+344
```

Arbitrary precision through mpfr_eint

```mathematica
In[11]:= N[ExpIntegralEi[2], 50]
Out[11]= 4.95423435600189016337950513022703527551805356242004
```

Arbitrary-precision complex

```mathematica
In[12]:= N[ExpIntegralEi[I], 30]
Out[12]= 0.3374039229009681346626462038893 + 2.516879397162079634172675005462*I
```

The derivative E^z/z with the chain rule

```mathematica
In[13]:= D[ExpIntegralEi[x^2], x]
Out[13]= (2 E^x^2)/x
```

## Algorithm

Mathilda -- the exponential integral Ei.

```text
  ExpIntegralEi[z]   Ei(z) = -PV Int_{-z}^Inf e^-t/t dt
```

Ei has a branch cut along the negative real axis (-Inf, 0); the principal value is taken on the cut. Evaluation is layered so each kind of argument takes the cheapest route:

```text
  exact special values    ->  -Infinity, Infinity, 0, +-I Pi, Indeterminate
  machine real x > 0       ->  MPFR mpfr_eint (correctly rounded)
  machine real x < 0       ->  real convergent series with ln|x|
  arbitrary real           ->  the same, at the input precision
  complex (any precision)  ->  the convergent series in MPFR with guard bits
  everything else          ->  stays symbolic (return NULL)
```

The convergent series (DLMF 6.6.2) is

```text
  Ei(z) = gamma + Log(z) + Sum_{k>=1} z^k / (k k!),
```

valid for all z != 0. On the real negative axis the real part ln|x| is used (principal value on the cut); for genuinely complex z the principal Log(z) supplies the +-i Pi jump across the cut. The partial sums of the series can reach magnitude ~e^|z| before the answer (~e^Re z) emerges, so the MPFR path adds (|z| + |Re z|)/ln2 guard bits to absorb that cancellation exactly. For x > 0 MPFR's native mpfr_eint is used directly (fast and correctly rounded, which also covers very-high-precision arguments).

Attributes: Listable, NumericFunction, Protected.

## Performance

Against other systems, from the benchmark suite (same input, results cross-checked for agreement):

| case | Mathilda | Wolfram | Python |
|---|---:|---:|---:|
| integrate Sin[x] Exp[x] | 6.62 s | 0.726 s | 5.47 s |
| integrate Exp[x]/x (non-elementary) | 5.38 s | 0.318 s | 33.5 s |
| integrate Tan[x]^3 | 4.07 s | 0.196 s | 3.62 s |
| integrate 1/(1+Exp[x]) | 2.43 s | 0.109 s | 3.31 s |
| integrate Log[x]^3 | 1.59 s | 1.85 s | 10.6 s |
| integrate x Exp[x^2] | 0.685 s | 0.087 s | 4.42 s |

## Implementation notes

**Algorithm.** `builtin_expintegralei` evaluates `Ei(z)`, which has a branch cut along `(-∞, 0)`; on the cut it returns the principal value. A machine call first passes the `sf_machine_scalar_args` gate and goes to the double kernels in `sf_machine_expint.c`.

**Real machine kernel** (`sf_machine_ei`), following Cody & Thacher's scheme with coefficients refitted from mpmath by `tools/sf_accuracy/fit_expint.py`:

- `0 < x <= 16`: `γ + log x + x·R(x)`, where `R` is split into six pieces, each a polynomial with positive coefficients (so Horner never cancels).
- `16 < x <= 64`: `e^x/x·(1 + u·H(u))` with `u = 1/x`, using rationals on two octaves.
- `x > 64`: the asymptotic series.
- `x < 0`: `-E1(-x)`, where `E1` is a polynomial for `x <= 0.6` and `e^{-x}/x·(1 - u·J(u))` beyond.

`e^x` is formed as `e^{x/2}·e^{x/2}`, so it does not overflow early.

**Complex machine kernel** (`sf_machine_ei_c`):

- The power series by Horner wherever its terms do not cancel (`|z| - Re z <= 2`).
- Otherwise `-E1(-z) ± iπ`, with `E1` from the continued fraction (DLMF 6.9.1). The fraction is evaluated backward from a fixed depth `n = 240/s + 8`, where `s = |w| + Re w`; forward Lentz evaluation measured 55 ulp.
- The asymptotic series for `|z| > 45`.

A decline (the pole at 0, overflow, underflow) falls through to `ei_one_arg`:

- **Exact special values:** `Ei[0] = -Infinity`, `Infinity → Infinity`, `-Infinity → 0`, `±I Infinity → ±I π`, and `ComplexInfinity`/`Indeterminate → Indeterminate`.
- **Real, `x > 0`, any precision:** `mpfr_eint`, which is correctly rounded. This also gives the extended-exponent answer for `ExpIntegralEi[800.]`.
- **Real, `x < 0`:** the series `γ + ln|x| + Σ x^k/(k·k!)` with `|x|/ln2` guard bits, or the asymptotic form when `|x|` is large.
- **Complex:** the same series with the principal `Log z` in the `ecx` toolkit, using `(|z| + |Re z|)/ln2` guard bits, or for large `|z|` the asymptotic form plus the `iπ sign(Im z)` jump.

An explicit low-precision request (an MPFR leaf below 53 bits) gets an MPFR result, so the requested digits print. Every result passes through `numeric_promote_result_if_degenerate`.

**Data structures.** `Expr` trees and the `ecx` toolkit (pairs of `mpfr_t`) for the MPFR paths. The machine side keeps coefficient tables in `sf_machine_expint_coef.h` and does complex arithmetic on explicit `(re, im)` doubles. The ND kernel is `NDKU_ExpIntegralEi = { expintegralei_machine_complex, ndk_ExpIntegralEi_r }`, registered with `REG_U`, so it is packed-aware. Its real arm declines the pole, and its complex arm is `sf_machine_ei_c`. `Compile[]` lowers `ExpIntegralEi` through this kernel for real and complex scalars and for rank-1 arrays. The scalar call, packed arrays, `NDArray` and `Compile[]` all use the same helpers and agree bit for bit.

**Complexity / limits.** The machine kernels are `O(1)`. Against mpmath they show a median error of 0.4 ulp (real) and 0.6 ulp (complex), with a maximum of 2.6 ulp away from the zero `x₀ = 0.3725`. The MPFR series costs about `2|z|` terms plus guard bits until the asymptotic regime takes over. Exact non-zero arguments (`ExpIntegralEi[2]`) stay unevaluated. Attributes: `Listable`, `NumericFunction`, `Protected`.

- Exact special values: `ExpIntegralEi[0] = -Infinity`,
  `ExpIntegralEi[Infinity] = Infinity`, `ExpIntegralEi[-Infinity] = 0`,
  `ExpIntegralEi[I Infinity] = I Pi`, `ExpIntegralEi[-I Infinity] = -I Pi`;
  `ComplexInfinity` and `Indeterminate` map to `Indeterminate`.
- Exact non-zero arguments stay symbolic (`ExpIntegralEi[2]`, `ExpIntegralEi[1/2]`);
  numeric values follow from a `Real`/MPFR argument or from `N`.
- Numeric evaluation (machine *and* arbitrary precision):
  - Real x > 0 → MPFR `mpfr_eint` (correctly rounded, fast even at very high
    precision: `ExpIntegralEi[1.8] = 4.24987`, `ExpIntegralEi[2.] = 4.95423`,
    `N[ExpIntegralEi[2], 50] = 4.9542343560018901633795051302270352755180535624200`).
  - Real x < 0 → the on-cut convergent series Ei(x) = γ + ln|x| + Σ xᵏ/(k·k!),
    in MPFR with `|x|/ln2` guard bits to absorb the partial-sum cancellation, and
    returns a **real** principal value (`ExpIntegralEi[-1.] = -0.219384`,
    `ExpIntegralEi[-5.] = -0.00114830`).
  - **Complex** → the same series with the principal `Log(z)`, evaluated in MPFR
    with `(|z| + |Re z|)/ln2` guard bits, so machine-precision complex results are
    fully accurate, e.g. `ExpIntegralEi[2. + I] = 4.06998 + 3.40094 I`,
    `N[ExpIntegralEi[2 + I], 30] = 4.06998094789392774228769025521 + 3.40094396980012162163040462603 I`.
    Approaching the cut from above gives +I Pi, from below −I Pi
    (`ExpIntegralEi[-1. + 10^-10 I] ≈ -0.219384 + 3.14159 I`). A `double complex`
    series is the `USE_MPFR=0` fallback.
  - **Large |z|** (real or complex) → the convergent series is infeasible
    (it would need ~`|z|/ln2` guard bits and ~`2|z|` terms), so once `|z|`
    exceeds roughly `prec·ln2` the **asymptotic expansion** takes over:
    `Ei(z) ~ (e^z/z) Σ_{k≥0} k!/z^k + i π sign(Im z)`, summed to its smallest
    term (DLMF 6.12.2; the `i π sign(Im z)` constant is the branch jump, so
    `ExpIntegralEi[±I Infinity] = ±I Pi` is recovered in the limit). The two
    regimes overlap and agree (`ExpIntegralEi[-50 I]` from either route gives
    `-0.00562839 − 3.12241 I`). This also fixes an `mpfr_init2` abort: the old
    guard-bit count `(long)(|z|/ln2)` overflowed for astronomically large `|z|`
    (e.g. `N[ExpIntegralEi[-10^60 I] + I Pi, 20]`).
- **Machine precision** (real and complex) runs on double kernels shared bit-for-bit
  with the packed/NDArray kernels and `Compile[]`
  (`src/special_functions/sf_machine_expint.c`; the scalar builtin calls the same
  function, so scalar == List == packed == NDArray == Compile). A decline (pole,
  overflow, an underflow to a subnormal) falls through to the MPFR path unchanged,
  and arbitrary precision is untouched. Real: Ei = γ + log x + x·R(x) for 0 < x ≤ 16
  (R on six pieces, each a positive-coefficient polynomial, so Horner never
  cancels), e^x/x·(1 + H(1/x)) rational on [16, 64], the asymptotic series beyond;
  Ei(−x) = −E1(x) with E1 from a polynomial for x ≤ 0.6 and e^−x/x·(1 − J(1/x))
  rationals beyond (Cody & Thacher's scheme, coefficients refitted with mpmath by
  `tools/sf_accuracy/fit_expint.py`). Complex: the power series by Horner where its
  terms do not cancel (|z| − Re z ≤ 2), else −E1(−z) ± iπ with E1 from the continued
  fraction evaluated backward from a calibrated depth, and the asymptotic series for
  |z| > 45. Accuracy vs mpmath: median 0.4 ulp (real) / 0.6 (complex), max ≤ 2.6 ulp
  off the zero x₀ = 0.3725 (scaled by |x f′| there).
- Derivative: `D[ExpIntegralEi[z], z] = E^z/z` (chain rule applies, e.g.
  `D[ExpIntegralEi[x^2], x] = (2 E^x^2)/x`); the origin Taylor series at a regular
  point follows from the generic `D`-based fallback.
- Wrong arity emits `ExpIntegralEi::argx` and stays unevaluated.

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## References

**See also:** [Real](../../other-advanced/Real/), [N](../../arithmetic/N/), [D](../../calculus/D/)

- DLMF §6.2.5, §6.6.2 — the exponential integral Ei(z) and its convergent series γ + Log(z) + Sum z^k/(k k!).
- DLMF §6.9.1 — the continued fraction for E1; DLMF §6.12.2 — the asymptotic expansion Ei(z) ~ (e^z/z) Sum k!/z^k.
- W. J. Cody and H. C. Thacher, *Rational Chebyshev approximations for the exponential integral E1(x)*, Math. Comp. 22 (1968) 641–649.
- W. J. Cody and H. C. Thacher, *Chebyshev approximations for the exponential integral Ei(x)*, Math. Comp. 23 (1969) 289–303.
- Source: [`src/special_functions/expintegralei.c`](https://github.com/stblake/mathilda/blob/main/src/special_functions/expintegralei.c)
- Specification: [`docs/spec/builtins/special-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/special-functions.md)
- Tests: [`tests/test_cherry_ei.c`](https://github.com/stblake/mathilda/blob/main/tests/test_cherry_ei.c)
- Tests: [`tests/test_cherry_li.c`](https://github.com/stblake/mathilda/blob/main/tests/test_cherry_li.c)
- Tests: [`tests/test_compile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile.c)
- Tests: [`tests/test_expint_machine.c`](https://github.com/stblake/mathilda/blob/main/tests/test_expint_machine.c)

## Notes & additional examples

### Notes

- `ExpIntegralEi[z] = -PV ∫_{-z}^∞ e^(-t)/t dt`, with a branch cut on the negative real axis where the principal (real) value is returned. Exact non-zero arguments stay unevaluated.
- Machine real and complex arguments use the `sf_machine_expint.c` kernels (Cody–Thacher pieces on the real line; Horner series, a backward-evaluated E1 continued fraction and the asymptotic series in the complex plane), shared bit for bit with packed lists, `NDArray` and `Compile[]`.
- Arbitrary precision: `mpfr_eint` for `x > 0`, the guarded convergent series otherwise, switching to the asymptotic expansion for large `|z|`.
- Attributes: `Listable`, `NumericFunction`, `Protected`. See also `ExpIntegralE`, `LogIntegral`, `CosIntegral`, `SinIntegral`, `CoshIntegral`, `SinhIntegral`.
