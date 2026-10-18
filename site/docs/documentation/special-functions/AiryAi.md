# AiryAi

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`AiryAi[z]`**

gives the Airy function Ai(z), the solution of y'' = z y that decays as z -\> +Infinity.

**`AiryAi[0] = 1/(3^(2/3) Gamma[2/3]), AiryAi[+-Infinity] = 0. An entire`**

<details>
<summary>Notes</summary>

function of z. Real and complex inputs evaluate numerically at machine or arbitrary (MPFR) precision; D\[AiryAi\[z\], z\] = AiryAiPrime\[z\]. Listable.

</details>

## Examples (17)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= AiryAi[0]
Out[1]= 1/(3^(2/3) Gamma[2/3])

In[2]:= AiryAi[1.8]
Out[2]= 0.0470362

In[3]:= AiryAi[2.5 + I]
Out[3]= -0.00191209 - 0.0180329*I
```

### Applications (14)

Closed-form value at the origin

```mathematica
In[4]:= AiryAi[0]
Out[4]= 1/(3^(2/3) Gamma[2/3])
```

Other exact arguments stay symbolic

```mathematica
In[5]:= AiryAi[2]
Out[5]= AiryAi[2]
```

Oscillatory on the left, decaying on the right

```mathematica
In[6]:= AiryAi[{-2., 0., 2.}]
Out[6]= {0.227407, 0.355028, 0.0349241}
```

The Airy equation y'' = z y

```mathematica
In[7]:= D[AiryAi[z], {z, 2}]
Out[7]= z AiryAi[z]
```

Asymptotic expansion with the essential factor kept symbolic

```mathematica
In[8]:= Series[AiryAi[x], {x, Infinity, 2}]
Out[8]= E^(-2/3 x^(3/2)) ((1/x)^(1/4)/(2 Sqrt[Pi]) - (5 (1/x)^(7/4))/(96 Sqrt[Pi]) + O[1/x]^(13/4))
```

Complex argument, machine kernel

```mathematica
In[9]:= AiryAi[2.5 + I]
Out[9]= -0.00191209 - 0.0180329*I
```

Far into the oscillatory region

```mathematica
In[10]:= AiryAi[-50.]
Out[10]= -0.161881
```

Below the double range; promoted instead of flushed to zero

```mathematica
In[11]:= AiryAi[120.]
Out[11]= 2.157498091032858e-382
```

Overflow on the imaginary axis returns an arbitrary-precision value

```mathematica
In[12]:= AiryAi[150. I]
Out[12]= 8.1951500071418076667e+374 + 6.3801431687816647218e+374*I
```

Fifty digits through arb

```mathematica
In[13]:= N[AiryAi[2], 50]
Out[13]= 0.0349241304232743791353220807918076097610602138975833
```

Precision follows the input

```mathematica
In[14]:= AiryAi[4.8`40]
Out[14]= 0.00017032552328643494849005822667939733412419
```

The Wronskian Ai Bi' - Ai' Bi = 1/Pi

```mathematica
In[15]:= N[AiryAi[-3] AiryBiPrime[-3] - AiryAiPrime[-3] AiryBi[-3], 40]
Out[15]= 0.3183098861837906715377675267450287240689
```

NDArray input runs element-wise

```mathematica
In[16]:= AiryAi[NDArray[{1., 2., 3.}]]
Out[16]= NDArray[{0.135292, 0.0349241, 0.00659114}]
```

Complex Compile matches the scalar value

```mathematica
In[17]:= cf = Compile[{{z, _Complex}}, AiryAi[z]]; cf[2.5 + I]
Out[17]= -0.00191209 - 0.0180329*I
```

## Algorithm

Mathilda -- the Airy function Ai.

```text
  AiryAi[z]   Airy function Ai(z), the solution of  y'' = z y  that tends to
              zero as z -> +Infinity along the real axis. Ai is an *entire*
              function of z (no branch cuts).
```

Evaluation is layered so each kind of argument takes the most accurate and cheapest route:

```text
  exact special values   ->  AiryAi[0] = 1/(3^(2/3) Gamma[2/3]),
                             AiryAi[+-Infinity] = 0
  machine real           ->  unified complex-MPFR core at 53 bits, real part
  arbitrary real         ->  unified complex-MPFR core at mpfr_get_prec bits
  complex (any precision) ->  unified complex-MPFR core, Complex[..] result
  everything else        ->  stays symbolic (return NULL)
```

The unified core `airy_ai_core` evaluates Ai(z) and Ai'(z) together in a file-local complex-MPFR toolkit (`acx`, pairs of mpfr_t -- no MPC library is available; this mirrors the `ecx`/`pcx`/`gcx` toolkits in erf.c/polylog.c/ gamma.c). It routes between two algorithms on r = |z| and the requested output precision P:

```text
  - Maclaurin series (small/moderate |z|). From Ai'' = z Ai the Taylor
    coefficients satisfy a_0 = Ai(0), a_1 = Ai'(0), a_2 = 0 and
    a_n = a_{n-3} / (n (n-1)) for n >= 3. The partial sums reach magnitude
    ~exp((2/3) r^{3/2}) before cancelling for complex / negative arguments,
    so the core adds  (2/3) r^{3/2} / ln2  guard bits to absorb that exactly.

  - Asymptotic series (large |z|), DLMF 9.7.5/9.7.6. With zeta = (2/3) z^{3/2}
        Ai(z)  ~ exp(-zeta)/(2 sqrt(pi) z^{1/4}) Sum (-1)^k u_k / zeta^k,
        Ai'(z) ~ -z^{1/4} exp(-zeta)/(2 sqrt(pi)) Sum (-1)^k v_k / zeta^k,
    summed to the optimal (smallest-term) truncation. The single series is
    accurate for |arg z| <= 2 pi / 3; closer to the negative real axis
    (2 pi / 3 < |arg z| <= pi) the core uses the connection relation
    (DLMF 9.2.12)
        Ai(z) = -[ w Ai(w z) + conj(w) Ai(conj(w) z)],  w = e^{2 pi i / 3},
    which maps the argument into two points with |arg| <= 2 pi / 3 where the
    direct series is accurate; the oscillation on the negative real axis then
    emerges naturally from the sum of the two rotated evaluations.
```

D[AiryAi[z], z] = AiryAiPrime[z] (see calculus/deriv.c); the Maclaurin series at 0 is produced by the generic Taylor-via-D path once AiryAi[0] / AiryAiPrime[0] have closed-form values.

AiryAiPrime[z] = Ai'(z) is a full numeric evaluator in its own right: because

```text
`airy_ai_core` returns Ai(z) and Ai'(z) together, AiryAiPrime reuses the very
```

same Maclaurin / asymptotic / connection machinery and simply selects the derivative component. Its exact values are AiryAiPrime[0] = -1/(3^(1/3) Gamma[1/3]) and AiryAiPrime[+Infinity] = 0; at -Infinity Ai' has no limit (oscillation with growing ~|z|^(1/4) amplitude) and is left unevaluated.

Attributes (both heads): Listable, NumericFunction, Protected.

## Implementation notes

**Algorithm.** `builtin_airyai` tries four layers in order. (1) **Machine scalar**: `airy_machine` uses the shared gate `sf_machine_scalar_args` (`sf_scalar.c`) and calls `sf_airy_scalar(SF_AIRY_AI, …)` in `src/special_functions/sf_machine_bessel.c`. For a **real** `x` with `-10 <= x <= 10`, the kernel evaluates piecewise polynomials from `sf_machine_airy_tab.c` (Chebyshev interpolation at 50 digits by `tools/sf_accuracy/gen_airy_tables.py`, converted to monomials per interval). On `[-10, 2]` they fit `Ai` itself; on `[2, 10]` they fit the scaled `Ai·e^zeta`, and `exp()` carries the exponential. Beyond that range the Poincaré expansions DLMF 9.7.5–9.7.12 apply, with `zeta = (2/3) x^{3/2}` formed in double-double: `Ai` is ill-conditioned in `zeta`, so a plain double `zeta` at `x = 100` would cost about 600 ulp. For a **complex** `z` the kernel uses the Maclaurin series for `|z| <= 1`. Beyond that it uses `Ai(z) = (1/π) sqrt(z/3) K_{1/3}(zeta)` (DLMF 9.6.1), taking `K` from the same Amos-style right-half-plane Bessel core as `BesselK` and continuing across `Re zeta < 0`. One core call yields all four of `I_{∓1/3}`, `K_{1/3}`, `I_{2/3}` and `K_{2/3}`. The low part of the double-double `zeta` is folded back by a first-order Taylor correction. A machine value past the double range declines. (2) **Exact values**: `AiryAi[0] = 1/(3^(2/3) Gamma[2/3])`, `AiryAi[±Infinity] = 0` and `Indeterminate → Indeterminate`. (3) **Arbitrary precision**: an MPFR real, or a complex argument above 53 bits, goes to arb's `acb_hypgeom_airy` via `flint_num_airy_ai` (`src/flint_num_bridge.c`). `Ai` is entire, so arb's value needs no branch adjustment. (4) **Fallback**: when arb declines or a machine call overflows or underflows, the hand-written complex-MPFR core `airy_ai_core` runs. It sums the Maclaurin series `a_n = a_{n-3}/(n(n-1))` with `(2/3) r^{3/2}/ln2` guard bits, or the asymptotic series DLMF 9.7.5 to optimal truncation for `|arg z| <= 2π/3`, or applies the connection DLMF 9.2.12 nearer the negative axis. `builtin_airyai` wraps the result in `numeric_promote_result_if_degenerate`, so a 53-bit underflow is re-evaluated with an extended exponent (`AiryAi[120.] ≈ 2.1575×10^-382`). `AiryAiPrime` shares every layer and selects the derivative component. `D[AiryAi[z], z] = AiryAiPrime[z]` lives in `calculus/deriv.c`.

**Data structures.** Scalars use `Expr` leaves. The machine kernel works on double `(re, im)` pairs with static `u_k`/`v_k` and Maclaurin coefficient tables, double-double helpers (`two_sum`, `two_prod`), and the generated `sf_machine_airy_tab.c`. The fallback core uses a file-local complex-MPFR toolkit `acx`. The ND kernel `NDKU_AiryAi` (`REG_U` in `ndkernels.c`) has both a real arm (`sf_airy_ai_r`) and a complex arm (`sf_airy_nd`), so packed and visible `NDArray` buffers of either kind run element-wise. `Compile[]` lowers `AiryAi` at scalar real, scalar complex and rank-1 shapes through the same helpers. Scalar, List, packed, `NDArray` and compiled values agree bit for bit.

**Complexity / limits.** At machine precision the cost is O(1) per element: one table lookup or a short asymptotic sum, or one complex Bessel core call. The documented accuracy is ≤ 4 ulp of the local magnitude on the real line, median ≈ 1 and max ≤ 10 for complex arguments. Arbitrary precision costs arb's rigorous evaluation, about 9× faster than the old MPFR core. `Ai` is entire. An exact non-zero argument (`AiryAi[2]`) stays unevaluated until `N`. `AiryAiPrime[-Infinity]` is left unevaluated (growing-amplitude oscillation).

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## References

**See also:** [Erf](../../special-functions/Erf/), [N](../../arithmetic/N/), [AiryAiPrime](../../other-advanced/AiryAiPrime/)

- DLMF §9.2, §9.7 — Airy functions: power series, asymptotic expansions, and the connection formula Ai(z) = -[w Ai(wz) + conj(w) Ai(conj(w) z)].
- DLMF §9.6.1-9.6.3 — Ai(z) = (1/pi) sqrt(z/3) K_{1/3}(zeta), the Bessel representations used by the complex machine kernel.
- D. E. Amos, A portable package for Bessel functions of a complex argument and nonnegative order, ACM TOMS 12 (1986) 265-273 (Algorithm 644; CAIRY).
- F. Johansson, Arb: efficient arbitrary-precision midpoint-radius interval arithmetic, IEEE Trans. Comput. 66 (2017) 1281-1292.
- Source: [`src/special_functions/airyai.c`](https://github.com/stblake/mathilda/blob/main/src/special_functions/airyai.c)
- Specification: [`docs/spec/builtins/special-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/special-functions.md)
- Tests: [`tests/test_airyai.c`](https://github.com/stblake/mathilda/blob/main/tests/test_airyai.c)
- Tests: [`tests/test_airybi.c`](https://github.com/stblake/mathilda/blob/main/tests/test_airybi.c)
- Tests: [`tests/test_compile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile.c)
- Tests: [`tests/test_dsolve.c`](https://github.com/stblake/mathilda/blob/main/tests/test_dsolve.c)

## Notes & additional examples

### Notes

- `AiryAi[z]` is the solution of `y'' == z y` that decays as `z -> +Infinity`. It is entire, `AiryAi[±Infinity] = 0`, and the origin value is `1/(3^(2/3) Gamma[2/3])`.
- Machine real and complex arguments use the raw-double kernel in `sf_machine_bessel.c`. On `[-10, 10]` it evaluates piecewise Chebyshev tables; beyond that the Poincaré expansions apply, with `zeta` in double-double. For complex `z` it uses `K_{1/3}`. The same helper serves packed lists, `NDArray` and `Compile[]` (scalar real, scalar complex and rank-1). Arbitrary precision goes through arb's `acb_hypgeom_airy`, and an overflow or underflow is promoted to an MPFR value.
- `D[AiryAi[z], z] = AiryAiPrime[z]`, which shares the same evaluators.
- Attributes: `Listable`, `NumericFunction`, `Protected`.
- See also: `AiryBi`, `AiryAiPrime`, `BesselK`.
