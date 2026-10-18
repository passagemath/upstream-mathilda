# AiryBi

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`AiryBi[z]`**

gives the Airy function Bi(z), the solution of y'' = z y that grows exponentially as z -\> +Infinity.

**`AiryBi[0] = 1/(3^(1/6) Gamma[2/3]), AiryBi[Infinity] = Infinity,`**

**`AiryBi[-Infinity] = 0. An entire function of z. Real and complex inputs`**

**`D[AiryBi[z], z] = AiryBiPrime[z]. Listable.`**

<details>
<summary>Notes</summary>

evaluate numerically at machine or arbitrary (MPFR) precision;

</details>

## Examples (15)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= AiryBi[0]
Out[1]= 1/(3^(1/6) Gamma[2/3])

In[2]:= AiryBi[1.8]
Out[2]= 2.59587

In[3]:= AiryBi[2.5 + I]
Out[3]= 0.512544 + 5.335*I
```

### Applications (12)

Closed-form value at the origin

```mathematica
In[4]:= AiryBi[0]
Out[4]= 1/(3^(1/6) Gamma[2/3])
```

The dominant solution grows without bound

```mathematica
In[5]:= AiryBi[Infinity]
Out[5]= Infinity
```

Oscillatory on the left, growing on the right

```mathematica
In[6]:= AiryBi[{-2., 0., 2.}]
Out[6]= {-0.412303, 0.614927, 3.29809}
```

Derivative head

```mathematica
In[7]:= D[AiryBi[z], z]
Out[7]= AiryBiPrime[z]
```

Maclaurin series; every third coefficient vanishes

```mathematica
In[8]:= Series[AiryBi[x], {x, 0, 4}]
Out[8]= 1/(3^(1/6) Gamma[2/3]) + (3^(1/6) x)/Gamma[1/3] + x^3/(6 3^(1/6) Gamma[2/3]) + x^4/(4 3^(5/6) Gamma[1/3]) + O[x]^5
```

Complex argument, machine kernel

```mathematica
In[9]:= AiryBi[2.5 + I]
Out[9]= 0.512544 + 5.335*I
```

Far into the oscillatory region

```mathematica
In[10]:= AiryBi[-50.]
Out[10]= -0.13715
```

Overflow past the double range returns an arbitrary-precision value

```mathematica
In[11]:= AiryBi[1000.]
Out[11]= 5.4077118391949465475e+9154
```

Fifty digits through arb

```mathematica
In[12]:= N[AiryBi[2], 50]
Out[12]= 3.29809499997821471028060442522345242200397596340362
```

Precision follows the input

```mathematica
In[13]:= AiryBi[1.8`30]
Out[13]= 2.595869356743906290059891254823
```

Complex arbitrary precision

```mathematica
In[14]:= N[AiryBi[-4 + 3 I], 30]
Out[14]= 39.53486809004638450919341130274 - 77.72465709666247736654596806913*I
```

Compiled, bit-identical to the scalar values

```mathematica
In[15]:= cf = Compile[{{v, _Real, 1}}, AiryBi[v]]; cf[{-2., 0., 2.}]
Out[15]= {-0.412303, 0.614927, 3.29809}
```

## Algorithm

Mathilda -- the Airy function Bi.

```text
  AiryBi[z]   Airy function Bi(z), the solution of  y'' = z y  that grows
              exponentially as z -> +Infinity along the real axis. Bi is an
              *entire* function of z (no branch cuts), the companion of Ai.
```

Evaluation is layered so each kind of argument takes the most accurate and cheapest route:

```text
  exact special values   ->  AiryBi[0] = 1/(3^(1/6) Gamma[2/3]),
                             AiryBi[+Infinity] = Infinity, AiryBi[-Infinity] = 0
  machine real           ->  unified complex-MPFR core at 53 bits, real part
  arbitrary real         ->  unified complex-MPFR core at mpfr_get_prec bits
  complex (any precision) ->  unified complex-MPFR core, Complex[..] result
  everything else        ->  stays symbolic (return NULL)
```

The unified core `airy_bi_core` evaluates Bi(z) and Bi'(z) together in a file-local complex-MPFR toolkit (`acx`, pairs of mpfr_t -- no MPC library is available; this mirrors the `acx`/`ecx`/`pcx`/`gcx` toolkits in airyai.c/ erf.c/polylog.c/gamma.c). It routes between three algorithms on r = |z|, theta = arg z, and the requested output precision P:

```text
  - Maclaurin series (small/moderate |z|, accurate everywhere). From
    Bi'' = z Bi the Taylor coefficients satisfy b_0 = Bi(0), b_1 = Bi'(0),
    b_2 = 0 and b_n = b_{n-3} / (n (n-1)) for n >= 3 -- identical recurrence
    to Ai, different seed constants. The partial sums reach magnitude
    ~exp((2/3) r^{3/2}) before cancelling for complex / negative arguments,
    so the core adds  (2/3) r^{3/2} / ln2  guard bits to absorb that exactly.

  - Dominant asymptotic series (large |z|, central sector), DLMF 9.7.7/9.7.8.
    With zeta = (2/3) z^{3/2}
        Bi(z)  ~ exp(zeta)/(sqrt(pi) z^{1/4}) Sum u_k / zeta^k,
        Bi'(z) ~ z^{1/4} exp(zeta)/sqrt(pi) Sum v_k / zeta^k,
    summed to the optimal (smallest-term) truncation. The u_k, v_k are the
    SAME coefficients as Ai's asymptotic series, but with no (-1)^k sign and
    prefactor 1/sqrt(pi) (not 1/(2 sqrt(pi))). The single dominant series is
    accurate only where the neglected recessive companion ~exp(-2 Re zeta)
    is below 2^-P, i.e. Re zeta = (2/3) r^{3/2} cos(3 theta/2) > (P ln2)/2.
    Bi's anti-Stokes line is |arg z| = pi/3, so this keeps the whole positive
    half-plane (including the exponentially large positive axis) at full
    precision.

  - Connection to Ai (large |z|, otherwise -- near and left of |arg z| = pi/3,
    covering the oscillatory negative real axis). DLMF 9.2.10:
        Bi(z)  = e^{ i pi/6} Ai(z e^{ 2 pi i/3}) + e^{-i pi/6} Ai(z e^{-2 pi i/3}),
        Bi'(z) = e^{i5pi/6} Ai'(z e^{ 2 pi i/3}) + e^{-i5pi/6} Ai'(z e^{-2 pi i/3}).
    The two rotated points have |arg| <= pi and |w| = |z| (large), so they are
    evaluated by a file-local Ai asymptotic kernel (direct series + Ai's own
    2 pi/3 connection wrapper, DLMF 9.2.12). The Bi oscillation on the
    negative axis emerges naturally from the two rotated Ai evaluations.
```

D[AiryBi[z], z] = AiryBiPrime[z] (see calculus/deriv.c); the Maclaurin series at 0 is produced by the generic Taylor-via-D path once AiryBi[0] / AiryBiPrime[0] have closed-form values.

AiryBiPrime[z] = Bi'(z) is a full numeric evaluator in its own right: because

```text
`airy_bi_core` returns Bi(z) and Bi'(z) together, AiryBiPrime reuses the very
```

same Maclaurin / asymptotic / connection machinery and simply selects the derivative component. Its exact values are AiryBiPrime[0] = 3^(1/6)/Gamma[1/3] and AiryBiPrime[+Infinity] = Infinity (Bi' is the dominant, growing solution); at -Infinity Bi' has no limit (oscillation with growing ~|z|^(1/4) amplitude) and is left unevaluated.

Attributes (both heads): Listable, NumericFunction, Protected.

## Implementation notes

**Algorithm.** `builtin_airybi` mirrors `AiryAi`'s four layers. (1) **Machine scalar**: `airy_machine` uses the shared gate `sf_machine_scalar_args` and calls `sf_airy_scalar(SF_AIRY_BI, …)` in `src/special_functions/sf_machine_bessel.c`. For **real** `x` on `[-10, 10]` it evaluates the generated piecewise polynomials in `sf_machine_airy_tab.c`: `Bi` itself on `[-10, 2]`, and the scaled `Bi·e^{-zeta}` on `[2, 10]`, with `exp()` carrying the growth. Beyond that range the Poincaré expansions DLMF 9.7 apply, with `zeta = (2/3) x^{3/2}` in double-double. For a **complex** `z` it uses the Maclaurin series for `|z| <= 1`, and otherwise DLMF 9.6.3, `Bi(z) = sqrt(z/3)(I_{-1/3}(zeta) + I_{1/3}(zeta))`. One call to the right-half-plane Bessel core (shared with `BesselI`/`BesselK`) supplies every order needed, and the continuation formulas handle `Re zeta < 0`. A value past the double range declines. (2) **Exact values**: `AiryBi[0] = 1/(3^(1/6) Gamma[2/3])`, `AiryBi[Infinity] = Infinity`, `AiryBi[-Infinity] = 0` and `Indeterminate → Indeterminate`. (3) **Arbitrary precision**: an MPFR real, or a complex argument above 53 bits, goes to arb's `acb_hypgeom_airy` via `flint_num_airy_bi` (`src/flint_num_bridge.c`). (4) **Fallback**: the file-local complex-MPFR core `airy_bi_core` has three methods. It sums the Maclaurin series (`b_n = b_{n-3}/(n(n-1))`, `(2/3) r^{3/2}/ln2` guard bits) for small `|z|`. In the central sector it uses the dominant asymptotic series DLMF 9.7.7/9.7.8, but only where the recessive companion `~exp(-2 Re zeta)` is below `2^-P`. Near and left of the anti-Stokes line `|arg z| = π/3`, including the oscillatory negative axis, it uses the connection to `Ai` (DLMF 9.2.10). A machine input whose value overflows (`AiryBi[1000.] ≈ 5.41×10^9154`) is returned as an arbitrary-precision number. `AiryBiPrime` shares every layer and selects the derivative component. `D[AiryBi[z], z] = AiryBiPrime[z]` lives in `calculus/deriv.c`.

**Data structures.** Scalars use `Expr` leaves. The machine kernel works on double `(re, im)` pairs with static coefficient tables, double-double `zeta`, and `sf_machine_airy_tab.c`. The fallback core uses the same file-local `acx` complex-MPFR toolkit as `airyai.c`. The ND kernel `NDKU_AiryBi` (`REG_U` in `ndkernels.c`) has a real arm `sf_airy_bi_r` and a complex arm `sf_airy_nd`, so real and complex packed or `NDArray` buffers run element-wise. `Compile[]` lowers `AiryBi` at scalar real, scalar complex and rank-1 shapes through the same helpers, bit-identical to the scalar.

**Complexity / limits.** At machine precision the cost is O(1) per element. The documented accuracy is ≤ 4 ulp of the local magnitude on the real line, median ≈ 1 and max ≤ 10 for complex arguments. Arbitrary precision is arb's rigorous evaluation. `Bi` is entire. An exact non-zero argument (`AiryBi[2]`) stays unevaluated until `N`. `AiryBiPrime[-Infinity]` is left unevaluated (growing-amplitude oscillation).

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## References

**See also:** [AiryAi](../../special-functions/AiryAi/), [N](../../arithmetic/N/), [AiryBiPrime](../../other-advanced/AiryBiPrime/)

- DLMF §9.2, §9.7 — Airy functions: dominant asymptotic series Bi(z) ~ exp(zeta)/(sqrt(pi) z^{1/4}) Sum u_k/zeta^k, and the Bi–Ai connection Bi(z) = e^{iπ/6} Ai(z e^{2πi/3}) + e^{-iπ/6} Ai(z e^{-2πi/3}).
- DLMF §9.6.3 — Bi(z) = sqrt(z/3) (I_{-1/3}(zeta) + I_{1/3}(zeta)), used by the complex machine kernel.
- D. E. Amos, A portable package for Bessel functions of a complex argument and nonnegative order, ACM TOMS 12 (1986) 265-273 (Algorithm 644; CBIRY).
- F. Johansson, Arb: efficient arbitrary-precision midpoint-radius interval arithmetic, IEEE Trans. Comput. 66 (2017) 1281-1292.
- Source: [`src/special_functions/airybi.c`](https://github.com/stblake/mathilda/blob/main/src/special_functions/airybi.c)
- Specification: [`docs/spec/builtins/special-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/special-functions.md)
- Tests: [`tests/test_airybi.c`](https://github.com/stblake/mathilda/blob/main/tests/test_airybi.c)
- Tests: [`tests/test_compile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile.c)
- Tests: [`tests/test_numeric_stress.c`](https://github.com/stblake/mathilda/blob/main/tests/test_numeric_stress.c)
- Tests: [`tests/test_sf_machine_bessel.c`](https://github.com/stblake/mathilda/blob/main/tests/test_sf_machine_bessel.c)

## Notes & additional examples

### Notes

- `AiryBi[z]` is the solution of `y'' == z y` that grows as `z -> +Infinity`. It is entire, `AiryBi[-Infinity] = 0`, and the origin value is `1/(3^(1/6) Gamma[2/3])`.
- Machine real and complex arguments use the raw-double kernel in `sf_machine_bessel.c`. On `[-10, 10]` it evaluates piecewise Chebyshev tables; beyond that the Poincaré expansions apply. For complex `z` it uses `I_{±1/3}`. That helper also serves packed lists, `NDArray` and `Compile[]`. Arbitrary precision goes through arb's `acb_hypgeom_airy`, and the MPFR core is the fallback.
- `D[AiryBi[z], z] = AiryBiPrime[z]`. `AiryBiPrime[-Infinity]` is left unevaluated.
- Attributes: `Listable`, `NumericFunction`, `Protected`.
- See also: `AiryAi`, `AiryBiPrime`, `BesselI`.
