# EllipticE

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`EllipticE[m]`**

is the complete elliptic integral of the second kind, Integrate\[Sqrt\[1 - m Sin\[t\]^2\], {t, 0, Pi/2}\], and EllipticE\[phi, m\] the incomplete one, with upper limit phi. The parameter argument is m = k^2, not the modulus k. EllipticE\[0\] is Pi/2 and EllipticE\[1\] is 1. EllipticE\[phi, 1\] is Sin\[phi\] only for |phi| \<= Pi/2, since E(phi|1) is the integral of Abs\[Cos\[t\]\]: EllipticE\[2, 1\] is 2 - Sin\[2\], not Sin\[2\].

## Examples (24)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (6)

```mathematica
In[1]:= EllipticK[0]
Out[1]= 1/2 Pi

In[2]:= N[EllipticK[1/2], 25]
Out[2]= 1.8540746773013719184338503

In[3]:= D[EllipticF[phi, m], phi]
Out[3]= 1/Sqrt[1 - m Sin[phi]^2]

In[4]:= N[EllipticPi[3/2, ArcSin[Sqrt[2] Sqrt[1/(1 + 7/5)]], 1/2], 20]
Out[4]= 0.987739699728522080208 - 2.7206990463513267759*I

In[5]:= EllipticK[2.]
Out[5]= 1.31103 - 1.31103*I

In[6]:= EllipticE[0.3 + 0.4 I]
Out[6]= 1.46251 - 0.175161*I
```

### Applications (18)

Complete, m = 0: the quarter period again

```mathematica
In[7]:= EllipticE[0]
Out[7]= 1/2 Pi
```

The integrand is Abs[Cos[t]], whose quarter-period integral is 1

```mathematica
In[8]:= EllipticE[1]
Out[8]= 1
```

One argument is complete, two incomplete -- and a full quarter period closes the gap

```mathematica
In[9]:= EllipticE[Pi/2, 1/2]
Out[9]= EllipticE[1/2]
```

Unlike K, E has no closed form at 1/2 or -1

```mathematica
In[10]:= N[EllipticE[1/2], 30]
Out[10]= 1.350643881047675502520174735339
```

E(phi|1) = Integrate[Abs[Cos[t]]], which is Sin[phi] on |phi| <= Pi/2

```mathematica
In[11]:= EllipticE[1/2, 1]
Out[11]= Sin[1/2]
```

Past Pi/2 that identity is false, so the call stays symbolic rather than answering Sin[2]

```mathematica
In[12]:= EllipticE[2, 1]
Out[12]= EllipticE[2, 1]
```

```mathematica
In[13]:= N[EllipticE[2, 1], 20]
Out[13]= 1.09070257317431830461

In[14]:= N[2 - Sin[2], 20]
Out[14]= 1.09070257317431830461
```

The parameter derivative of the complete integral

```mathematica
In[15]:= D[EllipticE[m], m]
Out[15]= (EllipticE[m] - EllipticK[m])/(2 m)
```

The amplitude derivative is the integrand

```mathematica
In[16]:= D[EllipticE[phi, m], phi]
Out[16]= Sqrt[1 - m Sin[phi]^2]
```

Same dedicated kernel as K, with the 1/(1-2k) factor

```mathematica
In[17]:= Series[EllipticE[m], {m, 0, 3}]
Out[17]= 1/2 Pi - 1/8 Pi m - 3/128 Pi m^2 - 5/512 Pi m^3 + O[m]^4
```

Certified decreasing in m, so the endpoints come back swapped

```mathematica
In[18]:= EllipticE[Interval[{1/4, 1/2}]]
Out[18]= Interval[{EllipticE[1/2], EllipticE[1/4]}]
```

Past m = 1 a machine argument gets a machine Complex, the limit from below the cut

```mathematica
In[19]:= EllipticE[2.]
Out[19]= 0.59907 + 0.59907*I
```

The arbitrary-precision Arb path agrees

```mathematica
In[20]:= N[EllipticE[2], 20]
Out[20]= 0.59907011736779610372 + 0.59907011736779610372*I
```

Complex parameter at machine precision, by the AGM

```mathematica
In[21]:= EllipticE[0.3 + 0.4 I]
Out[21]= 1.46251 - 0.175161*I
```

The defining integral by quadrature

```mathematica
In[22]:= NIntegrate[Sqrt[1 - (0.3 + 0.4 I) Sin[t]^2], {t, 0, Pi/2}]
Out[22]= 1.46251 - 0.175161*I
```

Compiled and interpreted calls share one kernel

```mathematica
In[23]:= f = Compile[{{z, _Complex}}, EllipticE[z]]; f[0.3 + 0.4 I] == EllipticE[0.3 + 0.4 I]
Out[23]= True
```

Listable; the unary kernel runs element-wise on the buffer

```mathematica
In[24]:= EllipticE[{0.1, 0.2, 0.3}]
Out[24]= {1.53076, 1.48904, 1.44536}
```

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

**Algorithm.** The head is arity-overloaded and `builtin_elliptice` dispatches on argument
count; a wrong count emits `EllipticE::argt`.

- **Complete, `E[m]`.** The machine fast path `ell_complete_machine` runs first, gated by
  the shared `sf_machine_scalar_args` (machine `Real` or machine `Complex`; exact and
  `EXPR_MPFR` arguments do not pass). A real `m ≤ 1` goes to `sf_machine_ellk_e` in
  `src/special_functions/sf_machine_ellk.c`: Cephes' structure `E = P(t) − t log(t) Q(t)`,
  `t = 1 − m`, with degree-10 `P`, `Q` re-fitted by `tools/sf_accuracy/fit_ellk.py`
  (relative minimax, fit error 5.5·10⁻¹⁹), `t = 0` answered as `E(1) = 1`, and `m < 0`
  folded onto `t' = 1/t` by DLMF 19.7.5. A complex `m`, or a real `m > 1`, goes to
  `sf_machine_ellk_e_c`: the AGM `K = π/(2 M(1, √(1−m)))` with Cox's "right choice" of
  square root, plus the running sum `E = K (1 − Σ 2ⁿ⁻¹ cₙ²)` (DLMF 19.8.6), all in raw
  `double` re/im. The branch cut is `[1, ∞)` and a real `m` on it takes the limit from
  below, so `EllipticE[2.]` is `0.59907 + 0.59907 I`. Only if that kernel declines does
  `elliptice_complete` run: Arb's `acb_elliptic_e` for an inexact argument, then the exact
  reductions `E[0] = π/2`, `E[1] = 1`, `E[∞] = ComplexInfinity`.
- **Incomplete, `E[φ, m]`.** `elliptice_incomplete` runs the numeric path first —
  `elliptic_machine_e_inc` for a real pair, else Arb's `acb_elliptic_e_inc` — then the
  reductions `E[0, m] = 0`, `E[φ, 0] = φ`, `E[π/2, m] = E[m]` and the oddness fold. The
  machine kernel is Carlson's: `elliptic_inc_real(φ, m, want_E = true)` reduces the
  amplitude to the principal strip, evaluates `s R_F(c², 1−ms², 1) − (m/3) s³ R_D(…)`, and
  closes the quasi-period with the complete `E(m)` above.

**`R_F` and `R_D` are one loop.** Every caller that wants `R_D` wants `R_F` at the same
arguments, and the two recurrences walk an identical `xₘ, yₘ, zₘ` sequence, so
`carlson_rf_rd` computes both from one duplication — three square roots per step instead
of six, which pays for the tighter `EC_ERRTOL_RD = 0.0015`.

**`E[φ, 1]` is `Sin[φ]` only on the principal strip.** `E(φ|1) = ∫₀^φ |cos t| dt`, which is
`Sin[φ]` for `|φ| ≤ π/2` and not beyond it. Applied unconditionally the rule made
`EllipticE[2, 1]` answer `Sin[2] = 0.909297` where the value is `2 − Sin[2] = 1.090703`.
The reduction is gated on `ell_in_principal_strip`, so an amplitude off the strip stays
symbolic rather than wrong, and `N[]` routes it to Arb.

**Data structures.** `static const double` coefficient tables and `double` scalars for the
machine kernels, `acb_t` with the bridge's accuracy ladder above machine precision. Both
arities are registered ND kernels: the complete form is `NDKU_EllipticE` (`REG_U`), whose
real arm is `elliptic_machine_e_complete` (delegating to `sf_machine_ellk_e`) and whose
complex arm is `sf_machine_ellk_e_c`, so a complex `NDArray` and `Compile[{{z, _Complex}}, …]`
run the same code as the scalar; the incomplete form is `NDKB_EllipticE` (`REG_B`), real
only. Scalar, `List`, packed, `NDArray` and `Compile[]` agree bit for bit
(`tests/test_elliptic_machine.c`). A real buffer containing an `m > 1` element is
abandoned to the `List` path, whose scalar call then takes the complex kernel.

**Complexity / limits.** `O(1)`, no allocation: about 0.12 µs per real complete call and
0.36 µs per complex one; 1.3 ulp on `[0, 0.95]` and 2.4 ulp on `[−50, 1)` real, 3.9 ulp
complex on `|Re m|, |Im m| ≤ 0.9`. The complete kernels decline only at overflow
(`|Re m|` or `|Im m| > 10³⁰⁰`) or a non-finite result. The incomplete machine kernel is
real-principal-domain only (`1 − m sin²φ ≥ 0`, real `φ`); a complex amplitude or parameter
goes to Arb even at machine precision, and comes back as a 53-bit arbitrary-precision
number rather than a machine `Complex`. `Interval[]` threads the complete form by
certified monotonicity (`E` decreasing in `m` below 1); `Series` at `m = 0` uses the
dedicated kernel `(π/2) Σ aₖ mᵏ/(1−2k)`.

- Exact reductions: `EllipticK[0] = π/2`, `EllipticK[1] = ComplexInfinity`,
  `EllipticE[0] = π/2`, `EllipticE[1] = 1`; `EllipticF[0, m] = 0`,
  `EllipticF[φ, 0] = φ`, `EllipticF[π/2, m] = EllipticK[m]`;
  `EllipticE[π/2, m] = EllipticE[m]`;
  `EllipticPi[0, m] = EllipticK[m]`, `EllipticPi[0, φ, m] = EllipticF[φ, m]`,
  `EllipticPi[n, 0, m] = 0`, `EllipticPi[n, π/2, m] = EllipticPi[n, m]`.
- Exact non-special arguments stay symbolic (`EllipticF[1/3, 1/2]`); a numeric
  value follows from an inexact argument or from `N`.
- **Numeric evaluation** routes through FLINT/Arb's `acb_elliptic_*`, which is
  rigorous, arbitrary-precision, defined on the whole complex plane, and already
  uses the parameter convention and branch placement above. That buys the three
  things a hand-rolled kernel most easily gets wrong:
  - the **quasi-periodic extension** off the principal strip
    (`F(φ + kπ | m) = F(φ | m) + 2k K(m)`, and likewise for `E` and `Π`), so φ of
    any size is handled;
  - **complex φ**, which `EllipticF[ArcSin[z], m]` produces routinely as soon as
    `|z| > 1` — the normal case for an elliptic pencil written in `ArcSin` form;
  - the **value past the pole** for `EllipticPi` with `n > 1`, where the path
    crosses `sin²t = 1/n`. It is not a real principal value — `Π(3/2 | 1/2)` is
    `−0.456720313453 − 2.72069904635 i`, mpmath and Arb agreeing — so the machine
    kernel declines there and Arb answers.
- **Complete integrals at machine precision** (`src/special_functions/sf_machine_ellk.c`).
  One kernel per (function, real | complex), called by the scalar builtin, the
  packed/NDArray element-wise path and `Compile[]` alike, so all five surfaces
  (scalar, List, packed, NDArray, Compile) agree bit for bit.
  - *Real `m < 1`*: Cephes' form (Moshier's `ellpk`/`ellpe`; the structure of
    A&S 17.3.34/17.3.36) `K = P(t) − log(t) Q(t)`, `E = P(t) − t log(t) Q(t)`,
    `t = 1 − m`, with degree-10 `P`, `Q` **re-fitted** for Mathilda by
    `tools/sf_accuracy/fit_ellk.py` (linear relative-error minimax by Lawson
    IRLS in mpmath, points down to `t = 10⁻³⁰⁰`; fit error ≤ 1.2·10⁻¹⁸).
    `m < 0` by `K(m) = K(m/(m−1))/√(1−m)`, `E(m) = √(1−m) E(m/(m−1))`
    (DLMF 19.7.5), i.e. the same polynomials at `t' = 1/t`. Accuracy vs mpmath:
    max 1.7 ulp (`K`), 1.3 ulp (`E`) on `[0, 0.95]`; ≤ 2.4 ulp on `[−50, 1)`.
  - *Complex `m`*: the AGM, `K = π/(2 M(1, √(1−m)))`,
    `E = K (1 − Σ 2ⁿ⁻¹ cₙ²)` (DLMF 19.8.5–6), in raw `double` re/im with the
    "right choice" of square root at each step (Cox 1984), which gives the
    principal branch with the cut on `[1, ∞)`; on the cut the value is the limit
    from below (`EllipticK[2.] = 1.31103 − 1.31103 I`, as Mathematica and
    mpmath). Max 3.4 ulp (`K`), 3.9 ulp (`E`) on `|Re m|, |Im m| ≤ 0.9`; ≤ 4.7
    ulp on `|Re m|, |Im m| ≤ 20`.
  - A machine **complex** argument, and a machine **real** `m > 1`, now return
    a machine `Complex`; they used to go through Arb and come back as 53-bit
    `EXPR_MPFR` numbers. The pole `K(1)` and `|m| > 10³⁰⁰` decline to the
    exact/Arb path.
  - Cost (net, per call): `EllipticK[0.6]` 0.19 → 0.12 µs, `EllipticE[0.6]`
    0.24 → 0.12 µs; complex 6.7 → 0.37 µs (`K`) and 19.9 → 0.36 µs (`E`).
    10⁶-element buffer: `K` 12.5 → 4.0 ms, `E` 18.5 → 3.3 ms (SciPy 16.7 / 16.5).
- **Machine kernels (incomplete).** `EllipticE[φ, m]` and `EllipticF`
  carry `double` kernels built on Carlson's symmetric forms `R_F` and `R_D`, used
  by the packed/NDArray element-wise paths. They cover the real principal domain
  and **decline** outside it (`m > 1`, or `1 − m sin²φ < 0`), which abandons the
  buffer so the List path answers exactly through Arb — slower, never wrong.
  Measured against mpmath on identical machine inputs, the buffer now agrees with
  the scalar path to 3 ulp (`K`, `F`) and 6 ulp (`E`); `R_F` and `R_D` share one
  duplication loop, since every caller that wants `R_D` wants `R_F` at the same
  arguments and the two recurrences walk an identical sequence.
  `EllipticPi` carries one too, on Carlson's `R_J` (and `R_C`, which `R_J`
  needs): the **complete** `EllipticPi[n, m]` is a binary (`REG_B`) kernel at
  41 ns/element and 2 ulp, and the **incomplete** `EllipticPi[n, φ, m]` is an
  n-ary (`REG_N`) one at 1 ulp, carrying the quasi-period
  `Π(n; φ+kπ|m) = Π(n; φ|m) + 2k Π(n|m)`. The element-wise NDArray layer tops
  out at arity 2, so only the complete form rides that path; `Compile[]` lowers
  both, the three-argument shape through `OP_KERNN`. Like the others the kernel
  **declines** outside the real principal domain: `n > 1` is not a real
  principal value at all (`Π(3/2 | 1/2)` is `−0.45672 − 2.72070 i`), so the
  honest machine answer there is to hand back to Arb.
  Because `packed_aware` is a property of the symbol and not of one arity,
  registering the two-argument kernel also stopped the transparency gate
  materialising packed `List`s for the three-argument form — which is why the
  `ndarray_delist_and_reeval` fallback in all three argument positions had to
  land first, and must stay.
- **Exact values.** `EllipticK[-1] = Γ(1/4)²/(4√(2π))` and
  `EllipticK[1/2] = 8π^(3/2)/Γ(-1/4)²` (the lemniscatic singular values; `E` has
  no closed form at either and keeps none). `EllipticPi[n, 0] = π/(2√(1-n))`,
  `EllipticPi[n, φ, 0] = ArcTanh[√(n−1) tan φ]/√(n−1)`. At infinity:
  `EllipticK[∞] = 0`, `EllipticE[∞] = EllipticE[ComplexInfinity] =
  ComplexInfinity`, `EllipticF[φ, ∞] = 0`,
  `EllipticPi[∞, m] = EllipticPi[n, ∞] = 0`. All three incomplete forms are
  **odd in the amplitude** — `EllipticF[−φ, m] = −EllipticF[φ, m]`, likewise `E`
  and `Π` — because each integrand is even in `t`; the fold uses the same
  superficial-negativity test the trig heads use, so `−2x` folds and `−x−y` does
  not.

- **`Series`.** The complete integrals expand at `m = 0` from a dedicated
  kernel: `K = (π/2) Σ aₖ mᵏ` and `E = (π/2) Σ aₖ mᵏ/(1−2k)` with `a₀ = 1`,
  `aₖ = aₖ₋₁((2k−1)/(2k))²`. It has to be dedicated — the generic
  Taylor-via-`D` path evaluates `d/dm K` at `m = 0`, which is `(π/2 − π/2)/0`
  → `Indeterminate`, and the expansion is then abandoned, so
  `Series[EllipticK[m], {m, 0, 1}]` did not previously emit even the leading
  `π/2`. Composition through an inner series works; an expansion about a regular
  point still takes Taylor-via-`D`. The incomplete forms expand in `φ` by the
  generic path, as before.

- **`Interval`.** The complete integrals thread by certified monotonicity —
  `K` increasing and `E` decreasing, both only below `m = 1`, beyond which the
  value is complex. The amplitude slot of `EllipticF`, `EllipticE[φ,m]` and
  `EllipticPi` threads through the general derivative certifier. The `m` slot
  declines: its derivative reproduces the head, so the certifier can never
  bottom out there.

- **Machine numbers.** An inexact machine scalar is answered by the `double`
  kernel and comes back as a machine `Real`, not a 53-bit `EXPR_MPFR` — so
  `Precision[EllipticK[0.5]]` is `MachinePrecision`, `MachineNumberQ` is `True`,
  and an elliptic-produced list packs for its consumers. An `EXPR_MPFR`
  argument still routes to Arb, so `SetPrecision[m, 16]` is never silently
  answered in `double`.

- **The `m = 1` continuation.** `E(φ | 1) = ∫₀^φ |cos t| dt`, which is `Sin[φ]`
  only on the principal strip `|φ| ≤ π/2`; off it the value continues as
  `Sin[φ − kπ] + 2k`. The reduction therefore fires only where it is *provably*
  valid — decided by evaluating `Abs[φ] ≤ π/2`, so an exact rational, a `Real`,
  an `MPFR` and a symbol carrying assumptions all get an answer. A symbolic
  amplitude is undecidable and stays symbolic rather than becoming wrong:
  `EllipticE[2, 1]` is left alone, and `N[EllipticE[2, 1]]` is `1.0907025731743`
  (`= 2 − Sin[2]`), agreeing with its neighbour at `m = 1 − 10⁻²⁵`.

- **Inexactness and poles.** The numeric path runs *before* the remaining exact
  reductions, so an inexact argument gets an inexact answer — `EllipticK[0.]` is
  `1.5707963267948966`, not the exact `π/2`, and `EllipticPi[0., 1/2]` is
  `1.8540746773013719`, not the symbolic `EllipticK[1/2]`. The exact reductions
  remain underneath as the fallback, which is what still reduces
  `EllipticE[π/2, 0.5]` (Arb reads numbers, not a `Times`). Poles are checked
  *ahead* of the numeric path, in every spelling of the argument, because Arb
  returns a non-finite ball there and the bridge renders that as "unevaluated":
  `EllipticK[1]`, `EllipticK[1.]` and `EllipticK[SetPrecision[1, 30]]` all give
  `ComplexInfinity`, and so does the complete `EllipticPi[1, m]` — whose
  divergence is real (`Π(1−ε | 1/2)` grows as `1/√ε`). The *incomplete* third
  kind has no pole at `n = 1`: `EllipticPi[1, 1, 1/2]` is `1.73199154202`.

- **Derivatives.** The φ-derivatives are the integrands, which is what makes a
  numeric verification of an antiderivative built from these kernels close:
  `D[EllipticF[φ, m], φ] = 1/√(1 − m sin²φ)`,
  `D[EllipticE[φ, m], φ] = √(1 − m sin²φ)`,
  `D[EllipticPi[n, φ, m], φ] = 1/((1 − n sin²φ) √(1 − m sin²φ))`. Also
  `D[EllipticK[m], m] = (E(m) − (1−m) K(m)) / (2m(1−m))`,
  `D[EllipticE[m], m] = (E(m) − K(m)) / (2m)` and
  `D[EllipticE[φ, m], m] = (E(φ,m) − F(φ,m)) / (2m)`.
  `D[EllipticF[φ, m], m] = −E(φ,m)/(2(m−1)m) − F(φ,m)/(2m) +
  sin(2φ)/(4(m−1)√(1−m sin²φ))`, and for the **complete** third kind
  `D[EllipticPi[n, m], n] = (nE + (m−n)K + (n²−m)Π)/(2(m−n)(n−1)n)` and
  `D[EllipticPi[n, m], m] = (E/(m−1) + Π)/(2(n−m))`. Each of those three was an
  inert placeholder until it had been checked against a central difference at 30
  digits.

  The `n`- and `m`-derivatives of the **incomplete** `EllipticPi` remain
  deliberately inert, and not for want of trying: a least-squares fit over the
  natural candidate basis (E, F, Π over `(n−m)`, `(m−1)`, `(n−1)`, plus the
  `sin(2φ)/√(…)` boundary term) does not recover them — residual 1.5 relative,
  with no simple rational coefficients — so the closed form involves terms that
  basis does not span. An inert derivative is honest where a guessed one corrupts
  every caller silently; `BesselJ` treats its order the same way. One visible
  consequence: `Interval[]` threading for the incomplete `EllipticPi` in its `m`
  slot stays symbolic, since the certifier differentiates to get a sign.
- Wrong arity emits `EllipticK::argx` / `EllipticF::argrx` (fixed arity) or
  `EllipticE::argt` / `EllipticPi::argt` (the overloaded pair) and stays
  unevaluated.
- Without FLINT (`USE_FLINT` undefined) the numeric path falls back to the
  machine kernels (Cephes-form / AGM for the complete integrals, Carlson for
  the incomplete ones) and otherwise declines, leaving the
  call symbolic.

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## References

**See also:** [EllipticK](../../special-functions/EllipticK/), [EllipticF](../../special-functions/EllipticF/), [EllipticPi](../../special-functions/EllipticPi/), [E](../../mathematical-constants/E/), [N](../../arithmetic/N/), [ArcSin](../../elementary-functions/ArcSin/), [Complex](../../arithmetic/Complex/), [List](../../other-advanced/List/)

- M. Abramowitz and I. A. Stegun, *Handbook of Mathematical Functions* (Dover, 1964), §17.3.36 — the `P(t) − t log(t) Q(t)` form of `E`.
- S. L. Moshier, *Cephes Math Library* (1984–2000), `ellpe` — the polynomial-with-logarithm structure (coefficients re-fitted for Mathilda).
- DLMF §19.7.5 — `E(m) = √(1−m) E(m/(m−1))`; §19.8.6 — `E` from the AGM.
- D. A. Cox, *The arithmetic-geometric mean of Gauss*, L'Enseignement Math. **30** (1984) 275-330.
- B. C. Carlson, *Computing elliptic integrals by duplication*, Numer. Math. **33** (1979) 1-16.
- B. C. Carlson, *Numerical computation of real or complex elliptic integrals*, Numer. Algorithms **10** (1995) 13-26.
- DLMF §19.25.7 — `E(φ|m) = s R_F(...) − (m/3) s³ R_D(...)`.
- Source: [`src/special_functions/elliptic.c`](https://github.com/stblake/mathilda/blob/main/src/special_functions/elliptic.c)
- Specification: [`docs/spec/builtins/special-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/special-functions.md)
- Tests: [`tests/test_elliptic.c`](https://github.com/stblake/mathilda/blob/main/tests/test_elliptic.c)
- Tests: [`tests/test_elliptic_machine.c`](https://github.com/stblake/mathilda/blob/main/tests/test_elliptic_machine.c)
- Tests: [`tests/test_meijerg.c`](https://github.com/stblake/mathilda/blob/main/tests/test_meijerg.c)

## Notes & additional examples

### Notes

`EllipticE` is arity-overloaded exactly as in the Wolfram Language: **one** argument is the
complete integral, **two** the incomplete one. A wrong count emits `EllipticE::argt` and
leaves the call unevaluated.

The `m = 1` pair of examples is the interesting one. `E(φ|1) = ∫₀^φ |cos t| dt` equals
`Sin[φ]` only for `|φ| ≤ π/2`; applied unconditionally it made `EllipticE[2, 1]` answer
`0.909297` where the value is `2 − Sin[2] = 1.090703`, and jump 0.18 against its own
neighbour at `m = 1 − 10⁻¹⁸`. The reduction is now gated on the principal strip, so the
exact call stays symbolic — correct but not closed-form — and `N[]` routes it to Arb, where
it agrees with `2 − Sin[2]` to every digit shown.

Both arities have `double` kernels, but different ones. The **complete** `E[m]` uses
`src/special_functions/sf_machine_ellk.c`: a Cephes-form `P(t) − t log(t) Q(t)` fit for
real `m ≤ 1` (about 1 ulp), and the arithmetic-geometric mean for complex `m` or real
`m > 1`, which returns a machine `Complex` on the principal branch (cut `[1, ∞)`, value on
the cut from below). The same kernels serve packed `List`s, real and complex `NDArray`s and
`Compile[]`, bit for bit. The **incomplete** `E[φ, m]` uses Carlson's `R_F` and `R_D`,
computed from one shared duplication loop, on the real principal domain only; a complex
amplitude or parameter is answered by Arb even at machine precision.

- Attributes: `Listable`, `NumericFunction`, `Protected`.
- See also: `EllipticK`, `EllipticF`, `EllipticPi`, `Hypergeometric2F1`
  (`E(m) = (π/2) ₂F₁(−½, ½; 1; m)`).
