# EllipticK

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`EllipticK[m]`**

is the complete elliptic integral of the first kind, Integrate\[1/Sqrt\[1 - m Sin\[t\]^2\], {t, 0, Pi/2}\]. The argument is the PARAMETER m = k^2, not the modulus k. EllipticK\[0\] is Pi/2 and EllipticK\[1\] is ComplexInfinity; exact arguments otherwise stay symbolic and inexact ones evaluate numerically at their precision.

## Examples (23)

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

### Applications (17)

At m = 0 the integrand is 1, so K is the quarter period Pi/2

```mathematica
In[7]:= EllipticK[0]
Out[7]= 1/2 Pi
```

The integrand loses its last denominator at t = Pi/2

```mathematica
In[8]:= EllipticK[1]
Out[8]= ComplexInfinity
```

A lemniscatic singular value: closed form in Gamma, not a decimal

```mathematica
In[9]:= EllipticK[-1]
Out[9]= Gamma[1/4]^2/(4 Sqrt[2 Pi])
```

The second singular value; EllipticE has no closed form at either

```mathematica
In[10]:= EllipticK[1/2]
Out[10]= (8 Pi^(3/2))/Gamma[-1/4]^2
```

The same number to 30 digits, through Arb

```mathematica
In[11]:= N[EllipticK[1/2], 30]
Out[11]= 1.854074677301371918433850347195
```

A dedicated kernel: the generic Taylor path hits 0/0 here

```mathematica
In[12]:= Series[EllipticK[m], {m, 0, 3}]
Out[12]= 1/2 Pi + 1/8 Pi m + 9/128 Pi m^2 + 25/512 Pi m^3 + O[m]^4
```

Past the branch point at m = 1 the value is complex

```mathematica
In[13]:= N[EllipticK[9/5], 20]
Out[13]= 1.41933775128651291929 - 1.348846512193268578*I
```

A machine m past the cut gets a machine Complex: the limit from below, from the complex AGM kernel

```mathematica
In[14]:= EllipticK[2.]
Out[14]= 1.31103 - 1.31103*I
```

Cross-check by the reciprocal-parameter identity K(m) = (K(1/m) - I K(1 - 1/m))/Sqrt[m] at m = 2

```mathematica
In[15]:= N[(EllipticK[1/2] - I EllipticK[1/2])/Sqrt[2], 20]
Out[15]= 1.31102877714605990523 - 1.31102877714605990523*I
```

Complex parameter at machine precision: AGM in raw double re/im

```mathematica
In[16]:= EllipticK[0.3 + 0.4 I]
Out[16]= 1.65024 + 0.209511*I
```

The defining integral, by quadrature, agrees

```mathematica
In[17]:= NIntegrate[1/Sqrt[1 - (0.3 + 0.4 I) Sin[t]^2], {t, 0, Pi/2}]
Out[17]= 1.65024 + 0.209511*I
```

A complex NDArray runs the same kernel element-wise

```mathematica
In[18]:= EllipticK[NDArray[{0.5 + 0.5 I, 2. + 0. I, -3. + I}]]
Out[18]= NDArray[{1.69595 + 0.322277*I, 1.31103 - 1.31103*I, 1.0669 + 0.0767626*I}]
```

And it lowers in Compile[] at a complex scalar

```mathematica
In[19]:= CompileDiagnostics[{{z, _Complex}}, EllipticK[z]]
Out[19]= {"Compiled" -> True, "ResultType" -> "Complex", "Instructions" -> 2, "CommonSubexpressions" -> 0, "InstructionsUnoptimized" -> 2}
```

Here dK/dm ~ 2^56, and all 20 digits are still right

```mathematica
In[20]:= N[EllipticK[1 - 10^-17], 20]
Out[20]= 20.9582676515692789829
```

The parameter derivative, in terms of E and K themselves

```mathematica
In[21]:= D[EllipticK[m], m]
Out[21]= (EllipticE[m] - EllipticK[m] (1 - m))/(2 m (1 - m))
```

Certified: K is increasing in m below 1

```mathematica
In[22]:= EllipticK[Interval[{1/4, 1/2}]]
Out[22]= Interval[{EllipticK[1/4], (8 Pi^(3/2))/Gamma[-1/4]^2}]
```

Real in, Real out -- so the result packs

```mathematica
In[23]:= Precision[EllipticK[0.5]]
Out[23]= MachinePrecision
```

## Implementation notes

**Algorithm.** `builtin_elliptick` tries, in this order: the pole (`ell_is_one`, which
recognises `1`, `1.` and an `EXPR_MPFR` one alike → `ComplexInfinity`); then the
**machine fast path** `ell_complete_machine`, gated by the shared
`sf_machine_scalar_args` (a machine `Real` or machine `Complex` argument; exact and
`EXPR_MPFR` arguments do not pass); then the exact reductions — `K[0] = π/2`,
`K[∞] = 0` (since `K(m) ~ (π − i log 16m)/(2√m)`), and the two lemniscatic singular values
`m = 1/2` and `m = −1` as closed forms in `Γ` (`ell_k_closed_half`,
`ell_k_closed_minus_one`). Any remaining inexact argument — in practice an arbitrary-precision
one, or a machine one the kernels declined — goes to `flint_num_elliptic_k`, i.e. Arb's
`acb_elliptic_k`, under the bridge's accuracy ladder (`nb_eval1` doubles the working
precision while `acb_rel_accuracy_bits` falls short of the digits requested). An exact
non-special argument stays symbolic.

The machine kernels live in `src/special_functions/sf_machine_ellk.c` and are chosen by the
argument's type and value:

- *real `m < 1`* — `sf_machine_ellk_k`: Cephes' structure `K = P(t) − log(t) Q(t)`,
  `t = 1 − m`, two degree-10 Horner chains and one `log`. The coefficients are not
  Cephes': `tools/sf_accuracy/fit_ellk.py` re-fitted them as a relative-error minimax
  (Lawson IRLS in mpmath, points down to `t = 10⁻³⁰⁰` so that `Q(0)` is pinned against
  `|log t| ≈ 690`), fit error 1.2·10⁻¹⁸. For `m < 0` the DLMF 19.7.5 transformation maps
  the argument to `t' = 1/t`, so the same polynomials serve `(−∞, 1)`;
- *complex `m`, or a real `m > 1`* — `sf_machine_ellk_k_c`: the AGM `π/(2 M(1, √(1−m)))` in
  raw `double` re/im (no `<complex.h>`), taking at every step the square root with
  `Re(aₙ₊₁ conj bₙ₊₁) ≥ 0` (Cox's "right choice"), which yields the principal branch with
  the cut on `[1, ∞)`. A zero imaginary part of either sign is read as `+0`, so on the cut
  the value is the limit from below: `EllipticK[2.]` is `1.31103 − 1.31103 I`. The final
  `(π/2)/M` uses Smith's division. A complex slot holding a real `m < 1` is sent back to the
  real kernel, so `0.9 + 0. I` in a complex buffer and the scalar `0.9` agree bit for bit.

**Data structures.** `static const double` coefficient tables and plain `double`
scalars; nothing is allocated. A complex result is built by `sf_machine_complex_result`,
which hands out one shared `Complex` head node by refcount. The same two kernels are the
`NDKU_EllipticK` unary ND kernel in `src/ndkernels.c` (real arm `elliptic_machine_k`, which
now delegates to `sf_machine_ellk_k`; complex arm `sf_machine_ellk_k_c`), registered with
`REG_U`, which also makes the head `packed_aware`, and `Compile[]` calls those registered
kernels at `_Real`, `_Complex` and rank-1 shapes. So scalar, `List`, packed `List`,
`NDArray` and `Compile[]` agree bit for bit (`tests/test_elliptic_machine.c`). A *real*
buffer stays real-closed: an element with `m > 1` makes the real arm decline, the buffer
is abandoned to the `List` path, and that element's scalar call takes the complex kernel.
The incomplete integrals keep Carlson's `R_F`/`R_D` duplication in `elliptic.c`, and the
quasi-period `2k K(m)` they add off the principal strip goes through the same
`elliptic_machine_k`.

**Complexity / limits.** `O(1)`: about 0.12 µs per real call and 0.37 µs per complex one
(the AGM converges quadratically, 4–6 square roots for `|m| ≤ 1`, at most 40 iterations).
Measured against mpmath the real kernel is within 1.7 ulp on `[0, 0.95]` and 2.4 ulp on
`[−50, 1)`; the complex one within 3.4 ulp on `|Re m|, |Im m| ≤ 0.9` and 4.7 ulp out to 20.
The kernels decline the pole, a non-finite result, and `|Re m|` or `|Im m| > 10³⁰⁰` (where
the AGM products would overflow); those cases, and every request above machine precision,
are Arb's, which is rigorous but slower and single-threaded. `Series` at `m = 0` comes from
a dedicated kernel rather than Taylor-via-`D`, because `d/dm K` at `m = 0` is
`(π/2 − π/2)/0`, and `Interval[]` threads by certified monotonicity below `m = 1`.

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

**See also:** [EllipticF](../../special-functions/EllipticF/), [EllipticE](../../special-functions/EllipticE/), [EllipticPi](../../special-functions/EllipticPi/), [E](../../mathematical-constants/E/), [N](../../arithmetic/N/), [ArcSin](../../elementary-functions/ArcSin/), [Complex](../../arithmetic/Complex/), [List](../../other-advanced/List/)

- M. Abramowitz and I. A. Stegun, *Handbook of Mathematical Functions* (Dover, 1964), §17.3.34 — the `P(t) − log(t) Q(t)` form of `K`.
- S. L. Moshier, *Cephes Math Library* (1984–2000), `ellpk` — the polynomial-with-logarithm structure (coefficients re-fitted for Mathilda).
- DLMF §19.7.5 — the imaginary-modulus transformation `K(m) = K(m/(m−1))/√(1−m)`.
- DLMF §19.8.5 — `K(m) = π/(2 M(1, √(1−m)))`, the arithmetic-geometric mean.
- D. A. Cox, *The arithmetic-geometric mean of Gauss*, L'Enseignement Math. **30** (1984) 275-330 — the "right choice" of square root for the complex AGM.
- R. L. Smith, *Algorithm 116: Complex division*, Comm. ACM **5** (1962) 435.
- B. C. Carlson, *Numerical computation of real or complex elliptic integrals*, Numer. Algorithms **10** (1995) 13-26.
- Source: [`src/special_functions/elliptic.c`](https://github.com/stblake/mathilda/blob/main/src/special_functions/elliptic.c)
- Specification: [`docs/spec/builtins/special-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/special-functions.md)
- Tests: [`tests/test_elliptic.c`](https://github.com/stblake/mathilda/blob/main/tests/test_elliptic.c)
- Tests: [`tests/test_elliptic_machine.c`](https://github.com/stblake/mathilda/blob/main/tests/test_elliptic_machine.c)
- Tests: [`tests/test_integrate_ramanujan.c`](https://github.com/stblake/mathilda/blob/main/tests/test_integrate_ramanujan.c)

## Notes & additional examples

### Notes

`EllipticK[m]` takes the **parameter** `m = k²`, not the modulus `k` — the Wolfram
convention, and the one place a reader of these signatures goes silently wrong.

Three argument classes are answered by three different routes. An exact `0`, `1`, `∞`,
`1/2` or `−1` reduces symbolically; any other exact argument stays symbolic, since a
decimal would be a loss of information rather than an evaluation. A machine argument
is answered in `double` (`src/special_functions/sf_machine_ellk.c`): a real `m < 1` by a
Cephes-form `P(t) − log(t) Q(t)` fit, about 1 ulp, and a complex `m` or a real `m > 1`
by the arithmetic-geometric mean, which returns a machine `Complex` on the principal
branch (cut on `[1, ∞)`, value on the cut taken from below, as in Mathematica). Any
request above machine precision goes to Arb's `acb_elliptic_k`, rigorous to the digits
asked for.

- Attributes: `Listable`, `NumericFunction`, `Protected`. The same machine kernels serve
  the scalar call, packed `List`s, real and complex `NDArray`s and `Compile[]` (`_Real`,
  `_Complex`, rank-1), bit for bit.
- `EllipticK[1.]` and `EllipticK[SetPrecision[1, 30]]` are `ComplexInfinity`, like the
  exact `1`: the pole is tested before any numeric path.
- See also: `EllipticE`, `EllipticF`, `EllipticPi`, `Hypergeometric2F1`
  (`K(m) = (π/2) ₂F₁(½, ½; 1; m)`).

The conditioning example is worth reading twice: `dK/dm ≈ 1/(2(1−m))`, which at
`m = 1 − 10⁻¹⁷` is about `2⁵⁶`, so an arbitrary-precision request has to be *planned* for
the loss rather than merely issued at the precision asked for. How much a function loses
is set by its conditioning, not by how its argument is spelled.
