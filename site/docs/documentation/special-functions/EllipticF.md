# EllipticF

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`EllipticF[phi, m]`**

is the incomplete elliptic integral of the first kind, Integrate\[1/Sqrt\[1 - m Sin\[t\]^2\], {t, 0, phi}\]. The second argument is the PARAMETER m = k^2, not the modulus k. EllipticF\[phi, 0\] is phi and EllipticF\[Pi/2, m\] is EllipticK\[m\]; phi may be complex and of any size (the quasi-period is applied).

## Examples (19)

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

### Applications (13)

A full quarter period is the complete integral

```mathematica
In[7]:= EllipticF[Pi/2, 1/2]
Out[7]= (8 Pi^(3/2))/Gamma[-1/4]^2
```

With m = 0 the integrand is 1

```mathematica
In[8]:= EllipticF[phi, 0]
Out[8]= phi
```

Odd in the amplitude, by the same superficial-negativity test the trig heads use

```mathematica
In[9]:= EllipticF[-2 x, m]
Out[9]= -EllipticF[2 x, m]
```

The amplitude derivative is the integrand itself

```mathematica
In[10]:= D[EllipticF[phi, m], phi]
Out[10]= 1/Sqrt[1 - m Sin[phi]^2]
```

The parameter derivative, checked against a 30-digit central difference before it was allowed in

```mathematica
In[11]:= D[EllipticF[phi, m], m]
Out[11]= -EllipticF[phi, m]/(2 m) - EllipticE[phi, m]/(2 m (-1 + m)) + Sin[2 phi]/(4 (-1 + m) Sqrt[1 - m Sin[phi]^2])
```

```mathematica
In[12]:= N[EllipticF[1/5, 1/2], 25]
Out[12]= 0.2006673105648029970320719

In[13]:= NIntegrate[1/Sqrt[1 - 1/2 Sin[t]^2], {t, 0, 1/5}]
Out[13]= 0.200667
```

The quasi-period F(phi + k Pi) = F(phi) + 2k K, residual at 20 digits

```mathematica
In[14]:= N[EllipticF[1/5 + 3 Pi, 1/2] - (EllipticF[1/5, 1/2] + 6 EllipticK[1/2]), 20]
Out[14]= 1.17549435082228750797e-38
```

|z| > 1 makes ArcSin[z] complex -- the normal case for an elliptic pencil

```mathematica
In[15]:= N[EllipticF[ArcSin[3/2], 1/2], 20]
Out[15]= 1.3820851603910208068 - 1.85407467730137191843*I
```

The amplitude that used to cost eight digits

```mathematica
In[16]:= N[EllipticF[Pi/2 - 10^-8, 99/100], 20]
Out[16]= 3.69563726298987467781
```

Listable, and the list rides the packed buffer

```mathematica
In[17]:= EllipticF[{0.3, 0.6, 0.9}, 0.5]
Out[17]= {0.302255, 0.618108, 0.960966}
```

Honest limitation: no machine-complex kernel, so Arb answers and the result is a 53-bit arbitrary-precision number, not a machine Complex

```mathematica
In[18]:= EllipticF[0.5 + 0.2 I, 0.3]
Out[18]= 0.5033102181386211 + 0.2069544103393203*I
```

The defining integral along the straight path to the complex amplitude

```mathematica
In[19]:= NIntegrate[(0.5 + 0.2 I)/Sqrt[1 - 0.3 Sin[(0.5 + 0.2 I) s]^2], {s, 0, 1}]
Out[19]= 0.50331 + 0.206954*I
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

**Algorithm.** `builtin_ellipticf` runs the numeric path first (inexact in, inexact out),
then the exact reductions `F(0, m) = 0`, `F(φ, 0) = φ`, `F(φ, ∞) = 0`,
`F(π/2, m) = K(m)`, and finally the oddness fold in the amplitude. The machine kernel is
`elliptic_inc_real(φ, m, want_E = false)`:

1. reduce the amplitude to the principal strip — `k = ⌊φ/π + 1/2⌋`, `r = φ − kπ`, so
   `|r| ≤ π/2`;
2. evaluate `F(r|m) = s · R_F(cos²r, 1 − m s², 1)` with `s = sin r` (DLMF 19.25.5, after
   using `R_F`'s homogeneity to clear the `csc²` scaling);
3. restore the shift with the quasi-period `F(φ + kπ | m) = F(φ|m) + 2k K(m)`, which calls
   `elliptic_machine_k` — since v0.438 a thin wrapper over the Cephes-form kernel
   `sf_machine_ellk_k` in `sf_machine_ellk.c` that `EllipticK` itself uses, so the
   shifted value is consistent with `EllipticK[m]` and declines whenever it does.

Outside the real principal domain (`1 − m sin²φ < 0`, or a complex amplitude, which
`EllipticF[ArcSin[z], m]` produces as soon as `|z| > 1`) the kernel declines and
`flint_num_elliptic_f` answers through Arb's `acb_elliptic_f`, which carries the
quasi-period and the branch placement itself. Unlike the complete `EllipticK`/`EllipticE`,
which have machine-complex kernels, `EllipticF` has none: a machine-complex argument is
answered by Arb and comes back as a 53-bit arbitrary-precision number, not a machine
`Complex`.

**One spelling carries eight digits.** The first argument of `R_F` is `Cos[r]^2` and must
be computed as `c*c`, never as `1 - s*s`: near `r = π/2` that subtraction's absolute error
is the size of the true value, so `R_F`'s first argument arrives with ~100% relative
error. Measured against the 30-digit path at the same `double`:
`EllipticF[π/2 − 1e−8, 0.99]` was wrong by **2.7e-08** relative (1.2e8 ulp) against
6.7e-16 at a generic amplitude; it is now 1.7e-15 worst case over a π/2 approach ladder.
A single `fmax`-style clamp (`if (a < 0) a = 0`) covers the one case squaring cannot — an
`r` at which `cos r` underflows to zero, where a rounding-negative `a` would make
`carlson_rf` decline at exactly the amplitude where `F` is simply `K`.

**Data structures.** `double` scalars and the shared `carlson_rf` duplication loop; `acb_t`
on the Arb path. `EllipticF` is registered as a binary (`REG_B`) ND kernel, so a packed or
visible `NDArray` pair runs element-wise through the same `double` code.

**Complexity / limits.** `O(1)` per element, ~3 ulp; `Compile[]` lowers it at scalar and
rank-1 shapes. `D[EllipticF[φ,m], m]` is a closed form (checked against a central
difference at 30 digits before it was allowed in), unlike the `EllipticPi` parameter
derivatives, which stay inert.

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

**See also:** [EllipticK](../../special-functions/EllipticK/), [EllipticE](../../special-functions/EllipticE/), [EllipticPi](../../special-functions/EllipticPi/), [E](../../mathematical-constants/E/), [N](../../arithmetic/N/), [ArcSin](../../elementary-functions/ArcSin/), [Complex](../../arithmetic/Complex/), [List](../../other-advanced/List/)

- B. C. Carlson, *Computing elliptic integrals by duplication*, Numer. Math. **33** (1979) 1-16.
- B. C. Carlson, *Numerical computation of real or complex elliptic integrals*, Numer. Algorithms **10** (1995) 13-26.
- W. H. Press et al., *Numerical Recipes in C*, 2nd ed. (Cambridge, 1992), §6.11.
- DLMF §19.25.5 — `F(φ|m) = s R_F(1−s², 1−m s², 1)`, `s = sin φ`.
- Source: [`src/special_functions/elliptic.c`](https://github.com/stblake/mathilda/blob/main/src/special_functions/elliptic.c)
- Specification: [`docs/spec/builtins/special-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/special-functions.md)
- Tests: [`tests/test_elliptic.c`](https://github.com/stblake/mathilda/blob/main/tests/test_elliptic.c)
- Tests: [`tests/test_parallelmixedspecial.c`](https://github.com/stblake/mathilda/blob/main/tests/test_parallelmixedspecial.c)

## Notes & additional examples

### Notes

The second argument is the **parameter** `m = k²`, not the modulus `k`.

The amplitude may be of any size and may be complex, and both cases matter in practice.
Size is handled by the quasi-period `F(φ + kπ | m) = F(φ | m) + 2k K(m)`: the machine
kernel reduces `φ` into `|r| ≤ π/2`, evaluates there, and adds back `2k K(m)` — so a
shifted amplitude declines exactly when the complete integral does. Complex amplitudes
arrive unbidden, since `EllipticF[ArcSin[z], m]` is complex as soon as `|z| > 1`, which is
the normal spelling of an elliptic pencil; those go to Arb, which is defined on the whole
plane.

The `π/2 − 10⁻⁸` example is a regression guard. `R_F`'s first argument is `Cos[φ]²`, and
computing it as `1 − Sin[φ]²` gave ~100% relative error there: the result was wrong by
2.7e-08 relative, an eight-digit loss, against 6.7e-16 at a generic amplitude. Squaring
`cos` instead brings it to 1.7e-15 worst case.

- Attributes: `Listable`, `NumericFunction`, `Protected`. The real machine kernel serves the
  scalar call, packed `List`s, real `NDArray`s and `Compile[]` (scalar and rank-1); the
  quasi-period term `2k K(m)` uses the same Cephes-form `K` kernel as `EllipticK`.
- Unlike `EllipticK`/`EllipticE[m]`, there is no machine-complex kernel: a complex machine
  argument is answered correctly by Arb but as a 53-bit arbitrary-precision number.
- See also: `EllipticK`, `EllipticE`, `EllipticPi`.
