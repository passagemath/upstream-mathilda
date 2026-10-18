# Simplify

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Simplify[expr]`**

performs a sequence of algebraic and other transformations on expr and returns the simplest form it finds.

**`Simplify[expr, assum]`**

does simplification using assumptions assum.

<details>
<summary>Notes</summary>

Options: Assumptions (default $Assumptions) -- facts assumed while simplifying. ComplexityFunction (default: leaf count plus integer-digit count, matching Mathematica) -- ranks candidate forms; the lowest-scoring form is returned. TransformationFunctions (default Automatic) -- the functions applied to try to transform parts of expr. Automatic uses the built-in collection; {f1, f2, ...} uses only the fi; {Automatic, f1, ...} uses the built-in functions together with the fi. TimeConstraint (default Infinity) -- a per-sub-expression wall-clock budget in seconds; when a sub-expression's search exceeds it, the best form found so far is returned instead of continuing. A list {tLoc, ...} uses tLoc as the per-sub-expression budget. The built-in collection tries Together, Cancel, Expand, Factor, FactorSquareFree, Apart, TrigExpand, TrigFactor, and a TrigToExp/ExpToTrig roundtrip, keeping the smallest result. Under positivity / reality assumptions Simplify also applies Log/Power identities -- Log\[a b\] -\> Log\[a\] + Log\[b\], (a b)^c -\> a^c b^c, (a^p)^q -\> a^(p q), Log\[a^p\] -\> p Log\[a\] and the like -- whenever the operand-domain conditions are provable from the assumption set. Assumptions can be equations, inequalities, domain specifications such as Element\[x, Integers\], or logical combinations of these. Lists of assumptions are converted to conjunctions. Simplify automatically threads over lists, equations, inequalities, and logic functions.

</details>

## Examples (44)

Every input below was run against the current Mathilda build and its output recorded.

### Scope (28)

```mathematica
In[1]:= Simplify[(x - 1)(x + 1)(x^2 + 1) + 1]
Out[1]= x^4

In[2]:= Simplify[3/(x + 3) + x/(x + 3)]
Out[2]= 1

In[3]:= Simplify[a x + b x + c]
Out[3]= c + (a + b) x

In[4]:= Simplify[Sin[x]^2 + Cos[x]^2]
Out[4]= 1

In[5]:= Simplify[2 Tan[x]/(1 + Tan[x]^2)]
Out[5]= Sin[2 x]

In[6]:= Simplify[(E^x - E^(-x))/Sinh[x]]
Out[6]= 2

In[7]:= Simplify[{Sin[x]^2 + Cos[x]^2, 3/(x + 3) + x/(x + 3)}]
Out[7]= {1, 1}

In[8]:= Simplify[Sqrt[2] Sqrt[3]]
Out[8]= Sqrt[6]

In[9]:= Simplify[Sqrt[6] - Sqrt[2] Sqrt[3]]
Out[9]= 0

In[10]:= Simplify[1 - (-1)^(1/3) + (-1)^(2/3)]
Out[10]= 0

In[11]:= Simplify[1 - (-1)^(1/5) + (-1)^(2/5) - (-1)^(3/5) + (-1)^(4/5)]
Out[11]= 0

In[12]:= Simplify[Sqrt[3 + 2 Sqrt[2]] - (1 + Sqrt[2])]
Out[12]= 0

In[13]:= Simplify[Sin[ArcCos[x]] - Sqrt[1 - x^2]]
Out[13]= 0

In[14]:= Simplify[ArcSin[x] + ArcCos[x] - Pi/2]
Out[14]= 0

In[15]:= Simplify[Log[72] - 3 Log[2] - 2 Log[3]]
Out[15]= 0

In[16]:= Simplify[4 Sin[x]^2 Cos[x]^2 + 4 Sin[x] Cos[x] + 1]
Out[16]= (Cos[x] + Sin[x])^4

In[17]:= Simplify[2 x - 4 y + 6 z - 10 == -8]
Out[17]= x + 3 z == 1 + 2 y

In[18]:= Simplify[-2 x < 4]
Out[18]= x > -2

In[19]:= Simplify[Abs[Sin[x]], 0 < x < Pi/2]
Out[19]= Sin[x]

In[20]:= Simplify[Abs[Tan[x]] + Tan[x], Pi/2 < x < Pi]
Out[20]= 0

In[21]:= Simplify[Sqrt[x^2], x > 0]
Out[21]= x

In[22]:= Simplify[Sqrt[x^2], Element[x, Reals]]
Out[22]= Abs[x]

In[23]:= Simplify[Log[a b], a > 0 && b > 0]
Out[23]= Log[a] + Log[b]

In[24]:= Simplify[(a^p)^q, a > 0 && Element[p, Reals]]
Out[24]= a^(p q)

In[25]:= Simplify[Sqrt[x^2 y^2], x > 0 && y < 0]
Out[25]= -x y

In[26]:= Simplify[Cos[k Pi], Element[k, Integers]]
Out[26]= (-1)^k

In[27]:= Simplify[Log[E^(x + y)], {Element[x, Reals], Element[y, Reals]}]
Out[27]= x + y

In[28]:= Simplify[x > 0, x > 0]
Out[28]= True
```

### Options (4)

```mathematica
In[29]:= Simplify[1/(x - 1) + 1/(1 - x), TransformationFunctions -> {Together}]
Out[29]= 0

In[30]:= Simplify[Sin[x]^2 + Cos[x]^2, TransformationFunctions -> {Together}]
Out[30]= Cos[x]^2 + Sin[x]^2

In[31]:= Simplify[Sin[x]^2 + Cos[x]^2, TransformationFunctions -> {Automatic, Together}]
Out[31]= 1

In[32]:= Simplify[a + b, TransformationFunctions -> {(# /. a -> 0 &)}]
Out[32]= b
```

### Worked examples (2)

```mathematica
In[33]:= Simplify[Root[1+#^4&,2]^3 + Root[1+#^4&,2]]
Out[33]= I Sqrt[2]

In[34]:= Simplify[D[Integrate[Sqrt[Tan[x]], x], x] - Sqrt[Tan[x]]]
Out[34]= 0
```

### Applications (10)

Rational cancellation through the SHAPE_RATIONAL fast path

```mathematica
In[35]:= Simplify[(x^2 - 1)/(x - 1)]
Out[35]= 1 + x
```

Constant-phase angle expansion proves the triple-tangent identity

```mathematica
In[36]:= Simplify[Tan[x] + Tan[x + Pi/3] + Tan[x + 2 Pi/3] - 3 Tan[3 x]]
Out[36]= 0
```

Trig at rational multiples of Pi is folded as an exact algebraic number

```mathematica
In[37]:= Simplify[Tan[Pi/7] Tan[2 Pi/7] Tan[3 Pi/7]]
Out[37]= Sqrt[7]
```

Machin's formula, certified by exp(i e) == 1 in qqbar

```mathematica
In[38]:= Simplify[4 ArcTan[1/5] - ArcTan[1/239] - Pi/4]
Out[38]= 0
```

Integer-assumption rules: parity, roots of unity, floor division, Gamma ratios, coprimality

```mathematica
In[39]:= Simplify[{(-1)^(n^2), I^(4 n + 2), Floor[n/2] + Ceiling[n/2], Gamma[n + 2]/Gamma[n], GCD[n, n + 1]}, Element[n, Integers]]
Out[39]= {(-1)^n, -1, n, n (1 + n), 1}
```

The Reduce-backed sign oracle resolves both absolute values

```mathematica
In[40]:= Simplify[Abs[x - 2] + Abs[x + 3], 0 < x < 1]
Out[40]= 5
```

Multiple-angle collapse, emitted only after Reduce proves it on the region

```mathematica
In[41]:= Simplify[ArcCos[1 - 2 x^2], 0 < x < 1]
Out[41]= 2 ArcSin[x]
```

Branch-gated power fold: Im of the exponent stays in the principal strip

```mathematica
In[42]:= Simplify[Sqrt[E^(2 I x)], -Pi/2 < x < Pi/2]
Out[42]= E^(I x)
```

Honest decline: Cos changes sign at Pi/2

```mathematica
In[43]:= Simplify[Abs[Cos[x]], 0 < x < Pi]
Out[43]= Abs[Cos[x]]
```

A scalar TimeConstraint bounds the whole call

```mathematica
In[44]:= Simplify[Expand[(x + y)^20]/(x + y)^18, TimeConstraint -> 5]
Out[44]= (x + y)^2
```

## Options & behaviour

### Specialised transformations

Beyond the generic algebraic transforms, `Simplify` reaches a number of
canonical forms that the generic pipeline alone does not:

- **Cross-base radical fusion** — distinct positive-integer radicals sharing an
  exponent are combined inside a `Times`.
- **Roots of unity** — `(-1)^(p/q)` and `E^(I p Pi/q)` atoms are reduced modulo
  the relevant cyclotomic polynomial.
- **Exact trig at rational multiples of Pi** — `Sin`/`Cos`/`Tan`/`Cot`/`Sec`/`Csc`
  of `r Pi` (`r` rational) is an algebraic number (built from the root of unity
  `exp(I pi r)`), extending exact closed forms past the small-denominator tables.
  A variable-free product or combination folds to its exact value:
  `Cos[Pi/7] Cos[2 Pi/7] Cos[3 Pi/7] -> 1/8`,
  `Tan[Pi/7] Tan[2 Pi/7] Tan[3 Pi/7] -> Sqrt[7]`, and `RootReduce[Sin[Pi/7]]`
  is an exact `Root[...]`.
- **Radical denesting** — `Sqrt[A + Sqrt[B]]` and cube-root towers collapse via
  the half-sum identity when the result is cleaner.
- **Inverse trig / hyperbolic identities** — standard relations such as
  `Sin[ArcCos[x]] == Sqrt[1 - x^2]` reduce, and the complementary-angle pairs
  `ArcSin[x]+ArcCos[x]`, `ArcTan[x]+ArcCot[x]`, `ArcSec[x]+ArcCsc[x]` all fold to
  `Pi/2`.
- **Inverse-trig rational-angle combinations** — a `Z`-linear combination of
  `ArcTan`/`ArcSin`/`ArcCos` of constants plus a rational multiple of `Pi` that
  is identically zero is decided exactly: `exp(I e)` is built from the Euler
  closed forms and proven to equal `1` over the algebraic numbers (so
  `e == 0 mod 2 Pi`), with a numeric screen only selecting the `2 Pi` branch.
  Machin-type identities reduce — `4 ArcTan[1/5] - ArcTan[1/239] - Pi/4 -> 0`,
  `ArcSin[3/5] + ArcSin[5/13] - ArcSin[56/65] -> 0`. The hyperbolic constant
  combinations (`ArcSinh`/`ArcCosh`/`ArcTanh`) reduce likewise via `e^e == 1`.
- **Conditional identities under `Assumptions`** — identities true only on a
  region/domain reduce when the assumptions establish it, via general machinery
  (never ad-hoc), with a `Reduce`/CAD entailment bridge for the sign/domain
  conditions the lightweight provers cannot settle:
  - radical-product combine `Sqrt[a] Sqrt[b] -> Sqrt[a b]` when a base is
    provably non-negative (`Simplify[Sqrt[x-1] Sqrt[x+1] - Sqrt[x^2-1], x>1] -> 0`)
  - inverse-of-forward `ArcTanh[Tanh[x]] -> x` (`x in Reals`), the hyperbolic
    analogue of `ArcTan[Tan[x]] -> x` on `-Pi/2 < x < Pi/2`; also `ArcSinh[Sinh]`,
    `ArcCsch[Csch]`, `ArcCoth[Coth]` (real argument) and `ArcCosh[Cosh[x]] -> |x|`,
    `ArcSech[Sech[x]] -> |x|` (sign from the assumptions)
  - **hyperbolic sign oracle** — the odd hyperbolics are sign-preserving, so
    `sign(Sinh[u]) = sign(Tanh[u]) = sign(Coth[u]) = sign(Csch[u]) = sign(u)` and
    `Sech[u] > 0`, feeding the principal-root / `Abs` collapse: `Sqrt[Cosh[x]^2-1]
    -> Sinh[x]` and `Sqrt[1-Sech[x]^2] -> Tanh[x]` on `x>0`, with the mirror sign
    flips on `x<0` (`Sqrt[Sinh[x]^2] -> -Sinh[x]`); and the half-angle
    `Sqrt[(Cosh[x]+1)/2] -> Cosh[x/2]`
  - **imaginary-period phase** — the hyperbolic functions have an imaginary
    period, so `Sinh[x + 2 n Pi I] -> Sinh[x]`, `Tanh[x + n Pi I] -> Tanh[x]`,
    `Cosh[n Pi I] -> (-1)^n` for a symbolic integer `n` (a real multiple of `Pi`
    is correctly *not* a period). These, with the complex-component and
    inverse-hyperbolic cases above, form the assumed-hyperbolic stress corpus
    (`tests/test_hyp_assume_simplify.c`), the sibling of the assumed-trig corpus.
  - `Abs[z]^2` / `Re` / `Im` via `ComplexExpand` when the variables are real
  - inverse-function **addition** identities on a region, e.g.
    `ArcTan[x]+ArcTan[y]-ArcTan[(x+y)/(1-x y)] -> 0` on `-1<x<1 && -1<y<1`, and the
    `ArcTanh`/`ArcSinh`/`ArcCosh` analogues — decided by a derivative-constancy
    engine (gradient `== 0`, branch-cut-free on the connected region by `Reduce`,
    value `0` at a sample point) or a pointwise forward-function recognizer.
    These are sound: nothing reduces without the region established (a wrong or
    too-large region, e.g. crossing `x y = 1`, is left unchanged).
- **Logarithm simplification** — `Log` of a positive rational is decomposed over
  its prime factors, and linear combinations of logs are fused
  (`Sum c_i Log[a_i] -> Log[Prod a_i^c_i]`).
- **Log-power symmetry under positivity** — `base^exp -> Exp[exp Log[base]]`
  when the base and every single-arg `Log` argument in the exponent are provably
  positive, exposing the symmetry `x^Log[y] = y^Log[x]` so it cancels
  (`Simplify[x^Log[y] - y^Log[x], x>0 && y>0] -> 0`). Kept only on a strict
  complexity win, so a standalone `x^Log[x]` is left unchanged; the two-sided
  positivity gate is required — a one-sided assumption (only `y>0`) correctly
  does **not** reduce, since the identity fails at `x<=0`.
- **Constant complex powers** — a variable-free `Power[c1, c2]` whose base is not
  a positive real folds to its principal value `Exp[c2 Log[c1]]`:
  `Simplify[I^I - E^(-Pi/2)] -> 0`, `Simplify[(-1)^I E^Pi] -> 1`. Taken only on a
  whole-expression complexity win, so `2^I`, standalone `I^I`, and surds such as
  `(-1)^(1/3)` are preserved.
- **Inverse-function logarithmic forms** — `ArcCos` joins `ArcSin`/`ArcSinh`/
  `ArcTanh` in the `TrigToExp`/`ExpToTrig` table, so
  `Log[x + I Sqrt[1-x^2]] = I ArcCos[x]` is recognised (principal-value general,
  unconditional). `ArcCosh[u] = Log[u + Sqrt[u^2-1]]` is recognised **under
  `u >= 1`** — the combined radical is branch-correct only there, so the gate is
  mandatory (it is left unreduced for `u <= -1`).
- **Range-gated inverse-of-direct trig** — `ArcSin[Sin[x]] -> x` on
  `-Pi/2 <= x <= Pi/2`, `ArcCos[Cos[x]] -> x` on `0 <= x <= Pi`,
  `ArcCot[Cot[x]] -> x` on `0 < x < Pi` (each inverse is a bijection from exactly
  that strip; mirrors the existing `ArcTan[Tan[x]]`). Off-range inputs are left
  unchanged.
- **Root-of-unity powers under integer assumptions** — a root of unity raised to
  an integer-valued exponent is reduced modulo the base's order. The exponent is
  reduced by a Boolean-ring mod-2 reduction (`F2[g]/(g^2-g)`: coefficients mod 2,
  `g^e -> g`, like-monomials XOR — sound because `g^e == g (mod 2)` for integer
  `g`), with `I = (-1)^(1/2)`, `-I = (-1)^(-1/2)` and `E^(I Pi X) = (-1)^X` routed
  to the same engine. So for integer `n, m`: `(-1)^(2n+1) -> -1`, `(-1)^(3n) ->
  (-1)^n`, `(-1)^(n^2) -> (-1)^n`, `(-1)^(n(n+1)) -> 1`, `I^(4n) -> 1`,
  `I^(2n) -> (-1)^n`, `(-I)^(2n) -> (-1)^n`, `(-1)^Mod[n,2] -> (-1)^n`, and
  `(-1)^(n-m) -> (-1)^(n+m)`. It declines (leaves the power unevaluated) on any
  exponent it cannot express as a polynomial in provably-integer generators, and
  off the integer assumption (e.g. `(-1)^(2n)` does **not** fold for real `n`).
- **Pi-periodicity with a compound integer coefficient** — `Sin[K Pi] -> 0`,
  `Cos[K Pi] -> (-1)^K`, and the shift `Sin[t + K Pi] -> (-1)^K Sin[t]` fire for
  any provably-integer `K` (not only a bare symbol), so e.g.
  `Sin[(n+m) Pi] -> 0` and `Cos[(n+m) Pi] - (-1)^(n+m) -> 0`. The emitted
  `(-1)^K` is then canonicalised by the root-of-unity reducer above. Also
  `E^(2 I Pi n + z) -> E^z` and `E^(I Pi n + z) -> (-1)^n E^z`.
- **Symbolic GCD/LCM and orthogonal-polynomial values under integer
  assumptions** — for integer `a, b, n`: `GCD[a,a] = Abs[a]`, `GCD[a,1] = 1`,
  `LCM[a,1] = Abs[a]`, `GCD[a,b] = 1` when `|a-b| = 1` (so `GCD[n,n+1] = 1`), and
  `LegendreP[n,-1] = (-1)^n` (`LegendreP[n,1] = 1` is unconditional). The core
  GCD/LCM evaluators are unchanged, so these fire only under the assumption.
- **Gamma / factorial / combinatorics ratios under integer assumptions** — the
  Factorial-ratio engine first converts `Gamma[z] -> Factorial[z-1]`,
  `Pochhammer[a,m] -> Factorial[a+m-1]/Factorial[a-1]`, and
  `Binomial[p,q] -> Factorial[p]/(Factorial[q] Factorial[p-q])` (restoring the
  Gamma form on output), so for integer `n, k`: `Gamma[n+1] = n!`,
  `Gamma[n+2]/Gamma[n] = n(n+1)`, `Binomial[n,k] = Binomial[n,n-k]`, Pascal's
  rule, `Pochhammer[1,n] = n!`, and `Factorial2[2n] = 2^n n!`. The conversion is
  gated on an active integer-assumption context so unassumed `Simplify` and the
  Sum engine are unaffected.
- **Rounding and modular reductions under integer assumptions** — for integer
  `n`: `Floor[n] -> n`, `Ceiling[n] -> n`, `Round[n] -> n`,
  `FractionalPart[n] -> 0`, `Mod[n, 1] -> 0`; integer-shift extraction
  `Floor[n + x] -> n + Floor[x]` (and `Ceiling`/`Round`/`FractionalPart`)
  floor-division `Floor[P/d] -> (P - Mod[P,d])/d` so
  `Floor[n/2] + Ceiling[n/2] -> n`; divisibility `Mod[n m, m] -> 0`; and
  parity-preserving `Mod[n^2, 2] -> Mod[n, 2]`. Each declines off the integer
  assumption (`Floor[n]` and `Mod[n,1]` are left alone for real `n`).
- **Branch-gated power of an exponential** — `(E^w)^r -> E^(r w)` when `Im[w]` is
  provably in `(-Pi, Pi]` (so the principal `Log` recovers `w`), e.g.
  `Sqrt[E^(2 I x)] -> E^(I x)` under `-Pi/2 < x < Pi/2`. This is the sound,
  assumption-gated counterpart of `PowerExpand`'s unconditional collapse, so it
  does **not** reduce off the strip.
- **Pythagorean completion and reduction** for trig and hyperbolic squares.
- **Exact trig/exp zero-recognition** — a `Plus` that is a rational function of a
  single exponential kernel `t = E^(I x)` and is identically zero (canonically a
  Risch antiderivative diff-back `D[G] - f`) is proven `0` by exact rational
  point-evaluation on a Nullstellensatz grid — no numeric sampling, no slow
  `Together`. When that rigorous test declines — because of bare polynomial
  dependence on the kernel variable (e.g. the `x` in the `x E^x Sin[x]`
  diff-back) or mixed real+imaginary exponential kernels
  (`E^((1+I) x) = E^x E^(I x)`) — an exact `TrigToExp`-collapse fallback catches
  the identity: `Simplify[D[Integrate[x E^x Sin[x], x], x] - x E^x Sin[x]] -> 0`,
  and angle-addition identities such as
  `Sin[x] Cos[y] + Cos[x] Sin[y] - Sin[x + y] -> 0` collapse too. A circular trig
  head with an affine argument `k x + c Pi` (`c` rational) is handled by
  expanding the constant phase out of the kernel (the `c Pi` part contributes
  exact root-of-unity coefficients), so `Tan[Pi/2 - x] - Cot[x] -> 0`,
  `Tan[x] + Tan[x + Pi/3] + Tan[x + 2 Pi/3] - 3 Tan[3 x] -> 0`, and
  `Sin[x] Sin[Pi/3 - x] Sin[Pi/3 + x] - 1/4 Sin[3 x] -> 0` (and `PossibleZeroQ`
  agrees). The same applies to a hyperbolic head with an affine *imaginary*
  phase `k x + i c Pi` (`Cosh[i c Pi] = Cos[c Pi]`, `Sinh[i c Pi] = i Sin[c Pi]`
  fold to the same root-of-unity coefficients), so `Tanh[x + I Pi] - Tanh[x]
  -> 0`, `Tanh[I Pi/2 - x] + Coth[x] -> 0`, and `Tanh[x] + Tanh[x + I Pi/3] +
  Tanh[x + 2 I Pi/3] - 3 Tanh[3 x] -> 0`.
- **Trig / radical-trig rational normal form** — rational functions of trig and
  hyperbolic kernels are reduced to a canonical fraction modulo the Pythagorean
  ideal. A quadratic radical of a kernel (e.g. `Sqrt[Tan[x]]`, `Tan[x]^(3/2)`) is
  carried as an algebraic generator `l` with `l^2 = g`, so rational functions of
  `Sqrt[Tan[x]]` reduce too — `Simplify[Tan[x]/Sqrt[Tan[x]]] -> Sqrt[Tan[x]]`,
  and `D[Integrate[Sqrt[Tan[x]], x], x] // Simplify -> Sqrt[Tan[x]]`. Radicands
  that are rational with an *odd* generator in the denominator (`Cot = Cos/Sin`,
  `Csc = 1/Sin`) are handled too: the inverse odd-generator powers the relation
  injects are cleared before the denominator is rationalised, so
  `Simplify[Cot[x]/Sqrt[Cot[x]]] -> Sqrt[Cot[x]]`.
- **Multi-generator radical rational normal form** — a rational function of two
  or more distinct radical bases (e.g. `a^(1/3)` and `(a+b x)^(1/3)`) is reduced
  in the quotient ring `K[g_1, ..., g_n] / <g_k^{q_k} - base_k>`: each base is
  carried as an algebraic generator, the terms are combined over a common
  denominator, reduced modulo the generator relations, and the denominator is
  rationalised, before the radicals are substituted back. This recovers cross-base
  cancellations the single-generator `Together`/`Cancel` path cannot — e.g.
  `D[Integrate[1/(x^3 (a+b x)^(1/3)), x], x] // Simplify -> 1/(x^3 (a+b x)^(1/3))`.
  **Algebraic *constant* radicals** such as `Sqrt[2]` are collected as generators
  too (with the relation `s^q - c`, e.g. `s^2 - 2`), so a rational function of `x`,
  one `x`-dependent radical, and a constant radical in its coefficients — the shape
  of an antiderivative re-differentiated against its integrand — reduces too:
  `D[Integrate[(x^2+1)/(x^3 Sqrt[2x^4-2x^2+1]), x], x] - integrand // Simplify ->
  0` (previously this fell into the algebraic-field `Together` over
  `Q(x)[R]/(R^2-q)` and did not terminate). The pass engages only when at least
  one base carries a free symbol — a purely numeric radical identity is left to the
  `RootReduce`/qqbar pass — and the result is adopted only when its `SimplifyCount`
  strictly improves, so it never regresses a case.
- **Equation / inequality rebalancing** — a binary relation is normalised by
  dividing through the GCD of integer coefficients and partitioning terms across
  the relation; the rebalanced form is kept when its `SimplifyCount` is lower.
  Strict inequalities flip when divided by a negative.

### Assumptions

`Simplify[expr, assum]` simplifies under `assum`, which may be equations,
inequalities, domain specifications such as `Element[x, Integers]`, or logical
combinations of these. A list `{a1, a2, ...}` is treated as the conjunction
`And[a1, a2, ...]`. Under provable positivity / reality, `Simplify` applies
`Log`/`Power` identities — `Log[a b] -> Log[a] + Log[b]`, `(a b)^c -> a^c b^c`,
`(a^p)^q -> a^(p q)`, `Log[a^p] -> p Log[a]` and the like — whenever the
operand-domain conditions follow from the assumption set. Per-symbol sign facts
drive `Sqrt[x^2] -> x` / `-x` / `Abs[x]`, and integer facts drive the
`Sin[n Pi] -> 0`, `Cos[n Pi] -> (-1)^n` family.

A **deep sign oracle** extends `Abs[g] -> ±g` and `Sqrt[g^2] -> ±g` to an
*arbitrary* real `g` whose sign is fixed on the assumed region — not just bare
symbols. It layers the structural provers, a decomposition of the pole-bearing
trig family (`sign[Tan] = sign[Sin] sign[Cos]`, `sign[Sec] = sign[Cos]`, …), and
a sound Reduce/CAD base case (a sign is proved only when the opposite strict
inequality is unsatisfiable on the region). It declines — leaving `Abs`/`Sqrt`
intact — wherever the sign is not constant (e.g. `Abs[Cos[x]]` on `0 < x < Pi`,
which flips at `Pi/2`, or an unbounded region).

When no positional assumption and no `Assumptions` option are given, `Simplify`
reads the current value of `$Assumptions`.

A predicate that appears literally among the assumed facts folds to `True`:

### Options

- **`Assumptions`** (default `$Assumptions`) — the facts assumed while
  simplifying. An explicit `Assumptions -> X` overrides `$Assumptions`; a
  positional assumption is conjoined with `$Assumptions`.
- **`ComplexityFunction`** (default the built-in `SimplifyCount` measure) — ranks
  candidate forms; `Simplify` returns the lowest-scoring one. A custom function
  `f` must return an integer or bigint for `f[candidate]`; otherwise the default
  is used. `ComplexityFunction -> Automatic` is a synonym for the default and
  takes the fast native scoring path. Compared with the default, `LeafCount`
  drops the integer-digit penalty.
- **`TransformationFunctions`** (default `Automatic`) — the functions applied to
  try to transform parts of `expr` (see [TransformationFunctions](#transformationfunctions)).
- **`TimeConstraint`** (default `Infinity`) — a wall-clock budget in seconds.
  A scalar `t` bounds the **whole call**: the budget is honoured on every path —
  the heuristic search, the specialised rational/polynomial/log-exp pipelines
  (including the SHAPE_RATIONAL input that dispatches straight to a pipeline),
  the bottom-up descent, and the seed phase — so the entire `Simplify[expr, …]`
  returns the best form found so far once `t` is exhausted. A list
  `{tLoc, tTot}` additionally caps each individual sub-expression's search at
  `tLoc` while the whole call is capped at `tTot` (matching `FullSimplify`); a
  one-element `{tLoc}` sets the per-sub-expression budget only. The check is
  synchronous between search steps, so it fails gracefully with no memory leak,
  but it does **not** interrupt a single long-running kernel call mid-flight (a
  pathological single `Together`/`Factor` is instead bounded by the poly engine's
  own degree/size guards). `Infinity` (the default) imposes no limit.

## Performance

Against other systems, from the benchmark suite (same input, results cross-checked for agreement):

| case | Mathilda | Wolfram | Python |
|---|---:|---:|---:|
| Discriminant of deg 20 | 2.51 s | 0.068 s | 0.182 s |
| Simplify quartic-to-Cos[4x] | 2.35 s | 0.003 s | 7.45 s |
| TrigReduce product of 4 sines | 0.894 s | 0.14 s | 16.5 s |
| Expand (1+x)^400 | 0.434 s | 0.107 s | 0.003 s |
| Cancel deg-60 over deg-58 | 0.337 s | 0.569 s | 7.37 s |
| PolynomialGCD, coprime deg 40 | 0.252 s | 0.087 s | 0.334 s |

## Implementation notes

**Algorithm.** Simplify is a complexity-weighted, memoized candidate-set search
over the existing battery of algebraic transforms. `builtin_simplify`
(`simp_builtins.c`) parses options — a positional or `Assumptions ->`
assumption (combined into the `$Assumptions` default with `And`), a
`ComplexityFunction`, and `TransformationFunctions` — builds an `AssumeCtx` of
normalised facts, then drives the search. Inexact inputs are rationalised on
entry and numericalised on exit (`internal_rationalize_then_numericalize`).
`Equal`/`Less`/`And`/... heads are threaded manually (with a relational
rebalancing candidate). The head carries only `Protected`. A `List`
argument is simplified as one expression, except inexact input, which
`simplify_thread_list` threads element by element so that numeric
contagion does not cross list boundaries.

The core driver is `simp_bottomup` (`simp_bottomup.c`), which descends into
every `Plus`/`Times`/`Power` child (memoizing each result in a `SimpMemo` hash
table) and dispatches each node through `simp_dispatch` → `simp_search`.
`simp_classify` routes pure-polynomial/rational shapes to dedicated pipelines
(`simp_pipeline_polynomial`/`_rational`/`_logexp`) and there are top-level fast
paths: a `SHAPE_RATIONAL` shortcut that runs `Together`/`Cancel`/`Factor` once
at the top, an algebraic-tower `Together[expr, Extension -> Automatic]` collapse,
and a `simp_trig_rational` substitution that maps trig/opaque subtrees to ground
symbols and works in the quotient ring.

`simp_search` is the heuristic engine. It seeds a `CandSet` with the input plus
the output of a long list of correctness-preserving rewriters (assumption rules,
log/exp identities, `SimpLogRules`, trig roundtrip, ExpToTrig, Pythagorean
square-completion / reduction / canonicalisation, trig-at-rational-Pi, tan-
addition, half-angle, radical product combine, sqrt/cube-root/algebraic
denesting, roots-of-unity, factorial decomposition, per-variable Collect). It
then runs `SIMP_ROUNDS` (= 2) rounds in which every seed is fed through the
`SIMP_TRANSFORMS` table — `Together, Cancel, Expand, ExpandNumerator,
ExpandDenominator, Factor, FactorSquareFree, FactorTerms, Apart, TrigExpand,
TrigFactor, TrigReduce, TrigToExp` — plus chained Pythagorean/radical/trig-Pi
passes. Each transform call is gated by `transform_can_fire` (cheap precondition
check) and the running best is updated by `update_best` against the complexity
score; candidates strictly worse than their parent are dropped from the next
round (with a loosened `2*parent + 8` bound for TrigExpand and the seed-phase
blow-up guard). A short-circuit (`simp_best_is_zero`) exits as soon as a literal
`0` is reached. Final polish passes apply `simp_lift_common_factor`,
`transform_pythag_reduce`, and `canon_negate_pairs`.

After the search, `builtin_simplify` runs a fixed tail of exact folds, each kept
only when it does not raise the score: a constant-algebraic fold
(`flint_qqbar_is_constant_algebraic` / `flint_qqbar_reduce_coeffs`, which treats
trig at any rational multiple of `Pi` as an algebraic number through the root of
unity `exp(i pi r)`, so `Tan[Pi/7] Tan[2 Pi/7] Tan[3 Pi/7]` becomes `Sqrt[7]`);
`simp_const_cpow_fold` (principal value of a constant power with a
non-positive-real base, so `I^I - E^(-Pi/2)` goes to `0`); `simp_invtrig_combo_is_zero` /
`simp_invhyp_combo_is_zero` (a Z-combination of inverse trig/hyperbolic
functions of constants is proved a multiple of `2 Pi` by checking `exp(i e) == 1`
in qqbar, then a numeric screen picks the branch, as in Machin's formula); and, under
assumptions, `simp_invhyp_addition_is_zero` and `simp_region_identity_is_zero`,
which prove an inverse-function identity on a box region from
derivative constancy. Every partial derivative must be exactly 0, every
branch-sensitive subterm is shown by `Reduce` to stay off its cut, and the
value at a `FindInstance` point must be exactly 0. Last comes a confluence
post-pass. When the result still holds a rounding or `Mod` head, it runs
`apply_assumption_rules` to a fixpoint, so `Simplify[Floor[n], n ∈ Integers]`
agrees with `Refine`. A `transform_trigexp_vanish` fast path ahead of the
bottom-up descent proves rational functions of trig/exp kernels identically
zero (`E^(k I x) -> t^k`, exact numerator test), so they never reach the
per-node search.

**Assumption engine.** `AssumeCtx` facts feed two layers. The first is the fast
structural provers (`prov_pos`, `prov_int`, `prov_re`, …) and
`assume_structural_rewrite` (`simp_assume_rewrite.c`). That file holds the
integer-assumption rules: the Boolean-ring parity reduction `simp_reduce_mod2`
and `canon_neg1_pow` route every root of unity `(-1)^Y`, `I^Y` and
`E^(I Pi Y)` to a reduced sign. `trig_pi_shift_symbolic` gives
`Sin[t + k Pi] -> (-1)^k Sin[t]` for a provably-integer `k`. It also holds the
floor-division and `Mod` reductions, `GCD`/`LCM`/`LegendreP[n,-1]` folds,
Gamma/Binomial/Pochhammer ratios via `simp_factorial`, and inverse-of-forward
collapses on a known branch: `invtrig_branch_reduce` takes a two-sided numeric
interval and verifies it lies in a single branch. `invtrig_multiangle_reduce`
emits `ArcCos[1-2x^2] -> 2 ArcSin[x]` only on a `Reduce` proof. The same file
also holds the Schwarz reflection `Conjugate[H[z]] -> H[Conjugate[z]]` for
branch-cut-free `H`. The second layer
is a `Reduce`/CAD entailment oracle (`assume_reduce_entails`, as used by
`Refine`). A fact `P` holds if `Reduce[A && !P]` is `False`. Behind it,
`assume_sign_nonneg_deep` / `assume_sign_nonpos_deep` decide the sign of an
arbitrary subterm, splitting `Tan`/`Sec`/`Csc` into `Sin`/`Cos` signs first. That
is what resolves `Abs[g] -> ±g` and `Sqrt[g^2] -> ±g` (`simp_logexp_abs.c`,
with `sqrt_prepare_radicand` turning `1 + Cos[2 u]` into `2 Cos[u]^2` inside a
radical). Assumption-gated power rules (`transform_powbase_to_exp_assuming`,
`simp_radical_combine_assuming`, `(a^p)^q -> a^(p q)` only when `p` is
provably real) complete the set. Every gated rule declines off its
assumption, and no assumption path uses `PossibleZeroQ`.

**TimeConstraint.** A scalar `TimeConstraint -> t` (or `{tLoc, tTot}`) arms a
whole-call deadline, `g_simp_call_deadline` (`simp_arm_call_deadline`). It is
checked in `simp_dispatch`, the specialised pipelines, `simp_bottomup`, the seed
phase and the round loop. A nested `Simplify` can only tighten the deadline. The
check is synchronous: it cannot interrupt a single running kernel call. The poly
layer's degree-monotonicity guards (`pseudo_rem_standard`) cover that case.

**Complexity scorer.** Candidates are ranked by `score_with_func`: with no user
`ComplexityFunction` this is `simp_default_complexity` (the SimplifyCount metric
— LeafCount with integers contributing their decimal-digit count) plus a
`nested_radical_penalty` (+3 per truly-nested radical). A user function is
evaluated as `f[candidate]`; an Integer result is used directly, a BigInt scores
`SIMP_SCORE_INF`, anything else falls back to the default. Lower wins; ties
favour the form that reached the score plateau first (force-take semantics let
assumption/log/factorial rewrites win even on a tie).

**Data structures.** `CandSet` (dynamic `Expr*` array, dedup via `expr_eq`, capped
at `SIMP_CAND_CAP` = 12); `SimpMemo` (256-bucket chained hash of input→best for
the bottom-up driver); `AssumeCtx` (flat fact array from `assume_ctx_from_expr`);
a per-Simplify `FactorMemo` shared by Factor/Trig* so duplicate subexpression
work hits the cache. Each transform is invoked as `f[candidate]` through the real
evaluator (`traced_call_unary`), so Simplify composes the same builtins exposed
in the REPL.

**Limits.** Simplify is not a complete decision procedure. The `Reduce`
entailment is capped by variable count, and a decline leaves the subterm
alone, so wins under assumptions depend on what the provers and CAD can
settle. Because the search is scored, a correct rewrite that does not lower
`SimplifyCount` can be dropped. On `Pi/2 < x < Pi`, `ArcSin[Sin[x]]` stays
as it is, while `Refine` returns `Pi - x`. Simplify is symbolic and has no
NDArray or `Compile[]` path.

- `Protected`. **Not** `Listable`: a `List` in the assumption position is a
  conjunction of facts (see below), not a threading axis.
- The built-in transformation collection tries `Together`, `Cancel`, `Expand`,
  `Factor`, `FactorSquareFree`, `Apart`, `TrigExpand`, `TrigFactor`, a
  `TrigToExp`/`ExpToTrig` roundtrip, and per-variable `Collect`, keeping the
  smallest result.
- The default complexity measure is `SimplifyCount` — total subexpression count
  plus the decimal-digit count of integer leaves — so `100 Log[2]` is preferred
  over its expanded `Log[2^100]` form.
- **`Root[...]` objects** are treated as the constant algebraic numbers they are:
  a qqbar pre-pass canonicalises constant-algebraic subexpressions
  (`Simplify[Root[1+#^4&,2]^3 + Root[1+#^4&,2]]` → `I Sqrt[2]`), and an expression
  whose coefficients are `Root` objects is simplified without the multivariate
  blow-up that treating each `Root` as a polynomial generator once caused. A
  rational-function identity with `Root` coefficients collapses to `0` when the
  `Root`s are radical-expressible (degree ≤ 4, or a binomial); e.g.
  `Simplify[D[Integrate[Sqrt[Tan[x]], x], x] - Sqrt[Tan[x]]]` → `0`.
- Threads manually over `List`, `Equal`, `Unequal`, `Less`, `LessEqual`,
  `Greater`, `GreaterEqual`, `And`, `Or`, and `Not`, carrying any options through
  into each sub-call.

**Attributes:** `Protected`.

## References

**See also:** [List](../../other-advanced/List/), [Together](../../algebra/Together/), [Cancel](../../algebra/Cancel/), [Expand](../../algebra/Expand/), [Factor](../../algebra/Factor/), [FactorSquareFree](../../algebra/FactorSquareFree/), [Apart](../../algebra/Apart/), [TrigExpand](../../elementary-functions/TrigExpand/)

- Joel S. Cohen, *Computer Algebra and Symbolic Computation: Mathematical Methods* (A K Peters, 2003).
- Source: [`src/simp/simp_builtins.c`](https://github.com/stblake/mathilda/blob/main/src/simp/simp_builtins.c)
- Specification: [`docs/spec/builtins/simplification.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/simplification.md)
- Tests: [`tests/test_assuming.c`](https://github.com/stblake/mathilda/blob/main/tests/test_assuming.c)
- Tests: [`tests/test_characteristicpolynomial.c`](https://github.com/stblake/mathilda/blob/main/tests/test_characteristicpolynomial.c)
- Tests: [`tests/test_cherry_dilog.c`](https://github.com/stblake/mathilda/blob/main/tests/test_cherry_dilog.c)
- Tests: [`tests/test_cherry_dilog_exp.c`](https://github.com/stblake/mathilda/blob/main/tests/test_cherry_dilog_exp.c)

## Notes & additional examples

### Notes

- **Search, not rewriting.** `Simplify` seeds a candidate set with the input
  and the output of many exact rewriters. It then runs `Together`, `Cancel`,
  `Expand`, `Factor`, `Apart`, `TrigExpand`, `TrigReduce`, `TrigToExp` and
  related transforms on each candidate, and returns the lowest-scoring form
  (`SimplifyCount` by default, or any `ComplexityFunction`). A correct
  rewrite that does not lower the score can still lose: on
  `Pi/2 < x < Pi`, `Simplify[ArcSin[Sin[x]]]` stays as it is, while
  `Refine` returns `Pi - x`.
- **Assumptions.** The assumption argument is a positional second argument
  or `Assumptions -> …`; a positional one is conjoined with `$Assumptions`
  (and with `Assuming`). The rules that depend on assumptions are sound by
  construction. Each needs its assumption proved, by the structural provers
  or a `Reduce`/CAD entailment, and declines otherwise: `(-1)^(2 n)` does not
  reduce for real `n`, and `Abs[Cos[x]]` on `0 < x < Pi` is left alone.
  `(a^p)^q -> a^(p q)` needs both `a > 0` and `p` real. Integer-only folds such as `Gamma[n + 2]/Gamma[n] -> n (1 + n)` fire only
  under an integer assumption.
- **Not a decision procedure for predicates.** `Simplify[x > 0, x > 1]`
  stays `x > 0`. Use `Refine` (or `Reduce`) to decide a relation under
  assumptions.
- **TimeConstraint.** A scalar `t` caps the whole call. `{tLoc, tTot}` is
  also accepted, with `tTot` as the whole-call cap. The default is
  `Infinity`.
- **Threading.** `Simplify` threads over lists, equations and logical
  connectives itself. Its only attribute is `Protected`. Inexact input is rationalised, simplified, and
  numericalised again.
- **No numeric surfaces.** This head is symbolic, with no NDArray or
  `Compile[]` path.
- See also: `FullSimplify`, `Refine`, `PossibleZeroQ`, `Together`, `Factor`,
  `TrigReduce`, `Assuming`.
