# Integrate

!!! warning "Status: Partial"
    implemented with documented limitations or caveats; some argument forms fall through to symbolic/unevaluated output.

## Description

**`Integrate[f, x]`**

gives the indefinite integral of f with respect to x.

**`Integrate[f, {x, xmin, xmax}]`**

gives the definite integral by the

**`Integrate[f, {x, xmin, xmax}, {y, ymin, ymax}, ...]`**

gives the iterated

**`Integrate[f, x, Method -> "<name>"]`**

dispatches directly to a single

**`Integrate[f, x, TimeConstraint -> t]`**

gives the cascade a wall-clock budget

<details>
<summary>Notes</summary>

fundamental theorem of calculus (Method -\> "NewtonLeibniz"). multiple integral (innermost/last spec integrated first; inner bounds may depend on outer variables).  See also Integrate\`SingularPoints. subroutine, bypassing the default cascade.  Accepted method names:

```text
"Automatic"          — the full dispatch cascade (default)
```

"BronsteinRational"  — Integrate\`BronsteinRational (polynomial / rational) "DerivativeDivides"  — Integrate\`DerivativeDivides (substitution u(x); direct + Eliminate/Solve)

```text
"LinearRadicals"     — Integrate`LinearRadicals (rationalise radicals of a x + b)
```

"QuadraticRadicals"  — Integrate\`QuadraticRadicals (Euler substitution for Sqrt\[a x^2 + b x + c\]) "LinearRatioRadicals" — Integrate\`LinearRatioRadicals (rationalise radicals of (a x + b)/(c x + d)) "ChebychevAlgebraic" — Integrate\`ChebychevAlgebraic (binomial x^p (a x^r + b)^q via Chebychev's theorem)

```text
"GoursatAlgebraic"   — Integrate`GoursatAlgebraic (pseudo-elliptic F/R^p, p in {1/2,1/3,2/3,1/4,3/4}, via Mobius eigendescent)
"Weierstrass"        — Integrate`Weierstrass (continuous tan(x/2) / tanh(x/2) substitution)
"RischTranscendental"       — Integrate`RischTranscendental (recursive transcendental Risch; correct by construction)
"CRCTable"           — Integrate`CRCTable (lazy-loaded CRC integral table)
```

"ParallelMixedTower" — Integrate\`ParallelMixedTower (parallel Risch-Norman over a simple radical in a mixed transcendental tower; Blake II) "ParallelMixedSpecial" — Integrate\`ParallelMixedSpecial (the same over a mixed tower with special-function kernels: Ei, li, erf, incomplete Gamma, Si, Ci, elliptic F/E/Pi)

```text
"Undefined"          — Integrate`Undefined (unknown functions u[x], u'[x]; Roach §1.7)
"NewtonLeibniz"       — real definite integrals via F(b)-F(a) (implicit for the {x,a,b} form)
"LineIntegral"        — complex contour integrals (implicit for the {x,z0,...,zn} form)
"Residue"             — improper/periodic real definite integrals by the residue theorem
                      (rational/Fourier on (-Inf,Inf), rational-in-Sin/Cos over a period,
                      principal values, even half-lines); tried before NewtonLeibniz under Automatic
"DiffUnderInt"         — parameter-dependent definite integrals by differentiation under the
```

("DifferentiationUnderIntegral") integral sign (Feynman's trick): Integrate\`DiffUnderInt; Laplace/Fourier, sinc, and even-rational half-line families; tried after Residue and NewtonLeibniz in the definite cascade "RamanujanMasterTheorem" — half-line Int\_0^Inf x^(s-1) f(x) dx by the Mellin transform /

```text
("Mellin")              Ramanujan Master Theorem: Integrate`RamanujanMasterTheorem;
                      exp/Gaussian/algebraic/Cos/Sin/ArcTan/Log/BesselJ/pFq/PolyLog
                      kernels (monomial x^k substitution; Erf, incomplete Gamma, BesselJ^2
                      reduced to pFq); also the exp-geometric kernel 1/(E^(c x)+g)
                      (Bose-Einstein / Fermi-Dirac -> Gamma*PolyLog), a Frullani pre-pass
                      (f(a x)-f(b x))/x -> (f(0)-f(Inf)) Log[b/a], and a Log[x]^k weight;
                      strip-gated, yielding a ConditionalExpression when
                      Assumptions do not prove convergence; after NewtonLeibniz under Automatic
```

Method -\> {"DerivativeDivides", "Substitution" -\> u} pins the kernel u(x), trialing only that substitution. Named methods are strict: failure returns unevaluated, with no fallback. of t seconds (default 3; Infinity for no budget).  It is checked BETWEEN stages, at every recursion depth, and inside the substitution stage's kernel and branch loops and its Eliminate/Solve calls -- so it stops the cascade progressing once spent, but it is not yet a hard bound: a single long-running evaluation inside one stage can overrun it.  It matters only for integrands Integrate cannot close; measured from inside the cascade (MATHILDA\_INTEGRATE\_PROFILE=1), two stages account for 90.8% of that cost -- the linearity split (43.5%) and the substitution search (47.3%) -- while the eight cheap exact stages ahead of them are together under 0.1%.  An explicit Method -\> name is a deliberate request for the full search and carries no budget unless you pass TimeConstraint yourself. The CRCTable rules are loaded from disk on first use only. An applied 1-D InterpolatingFunction integrates to its antiderivative InterpolatingFunction (mirroring D).

</details>

## Examples (73)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (7)

```mathematica
In[1]:= Integrate[3 + 5 x + 2 x^2, x]
Out[1]= 3 x + 5/2 x^2 + 2/3 x^3

In[2]:= Integrate[2 x/(x^2 + 1), x]
Out[2]= Log[1 + x^2]

In[3]:= Integrate[1/(x - a)^2, x]
Out[3]= 1/(a - x)

In[4]:= Integrate[(2x+3)/(x^2+3x+5)^2, x]
Out[4]= -1/(5 + 3 x + x^2)
```

Phase 2 LRT closes this

```mathematica
In[5]:= Integrate[1/((x-1)(x-2)(x-3)), x]
Out[5]= -Log[-2 + x] + 1/2 Log[3 - 4 x + x^2]
```

Phase 4 LogToReal

```mathematica
In[6]:= Integrate[1/(x^2 + 1), x]
Out[6]= ArcTan[x]
```

Two quadratic factors

```mathematica
In[7]:= Integrate[1/(x^4 + x^2 + 1), x]
Out[7]= 1/4 Log[1 + x + x^2] + ArcTan[(-1 + 2 x)/Sqrt[3]]/(2 Sqrt[3]) + ArcTan[(1 + 2 x)/Sqrt[3]]/(2 Sqrt[3]) - 1/4 Log[1 - x + x^2]
```

### Scope (38)

```mathematica
In[8]:= Integrate[x^2, {x, 0, 1}]
Out[8]= 1/3

In[9]:= Integrate[1/(1 + x^2), {x, 0, Infinity}]
Out[9]= 1/2 Pi
```

Improper, convergent

```mathematica
In[10]:= Integrate[1/Sqrt[x], {x, 0, 1}]
Out[10]= 2
```

Continuous branch form

```mathematica
In[11]:= Integrate[1/(2 + Cos[x]), {x, 0, 2 Pi}]
Out[11]= (2 Pi)/Sqrt[3]
```

Iterated

```mathematica
In[12]:= Integrate[x y, {x, 0, 1}, {y, 0, 1}]
Out[12]= 1/4
```

```mathematica
In[13]:= Integrate`SingularPoints[1/((x - 1)(x - 2)), {x, 0, 3}]
Out[13]= {1, 2}

In[14]:= Integrate[1/x, {x, 1 - I, 2 + 3 I}]
Out[14]= Log[2 + 3*I] - Log[1 - I]

In[15]:= Integrate[z^2, {z, 0, 1 + I}]
Out[15]= -2/3 + 2/3*I
```

Branch-point endpoint

```mathematica
In[16]:= Integrate[1/Sqrt[z], {z, 0, 1 + I}]
Out[16]= 2 Sqrt[1 + I]
```

```mathematica
In[17]:= Integrate`PathSingularPoints[1/z, {z, -1 - I, 1 + I}]
Out[17]= {0}
```

CCW loop about 0

```mathematica
In[18]:= Chop[N[Integrate[1/z, {z, 1, I, -1, -I, 1}]]]
Out[18]= 0.0 + 6.28319*I
```

```mathematica
In[19]:= Integrate[1/(1 + x^4), {x, -Infinity, Infinity}]
Out[19]= Pi/Sqrt[2]
```

Order-2 pole

```mathematica
In[20]:= Integrate[1/(1 + x^2)^2, {x, -Infinity, Infinity}]
Out[20]= 1/2 Pi
```

```mathematica
In[21]:= Integrate[Cos[x]/(1 + x^2), {x, -Infinity, Infinity}]
Out[21]= Pi/E

In[22]:= Integrate[1/(2 + Cos[x]), {x, 0, 2 Pi}]
Out[22]= (2 Pi)/Sqrt[3]
```

Principal value

```mathematica
In[23]:= Integrate[Sin[x]/x, {x, -Infinity, Infinity}]
Out[23]= Pi
```

Even half-line

```mathematica
In[24]:= Integrate[1/(1 + x^4), {x, 0, Infinity}]
Out[24]= Pi/(2 Sqrt[2])
```

```mathematica
In[25]:= Integrate[x f'[x] + f[x], x]
Out[25]= x f[x]

In[26]:= Integrate[f'[x] g[x] + f[x] g'[x], x]
Out[26]= f[x] g[x]

In[27]:= Integrate[f'[x]/f[x], x]
Out[27]= Log[f[x]]

In[28]:= Integrate[(f'[x] g'[x] - f[x] g''[x])/(f[x]^2 + g'[x]^2), x]
Out[28]= -ArcTan[Derivative[1][g][x]/f[x]]
```

Composite argument

```mathematica
In[29]:= Integrate[2 x f'[x^2], x]
Out[29]= f[x^2]
```

```mathematica
In[30]:= Integrate[(f[x] - x f[x] + f[x] Log[x f[x]] + x f'[x])/f[x], x]
Out[30]= -1/2 x^2 + x Log[x f[x]]

In[31]:= Integrate[Log[f[x]] f'[x]/f[x], x]
Out[31]= 1/2 Log[f[x]]^2
```

Direct, correct branch

```mathematica
In[32]:= Integrate[Sin[x] Sqrt[1 - Cos[x]], x]
Out[32]= 2/3 (1 - Cos[x])^(3/2)
```

```mathematica
In[33]:= Integrate[1/(x Log[x]), x]
Out[33]= Log[Log[x]]

In[34]:= Integrate`DerivativeDivides[2 x Exp[x^2], x]
Out[34]= E^x^2
```

U = Sqrt[Tan[x]]

```mathematica
In[35]:= Integrate[Sqrt[Tan[x]], x] - Log[1 + Tan[x] + Sqrt[2] Sqrt[Tan[x]]]/(2 Sqrt[2]) + Log[1 + Tan[x] - Sqrt[2] Sqrt[Tan[x]]]/(2 Sqrt[2])
Out[35]= ArcTan[-1 + Sqrt[2] Sqrt[Tan[x]]]/Sqrt[2] + ArcTan[1 + Sqrt[2] Sqrt[Tan[x]]]/Sqrt[2] + Log[1 + Tan[x] - Sqrt[2] Sqrt[Tan[x]]]/Sqrt[2] - Log[1 + Tan[x] + Sqrt[2] Sqrt[Tan[x]]]/Sqrt[2]
```

U = Sqrt[Cot[x]]

```mathematica
In[36]:= Integrate[Sqrt[Cot[x]], x] - Log[1 + Cot[x] - Sqrt[2] Sqrt[Cot[x]]]/(2 Sqrt[2]) + Log[1 + Cot[x] + Sqrt[2] Sqrt[Cot[x]]]/(2 Sqrt[2])
Out[36]= -ArcTan[-1 + Sqrt[2] Sqrt[Cot[x]]]/Sqrt[2] - ArcTan[1 + Sqrt[2] Sqrt[Cot[x]]]/Sqrt[2] + Log[1 + Cot[x] + Sqrt[2] Sqrt[Cot[x]]]/Sqrt[2] - Log[1 + Cot[x] - Sqrt[2] Sqrt[Cot[x]]]/Sqrt[2]
```

```mathematica
In[37]:= Integrate[1/Sqrt[x + 1], x]
Out[37]= 2 Sqrt[1 + x]

In[38]:= Integrate[Sqrt[x]/(1 + Sqrt[x]), x]
Out[38]= -2 Sqrt[x] + x + 2 Log[1 + Sqrt[x]]

In[39]:= Integrate[1/Sqrt[x^2 + 1], x]
Out[39]= ArcSinh[x]

In[40]:= Integrate[1/Sqrt[1 - x^2], x]
Out[40]= ArcSin[x]
```

Symbolic leading coeff

```mathematica
In[41]:= Integrate[1/Sqrt[a x^2 + 1], x]
Out[41]= -Log[-Sqrt[a] x + Sqrt[1 + a x^2]]/Sqrt[a]
```

```mathematica
In[42]:= Integrate[1/Sqrt[(x + 1)/(x - 1)], x]
Out[42]= (2 Sqrt[(1 + x)/(-1 + x)])/(-1 + (1 + x)/(-1 + x)) - 2 ArcTanh[Sqrt[(1 + x)/(-1 + x)]]
```

Paper eq. (10)

```mathematica
In[43]:= Integrate[3/(5 - 4 Cos[x]), x]
Out[43]= 2 Pi Floor[(-Pi + x)/(2 Pi)] + 2 ArcTan[3 Tan[1/2 x]]
```

Continuous form

```mathematica
In[44]:= Integrate[1/(2 + Cos[x]), x]
Out[44]= (2 Pi Floor[(-Pi + x)/(2 Pi)])/Sqrt[3] + (2 ArcTan[Tan[1/2 x]/Sqrt[3]])/Sqrt[3]
```

Hyperbolic: no Floor

```mathematica
In[45]:= Integrate[1/(2 + Cosh[x]), x]
Out[45]= (2 ArcTanh[Tanh[1/2 x]/Sqrt[3]])/Sqrt[3]
```

### Options (7)

Strict, no fallback

```mathematica
In[46]:= Integrate[Sin[x], x, Method -> "RischTranscendental"]
Out[46]= -Cos[x]
```

```mathematica
In[47]:= Integrate[x^3, x, Method -> "BronsteinRational"]
Out[47]= 1/4 x^4

In[48]:= Integrate[Sin[x]/Cos[x]^2, x, Method -> "DerivativeDivides"]
Out[48]= Sec[x]

In[49]:= Integrate[1/(1 + x^(1/3)), x, Method -> "LinearRadicals"]
Out[49]= -3 x^(1/3) + 3/2 x^(2/3) + 3 Log[1 + x^(1/3)]

In[50]:= Integrate[Sqrt[2 x + 3]/x, x, Method -> "LinearRadicals"]
Out[50]= 2 Sqrt[3 + 2 x] - 2 Sqrt[3] ArcCoth[Sqrt[3 + 2 x]/Sqrt[3]]

In[51]:= Integrate[1/Sqrt[x^2 - 1], x, Method -> "QuadraticRadicals"]
Out[51]= Log[2 x + 2 Sqrt[-1 + x^2]]

In[52]:= Integrate[1/Sqrt[(2 x + 1)/(x + 3)], x, Method -> "LinearRatioRadicals"]
Out[52]= -(5 Sqrt[(1 + 2 x)/(3 + x)])/(-4 + (2 (1 + 2 x))/(3 + x)) + (5 ArcTanh[Sqrt[(1 + 2 x)/(3 + x)]/Sqrt[2]])/(2 Sqrt[2])
```

### Applications (21)

Rational: Hermite reduction + Lazard–Rioboo–Trager, real Log/ArcTan form

```mathematica
In[53]:= Integrate[1/(x^3 + 1), x]
Out[53]= 1/3 Log[1 + x] + ArcTan[(-1 + 2 x)/Sqrt[3]]/Sqrt[3] - 1/6 Log[1 - x + x^2]
```

The PolyLogLadder stage: R(x) Log[linear]^k to PolyLog of rising order

```mathematica
In[54]:= Integrate[Log[x]^2/(1 + x), x]
Out[54]= Log[x]^2 Log[1 + x] + 2 (-PolyLog[3, -x] + Log[x] PolyLog[2, -x])
```

Polynomial times trig rational: TrigToExp, then the exponential-tower polylog ladder

```mathematica
In[55]:= Integrate[x/Sin[x], x]
Out[55]= I (-I x Log[1 - E^(-I x)] + PolyLog[2, E^(-I x)]) + I (I x Log[1 + E^(-I x)] - PolyLog[2, -E^(-I x)])
```

Symbolic exponent: the GammaPower recogniser gives the incomplete Gamma at once

```mathematica
In[56]:= Integrate[x^n Exp[-x], x]
Out[56]= -Gamma[1 + n, x]
```

One integration by parts closes to CosIntegral

```mathematica
In[57]:= Integrate[Log[x] Sin[x], x]
Out[57]= CosIntegral[x] - Cos[x] Log[x]
```

Non-elementary: the last stage, ParallelMixedSpecial, answers with an elliptic integral

```mathematica
In[58]:= Integrate[(1 + x^2)/Sqrt[1 - x^4], x]
Out[58]= EllipticE[ArcSin[x], -1]
```

Residue theorem, Fourier (Jordan) family: upper half-plane pole at I

```mathematica
In[59]:= Integrate[Cos[x]/(1 + x^2), {x, -Infinity, Infinity}]
Out[59]= Pi/E
```

Residue theorem on the unit circle for a periodic integrand

```mathematica
In[60]:= Integrate[1/(2 + Cos[x]), {x, 0, 2 Pi}]
Out[60]= (2 Pi)/Sqrt[3]
```

Symbolic exponent: convergence for 0 < s < 1 is proved from the Assumptions

```mathematica
In[61]:= Integrate[x^(s - 1)/(1 + x), {x, 0, Infinity}, Assumptions -> 0 < s < 1]
Out[61]= Pi Csc[Pi s]
```

The Mellin table gives the same value as Beta[s, 1 - s] = Pi Csc[Pi s]

```mathematica
In[62]:= Integrate[x^(s - 1)/(1 + x), {x, 0, Infinity}, Assumptions -> 0 < s < 1, Method -> "Mellin"]
Out[62]= Beta[s, 1 - s]
```

Laplace transform of a Bessel function

```mathematica
In[63]:= Integrate[Exp[-a x] BesselJ[0, b x], {x, 0, Infinity}, Assumptions -> a > 0 && b > 0]
Out[63]= 1/Sqrt[a^2 + b^2]
```

Recognised as the integral representation of BesselK

```mathematica
In[64]:= Integrate[Exp[-z Cosh[t]], {t, 0, Infinity}, Assumptions -> z > 0]
Out[64]= BesselK[0, z]
```

Expand the geometric series, integrate term by term, resum to Catalan's constant

```mathematica
In[65]:= Integrate[Log[x]/(1 + x^2), {x, 0, 1}, Method -> "SeriesExpansion"]
Out[65]= -Catalan
```

Frullani-type logarithm, by differentiating under the integral sign in the exponent

```mathematica
In[66]:= Integrate[(x^a - x^b)/Log[x], {x, 0, 1}, Assumptions -> a > 0 && b > 0]
Out[66]= Log[1 + a] - Log[1 + b]
```

King reflection: f(x) + f(Pi/2 - x) = 1, with no elementary antiderivative needed

```mathematica
In[67]:= Integrate[1/(1 + Tan[x]^Sqrt[2]), {x, 0, Pi/2}]
Out[67]= 1/4 Pi
```

Reflection x -> Pi - x removes the factor x

```mathematica
In[68]:= Integrate[x Sin[x]/(1 + Cos[x]^2), {x, 0, Pi}]
Out[68]= 1/4 Pi^2
```

Piecewise: split at the breakpoint x = 1 and integrate each branch

```mathematica
In[69]:= Integrate[Max[x, x^2], {x, 0, 2}]
Out[69]= 17/6
```

Cauchy principal value across the pole at x = 1

```mathematica
In[70]:= Integrate[1/(x^3 - 1), {x, 0, 2}, PrincipalValue -> True]
Out[70]= 1/18 (-3 Log[7] + Pi Sqrt[3] - 6 Sqrt[3] ArcTan[5/Sqrt[3]])
```

Differentiation under the integral sign in a, then integration back

```mathematica
In[71]:= Integrate[Log[1 + a^2 x^2]/(1 + x^2), {x, 0, Infinity}, Assumptions -> a > 0]
Out[71]= Pi Log[1 + a]
```

Two kernels, each a Meijer G of x^2, combined by the convolution theorem

```mathematica
In[72]:= Integrate[BesselJ[0, x] E^(-x^2), {x, 0, Infinity}, Method -> "MeijerConvolution"]
Out[72]= (BesselI[0, 1/8] Sqrt[Pi])/(2 E^(1/8))
```

A known gap: the value is Pi/4 but no stage closes it, so the integral is returned unevaluated

```mathematica
In[73]:= Integrate[Exp[-x] Sin[x]/x, {x, 0, Infinity}]
Out[73]= 1/2 MeijerG[{{1/2}, {1}}, {{0, 1/2}, {}}, 1]
```

## Options & behaviour

### Differentiation under the integral sign (`Integrate\`DiffUnderInt`)

For a parameter-dependent definite integral `I(p) = Integrate[f(x,p), {x,a,b}]`,
this method (Leibniz rule / "Feynman's trick") differentiates the integrand with
respect to a free parameter `p`, evaluates the resulting simpler definite
integral `J(p) = Integrate[D[f,p], {x,a,b}]`, integrates `J(p)` back over the
parameter, and fixes the constant of integration from an **exact** base value
`I(p0)` (a `p` where `f` vanishes identically, or reduces to a directly-
integrable form). The general case is the first-order linear ODE
`I'(p) = lambda(p) I(p) + M(p)` (Boulnois 2023). Verification is symbolic and
correct-by-construction (`Simplify[D[I,p] - (lambda I + M)] === 0` plus the exact
base) — there is **no** numeric crosscheck. Assumptions (`a > 0`, …) are honoured
and used to clean the closed forms.

Two stages solve that ODE. **Stage A (`lambda = 0`, pure quadrature)**: when the
inner integral `J` closes to an x-free form, `I = Integral J dp + const`. **Stage
B (`lambda != 0`, integrating factor)**: when `J` does *not* close but equals
`lambda(p) I + M(p)` — i.e. `D[f,p] - lambda(p) f = D[U,x]` by integration by
parts in `x`, for an x-free `lambda` and a `U` in the span of `f`'s transcendental
atoms — the integral is recovered as `I = (1/mu)(Integral mu M dp + C)` with
`mu = Exp[-Integral lambda dp]`. The ansatz `U = Sum c_k B_k` (the `B_k` the
atom-product shapes of `f` and `D[f,p]`) and `lambda` are found by matching
coefficients of the independent atoms via `Solve`, and the IBP identity is then
verified literally. This is what closes the classic Feynman Gaussian
`Integrate[Exp[-a^2 x^2] Cos[b x], {x,0,Infinity}]`: differentiating in `b` gives
`J = -(b/(2a^2)) I`, so `I = (Sqrt[Pi]/(2a)) Exp[-b^2/(4a^2)]`. All parameter
antiderivatives route through the bounded back-integration (never the unbounded
engine on a radical/trig/Gaussian-of-parameter form).

A **self-similar** recognizer handles the Laurent-Gaussian
`Exp[c0 + c2 x^2 + cm2/x^2]` (`c2, cm2 < 0`) that integration by parts in `x`
cannot reach: the reciprocal substitution `x -> Sqrt[cm2/c2]/x` leaves the
exponent invariant, so `J = lambda I` with `lambda` read off directly. This closes
`Integrate[Exp[-a^2 x^2 - b^2/x^2], {x,0,Infinity}]`: differentiating in `b` gives
`lambda = -2a`, so `I = (Sqrt[Pi]/(2a)) Exp[-2 a b]`.

When neither the inner integral nor Stage B closes for a parameter, selection
prefers a parameter with an **exact zero base** (`I(p0) = 0`) over one needing a
sign-sensitive computed base, since the `D[I,p]` check cannot catch a wrong
constant base.

Two further closers run after the quadrature stages. **Repeated differentiation**:
when a single differentiation of `trig^n/x^m` leaves a `1/x^{m-1}` the families
cannot close, the integrand is differentiated `k` times until the inner integral
closes, then back-integrated `k` times with the `k` lower derivatives pinned at an
exact base (`Integrate[Sin[a x]^3/x^3, {x,0,Inf}] = 3 Pi a^2/8`: `k = 2`, `I''(a) =
3 Pi/4`, `I(0) = I'(0) = 0`). **Reverse recognition**: when the integrand is itself
a parameter-derivative of a known integral -- `f = Log[x]^k h` with `d/da h =
Log[x] h` -- the value is `d^k/da^k Integrate[h, {x,a,b}]`
(`Integrate[Exp[-x] x^(a-1) Log[x], {x,0,Inf}] = d/da Gamma[a] = Gamma[a]
PolyGamma[0,a]`); the base integral `Integrate[h]` is taken from the families or,
for a decaying Laplace/Gamma shape `C Exp[alpha x] x^(s-1)`, from the engine (gated
to that shape so it cannot hang).

The output is finalized by two equality-gated transforms: 1-arg `PowerExpand`
(readability) and an assumption-aware `FullSimplify` (taken only when strictly
simpler), which contracts a split Log the per-family Simplify leaves standing
(`(Pi(-2Log[b]+Log[a b+b^2]))/b -> (Pi(Log[a+b]-Log[b]))/b`).

Because the general integrator is slow/hangs on the parameter-dependent inner
integrals Feynman's trick produces, `DiffUnderInt` evaluates the standard
families itself with closed-form formulas: the **Laplace/Fourier half-line**
`∫₀^∞ xⁿ e^{-p x}{1,cos,sin} dx`, the **sinc/Frullani** `∫₀^∞ …/x dx`, the
**even-rational half-line** `∫₀^∞ P(x)/Q(x²) dx`, the **general (non-even)
rational half-line** `∫₀^∞ R(s) ds` (real `ArcTan`/`Log` boundary values — this
is what closes a *decaying* sinc such as `∫₀^∞ e^{-p x} Sin[q x]/x dx = ArcTan[q/p]`
whose Laplace image is non-even), and the **Gaussian moment** family
`∫₀^∞ xⁿ e^{-p x²}{1,cos} dx` in `Sqrt[Pi]`/`e^{-q²/4p}`. The Gaussian
parameter back-integration `∫ c e^{-k p²} dp` is supplied directly as an `Erf`
(the engine does not produce it).

The **parameter back-integration** `∫ J dp` is bounded by construction: `J` is
first `PowerExpand`-ed (so a `Log[1/poly^rational]` image is pulled apart into
plain `Log`s — the reciprocal-radical spelling otherwise drives the indefinite
engine into an uninterruptible algebraic grind), then declined unless it is free
of radical/trig/Gaussian *of the parameter*, leaving only the elementary
`Log`/`ArcTan`/rational forms the families emit. A **product of trig-of-x**
(`Sin[b x] Sin[c x]`, …) is `TrigReduce`-normalised to a sum of single trig terms
before the families run, and parameter selection **prefers an exact zero base**
(an integrand that vanishes at the base point, `I(p₀)=0`) over a sign-sensitive
computed base — together these make `∫₀^∞ e^{-a x} Sin[b x] Sin[c x]/x dx` close
to the clean, all-signs-real `(1/4) Log[(a²+(b+c)²)/(a²+(b−c)²)]` rather than a
`b>c`-only form that goes non-real for `b<c`.

Forms outside these families (finite-period
trig, piecewise/`Min`-`Max` results, the Sin-Gaussian Dawson/Erfi moment) are
declined — the integral is returned unevaluated, fast, never a wrong value. The
two rational half-line families gate on the inner integrand being a **rational
function of `x`**: a differentiated exp-geometric/Mellin form (still carrying
`e^x` or `x^{s-1}`) is not, and feeding it to `Apart[·, x]` otherwise drives a
non-terminating rewrite — so such forms are declined up front and left for the
Ramanujan/Mellin method.

  (Stage B)
  (Stage B, self-similar)
  (reverse recognition)
  (repeated differentiation).

Still **deferred**: finite-period log-trig integrals (`Integrate[Log[a^2 - 2a Cos[x]
+ 1], {x,0,π}]`, `Integrate[Log[1 - a^2 Sin[x]^2], {x,0,π/2}]`). Their differentiated
inner integral (a rational function of Sin/Cos over a period) is reachable via the
residue unit-circle engine, but the Feynman back-integration produces a
complex-branch antiderivative singular at the natural base (`∝ 1/a` or
`1/√(1-a^2)`), and `Limit` cannot resolve the canceling log-singularity there; a
branch-correct back-integration is needed first.

Three **finite-domain** families need neither a pre-existing parameter (the first
two are purely numeric) nor an engine-safe inner integral (all three differentiate
to a finite-interval trig/radical integral the general engine cannot do — it hangs,
or for the radical returns a *wrong* value). Their closers canonicalise with a
change of variables, **introduce** an artificial Feynman parameter (or evaluate
directly), and supply the inner integral in closed form:

| Family | Integrand → value (change of variables) |
|---|---|
| power-log | `Log[1 + c x^p]/(x Sqrt[1 - x^(2p)])` on `{0,1}` → `(π²/8 − ArcCos[c]²/2)/p`  (`u = x^p`) |
| secant-radical | `Sec[2x] Log[1 + c Sqrt[1 - Tan[x]²]]` on `{0,π/4}` → `π²/8 − ArcCos[c]²/2`  (`t = Tan[x]`) |
| tangent-power | `Csc[2x]² Log[1 + Tan[x]^a]` on `{0,π/4}` → `(π Csc[π/a] − a)/4`  (`t = Tan[x]`) |

The `p` (power), `c` (coefficient), and `a` (exponent) are read from the integrand,
so the closed form is genuinely parametric. Families 1 and 2 share the inner
integral `∫₀¹ du/((1+q u)√(1−u²)) = ArcCos[q]/√(1−q²)` (emitted directly — the
engine leaves it unevaluated at a symbolic parameter, and mis-evaluates the
secant-radical variant). Family 3 is a direct Beta/digamma reflection, verified
correct-by-construction against the value-independent `Csc[2x]² = (1+Tan[x]²)²/(4
Tan[x]²)` identity and the exact rational anchor `a₀ = 3`. Examples:

.
(Family 3 requires the explicit `Method -> "DiffUnderInt"`; under the fully
automatic cascade an earlier method attempts the — separately hang-prone —
indefinite integral first.)

### Mellin / Ramanujan Master Theorem (`Integrate\`RamanujanMasterTheorem`)

The series/transform-based mechanism for half-line integrals
`∫₀^∞ x^{s-1} f(x) dx` of a *transcendental* `f` (the class residue and FTC do
not close). By Ramanujan's Master Theorem, if `f(x) = Σ (-1)^k φ(k) x^k / k!`
then the Mellin transform is `Γ(s) φ(-s)` on the fundamental strip. The
integrand is Expanded and, term by term, decomposed into `C · x^ρ · f(λx)`; the
kernel `f` is matched against a table of proven base Mellin transforms and the
power prefactor sets `s = ρ + 1`:

| kernel `f` | `∫₀^∞ x^{s-1} f dx` | strip |
|------------|--------------------|-------|
| `Exp[c x]`, `Re c<0` | `Γ(s) (-c)^{-s}` | `0<Re s` |
| `Exp[c x^2]`, `Re c<0` | `½ (-c)^{-s/2} Γ(s/2)` | `0<Re s` |
| `(p + q x^m)^{-a}`, `p,q>0` | `(1/m) p^{s/m-a} q^{-s/m} B(s/m, a-s/m)` | `0<Re s<m Re a` |
| `Cos[λ x]`, `λ>0` | `π / (2 Sin(πs/2) Γ(1-s)) λ^{-s}` | `0<Re s<1` |
| `Sin[λ x]`, `λ>0` | `π / (2 Cos(πs/2) Γ(1-s)) λ^{-s}` | `-1<Re s<1` |
| `Log[1 + λ x]`, `λ>0` | `π / (s Sin(π s)) λ^{-s}` | `-1<Re s<0` |
| `ArcTan[λ x]`, `λ>0` | `-π / (2 s Cos(πs/2)) λ^{-s}` | `-1<Re s<0` |
| `BesselJ[ν, λ x]`, `λ>0` | `2^{s-1} λ^{-s} Γ((ν+s)/2)/Γ((ν-s)/2+1)` | `-Re ν<Re s<3/2` |
| `pFq[{a}; {b}; -λ x]`, `λ>0` | `(∏Γ(b_j)/∏Γ(a_i)) Γ(s) (∏Γ(a_i-s)/∏Γ(b_j-s)) λ^{-s}` | `0<Re s<min Re a_i` |
| `PolyLog[ν, -λ x]`, `λ>0` | `π (-s)^{-ν} λ^{-s} / Sin(π s)` | `-1<Re s<0` |
| `1/(e^{c x}+γ)`, `c>0`, `γ≥-1` | `c^{-s} Γ(s) (-1/γ) PolyLog(s, -γ)` | `0<Re s` (`1<Re s` if `γ=-1`) |
| `BesselK[ν, λ x]`, `λ>0` | `2^{s-2} λ^{-s} Γ((s+ν)/2) Γ((s-ν)/2)` | `Re s>\|Re ν\|` |
| `AiryAi[λ x]`, `λ>0` | `λ^{-s} Γ(s) / (3^{(s+2)/3} Γ((s+2)/3))` | `0<Re s` |
| `Erfc[λ x]`, `λ>0` | `λ^{-s} Γ((s+1)/2) / (√π s)` | `0<Re s` |
| `SinIntegral[λ x]`, `λ>0` | `-λ^{-s} Γ(s) Sin(πs/2) / s` | `-1<Re s<0` |
| `CosIntegral[λ x]`, `λ>0` | `-λ^{-s} Γ(s) Cos(πs/2) / s` | `0<Re s<1` |
| `ExpIntegralE[n, λ x]` (general `n`), `λ>0` | `λ^{-s} Γ(s) / (s+n-1)` | `Re s>0, Re(s+n)>1` |
| `StruveH[0, λ x]`, `λ>0` | `π 2^{s-1} λ^{-s} / (Cos(πs/2) Γ(1-s/2)²)` | `-1<Re s<1` |
| `Gamma[ν, λ x]` (upper incomplete), `λ>0` | `λ^{-s} Γ(s+ν) / s` | `Re s>0, Re(s+ν)>0` |
| `Csch[λ x]`, `λ>0` | `2 λ^{-s} Γ(s) (1-2^{-s}) ζ(s)` | `1<Re s` |
| `Sech[λ x]`, `λ>0` | `2^{1-s} λ^{-s} Γ(s) LerchPhi(-1, s, 1/2)` | `0<Re s` |
| `Coth[λ x] - 1`, `λ>0` | `2^{1-s} λ^{-s} Γ(s) ζ(s)` | `1<Re s` |
| `Sin[λ x]^k` (`k≥2` int), `λ>0` | linearize (`TrigReduce`) → Σ harmonics via `Sin/Cos` rows | `-k<Re s<` (`0` even `k`, `1` odd) |
| `Cos[λ x]^k` (`k≥2` int), `λ>0` | linearize → Σ harmonics (mean dropped) | `0<Re s<` (`0` even → diverges, `1` odd) |
| `ArcTan[λ x]^2`, `λ>0` | `(π λ^{-s}/(2s)) Csc(πs/2) (ψ((1-s)/2) − ψ(1/2))` | `-2<Re s<0` |

The kernels below the exponential-geometric row are closed through three general
layers rather than one-off entries. `BesselK` and `AiryAi` are **dedicated base
transforms** (each a cancelling combination of individually-divergent series, so
no single `pFq` exists). `SinIntegral` and `StruveH` are **reduced to a single
entire `₁F₂`** (`Si(z)=z·₁F₂(½;3/2,3/2;−z²/4)`,
`H_ν(z)=(z/2)^{ν+1}(2/(√π Γ(ν+3/2)))·₁F₂(1;3/2,ν+3/2;−z²/4)`) and closed by the
`pFq` recognizer. `Erfc` and `CosIntegral` are closed by **Mellin operational
calculus**: `∫₀^∞ x^{s-1} f = −(1/s) ∫₀^∞ x^s f'` by parts whenever `f'` is itself
a recognized kernel (`Erfc'`=Gaussian, `Ci'`=cos/x), with the lower boundary
`Re s>0` added from `f`'s constant-or-logarithmic value at the origin.
`ExpIntegralE[n, λx]` is a **dedicated recognizer for general (symbolic) `n`**
giving `λ^{-s} Γ(s)/(s+n-1)` — the operational-calculus `E_n'=−E_{n-1}` chain
cannot bottom out for symbolic `n`. Its strip is `Re s>0` (the `E_n(0)` constant)
**and** `Re(s+n)>1` (the `z^{n-1}` branch term, binding for small `n`); the latter
is *not* implied by `{Re s>0, Re(s−n)<0}`, so under those assumptions the honest
result is a `ConditionalExpression` carrying `Re(s+n)>1`. The **upper incomplete gamma**
`Γ[ν, λx]` is a dedicated recognizer giving `λ^{-s} Γ(s+ν)/s` (the operational-
calculus identity `Γ'(ν,x)=−x^{ν-1}e^{-x}`, with `x^{ν-1}` folded into the power).
`Csch` and `Sech` are the two-exponential **hyperbolic Dirichlet** siblings of the
exponential-geometric kernel: expanding `e^{-λx}` geometrically lands on an
odd-argument Dirichlet series, `Csch[λx] → 2λ^{-s}Γ(s)(1-2^{-s})ζ(s)` (`Re s>1`)
and `Sech[λx] → 2^{1-s}λ^{-s}Γ(s)·LerchPhi(-1,s,1/2)` (`= 2λ^{-s}Γ(s)β(s)`,
`Re s>0`). The **Bose kernel** `Coth[λx] − 1 = 2/(e^{2λx}−1)` is reduced (before
Expand, so the divergent-on-its-own `Coth` and `−1` are never split) to the
exponential-geometric kernel and closed by it as `2^{1-s}λ^{-s}Γ(s)ζ(s)` (`Re s>1`)
bare `Coth[λx]` (no strip — it diverges at both ends) is left unreduced and declines.

The convergence strip is matched against the user `Assumptions` with `Simplify`
and then `Refine`, which reasons about real parts where `Simplify` does not — so a
`Re[s]>0` (or a sign `a>0`) assumption that pins the whole strip collapses the
`ConditionalExpression` to the bare transform. A strip only partly proven is
carried intact, so the reported convergence region is always complete. An
integrand carrying an *arbitrary function of the integration variable* (`f[x]`
with no definition) still skips this method; an x-free symbolic parameter does
not.

The last row is the **exponential-geometric** kernel of the statistical-mechanics
integrals: expanding `1/(e^{cx}+γ) = (-1/γ) Σ_{j≥1} (-γ)^j e^{-jcx}` and
integrating term by term lands on `PolyLog`. Its two headline specialisations are
`γ=-1` **Bose–Einstein** `∫₀^∞ x^{s-1}/(e^{cx}-1) = c^{-s} Γ(s) ζ(s)` (Planck /
Debye; the denominator zero at `x=0` tightens the strip to `Re s>1`) and `γ=+1`
**Fermi–Dirac** `∫₀^∞ x^{s-1}/(e^{cx}+1) = c^{-s} Γ(s) η(s)` (emitted as
`-Γ(s) PolyLog(s,-1)`, which stays finite at `s=1` where `(1-2^{1-s})ζ(s)` would
be `0·∞`). The admitted range is the exact **convergence** region `γ≥-1`, not the
series' `|γ|≤1`: for `γ>1` — the high-fugacity / positive-chemical-potential
**degenerate Fermi gas** (electrons in metals, white-dwarf matter) — the
denominator `e^{cx}+γ` has no interior zero, so the integral converges for
`Re s>0` and equals the same `-c^{-s} Γ(s) PolyLog(s,-γ)/γ` by the Fermi–Dirac
analytic continuation `F_{s-1}(η)=-Li_s(-e^{η})` (Dingle/Blakemore), even though
the geometric series diverges there; e.g. `∫₀^∞ 1/(e^x+2) dx = ½ ln 3`,
`∫₀^∞ x/(e^x+2) dx = -½ PolyLog(2,-2)`. Only the divergent Bose region `γ<-1`
(interior pole at `x=ln(-γ)>0`) stays declined. A **symbolic fugacity** is
admitted too — the general Bose integral `∫₀^∞ x^{s-1}/(z^{-1} e^x - 1) dx =
Γ(s) PolyLog(s, z)` closes for `0<z<1` and the degenerate Fermi
`∫₀^∞ x^{s-1}/(e^x + z) dx = -Γ(s) PolyLog(s, -z)/z` closes for `z>0`. The
built-in assumption engine only discharges syntactic matches (it proves neither
`1/z>0` nor `γ'≥-1` from the raw bound), so the `γ'≥-1` gate is decided by a
small **sound interval-bound prover** over the parameter box read off the
`Assumptions`: interval arithmetic yields a guaranteed enclosure, so the gate
never accepts an inadmissible fugacity (an unbounded-below or out-of-range `z`
simply declines).

Four operational layers extend the table:

- **Monomial substitution** `g(x^k)` (`k≠1`) via `y = x^k`:
  `∫ x^{s-1} g(x^k) = (1/|k|) ∫ y^{s/k-1} g(y)`, so `Sin[√x]`, `BesselJ[ν,2√x]`,
  `ArcTan[√x]`, `Cos[x²]` reduce to the linear table at `s/k`. The Jacobian is
  `1/|k|`, **not** `1/k`: for `k<0` the substitution reverses the `(0,∞)` limits,
  contributing the sign that keeps the factor positive — so a reciprocal argument
  like `ArcTan[a/x]` (`k=−1`) gives the correct `+(π/2s) Sec(πs/2) a^s`.
- **Trig powers** `Sin[λx]^k` / `Cos[λx]^k` (`k≥2` integer): `TrigReduce`
  linearizes into a sum of harmonics `c_j {Cos,Sin}(b_j x)` plus a possible mean
  each harmonic goes through the `Sin/Cos` base rows at symbolic `s` and the
  scaleless mean is dropped (its Mellin transform is `0` on the `Re s<0` half
  where the even-power mean lives). The strip is read off the whole function's
  asymptotics — left edge `−k` for `Sin^k` (`∼(λx)^k` at 0), `0` for `Cos^k`
  right edge `0` when a nonzero mean survives (even power), `1` when the mean is 0
  (odd). An empty strip declines (`Cos^{even}`: `0<Re s<0`, genuinely divergent).
- **`ArcTan[λx]^2`**: the Mellin convolution of `ArcTan` with itself, whose Barnes
  residue sum closes to a digamma difference,
  `(π λ^{-s}/(2s)) Csc(πs/2) (ψ((1−s)/2) − ψ(1/2))` on `−2<Re s<0` (Mathematica
  returns this unevaluated).
- **Hypergeometric reduction** (applied before Expand, so a cancellation kernel
  is never split): `Erf[u] → u·₁F₁`, `Γ[a]-Γ[a,x] → x^a/a·₁F₁` (lower incomplete
  gamma), the product `BesselJ[ν,·]² → ₁F₂`, and the complete elliptic integrals
  `EllipticK[z] → (π/2)·₂F₁(½,½;1;z)`, `EllipticE[z] → (π/2)·₂F₁(-½,½;1;z)`. The
  Weber function reduces to a `cos/sin` combination of two `BesselJ`,
  `BesselY[ν,z] → (cos(νπ) J_ν(z) − J_{-ν}(z))/sin(νπ)`, each half closed by the
  `BesselJ` base transform (so `∫₀^∞ x^{s-1} Y_ν(ax)` closes on
  `|Re ν|<Re s<3/2`).
- **Parametric differentiation** for `Log[1+λx]^n (1+λx)^{-w₀}`:
  `M = (-1)^n ∂ⁿ_w[λ^{-s} B(s, w-s)]|_{w=w₀}`, strip `-n<Re s<w₀`.
- **`Log[x]^k` weight** (a bare `Log[x]`, distinct from the `Log[1+λx]` kernel):
  since `∂_s x^{s-1} = x^{s-1} Log x`, a `Log[x]^k` factor is the `k`-th
  `s`-derivative of the base transform `M_R(s)`. The open strip carries unchanged
  (`Log^k` is dominated by `x^{±ε}`), so e.g. `∫₀^∞ Log[x]/(1+x²) dx = 0` and
  `∫₀^∞ x Log[x] e^{-x} dx = Γ'(2) = 1-γ`.
- The `pFq` transform is the master kernel — `1F1`, `2F1`, `3F2`, … close
  uniformly (`Hypergeometric1F1`/`2F1` are stored as `HypergeometricPFQ`).
- **Mellin convolution** — a product of **two** transcendental kernels
  `M[K₁K₂](s) = (1/2πi)∫ M[K₁](z) M[K₂](s-z) dz` is a Barnes integral of Gamma
  ratios; it closes in two regimes. At **equal internal scale** the Barnes
  integral is a pure Gamma ratio: two `BesselJ` give Weber–Schafheitlin (through
  the `J_μ J_ν → ₂F₃` identity and `pFq`), two `BesselK` give Barnes' first
  lemma, and `J_ν K_ν` (equal order) gives Kummer's theorem —

  | product (scale `a`) | `∫₀^∞ x^{s-1} · dx` | strip |
  |---------------------|---------------------|-------|
  | `J_μ(ax) J_ν(ax)` | `2^{s-1} a^{-s} Γ(1-s) Γ(\tfrac{μ+ν+s}{2}) / (Γ(\tfrac{μ-ν-s}{2}+1)Γ(\tfrac{ν-μ-s}{2}+1)Γ(\tfrac{μ+ν-s}{2}+1))` | `Re(μ+ν+s)>0, Re s<1` |
  | `K_μ(ax) K_ν(ax)` | `2^{s-3} a^{-s} Γ(\tfrac{s+μ+ν}{2})Γ(\tfrac{s+μ-ν}{2})Γ(\tfrac{s-μ+ν}{2})Γ(\tfrac{s-μ-ν}{2}) / Γ(s)` | `Re s>\|Re μ\|+\|Re ν\|` |
  | `J_ν(ax) K_ν(ax)` | `2^{s-2} a^{-s} Γ(\tfrac{s}{2}) Γ(1+\tfrac{ν}{2}+\tfrac{s}{4}) / ((ν+\tfrac{s}{2})Γ(1+\tfrac{ν}{2}-\tfrac{s}{4}))` | `Re s>0, Re(s+2ν)>0` |

  At **distinct scales** the convolution carries the scale ratio and lands on a
  Gauss/Kummer hypergeometric: `e^{-ax} J_ν(bx) → (b/2)^ν Γ(s+ν)/(a^{s+ν}Γ(ν+1))·₂F₁(\tfrac{s+ν}{2},\tfrac{s+ν+1}{2};ν+1;-b²/a²)`
  (`Re a>0`), and the Gaussian case (Weber's second exponential integral)
  `e^{-a²x²} J_ν(bx) → b^ν Γ(\tfrac{s+ν}{2})/(2^{ν+1}(a²)^{(s+ν)/2}Γ(ν+1))·₁F₁(\tfrac{s+ν}{2};ν+1;-b²/(4a²))`,
  both on `Re(s+ν)>0`. The **exponential × trigonometric** family
  `e^{-ax} Sin(bx)` / `e^{-ax} Cos(bx)` (`a,b>0`) closes in elementary form via the
  complex frequency `c = a − i b`, `M[e^{-cx}] = c^{-s}Γ(s)`:
  `e^{-ax} Sin(bx) → Γ(s)(a²+b²)^{-s/2} Sin(s·ArcTan(b/a))` (`Re s>-1`) and
  `e^{-ax} Cos(bx) → Γ(s)(a²+b²)^{-s/2} Cos(s·ArcTan(b/a))` (`Re s>0`) — a genuine
  generalization (`a→0` recovers the `Sin/Cos` rows, `b→0` the `Exp` row).

A **Frullani pre-pass** (run on the whole integrand before Expand, since each
half is individually divergent) recognises `(f(a x)-f(b x))/x` and returns
`(f(0⁺)-f(∞)) Log(b/a)` (`a,b>0`): the scale ratio is read structurally and the
pairing verified by the exact identity `(t₁ /. x→ρx)+t₂ = 0`; the boundary values
are the finite limits of `f`. So `∫₀^∞ (e^{-2x}-e^{-5x})/x dx = Log(5/2)` and
`∫₀^∞ (ArcTan(5x)-ArcTan(2x))/x dx = (π/2) Log(5/2)`.  For a **bounded oscillatory**
`f` whose ordinary limit at ∞ does not exist — `Cos[k x]`, `Sin[k x]` with `k`
proved real and nonzero — the boundary value Frullani uses at ∞ is the **Cesàro
mean, 0**, so `∫₀^∞ (Cos(a x)-Cos(b x))/x dx = Log(b/a)` (`a,b>0`).  This is gated
to a single trig factor of a linear real argument, so `Cos[k x]²` (mean ½) and
`Cosh` (`k` not provably real) are left to the ordinary cascade.

Each application is **gated on its convergence strip** — checked by `Simplify`
(numerically for a numeric `s`, or against the supplied `Assumptions` for a
symbolic `s`), so every result is correct by construction. **When the
assumptions do not prove the strip, the value is returned as a
`ConditionalExpression[value, strip]`** (matching Wolfram) — it collapses to the
bare value once the strip is proved and to `Undefined` if it is refuted. A
provably-violated strip declines. Verification is symbolic; there is **no**
numeric crosscheck (the trig/PolyLog transforms use reflection-formula forms
regular at `s=0`, so e.g. `∫₀^∞ Sin[x]/x dx = π/2` falls out with no limit). A
sum is integrated term by term (each term must converge on its own). Out of
scope — products of three or more transcendental kernels, two-kernel products
outside the Bessel / exp·Bessel / Gaussian·Bessel families above (e.g.
`AiryAi[x]²`, a genuine Meijer-G with 1/3-step Gamma coefficients), finite
intervals, and two-sided reductions — return unevaluated, never a wrong value.

 (Ramanujan's canonical example)

 (Debye)
 (Fermi–Dirac)

 (Frullani)

 (Weber–Schafheitlin / Barnes)
`Integrate[x^(s-1) Exp[-a x] BesselJ[ν,b x], {x,0,Infinity}, Assumptions→{s+ν>0,a>0,b>0}]` → a `₂F₁` in `-b²/a²`.

### Definite integration (Newton-Leibniz)

`Integrate[f, {x, xmin, xmax}]` gives the definite integral by the
**fundamental theorem of calculus** (`src/calculus/integrate_newton_leibniz.c`,
exposed explicitly as `Integrate\`NewtonLeibniz[f, {x, xmin, xmax}]` and as
`Method -> "NewtonLeibniz"`).  The mechanism:

1. Antidifferentiate `f` with the ordinary indefinite cascade to get `F`.  If
   no closed antiderivative exists, the definite integral is left unevaluated
   — never assigned a wrong value.
2. Locate the real poles of the **integrand** `f` strictly inside
   `(xmin, xmax)` via `Integrate\`SingularPoints` (roots of
   `Denominator[Together[f]]`, **plus** the `Cos[u] == 0` / `Sin[u] == 0` zeros
   under any `Tan`/`Cot`/`Sec`/`Csc` head — these are invisible to
   `Denominator[Together[·]]` because the evaluator keeps the atomic trig head and
   re-canonicalizes `Sin/Cos` back to it, so `Tan[x]` has denominator `1`). A
   pole family that `Solve` returns periodically (`Pi/2 + C[1] Pi`) is enumerated
   to its concrete interior members. A continuous `f` has an antiderivative that
   is pole-free where `f` is, so only `f`'s own poles bound the segments.
3. Split `[xmin, xmax]` at those poles and form the telescoping sum
   `Σ (F(p_{i+1}⁻) − F(p_i⁺))`, evaluating each boundary through the `Limit`
   engine: a plain limit at infinite endpoints, a one-sided limit
   (`Direction -> "FromBelow"/"FromAbove"`) at interior poles, and direct
   substitution at ordinary finite endpoints.  Improper integrals thereby
   acquire their correct finite value; a genuinely divergent integral emits
   `Integrate::idiv` and is left unevaluated (matching Mathematica).

`Integrate[f, {x, a, b}, {y, c, d}, …]` gives the iterated multiple integral,
reduced innermost-first (the last spec is the inner integral), so an inner
bound may depend on an outer variable.

`Integrate\`SingularPoints[expr, {x, a, b}]` returns the sorted list of the
interior real poles used for the split — exposed for inspection and reuse.

**Cauchy principal value.** `Integrate[f, {x, a, b}, PrincipalValue -> True]`
computes the Cauchy principal value across an interior pole. It is defined only
for **odd-order** poles (the integrand changes sign across the pole, so the
symmetric one-sided divergences cancel); the value is then `Re[F(b) − F(a)]`
with the real-branch antiderivative — the branch-cut crossing at each pole
contributes a purely imaginary `I Pi` (residue) that `Re` removes. An
**even-order** interior pole has no principal value and emits `Integrate::idiv`.
With no interior pole the option is a no-op. Examples: `∫₋₁¹ dx/x = 0`,
`∫₀³ dx/(x−2) = −Log 2`, `∫₀² x/(x²−1) dx = ½ Log 3`.

**Continuous branch-form antiderivatives.** Many continuous periodic integrands
(e.g. `1/(2 + Cos[x])`) antidifferentiate to a Weierstrass branch form that is
*already continuous* — an `ArcTan[Tan[x/2] …]` term whose jump is exactly
cancelled by an accompanying `Floor` step — so `F(b) − F(a)` by direct
substitution is the correct value with no interior split.  To stay correct when
an antiderivative could instead carry a genuine *uncorrected* branch jump, a
result built from a step head (`Floor`, `Sign`, …), an inverse-trig node over
a pole-bearing trig head, or a **radical whose radicand vanishes inside the
interval** (the real branch can kink there) is accepted only when a numeric
`NIntegrate` cross-check agrees; otherwise the integral is left unevaluated
rather than risking a wrong value. As a final net, a real-looking integrand
whose folded value has a **provably-nonzero imaginary part** is treated as a
branch-crossing artifact (`Integrate::idiv`), never returned as a complex number.

### Complex line / contour integration

When a `{x, a, b}` spec has a **non-real endpoint**, or when a spec lists more
than two points `{x, z0, z1, …, zn}`, `Integrate` evaluates the **contour
integral** of `f` along the straight segments `z0 → z1 → … → zn` in the complex
plane (`src/calculus/integrate_line.c`, exposed explicitly as
`Integrate\`LineIntegral[f, {x, z0, …, zn}]`).  Real-endpoint two-point specs
still go through the real-axis Newton-Leibniz path above.

Each segment `a → b` is parametrised by a **real** parameter,
`γ(t) = a + t (b − a)`, `t ∈ [0, 1]`, which reduces the complex problem to the
real machinery already in place:

1. Antidifferentiate `f` in `x` to get `F` (bail if unknown).
2. **On-path singularities** become real roots `t* ∈ (0, 1)` of
   `Denominator[Together[f(γ(t))]]` — a singularity strictly on the contour
   makes the integral divergent (`Integrate::idiv`, left unevaluated).
   `Integrate\`PathSingularPoints[f, {x, z0, …, zn}]` returns those points.
3. The segment value is the continuous change of `F` along the segment, with
   endpoint values taken as **real one-sided limits in `t`** when substitution
   is singular (so a complex-ray approach is a real one-sided limit the `Limit`
   engine can take).  For rational integrands whose antiderivative is a sum of
   logarithms / inverse-tangents of **affine** arguments, the branch-correct
   value is recovered by combining each into a single principal `Log` of a ratio
   `Log[(u(b))/(u(a))]` — exact because a straight segment subtends an angle
   `< π` at any point off it, so closed-contour residues come out exactly.
4. Every segment value is numerically cross-checked against a complex quadrature
   of `f(γ(t)) γ'(t)`; an uncorrectable branch crossing leaves the integral
   unevaluated rather than returning a wrong branch.

### Contour / residue-theorem definite integration

For the classical families of **improper** and **periodic** real integrals that
complex analysis dispatches by summing residues over the poles enclosed by a
standard contour, `Integrate[f, {x, a, b}]` runs a residue-theorem method
**before** Newton-Leibniz (also reachable as `Method -> "Residue"` or
`Integrate`ContourResidue[f, {x, a, b}]`).  Each answer is
**correct-by-construction**: once a family's structural gates hold, the residue
theorem gives the exact value, so there is *no* numeric quadrature crosscheck.
The only post-hoc gate is a self-consistency check that the residue sum closed to
a scalar (no surviving `x`/`Root`) and — for the real-valued families — that its
imaginary part vanishes (a residual `Im` would betray a mis-classified pole); a
failure returns unevaluated and the Newton-Leibniz path takes over.  Four
recognizers:

- **Rational on `(-∞, ∞)`** — `f = P/Q` with `deg Q ≥ deg P + 2` and **no real
  pole**: value `= 2 π i · Σ Res` over the poles in the upper half-plane.
- **Fourier / Jordan on `(-∞, ∞)`** — `f = R(x) · K` with
  `K ∈ {Cos[a x], Sin[a x], Exp[I a x]}`, `R` rational with `deg`-drop `≥ 1`,
  `a` a nonzero real: with `J = 2 π i · Σ_UHP Res[R Exp[I a x]]` (for `a > 0`),
  `∫ R Cos = Re[J]` and `∫ R Sin = Im[J]`. The bare complex-exponential kernel is
  recognised in both spellings — `Exp[I a x]` and the evaluator-normalised
  `Power[E, I a x]` — and returns `J` directly (`∫ R Exp[I a x] = J`); for
  negative `a` the lower half-plane is closed. Real-exponent `Exp[a x]` is *not*
  a Fourier kernel (it routes to the rectangular family). `Cos`/`Sin` are closed
  via `Re`/`Im` of the single decaying exponential, not split into two.
- **Rational-in-`{Sin, Cos}` over a full period** `(0, 2π)` or `(-π, π)` — via
  `z = Exp[I x]` on the unit circle: value `= 2 π i · Σ Res` over the poles
  inside the unit disk.  The parameters may be **symbolic under assumptions**,
  including a *coupled* constraint such as `a > b > 0`: the in-disk pole is
  selected at a representative point of the assumption region (from the
  independent bounds, or `FindInstance` when a parameter is bounded only by
  another), and a higher-order parametric pole's residue is taken by the
  analytic-part derivative `Res = (1/(m-1)!) d^{m-1}[N/(lead ∏_{r≠z₀}(z-r)^{m_r})]|_{z₀}`
  — fast where a Laurent series of a symbolic radical pole is not.  E.g.
  `Integrate[1/(a + b Cos[x])^3, {x,0,2π}, Assumptions -> a > b > 0]
  = π(2a²+b²)/(a²−b²)^(5/2)`, and `1/(a + b Cos[x]) → 2π/√(a²−b²)`.  A combined
  `b Cos[x] + c Sin[x]` is handled too: if the `FindInstance` representative is
  degenerate (e.g. `b = c = 0` for `a > √(b²+c²)`, where the pole quadratic loses
  degree), the family retries at a point with every parameter nonzero — sound
  because the integral is analytic wherever no pole sits on the unit circle (a
  circle pole still returns `Integrate::idiv`).  So
  `Integrate[1/(a + b Cos[x] + c Sin[x])^2, {x,0,2π}, Assumptions -> a > √(b²+c²)]
  = 2π a/(a²−b²−c²)^(3/2)`.
- **Parametrized contour over a full period** — a *transcendental* integrand of
  the form `g(Exp[I t]) · (I c Exp[I t])` (the Jacobian of `z = c Exp[I t]`):
  the same `t → −I Log[u]` substitution leaves `G(u) = integrand/(I u)`, whose
  period integral is the contour integral `∮_{|z|=c} g = 2 π i · Σ_{|u_k|<1}
  Res[G, u_k]`.  This is the non-rational sibling of the unit-circle family,
  covering an **essential singularity** at `0`: when the ordinary residue engine
  declines there, the residue is read as the `w^1` coefficient of `G(1/w)` (the
  Laurent `u^{-1}` coefficient), so `∮_{|z|=1} Exp[1/z] Sin[1/z] = 2 π i`,
  `∮_{|z|=2} z^5 Cos[1/z²] = 0`, and `∮_{|z|=1} Exp[2z]/(z⁴(z−½)) =
  2 π i (16 e − 128/3)` (order-4 pole at `0` plus a simple pole at `½`).  The
  value is genuinely complex, so only the finite-scalar (not the real) gate
  applies.  A pole on `|u| = 1` makes the period integral divergent.
- **Removable axis singularity (Fourier)** — a **simple** real-axis pole of `R`
  at which the kernel vanishes (so `f = R·K` is analytic there, e.g. `Sin[x]/x`
  at `0`) contributes a **half residue** `π i · Res`, the indented-contour value,
  giving `∫ Sin[x]/x = π`.  A *genuine* axis pole (kernel nonzero there, e.g.
  `Cos[x]/x`) makes the ordinary integral diverge and returns unevaluated —
  plain `Integrate` does not compute a principal value.  With **symbolic
  parameters** the axis half-residue and an enclosed pole may carry opposite
  signs of `I`; the `Re`/`Im` extraction uses `ComplexExpand` (parameters are
  real) so both are conjugated correctly, e.g.
  `∫ Sin[a x]/(x (x²+b²)) = (π/b²)(1 − e^{−a b})` for `a, b > 0`.
- **Multi-frequency Fourier (difference of kernels) on `(-∞, ∞)`** — `f = R(x) ·
  Σ_i C_i K[ω_i x]` with the `K` all one type (`Cos` or `Sin`), the `ω_i` of one
  sign, and `R` rational.  The sum lifts to `Σ_i C_i e^{i ω_i x}` and the single
  half-plane contour closes it; `∫ R·(Σ C_i Cos) = Re[π i Σ Res]`,
  `Σ Sin → Im`.  The point is a **difference whose pieces diverge separately** but
  whose sum is integrable because the numerator cancels the axis singularity: a
  real-axis pole of `R` is admitted only when `Limit[f, x → z₀]` is finite there
  (the cancellation has reduced the order), contributing the indented half
  residue — so `∫₀^∞ (Cos[a x] − Cos[b x])/x² = π(b − a)/2` (the `x = 0` double
  pole of `1/x²` is halved to a simple one by `Cos[a x] − Cos[b x] = O(x²)`),
  while a single `Cos[a x]/x²` (no cancellation, genuine double pole) correctly
  **declines** as divergent.
- **Log-over-quadratic on `(-∞, ∞)`** — `f = C Log[A + x²]/(E x² + D)`, `A > 0`,
  `D/E > 0`: the Feynman-parameter value `C (2π/√(D/E)) Log[√A + √(D/E)]`
  (differentiate in the log parameter `A`, integrate the resulting rational).
  The log argument is monic in `x²` with no linear term.  So
  `∫₋∞^∞ Log[x²+a²]/(x²+b²) = (2π/b) Log(a+b)` and, via the even-fold,
  `∫₀^∞ Log[1+x²]/(1+x²) = π Log 2`.  The companion **`Log[Abs[affine]]` over a
  quadratic** `∫₋∞^∞ C Log[Abs[α x+β]]/(E x²+D) = (C/E)(π/p)(Log|α| +
  ½ Log(r²+p²))`, `p=√(D/E)`, `r=−β/α` (`Log|x−r| = Re Log(x−r)`), is a dedicated
  recognizer that runs *before* the piecewise splitter — which would otherwise
  cut at `x=r` into two half-infinite pieces that do not close — so
  `∫₋∞^∞ Log|x−a|/(x²+b²) = π Log(a²+b²)/(2b)`.

On the **half-line** `(0, ∞)` a closed-form **hyperbolic-kernel table** runs
*first* (before the even-fold), covering
`∫₀^∞ C x^k {1, Sin[p x], Cos[p x]} {Sech, Csch}[q x]^n dx` for `q > 0` (and `p`
real on the trig arms): the five Gradshteyn–Ryzhik members
`Cos[p x]/Cosh[q x] → (π/2q) Sech[pπ/2q]`,
`Sin[p x]/Sinh[q x] → (π/2q) Tanh[pπ/2q]`,
`x Sin[p x]/Cosh[q x]² → (π/2q²) Csch[pπ/2q](pπ/2q Coth[pπ/2q] − 1)`,
`x/Sinh[q x] → π²/4q²`, `x²/Sinh[q x]² → π²/6q³`.  Only these exact tuples match
(other shapes decline and fall through); running first is what prevents the
Sinh/Csch members from otherwise hanging in the antiderivative search or the
symbolic-parameter keyhole path.

A one-line **half-line** add-on then covers even integrands:
`∫₀^∞ f = ½ ∫₋∞^∞ f`.  Three further recognizers handle branch-cut and
symbolic-exponent contours:

- **Keyhole / Mellin on `(0, ∞)`** — a branch power times a rational function,
  `f = x^p R(x)` with `p` non-integer (so `s = p + 1 ∉ ℤ`): value
  `= −π · Σ_k Res[ z^{s-1} R(z), z_k ] · e^{−iπs} / sin(π s)` over the poles
  `z_k` of `R`, with `z^{s-1}` on the branch `arg ∈ (0, 2π)`.  The residue is of
  the **full** integrand `z^{s-1} R(z)` — for a simple pole this equals
  `z_k^{s-1} Res(R, z_k)`, but for an order-≥2 pole it does not (a pure double
  pole has `Res(R)=0`), so poles of any order are handled by shifting `w=z−z_k`
  and expanding the analytic factor `(1+w/z_k)^{s-1}`.  The `e^{−iπs}/sin(π s)`
  prefactor is the exact reduction of the keyhole jump `1/(1 − e^{2πi s})`,
  landing a numeric `s` on an algebraic multiple of `π` (e.g.
  `∫₀^∞ x^{1/3}/(x²+1) = π/√3`, `∫₀^∞ √x/(1+x)² = π/2`).  Requires
  `0 < Re(s) < deg Q − deg P`.  The exponent may also be a **symbolic
  parameter** `a` under an interval assumption — the convergence interval is read
  off `s = a + 1` (e.g. `Integrate[x^a/(x+1)^3, {x,0,∞}, Assumptions -> -1 < a < 2]
  = π a(1−a)/(2 sin π a)`).  The denominator coefficients may be symbolic too,
  including under a *coupled* assumption that leaves them unbounded: the directly
  bounded exponent still fills the convergence gate while the (complex-conjugate)
  poles are classified at a `FindInstance` point, so `Integrate[x^a/(x²+2b x+c),
  {x,0,∞}, Assumptions -> -1 < a < 1 && c > b² && b > 0]` closes.  (A symbolic
  exponent is kept even when that representative point lands on an integer `s` —
  the generic `Csc[π s]` form is the right generic-parameter answer.)  Under
  `PrincipalValue -> True` a **simple pole on the positive axis** is admitted as a
  Cauchy PV half-residue — the indented contour averages the pole's two
  keyhole-branch values `z₀^p` (arg 0) and `z₀^p e^{2πi p}` (arg 2π) — so
  `∫₀^∞ √x/(x²−a²) dx = π/(2√a)` and `∫₀^∞ x^{1/3}/(x²−1) dx = π/(2√3)`; a
  higher-order axis pole (divergent PV) declines.
- **Keyhole with a logarithm on `(0, ∞)`** — `f = x^p (Log x)^m R(x)`, `R`
  rational, `m ≥ 1`.  Two complementary branches, by whether `x^p` jumps across
  the cut:
  - **`p` a non-negative integer** (folded in as `R' = x^p R`, so `a = 0`).  The
    keyhole contour of `(Log z)^{k+1} R(z)` (branch `arg ∈ [0,2π)`) gives, for
    each `k = 0..m`, a lower-triangular relation among `I_j = ∫₀^∞ (Log x)^j R dx`
    and the pole-residue sums `S_{k+1} = Σ_l Res[(Log z)^{k+1} R, z_l]`; solving
    it upward from `I_0` yields the requested `I_m` **entirely from residues**,
    no auxiliary quadrature.  E.g. `∫₀^∞ Log[x]/(1+x⁶) = −√3 π²/18`,
    `∫₀^∞ (Log x)²/(x²+x+1) = 16 π³/(81√3)`, `∫₀^∞ Log[x]/(1+x⁴) = −√2 π²/16`.
  - **`p` non-integer** (so `z^p` itself jumps, `1 − e^{2πi p} ≠ 0`).  The
    keyhole of `z^p (Log z)^m R(z)` gives the lower-triangular system
    `(1 − e^{2πi p}) I_k = 2πi T_k + e^{2πi p} Σ_{j<k} C(k,j) (2πi)^{k−j} I_j`,
    with `T_k = Σ_l Res[z^p (Log z)^k R, z_l]` on the branch; the `k=0` row is
    the plain keyhole/Mellin value.  E.g.
    `∫₀^∞ √x Log[x]/(x²+1)² = π(π−4)/(8√2)`.
  A simple pole uses the direct `z_l^p (Log z_l)^{k} Res(R, z_l)`; an order-≥2
  pole the shifted form.  Because the pole-residue sums carry the keyhole
  logarithm `Log z_l = Log|z_l| + i·arg z_l` evaluated at the (algebraic) poles,
  the closed form collapses only because `Arg` reduces exactly on roots of unity
  (`Arg[(-1)^{1/6}] = π/6`); a numeric `p`/`m` lands on an algebraic multiple of
  `π^{m+1}`.  A **simple pole on the branch cut `(0, ∞)`** (integer-`p` branch)
  enters as a Cauchy principal value: the contour indents around it above and
  below, so it contributes the *average* of its two keyhole-branch residues,
  `½ Res(R, z₀)[(Log z₀)^k + (Log z₀ + 2πi)^k]`.  A branch-point pole at `z₀ = 1`
  carrying the `Log` factor is **removable** (`Log 1 = 0` cancels it) and is
  admitted with no option — `∫₀^∞ Log[x]/(x³−1) = 4π²/27` (the integrand is in
  fact continuous at `x = 1`).  A genuine axis pole (`z₀ ≠ 1`) needs
  `PrincipalValue -> True` and then returns the principal value, e.g.
  `∫₀^∞ Log[x]/(x²−4) dx = π²/8`; a higher-order axis pole (divergent PV) declines.
  Under a symbolic parameter the keyhole close skips `ComplexExpand` (which
  balloons on the symbolic pole through the `(Log z)^m` powers and sends
  `Simplify` into exponential matcher backtracking) and takes the plain
  `Simplify` under the assumptions instead.
- **Pure rational principal value on `(0, ∞)`** — the same `m = 0` keyhole, under
  `PrincipalValue -> True`, closes a rational `R(x)` with a simple positive-axis
  pole: `Integrate[1/(x³−a³), {x,0,∞}, PrincipalValue -> True, Assumptions ->
  a>0] = −π/(3√3 a²)`, `∫₀^∞ x²/(x⁴−a⁴) = π/(4a)`.  It runs before the even-fold,
  so the even members route to the keyhole rather than the rational family's
  divergence preemption.  A scaled positive-axis pole (`z₀ = b ≠ 1`) is handled
  too: the `(1 + e^{2πip})` branch average can leave a spurious imaginary part that
  `Simplify` does not discharge for symbolic `b`, so a real-integrand on-axis PV is
  projected onto its (provably real) real part — `∫₀^∞ x^(a−1)/(x²−b²) PV =
  −(π/(2b^{2−a})) Cot[πa/2]`.
- **Whole-line oscillatory principal value** — `Integrate[R(x){Cos,Sin,Exp[I a x]},
  {x,−∞,∞}, PrincipalValue -> True]` with `R` rational and a **simple genuine**
  pole on the real axis.  Closing in the half-plane Jordan's lemma selects (`a > 0`
  upper, `a < 0` lower), the axis pole contributes the indented half-residue
  `±iπ·Res` oriented by the closure: `∫ Exp[I a x]/(x−b) dx = iπ e^{iab}`,
  `∫ Exp[I a x]/(x²−b²) = −π Sin[ab]/b`, `∫ x Exp[I a x]/(x²−b²) = iπ Cos[ab]`
  `R·Cos`/`R·Sin` take the real/imaginary part (`R` real on the axis), and the
  even members fold to the half-line — `∫₀^∞ Cos[a x]/(x²−b²) PV = −π Sin[ab]/(2b)`,
  `∫₀^∞ x Sin[a x]/(x²−b²) PV = π Cos[ab]/2`.  A removable pole mixed with a genuine
  PV pole declines (the Re/Im split is unreliable there), as does a double axis pole
  (divergent PV).
- **Log-weighted Mellin on `(0, ∞)`** — `Integrate[x^(c−1) Log[x]^k/(1 ∓ x^d), {x,
  0,∞}]`, `k ≥ 1`: the `k`-th derivative in the Mellin exponent of the base value
  (`(π/d)Cot[cπ/d]` for `1 − x^d`, taken as a principal value; `(π/d)Csc[cπ/d]` for
  the convergent `1 + x^d`), evaluated at `s = c` — `∫₀^∞ Log[x] x^(a−1)/(x+1) =
  −π² Csc[πa] Cot[πa]`, `∫₀^∞ Log[x]² x^(a−1)/(x+1) = π³ Csc[πa]³ (1 + Cos[πa]²)`,
  `∫₀^∞ Log[x] x^(a−1)/(x−1) PV = π² Csc[πa]²`.  The strip `0 < c < d` is proved
  against the assumptions.
- **Finite interval `(-1, 1)` with the Chebyshev weight** — `Integrate[R(x) /
  Sqrt[1 - x^2], {x, -1, 1}]`, `R` rational: `x = Cos[t]` maps it to
  `(1/2) ∫₀^{2π} R(Cos t) dt`, evaluated by the unit-circle family.  It runs in
  the residue stage, *before* Newton–Leibniz, so it also pre-empts the FTC
  branch's wrong-sign antiderivative — `∫₋₁¹ dx/((1+x²)√(1−x²)) = π/√2` (FTC
  alone returned `−π/√2`), `∫₋₁¹ dx/√(1−x²) = π`.
- **Mellin–Barnes / Bromwich vertical line** — `Integrate[F(s) x^{-s}, {s, c−i∞,
  c+i∞}]` with `F(s) = const · Γ(A + B s)` (`B > 0`): the inverse-Mellin contour,
  closed to the left, equals `2πi · Σ_k Res` over the Γ-pole ladder `s = −(A+k)/B`
  (`Res[Γ(A+Bs)] = (1/B)(−1)^k/k!`); the residue series closes to the transform's
  closed form — `∫ Γ(s) x^{−s} ds = 2πi e^{−x}` (Case 11, formerly a silent `0`),
  `∫ Γ(2s) x^{−s} ds = iπ e^{−√x}`.  The endpoint `c ± i∞` is recognised as a
  directed infinity (it no longer collapses to a real `∞`); the orientation sets
  the sign, and a non-Mellin–Barnes integrand or a non-closing series declines.
- **Sector on `(0, ∞)`** — `f = C x^m/(c + x^n)` with the exponents `m`, `n`
  possibly **symbolic parameters**: the wedge of angle `2π/n` gives
  `C (π/n) c^{s/n − 1} csc(π s/n)`, `s = m + 1`, under `c > 0`, `n > 0`,
  `0 < s < n` (proved over the assumed region).  This is the one family admitting
  a symbolic `n` (the keyhole cannot enumerate `n` poles), powering
  `Integrate[1/(1 + x^n), {x, 0, ∞}, Assumptions -> n > 1] = (π/n) csc(π/n)` and
  `Integrate[x^{2m}/(1 + x^{2n}), …] = (π/(2n)) csc(π(2m+1)/(2n))` (state the
  convergence as `n ≥ m + 1` — `Refine` does not take the integer step
  `n > m ∧ ℤ ⟹ n ≥ m + 1`).
- **Generalized Beta on `(0, ∞)`** — `C x^a (x + b)^{−c}` with a NON-integer `c`
  (a branch point at `−b` the rational Mellin engine cannot take): `x = b t`
  gives `C b^{a+1−c} B(a+1, c−a−1) = C b^{a+1−c} Γ(a+1) Γ(c−a−1)/Γ(c)`, under
  `b > 0`, `a + 1 > 0`, `c − a − 1 > 0`.  So
  `Integrate[x^a/(x+b)^c, {x, 0, ∞}, Assumptions -> b>0 && c>a+1 && a>-1] =
  b^{a+1−c} Γ(a+1) Γ(c−a−1)/Γ(c)`.
- **Mellin after a power substitution on `(0, ∞)`** — `C x^{μ−1} G(κ x^ν)` with
  `G ∈ {Exp, Sin, Cos}` and `ν` a positive integer.  `u = x^ν` turns it into
  `(C/ν)` times the Mellin transform of `G` at `s = μ/ν` — the
  sector/wedge contour: `Exp[−c u] → Γ(s) c^{−s}` (`Re s > 0`),
  `Sin[k u] → Γ(s) k^{−s} Sin[π s/2]`, `Cos[k u] → Γ(s) k^{−s} Cos[π s/2]`
  (`0 < Re s < 1`).  The Gamma / generalized-Fresnel family:
  `∫₀^∞ x^{s−1} e^{−a x²} = ½ a^{−s/2} Γ(s/2)` (state `s > 0`, not `Re[s] > 0`),
  `∫₀^∞ x^p Sin[x²] = ½ Γ((p+1)/2) Sin[π(p+1)/4]`, `∫₀^∞ Sin[x³] = Γ(4/3)/2`.
- **Rectangular / quasi-periodic on `(-∞, ∞)`** — reduced to the half-line core
  by the exact substitution `w = Exp[x]`
  (`∫_{-∞}^∞ f dx = ∫₀^∞ f(Log w)/w dw`).  Two shapes close: `f = Exp[c x] R(Exp[x])`
  (period `2πi`) becomes the branch-power Mellin integrand `w^{c-1} R(w)`, so
  `Integrate[Exp[a x]/(Exp[x]+1), {x, -∞, ∞}, Assumptions -> 0 < a < 1] = π csc(π a)`
  and `f = x^k/Sinh x` and its hyperbolic relatives — the `Sinh`/`Cosh`/… of
  `Log w` rationalise under `TrigToExp` (`Sinh[Log w] = (w − 1/w)/2`) while a
  polynomial `x^k` factor becomes `(Log w)^k`, handing the **keyhole-log** family
  the integrand `(Log w)^k R(w)` with its removable `w = 1` (i.e. `x = 0`) pole —
  `∫_{-∞}^∞ x/Sinh[x] dx = π²/2`.
- **Gaussian / shifted rectangle on `(-∞, ∞)`** — an *entire* integrand
  `E^(A x² + B x + C)` times a Fourier/exponential kernel, `Re A < 0`.  `e^{A z²}`
  is entire, so the line of integration shifts freely (a rectangle enclosing no
  poles) — the pole-free degenerate case of the rectangular contour.  `TrigToExp`
  turns any `Cos`/`Sin`/`Cosh`/`Sinh` kernel into a sum of `E^(linear)`, so the
  integrand becomes `Σ_i C_i E^(A x² + B_i x + D_i)` and each Gaussian closes by
  completing the square, `∫_{-∞}^∞ e^{A x² + B x + D} dx = √(−π/A)·e^{D − B²/(4A)}`:
  `∫_{-∞}^∞ e^{−x²} Cos[2 a x] dx = √π e^{−a²}`.  The leading coefficient's sign is
  verified over the whole assumption region — by `Refine[A < 0, assumptions]`, so
  a *symbolic* coefficient such as `A = −a` under `a > 0` is admitted (the plain
  interval test saw the compound `−a` as unbounded), while a positive one (a
  divergent integral) declines.  With a symbolic `A` and a trig kernel the result
  carries the parameters: `∫_{-∞}^∞ e^{−a x²} Cos[b x] dx = √(π/a) e^{−b²/(4a)}`
  (`a > 0`), and the even half-line form is half of it.  Not a residue
  computation (no poles), but the standard contour route for the family.
- **Periodic-strip (hyperbolic) rectangle on `(-∞, ∞)`** — `N(x)/Cosh[b x]^n`,
  integer `n ≥ 1`.  The quasi-period `Cosh[b(x + iπ/b)] = −Cosh[b x]` folds the
  infinite ladder of (order-`n`) poles into the single enclosed pole at
  `x = iπ/(2b)`: `∫_{-∞}^∞ e^{α x}/Cosh[b x] dx = (π/b) Sec[α π/(2b)]` for
  `|Re α| < Re b`.  `TrigToExp` turns a `Cosh`/`Sinh`/`Exp` numerator into a sum of
  such terms, so `∫ e^{a x}/Cosh[π x] = Sec[a/2]`,
  `∫ Cosh[a x]/Cosh[b x] = (π/b) Sec[πa/(2b)]`, and
  `∫ e^{i a x}/Cosh[b x] = (π/b) Sech[πa/(2b)]`.  For `n ≥ 2` the per-term value is
  the Gamma reflection of the Mellin-Barnes integral (equivalently the order-`n`
  residue of the fold), `2^{n−1}/(b (n−1)!)·R(w)·∏(w)`, `w = α/(2b)`, with
  `R = π/Cos[πw]`, `∏_{j=0}^{m−1}((j+½)²−w²)` for odd `n = 2m+1` and
  `R = πw/Sin[πw]`, `∏_{j=1}^{m−1}(j²−w²)` for even `n = 2m`; the convergence gate
  widens to `|Re α| < n Re b`.  So the double pole
  `∫_{-∞}^∞ Cos[a x]/Cosh[b x]² dx = π a Csch[πa/(2b)]/b²`, the triple pole
  `Cos[a x]/Cosh[b x]³ → π(a²+b²) Sech[πa/(2b)]/(2b³)`, and the `Sin`/`Sinh`
  numerators give a correct (unsimplified) `Csc`-sum.  Handles the symbolic
  non-integer `b` the `Exp[x]` rectangular reduction cannot rationalise.  A common
  positive scale on `x` in the hyperbolic arguments is normalised first
  (`Sinh[a x] → Sinh[x]`), so `∫₀^∞ x/Sinh[a x] = π²/(4a²)`.

**Assumptions and symbolic parameters.**  An `Integrate[f, {x, a, b},
Assumptions -> …]` option lets the residue families evaluate integrals whose
parameters are symbolic (`a > 0`, `0 < a < 1`, `n > 1`, …).  The recognizers
classify a parameter-dependent pole or kernel frequency by reading its sign at a
single generic point of the region the assumptions pin (a sign-consistent
instantiation), while the residue arithmetic stays fully symbolic — so the
closed form is still **correct by construction, with no numeric crosscheck**.
Convergence/applicability gates (`n > m + 1`, `0 < s < deg`-drop, `c > 0`) are
verified against the assumption-**guaranteed** interval bounds, not the sample
point, so an under-constrained problem (e.g. only `n > 0` for the sector family)
is refused rather than guessed; a parameter the assumptions leave two-sided
unbounded is likewise refused.  A sign gate that an interval cannot settle (a
compound coefficient like `−a`, or the open bound `a > 0` whose recorded lower
end is `0`) is proved instead by `Refine[… , assumptions]` — so the Gaussian's
`Re A < 0` admits `A = −a` under `a > 0`.  Radical pole locations are
`PowerExpand`-cleaned under all-positive parameters (`Sqrt[-4 a²] → 2 I a`) and
the closed form of a two-conjugate-pair rational is further collapsed with a
positivity-aware `PowerExpand`/`Refine` pass (`Sqrt[(a+b)²(a−b)²] → (a+b)(a−b)`),
so `∫₀^∞ 1/((x²+a²)(x²+b²)) = π/(2 a b (a+b))` rather than a nested-surd surface.
Real/imaginary parts of a Fourier contour value are taken with
`ComplexExpand[Re/Im[·]]` (parameters are real), which conjugates both `+I` and a
stored `−I` atom correctly.  *Not yet simplified* (a `Simplify`/`Refine` gap, not
a residue-method one): a keyhole-log answer with a symbolic positive parameter can
retain `Arg[I a]` (e.g. `∫₀^∞ Log[x]/(x²+a²)`, value `π Log[a]/(2a)`) or a
`(1 − E^(2 I π a))` factor.

**A deliberate decline.**  `∫₀^∞ ArcTan[a x]/(x (1 + b² x²)) dx = (π/2) Log[1 + a/b]`
returns unevaluated under `Method -> "Residue"`.  The integrand has a **branch
cut** (`ArcTan`), not isolated poles; its value is reached by parametric
differentiation (Feynman — `∂_a I = ∫₀^∞ dx/((1 + a²x²)(1 + b²x²)) = π/(2(a+b))`,
then integrate in `a`), which is not the residue theorem.  No Feynman engine is
wired onto the residue path, so this is an honest decline rather than a
misclassification.

The residue sums are closed with `RootReduce` (algebraic families) or the
`Conjugate` identities for `Re[J]`/`Im[J]` (the Fourier family); a value that
still contains an unreduced `Root`, or a surviving imaginary part, is treated as
a mis-fire and returned unevaluated.  Negative controls such as
`Integrate[1/(1 + x^3), {x, -Infinity, Infinity}]` and
`Integrate[1/(x^2 - 1), {x, -Infinity, Infinity}, Method -> "Residue"]`
(genuine real-axis poles),
`Integrate[Cos[x]/x, {x, -Infinity, Infinity}, Method -> "Residue"]`
(genuine axis pole, kernel nonzero),
`Integrate[1/Sqrt[1 + x^4], {x, -Infinity, Infinity}]` (branch point, not
rational), and `Integrate[1/(2 + Cos[x]), {x, 0, Pi}]` (not a full period) all
stay unevaluated.  The keyhole/Mellin, keyhole-with-log (integer *and*
non-integer power, including a **principal-value** pole on the positive axis),
Chebyshev-weight, Mellin–Barnes, sector, rectangular/hyperbolic and
Gaussian/shifted-rectangle families described above extend the reach to
branch-cut, log-weighted, principal-value, Gaussian and symbolic-exponent
contours.  A few integrands remain **deliberate, documented declines** (returned
unevaluated, never a forced or wrong value): `∫₀^∞ Log[x]/Cosh[x] dx` — its value
is the Mellin transform of `sech` differentiated at `s = 1`, a Dirichlet-beta
derivative (a `Γ(1/4)`-level constant), not a residue sum, and the `w = Exp[x]`
reduction fails because `Log x` becomes `Log[Log w]`; and the Hankel-type
`∫₁^{-∞} eˣ x^{-s} dx`, whose contour is ambiguous / divergent in the ordinary
sense.

**Output simplification.**  Every family's value is finally passed, once, through
a sound assumption-aware cleanup: `Refine` collapses `Sqrt[a²]→a`, `Abs[a]→a`,
`Arg[pos]→0`, `Arg[neg]→π` **only** where the assumptions prove the sign (never the
unsound `PowerExpand`), a plain `Simplify` combines the terms, and — on a small,
radical- and trig-free elementary rational-log form carrying a cancelling imaginary
unit — a bounded `FullSimplify` finisher lands the real closed form.  The result
replaces the raw form only when strictly smaller, so a legitimately-messy value
(an `Arg[I a]` that `Refine` cannot reduce; a mixed-sign radical whose sign is not
assumed; an already-clean `Sec` form) is never enlarged or corrupted.  This is what
turns the principal-value `∫₀^∞ Log[x]/(x²−a²) dx` (`a > 0`) from a large
`Arg`/`Abs`/`(a⁴)^{1/4}` expression into the clean `π²/(4a)`.

The `Integrate`` package also exposes the lower-level helpers
`Integrate`HermiteReduce`, `Integrate`IntegratePolynomial`,
`Integrate`BronsteinRational` (the explicit form),
`Integrate`IntRationalLogPart` (Phase 2's LRT computation),
`Integrate`RischTranscendental` (recursive transcendental Risch), `Integrate`LinearRadicals`
(linear-radical substitution), `Integrate`QuadraticRadicals`
(quadratic-radical Euler substitution), `Integrate`LinearRatioRadicals`
(linear-fractional / Möbius radical substitution), `Integrate`Weierstrass`
(continuous `Tan[x/2]` / `Tanh[x/2]` substitution), `Integrate`CRCTable`
(table lookup), `Integrate`Undefined` (unknown-function integration,
Roach §1.7), and the unit-test helpers `Integrate`Helpers`Content`,
`...`Primitive`, `...`Monic`, `...`LeadingCoefficient`,
`...`SquareFree`, `...`ExtractConstants`, `...`ApartList`.  All are
`Protected`.

### Integrate`CRCTable

`Integrate`CRCTable[f, x]` looks `f` up in the CRC Standard Mathematical
Tables (31st ed., 600+ formulas).  The rules live in
`src/internal/CRCMathTablesIntegrals.m`, internally on the head
`IntegrateTable`; `Integrate`CRCTable` is a thin wrapper around it.
The .m file is `Get`-loaded on the first invocation of the CRCTable
stage rather than at startup, so sessions that never call `Integrate`
pay nothing for the table.

Every recursive rule in the table carries an `IntegerQ[index] && index
> base` (or `< base`) guard sufficient to guarantee termination via
first-principles analysis of the reduction direction.  Without these
guards rules such as Formula 49 (the `1/(x^2 - c^2)^n` reduction)
would diverge on negative or non-integer `n`.  As defence-in-depth,
the C dispatcher caps CRC-rule recursion depth at 256 levels and
emits `Integrate`CRCTable::depth` rather than locking up on any rule
that escapes the audit.

The table covers the inverse-trig and inverse-hyperbolic families
(Formulas 427–464 and their `427h`–`464h` analogs) among many others.
Rules whose denominator carries a *squared* coefficient
(`Sqrt[1 ∓ a^2 x^2]`, `(1 ± a^2 x^2)`) bind that coefficient linearly as
`a_` via `c_ + a_. x_^2` and recover the linear coefficient as `Sqrt[±a]`
on the RHS (with a `Condition` linking it to the numerator coefficient) —
the matcher does not invert a square, so `a_^2` in a pattern will not
bind. Some more elaborate multi-argument `/;`-guarded rules still do not
fire; this is a separate issue tracked under the matcher work.

### Integrate`Undefined

`Integrate`Undefined[f, x]` integrates expressions that are rational in
**undefined functions** `u[x]` and their derivatives, following Kelly
Roach, "Indefinite and Definite Integration" (1992), §1.7 ("Undefined
Functions").  Each undefined function value `u[g]` and its derivative
tower `u'[g], u''[g], …` is treated as a differential-field generator
the integrand is reduced by recognising integration-by-parts /
total-derivative structure in the top generator.  A single inner call to
the rational integrator over a substituted generator symbol subsumes
Roach's polynomial, fraction, and log parts (so `1/u -> Log[u]`,
`1/u^2 -> -1/u`, and `a/(a^2+u'^2) -> ArcTan` all fall out).  Composite
arguments are handled via the chain rule (`g' = D[g, x]`).  A logarithm of
an unknown-function expression, `Log[eta]`, is itself recognised as a
transcendental generator and reduced by parts
(`Integrate[C L + D, x] = G L - Integrate[G L' - D, x]` with
`G = Integrate[C, x]`, `L' = eta'/eta`), with self-referential resolution
for perfect powers (`Integrate[Log[f] f'/f, x] = Log[f]^2/2`).  The stage
is gated to only run when `f` contains an undefined-function derivative or
such a logarithm; genuinely non-elementary integrands (e.g. `f'[x] g'[x]`,
`f'[x]^2`) are left unevaluated, with a cycle guard preventing the by-parts
recursion from looping.  `Protected`.

Known limitations: transcendental generators other than `Log` (e.g.
`ArcTan[eta]`, `Exp[eta]` with `eta` containing an unknown function) are
not yet recognised; nested unknown arguments `f[g[x]]` (with `g` itself
undefined) are deferred.

### Integrate`DerivativeDivides

`Integrate`DerivativeDivides[f, x]` integrates **by substitution** — the
classical "derivative-divides" method (Moses' SIN, Maxima's `diffdiv`).  It
recognises an integrand of the shape `f(x) = c · h(u(x)) · u'(x)`, reduces the
problem to `Integrate[h[u], u]`, and back-substitutes `u -> u(x)`.
Implemented in `src/calculus/integrate_derivdivides.c`.

Candidate kernels `u(x)` are every distinct subexpression of `f` that depends
on `x`, except `x` and `f` themselves (the C analogue of
`Cases[Union[Level[f, {0, Infinity}]], e_ /; !FreeQ[e, x]]`).  Two
complementary strategies are tried per kernel, each followed by a
**verification gate** `PossibleZeroQ[D[result, x] - f] === True` (numeric
sampling; as of 2026-06-09 the former rigorous `Simplify` confirm was dropped —
it cost ~1.1 s of `trigrat` normalization on the radical-trig cases for no gain
over the sampler):

1. **Direct quotient** — `q = Cancel[Together[f / D[u(x), x]]]`, then
   `q /. u(x) -> u`; accepted when free of `x`.  Cheap, emits no diagnostics,
   handles transcendental compositions (`x Exp[x^2]`), and **selects the
   correct radical branch inherently** (no squaring).  Tried first.

2. **Eliminate / Solve** — builds the differential relation
   `Eliminate[{Dt[y] == f Dt[x], u == u(x), Dt[u == u(x)]}, {x, Dt[x]}]`,
   solves for `Dt[y]`, and reduces each branch with `Factor //@ `,
   `PowerExpand` and `Cancel[. / Dt[u]]`.  It closes integrands the direct
   strategy cannot — including **radical substitutions** like
   `u = Sqrt[Tan[x]]` (reducing `Integrate[Sqrt[Tan[x]], x]` to the rational
   integral `2 u^2/(1 + u^4)`), and cases where `Cancel` canonicalises
   `1/Cos[x]` to `Sec[x]`, breaking the syntactic substitution.  Because the
   algebraisation can square radicals and invert functions, `Solve` returns
   several branches; the verification gate is what **selects the branch that
   differentiates back to `f`**.  This strategy runs in the Automatic cascade as
   well as under the explicit method (its Eliminate `::ifun` / `::alg`
   diagnostics are muted while the integrator drives it).  Because it is
   heavyweight (~0.1–1 s per kernel), as of 2026-06-09 it runs **only on the
   outermost integrand** — reduced sub-integrals are finished by the direct
   strategy and the rest of the cascade.

The reduced integral re-enters the full `Integrate`, so substitutions compose.
`Integrate` takes a **`TimeConstraint`** option (default **3**, `Infinity` for no
budget; `Options[Integrate]`). It gives the cascade a wall-clock budget, honoured
between stages at every recursion depth and inside this stage's kernel and branch
loops and its `Eliminate`/`Solve` calls. It is **not yet a hard bound** — a single
long-running evaluation inside one stage can overrun it — and it is deliberately
**not** allowed to abandon the rest of the cascade: doing so was measured to cost
eight closes in 125, because the stages after the expensive one are cheap
(ParallelMixedTower 0.27 s, ParallelMixedSpecial 0.34 s mean) and are often the
ones that actually close the integrand. An explicit `Method -> name` is a
deliberate request for the full search and carries no budget unless you pass
`TimeConstraint` yourself.

Where the time actually goes, measured from inside the cascade with
`MATHILDA_INTEGRATE_PROFILE=1` over 25 integrands `Integrate` cannot close: the
eight cheap exact stages ahead of this one cost **0.0%** (0.14 s in total), while
the linearity split is **43.5%** and this stage **47.3%** — 90.8% in two stages,
both because they recurse speculatively into the whole cascade (per `Plus` term and
per candidate kernel respectively, and the cascade runs at every depth). So the
cascade is ordered correctly; the cost is the fan-out.

Three guards keep the recursion finite and cheap: an **integrand memo** that
short-circuits any integrand (canonicalised by renaming the integration variable
to a fixed sentinel) already attempted in the current top-level descent — this
breaks circular substitution chains and collapses overlapping subproblems that
would otherwise fan out exponentially (e.g. `Integrate[x Sin[x^2], x]`); the
**outermost-only** restriction on the Eliminate/Solve search above; and a hard
depth backstop (8) with per-call fresh substitution symbols.  Strict: returns
unevaluated when no substitution closes the integral.  `Protected`.

Known limitations: kernels must appear **literally** in `f` (so `Tan[x]`,
which Mathilda keeps atomic rather than `Sin[x]/Cos[x]`, exposes no `Cos[x]`
kernel — such integrands are handled by the transcendental Risch stage instead); the reduced
integral must itself close under the other methods.

### Integrate`LinearRadicals

`Integrate`LinearRadicals[f, x]` integrates a **rational function of `x` and
radicals `(a x + b)^(m/n)` that share one linear argument** by the classical
rationalising substitution `u = (a x + b)^(1/n)`, where `n = LCM` of the radical
denominators.  Implemented in `src/calculus/integrate_linrad.c`.

It first scans `f` for every `Power[base, p/q]` with `q > 1`: each `base` must be
a degree-1 polynomial in `x` and all must be the **same** linear form `a x + b`
(distinct bases — e.g. `Sqrt[x] + Sqrt[x+1]` — are rejected, as are radicals of a
non-linear argument such as `Sqrt[x^2+1]`).  With `x = (u^n - b)/a` and
`dx = (n/a) u^(n-1) du` the integrand becomes

a **rational function of `u`** that `Integrate`BronsteinRational` closes exactly.
The radical substitution reuses `poly_subst_radical_to_gen` (shared with the
algebraic-factoring path); after the substituted integrand is confirmed rational
in `{x, u}`, the reduced integral re-enters the full `Integrate` and the result
is back-substituted `u -> (a x + b)^(1/n)`.  Unlike the heuristic methods, the
result is **not** put through a `Simplify[D[result, x] - f] === 0` gate: the
substitution is an exact bijection that introduces no branch issues, so the
antiderivative is correct by construction once the rational sub-integral closes.
(Skipping the gate also avoids a prohibitively expensive — and on integrands with
symbolic parameters, non-terminating — `PossibleZeroQ`/`Simplify` over the
symbolic-radical residue.)  A depth guard (8) and per-call fresh substitution
symbols keep the recursion finite and collision-free.  Strict: returns
unevaluated when `f` is not of this form or the reduced integral does not close.
`Protected`.

### Integrate`QuadraticRadicals

`Integrate`QuadraticRadicals[f, x]` integrates a **rational function of `x` and
square roots `(a x^2 + b x + c)^(m/2)` that share one quadratic argument** by a
classical **Euler substitution**.  Implemented in
`src/calculus/integrate_quadrad.c`.

It scans `f` for every `Power[base, p/q]` whose `x`-dependent `base` is a
degree-2 polynomial: each must be a **square root** (`q == 2`) and all must be
the **same** radicand `rad = a x^2 + b x + c` (a cube root such as
`(x^2+1)^(1/3)`, a radical of a cubic such as `Sqrt[x^3+1]`, or distinct
radicands are rejected; a purely linear radical belongs to
`Integrate`LinearRadicals`).  Substituting the radicals to a fresh symbol `y`
(`rad^(p/2) -> y^p`, via `poly_subst_radical_to_gen`) must leave a rational
function of `{x, y}`.

To keep the antiderivative **real-valued**, exactly **one** substitution is
chosen by the sign of the leading coefficient `a` (and, when `a < 0`, of the
discriminant `b^2 - 4 a c`) — the routine does *not* try all three Euler forms:

| Condition | Euler substitution | Image `y = Sqrt[rad]` |
|-----------|--------------------|-----------------------|
| `a > 0` | first | `Sqrt[a] x + u` |
| `a < 0` and `b^2 - 4 a c > 0` | third | `(x - alpha) u`, with `alpha` a real root |
| `a` symbolic | first (best-effort `a > 0` branch) | `Sqrt[a] x + u` |

For numeric `a < 0` a real radicand requires real roots, so the third
substitution subsumes the second (the `c > 0` form).  Each branch yields
`x = X(u)`, `dx = X'(u) du` and the radical image above; the rationalised
integrand `f dx /. {y -> ..., x -> X}` re-enters the full `Integrate` and the
result is back-substituted `u -> U(x)`.  As with `Integrate`LinearRadicals`, the
Euler substitution is an exact bijection on the relevant domain, so the result
is **not** put through a `Simplify[D[result, x] - f] === 0` gate — it is correct
by construction once the rational sub-integral closes.  A depth guard (8) and
fresh per-call substitution symbols keep the recursion finite.  Strict: returns
unevaluated when `f` is not of this form or the reduced integral does not close.
`Protected`.

### Integrate`LinearRatioRadicals

`Integrate`LinearRatioRadicals[f, x]` integrates a **rational function of `x`
and radicals `((a x + b)/(c x + d))^(m/n)` that share one linear-fractional
(Möbius) argument** by a rationalising substitution.  Implemented in
`src/calculus/integrate_linratiorad.c`.

It scans `f` for every `Power[base, p/q]` (`|q| > 1`) whose `base` depends on
`x`; all such bases must be **structurally identical** and `n = LCM` of the
radical denominators (distinct bases are rejected).  Unlike
`Integrate`LinearRadicals` / `Integrate`QuadraticRadicals`, the base is **not**
required to be a polynomial — it is the ratio
`Times[a x + b, Power[c x + d, -1]]`, which is exactly what partitions this
method from those two (their scans demand a polynomial base, so a ratio radical
falls through to here).  The shared base is canonicalised with
`Cancel[Together[.]]` and the Möbius coefficients read off its numerator
(`a, b`) and denominator (`c, d`); the denominator must be genuinely linear
(degree 1 — a constant denominator is `Integrate`LinearRadicals`' job) and the
map non-degenerate (`a d - b c != 0`).

With `m = n` the substitution `u = ((a x + b)/(c x + d))^(1/m)` inverts the
Möbius map,

and the radicals rewrite to `u` (`r^(p/q) -> u^(p m/q)`, via
`poly_subst_radical_to_gen`).  The result must be rational in `{x, u}`; the
rationalised integrand re-enters the full `Integrate` and the antiderivative is
back-substituted `u -> ((a x + b)/(c x + d))^(1/m)`.  As with the other radical
stages, the Möbius substitution is an exact bijection on the relevant domain, so
the result is **not** put through a `Simplify[D[result, x] - f] === 0` gate — it
is correct by construction once the rational sub-integral closes.  A depth guard
(8) and fresh per-call substitution symbols keep the recursion finite.  Strict:
returns unevaluated when `f` is not of this form or the reduced integral does
not close.  `Protected`.

### Integrate`Weierstrass

`Integrate`Weierstrass[f, x]` integrates a **rational function of the
trigonometric kernels** `Sin/Cos/Tan/Cot/Sec/Csc[x]` — or the **hyperbolic
kernels** `Sinh/Cosh/Tanh/Coth/Sech/Csch[x]` — by the continuous Weierstrass
substitution of Jeffrey & Rich (*The Evaluation of Trigonometric Integrals
Avoiding Spurious Discontinuities*, ACM TOMS 20(1), 1994).  Added 2026-06-09.

Algorithm: substitute `u = Tan[x/2]` (`u = Tanh[x/2]` for hyperbolic), turning
`f` into a rational function of `u`; integrate that (recursing through
`Integrate`BronsteinRational`); back-substitute; and — for the trigonometric
case — add the secular correction `K Floor[(x - b)/p]` (`b = Pi`, `p = 2 Pi`)
that removes the spurious jump discontinuities the classical substitution
introduces at the poles of `Tan[x/2]` (odd multiples of `Pi`).  The jump `K` is
the difference of the one-sided limits of the `u`-antiderivative at `±Infinity`
if that limit diverges (a *genuine* singularity of the integrand) no correction
is applied.  A `TrigExpand` pre-pass reduces multiple/sum-angle arguments
(`Cos[2 x]`, `Cosh[x] Cosh[2 x]`, ...) to kernels of the bare variable.

### Hyperbolic case needs no `Floor` correction

`Tanh[x/2]` is a smooth,
strictly monotone bijection `R -> (-1, 1)` with no poles, so the substitution
introduces no spurious discontinuity — the back-substituted antiderivative is
already continuous (genuine singularities such as `Cosh[x] = 2` are real poles
and are correctly left in place).

The substitution is an exact identity, the rational sub-integral is closed by a
verified integrator, and `Floor' = 0` almost everywhere, so the result is
**correct by construction** — no differentiate-back gate is applied (and the
`Floor` term would defeat symbolic `D` anyway).  In the Automatic cascade only
genuine rational integrands (a kernel in a denominator) are intercepted, so
polynomial trig such as `Integrate[Sin[x], x]` keeps its cleaner table form; the
explicit `Method -> "Weierstrass"` has no such gate.  Strict: returns unevaluated
when `f` is not a rational function of the trig/hyperbolic kernels of `x` (e.g.
`x` outside a kernel, a kernel of a nonlinear argument, mixed trig + hyperbolic,
or a radical of a kernel).  `Protected`.

### InterpolatingFunction integrands

`Integrate[InterpolatingFunction[...][x], x]` returns the antiderivative as a
fresh applied `InterpolatingFunction[...][x]`, mirroring how
`D[InterpolatingFunction[...][x], x]` differentiates such objects.
Differentiation only bumps the object's derivative-order annotation; the
per-window evaluation kernels evaluate derivatives of order `>= 0` only and
cannot produce an antiderivative, so integration instead builds a genuinely new
interpolant (`src/calculus/integrate_interp.c`):

1. Read the grid x-coordinates from the object's stored table.
2. Sample the original function values `y_i = ifun[x_i]`.
3. Accumulate the antiderivative node values `F_0 = 0`,
   `F_i = F_{i-1} + Integrate_{x_{i-1}}^{x_i} ifun` by 5-point Gauss-Legendre
   quadrature (exact through degree 9 — i.e. the default/Spline/Hermite
   piecewise-cubic interpolants; very high explicit `InterpolationOrder` panels
   incur a small quadrature error).
4. Build a Hermite `InterpolatingFunction` through `{{x_i}, F_i, y_i}` — because
   the antiderivative's exact derivative is the original function, supplying
   `F'(x_i) = y_i` makes `D[Integrate[ifun[x], x], x]` round-trip to `ifun[x]`.

Only the 1-D, direct case (the applied argument is the integration variable
itself) is reduced; `Integrate[ifun[g[x]], x]` is not generally expressible as
an `InterpolatingFunction` and is left to the cascade above. The construction
also handles the derivative-annotated objects produced by `D` (integrating the
sampled derivative recovers the lower-order antiderivative). Computations use
machine doubles.

### Soundness gates on symbolic and tower stages (v0.458)

- **Symbolic-parameter residue families.** A pole's half-plane is classified at a sample point of the assumed region. The closed form is then re-derived with the poles classified at up to three further generic points of the same region, and it is withdrawn when they disagree. Radical root expressions can swap identities inside the region: `1/(1 + 2 a x^2 + b x^4)` across `a^2 = b`. A withdrawn value leaves the integral to the later stages.
- **ParallelMixedTower** (cascade and `Method -> "ParallelMixedTower"`) declines an integrand in which a kernel outside its elementary class (anything other than `Plus`, `Times`, `Power`, `Exp`, `Log`, trig, hyperbolic, their inverses and radicals) is applied to an x-dependent argument. Both it and the ParallelMixedSpecial cascade stage reject a result that mentions a symbol absent from the integrand. ParallelMixedSpecial also rejects `$Aborted`.

Regression coverage: `tests/test_mathar_regressions.c`; corpus: `tests/MatharDefiniteIntegralsCorpus.m`, `tests/MatharIndefiniteIntegralsCorpus.m`.

## Algorithm

integrate.c

```text
`Integrate[f, x]` System` dispatcher.  Cascades through three
```

stages and supports an explicit `Method -> "..."` option:

```text
  1. Integrate`BronsteinRational      — polynomial / rational integrands
  2. Integrate`RischTranscendental    — recursive transcendental Risch
  3. Integrate`CRCTable               — CRC integral table (lazy-loaded)
  4. Integrate`ParallelMixedTower     — parallel Risch-Norman over a mixed
                                        radical tower (last resort; .m)
```

Method values: "Automatic" (default, full cascade), "BronsteinRational", "RischTranscendental", "CRCTable", "ParallelMixedTower" (strict passthrough,

```text
no fallback).  RischNorman / RischNormanBlake were removed in v0.163 --
```

ParallelMixedTower subsumes both.

The CRC table is large and most sessions never need it, so its .m file is Get-loaded on first invocation of try_crctable() rather

```text
than at startup.  See LAZY_LOAD_CRC below.
```

## Performance

Against other systems, from the benchmark suite (same input, results cross-checked for agreement):

| case | Mathilda | Wolfram | Python |
|---|---:|---:|---:|
| NI 50-digit Gaussian | 6.87 s | 2.11 s | 1.1 s |
| integrate Sin[x] Exp[x] | 6.62 s | 0.726 s | 5.47 s |
| integrate Exp[x]/x (non-elementary) | 5.38 s | 0.318 s | 33.5 s |
| integrate Tan[x]^3 | 4.07 s | 0.196 s | 3.62 s |
| NI 2-D ridge | 2.54 s | 43.5 s | 0.685 s |
| Discriminant of deg 20 | 2.51 s | 0.068 s | 0.182 s |

## Implementation notes

**Algorithm.** `builtin_integrate` (src/calculus/integrate.c) first threads over a `List`
integrand, then splits on the second argument. A range or contour spec goes to
`integrate_definite`; a symbol goes to the indefinite cascade. Before that cascade runs, an
applied `InterpolatingFunction` gives a new antiderivative `InterpolatingFunction`, and a
`SeriesData` is integrated term by term (`series_integrate`).

*Indefinite cascade.* Under `Method -> Automatic` each stage runs only if every earlier stage
declined (`IP_STAGE`). The order is:

1. Undefined-function integrands (`try_undefined`, Roach §1.7).
2. `BronsteinRational`: Hermite reduction, then the Lazard–Rioboo–Trager log part with
   Rioboo's `LogToReal`.
3. Radical substitutions: `LinearRadicals`, `QuadraticRadicals` (Euler), `LinearRatioRadicals`
   (Möbius), and `ChebychevAlgebraic` (binomial differentials).
4. Recognisers: `Fresnel` (Sin/Cos of a quadratic) and `GammaPower` (`x^p E^(a x^m)` with a
   symbolic `p`, giving an incomplete Gamma).
5. Linearity over a `Plus`, then linear-argument substitution.
6. `Weierstrass`: the Jeffrey–Rich continuous half-angle form, with a secular
   `K Floor[(x-b)/p]` term.
7. `DerivativeDivides`.
8. `RischTranscendental`: recursive transcendental Risch with Cherry's special-function
   extensions (`Erf`, `ExpIntegralEi`, `LogIntegral`, `PolyLog`).
9. `CRCTable`, loaded lazily.
10. `PolyLogLadder`, loaded lazily from src/internal/PolyLogLadder.m. It covers `R(x) Log[a x + b]^k`
    over partial fractions, and a polynomial times a trig rational via `TrigToExp`.
11. `LogByParts`: `Log[g]` times a trig kernel, integrated by parts once.
12. `ParallelMixedTower`: the parallel Risch–Norman method over a mixed radical tower, from
    src/internal/mixed/ParallelMixed.m. It is skipped for pseudo-elliptic `F/R^p`.
13. `GoursatAlgebraic`: pseudo-elliptic reductions, and their cube- and fourth-root versions,
    to genus-0 curves.
14. `ParallelMixedSpecial`: the only stage allowed to answer with non-elementary kernels
    (`Erf`, `Gamma[s, z]`, `SinIntegral`/`CosIntegral`, `EllipticF`/`E`/`Pi`). It runs last
    because an elementary answer is always preferred.

The heuristic stages are accepted only after a differentiate-back check.
`Integrate::nonelem` is held until the cascade ends, and is emitted only if no later stage
answers.

*Definite driver.* `integrate_definite` folds the specs from the innermost outward. A complex
or polyline spec goes to `integrate_line_from_spec`. A real interval runs the following chain;
apart from step 1, every step is limited to `Automatic` or its own pinned method:

1. Exact identities: an empty interval gives 0, and `DiracDelta` sifting gives a point value.
2. Whole-line `Log[Abs[affine]]/quadratic`.
3. `Piecewise`: split at the breakpoints of `Abs`, `Sign`, `UnitStep`, `Min`/`Max` and
   `Piecewise`. Each branch's sign is proved.
4. `Residue`: upper half-plane, unit circle and Fourier/Jordan families. A pole on the axis
   gives `Integrate::idiv`.
5. `Symmetry`: odd integrands give 0, and even ones halve the interval.
6. The Mellin `Log`-power family, then `RationalLog` on [0, 1], `SeriesExpansion` (term by
   term, resummed to Catalan / `DirichletBeta` / `Zeta`) and `FeynmanParameter` (Frullani-type
   `(Σ c_i x^e_i)/Log[x]`).
7. `Beta`: Euler Beta, Beta·₂F₁ and half-Beta.
8. `TrigPower`, log-sine and Poisson log integrals.
9. King reflection: `x g(x)`, and the general `f(x) + f(a+b-x) = const` case.
10. Raabe's integral.
11. `Mellin`: the Ramanujan Master Theorem table, with each convergence strip proved under
    `Assumptions`.
12. The principal-value `x^(c-1)/(1-x^d)` case.
13. `IntegralRepresentation`: Laplace–Bessel, Bessel-K and Airy.
    `MeijerConvolution` follows it: `C x^alpha K1(x) K2(x)` on (0, ∞) for two kernels from
    `BesselJ`/`Y`/`K`, `E^(-p x)`, `E^(-p x^2)`, `Cos`, `Sin` and `(a^2 + x^2)^-r`. Each is written
    as a Meijer G of x², and the convolution theorem (PBM 2.24.1.1) gives one `MeijerG`.
    Convergence is proved under the assumptions; an unproved condition returns a
    `ConditionalExpression` (`src/internal/MeijerConvolution.m`, v0.460).
14. `NewtonLeibniz`: an indefinite antiderivative, split at the interior poles, with the
    endpoints and poles evaluated by one-sided `Limit`. With `PrincipalValue -> True`, interior
    poles are combined symmetrically.
15. `DiffUnderInt`: differentiation under the integral sign as a first-order ODE in the
    parameter. It is skipped under `PrincipalValue`.

A result that still contains the integration variable is discarded. With an explicit
`Assumptions`, the folded value is `Refine`d, and the refinement is rejected if it makes a
real value complex.

**Data structures.** Everything is an `Expr*` tree. Stages talk to each other through
context-qualified heads (`call_stage` evaluating ``Integrate`Stage[f, x]``): a stage that returns its
own head unevaluated, or `$Failed`, has declined. `IntegrateMethod` holds the parsed
`Method`: `parse_method_option` handles the indefinite form, including
`{"DerivativeDivides", "Substitution" -> u}`, and `definite_parse_method` handles the
definite form. There are two kinds of per-evaluation state, both keyed on
`eval_toplevel_id()` so they reset with each user command:
- a 32-slot fail memo (`intg_fail_*`), so a fixed-point re-entry does not search the cascade
  a second time;
- one absolute `TimeConstraint` deadline, checked between stages at every depth.

Inexact integrands are rationalised by `common_rationalize_input`, integrated exactly, and
numericalised back at the input's minimum precision. Large `.m` stages (the CRC table,
PolyLogLadder, ParallelMixed*) load on first use.

**Complexity / limits.**
- The rational stage is complete: Hermite reduction plus LRT.
- Transcendental coverage is a decision procedure for the towers that
  `RischTranscendental` models, and a heuristic beyond them.
- The automatic indefinite cascade is bounded by `TimeConstraint` (default 3 s).
  `DerivativeDivides` is measured as nearly all of the cost on integrands that fail. An
  explicit `Method` is unbounded unless a `TimeConstraint` is given.
- Definite stages are narrow, shape-gated recognisers. They either hold by construction (the
  residue theorem, Beta and Mellin strips) or are checked with `NIntegrate` (King,
  series, Feynman). Anything undecidable is left unevaluated.
- `Integrate` is symbolic only. It has no packed, NDArray or `Compile[]` surface; the
  numerical counterpart is `NIntegrate`.
- Multivariate integration is supported only as iterated definite integrals.

- `Protected`, `Listable`.
- Eleven-stage dispatch cascade (`DerivativeDivides`, `LinearRadicals`,
  `QuadraticRadicals` and `LinearRatioRadicals` added 2026-06-06; `Weierstrass`
  added 2026-06-09; `ChebychevAlgebraic` and `GoursatAlgebraic` added 2026-06-29):
  `Integrate[f, x]` (Method -> Automatic, default) tries each subroutine in
  order and returns the first non-`NULL` result:

  **Execution order (as of v0.191):** the `GoursatAlgebraic` stage (item 7
  below) now runs **last**, after `ParallelMixedTower` (item 13) — Goursat is
  the deterministic specialist for pseudo-elliptic `F/R^p`, so every general
  method is given first crack at any integrand it also handles.  To keep that
  reorder cheap, the three *search* stages that otherwise grind for tens of
  seconds on such an integrand before declining — `try_linearity`'s
  product-over-sum split, `DerivativeDivides`'s Eliminate/Solve branch-search,
  and `ParallelMixedTower`'s tower search — are skipped when `f` carries a
  fractional power of a degree-≥3 polynomial in `x` (the pseudo-elliptic shape),
  so it reaches the Goursat stage promptly (e.g.
  `(x-1)/((x+2) Sqrt[x^3-1])` and `(x^4+2x^3-4)/(x^2 Sqrt[(x^2-1)(x^2-4)])` went
  from an effective hang / 45 s to ~0.07 s).  Derivative-divides still runs its
  cheap **direct-quotient** fold on those (so `x^2/(x^3-1)^(1/3)` is unaffected),
  and the explicit `Method -> "ParallelMixedTower"` is not gated.  The item
  numbers below group methods by kind and are unchanged.

  **Symbolic-exponent power times an exponential (as of v0.253).** A dedicated
  recogniser sits just after the Fresnel row, ahead of every general stage:

      Integrate[x^p E^(a x^m), x]  ==  -(1/m) (-a)^(-s) Gamma[s, -a x^m],  s = (p+1)/m

  gated to a **symbolic** exponent `p` (`m` a positive integer, `a` free of `x`) and
  emitted only behind an exact differentiate-back certificate.  It exists because the
  general stages do not merely decline on this shape, they *search*:
  `Integrate[x^n E^(-x), x]` cost a measured 12.9 s to come back unevaluated, and
  `Integrate[x^n E^(-x^2), x]` 12.4 s; both are now ~0.03 s, and
  `DSolve[y'' - y == x^n, y, x]` — whose variation-of-parameters particular is exactly
  that pair of integrals — went from an 8 s abort to a closed form in 0.09 s.  A
  **numeric** exponent is deliberately left to the stages below, which give the
  elementary answer for a non-negative integer `p` and the cleaner `Erf` form
  otherwise.  The certificate is checked twice, since the closed form is stated in the
  principal-branch convention `(-a x^m)^k == (-a)^k x^(m k)`: an exact `Simplify` pass
  first, then one under `PowerExpand`, which quotients out that convention and nothing
  else.

  **`Log` times a trig/hyperbolic kernel (as of v0.343).** A `LogByParts` stage sits
  after the cheap elementary stages (`DerivativeDivides` / `RischTranscendental` /
  `CRCTable`) and before the `ParallelMixedTower` / `ParallelMixedSpecial` tail:

      Integrate[Log[g] K, x]  ==  Log[g] V - Integrate[V (g'/g), x],   V = Integrate[K, x]

  for an integrand `c Log[g(x)] K(x)` with exactly one `Log` factor, `K` free of `Log`
  and carrying a trig/hyperbolic kernel of `x`.  On this family `V` is an elementary trig
  antiderivative and `V (g'/g) = (trig)/x` closes to `SinIntegral`/`CosIntegral`.  It
  exists because a `Log·trig` integrand is outside every cheap stage and the special
  stage *searches*: `Integrate[Log[x] Sin[x], x]` cost ~4.7 s to close (the sibling
  `x^3 Sin[x] Log[x]^2` ~12 s), so `DSolve[y'' + 4y == Log[x], y, x]` — whose VoP
  particular is exactly `Integrate[Log[x] Sin[2x], x]` and `Integrate[Log[x] Cos[2x], x]`
  — took ~9.4 s; both integrals are now ~0.09 s and the ODE 0.35 s, the answer unchanged.
  The cascade's "first non-`NULL` wins" rule means the stage only ever sees integrands the
  cheap stages declined (nothing already-fast is perturbed), and its recursive
  sub-integrals suppress the heavy tail (`g_integrate_no_special`) so a `Log·trig` case
  whose by-parts residual is itself non-elementary declines promptly rather than searching
  `ParallelMixedSpecial` twice.  Emitted only behind an exact `Simplify` differentiate-back
  certificate (`Sinc[z]` from `D[SinIntegral[z]]` rewritten to `Sin[z]/z` first).

  **Execution order (as of v0.238):** `ParallelMixedSpecial` (item 14) is now the
  **last** stage, after Goursat.  It is the only stage that may answer with a
  *non-elementary* function — `ExpIntegralEi`, `Erf`, the elliptic family — and an
  elementary antiderivative is always the better answer, so every stage that can
  produce one (Goursat's pseudo-elliptic reductions included) is given first
  crack.  Being a general search, it also costs least where it is reached only by
  integrands nothing cheaper could close.  Unlike `ParallelMixedTower` it is
  **not** gated off the pseudo-elliptic shape: that gate keeps the
  elementary-only stage away from a genus>0 `F/R^p` curve, on which it grinds and
  closes nothing, and an elliptic pencil is exactly what this stage *does* close
  (`1/Sqrt[x^3-x]` → `-Sqrt[2] EllipticF[ArcSin[Sqrt[2]/Sqrt[1+x]], 1/2]`).  The
  cascade calls only its complete-answer surface, so plain `Integrate[f, x]` can
  never return a partial `answer + Inactive[Integrate][remainder, x]`.
  1. `Integrate\`Undefined[f, x]` — when `f` contains an undefined-function
     derivative (e.g. `f'[x]`); see below.
  2. `Integrate\`BronsteinRational[f, x]` — when `PolynomialQ[f, x] ||
     rationalQ[f, x]`.
  3. `Integrate\`LinearRadicals[f, x]` — rational functions of `x` and radicals
     `(a x + b)^(m/n)` of one shared linear argument; rationalised by
     `u = (a x + b)^(1/n)`.
  4. `Integrate\`QuadraticRadicals[f, x]` — rational functions of `x` and square
     roots `(a x^2 + b x + c)^(m/2)` of one shared quadratic argument;
     rationalised by a single real-valued Euler substitution.
  5. `Integrate\`LinearRatioRadicals[f, x]` — rational functions of `x` and
     radicals `((a x + b)/(c x + d))^(m/n)` of one shared linear-fractional
     argument; rationalised by `u = ((a x + b)/(c x + d))^(1/n)`.
  6. `Integrate\`ChebychevAlgebraic[f, x]` — Chebychev binomial differentials
     `x^p (a x^r + b)^q` (`p, q, r` rational, `a, b` free of `x`).  Elementary
     iff one of `q`, `(p+1)/r`, `q+(p+1)/r` is an integer (Chebychev's theorem),
     with substitutions `x = u^N` (Type I), `u^s = a x^r + b` (Type II), or
     `u = x^r` then `t^s = (a u + b)/u` (Type III) that rationalise `f`.
     Recognition is a single structural scan, so it runs ahead of
     `DerivativeDivides`'s Eliminate/Solve search.  Non-elementary binomials
     return `NULL` (the cascade falls through to later methods).
  7. `Integrate\`GoursatAlgebraic[f, x]` — pseudo-elliptic integrands
     `F(x) R(x)^q` (`F` rational, `R` a polynomial) with `q` any rational of
     reduced denominator `2`, `3`, or `4` by Goursat's algorithm and its
     cube-/fourth-root generalisations (Blake 2026).  The exponent is split
     `R^q = R^k R^(-p)` with radical order `p in {1/2, 1/3, 2/3, 1/4, 3/4}` and
     the integer `R^k` absorbed into `F`, so positive-power radicals such as
     `(1-x^3)^(1/3)/x` are handled, not only radicals already in a denominator.
     A Mobius automorphism cycling the roots of `R` splits the integrand into
     eigencomponents that descend to genus-0 curves when the elementarity
     criterion holds (`p=1/2`: Klein four-group `V4`, trivial projection
     vanishes; `p=1/3,2/3`: order-3 cycle; `p=1/4,3/4`: order-4 cycle on
     harmonic roots).  For `p=1/2` with `R` a cubic carrying the `t^3-1` higher
     symmetry, when `V4` declines a Section-4 (Goursat 1887) period-3 reduction
     is tried: an order-3 Mobius `S` fixes one ramification point and cycles the
     other three, and the integral is elementary when `F` is a non-trivial
     period-3 character `F(S) = Exp[2 Pi I/3] F` (so `(x-1)/((x+2) Sqrt[x^3-1])`
     integrates).  For `p=1/3` with `R` a cubic, when the order-3 eigendescent is
     obstructed (`H1 != 0`) but `F` has a pole at a non-branch point, a
     constructive third-kind logarithmic-part reduction is tried instead: the
     antiderivative is
     `C Sum_{j=0..2} Exp[4 Pi I j/3] Log[R^(1/3) - Exp[2 Pi I j/3] (lead R)^(1/3) x]`,
     so the parametric family
     `(2-(k+1)x)/((1-(k+1)x) (x(1-x)(1-k x))^(1/3))` integrates over `Q(k)`.
     The rational reductions are integrated recursively and
     back-substituted.  Obstructed (genuinely elliptic) integrands,
     non-harmonic quartics, and the cross-character A4 cases (Section 5, e.g.
     `t/((t^3+8) Sqrt[t^3-1])`) return `NULL`.  A differentiate-back guard rejects the
     rare cases where the eigenspace zero-test misfires on deeply nested radical
     roots; the whole attempt runs under a CPU-time budget so a cyclotomic-root
     `R` with an unlucky cofactor (where algebraic-number `Together`/`Cancel`
     blows up) declines rather than hanging the cascade.  Uses
     `Solve[..., Cubics -> True, Quartics -> True]` (the Ferrari quartic solver,
     added 2026-06-29).  The Boolean global `Integrate\`GoursatDebug` (default
     `False`; added 2026-06-30) traces the descent to stderr when set `True`:
     whether the integrand matches the `F(x) R(x)^(-p)` form (and the recognised
     `F`, `R`, `p`), which involution / eigenspace criterion is tested and
     whether it holds (`V4` trivial projection, the order-3/order-4
     ω-eigencomponents, the period-3 trivial projection per fixed point), and the
     differentiate-back verdict.  The flag is latched once at the outermost call,
     so recursive genus-0 reductions share it and indent by depth.
     Before returning (square-root case), the antiderivative is normalised so its
     one radical is the single generator `Sqrt[R]`: any term left carrying `Sqrt[R]`
     *split* across a numerator/denominator as a product of factor roots
     (`Sqrt[x] Sqrt[(1-x)(1-k^2 x)] ...`) has its `x`-dependent half-power factors
     recombined into one radicand and reduced over `R` (to a rational when it is a
     perfect square, else `rational*Sqrt[R]`), removing the spurious branch point and
     keeping a downstream `Simplify` to a single radical generator.  Constant radicals
     are left intact, and the rewrite is kept only if the differentiate-back guard
     still passes.
     A graded battery of worked examples (every exponent `p`, both numerator and
     denominator radicals, and every involution equation, with the negative
     controls that decline) is collected in
     [`GOURSAT_EXERCISES.md`](../../../GOURSAT_EXERCISES.md) and mirrored as the
     `test_graded` ladder in `tests/test_integrate_goursat.c`.
     - **Fresnel stage** (Automatic cascade, after the algebraic-radical
       stages — where Goursat sat before its v0.191 move to the end): a
       Gaussian-phase trig integrand `K Sin[a x^2 + b x + c]` or `K Cos[...]`
       (`a != 0`, `K` free of `x`) closes to `FresnelS`/`FresnelC` by completing
       the square — the trigonometric sibling of the `K E^(a x^2+b x+c) -> Erf`
       recognizer, deterministic and diff-back verified
       (`src/calculus/integrate_fresnel.c`): `Integrate[Sin[x^2],x] =
       Sqrt[Pi/2] FresnelS[Sqrt[2/Pi] x]`, `Integrate[Cos[x^2],x] =
       Sqrt[Pi/2] FresnelC[Sqrt[2/Pi] x]`, and completing-square cases like
       `Sin[x^2+x+1]`.
  8. `Integrate\`Weierstrass[f, x]` — rational functions of the trig kernels
     `Sin/Cos/Tan/Cot/Sec/Csc[x]` (or hyperbolic `Sinh/Cosh/.../Csch[x]`) with a
     kernel in a denominator; continuous `Tan[x/2]` / `Tanh[x/2]` substitution
     (Jeffrey & Rich 1994).  Its entry first applies the **reduction formula** to a
     pure power `c Sec[x]^n` / `c Csc[x]^n` (`n ≥ 3`):
     `∫Secⁿ = Secⁿ⁻² Tan/(n−1) + (n−2)/(n−1) ∫Secⁿ⁻²` (Csc analogue) — the one
     trig-power family that previously spun past any budget in `DerivativeDivides`
     (`Integrate[Sec[x]^4, x] → 2/3 Tan[x] + 1/3 Sec[x]^2 Tan[x]`, instant).  Runs ahead of `DerivativeDivides`: it is
     domain-specific, deterministic, correct by construction, and yields a real,
     continuous antiderivative rather than a complex-logarithm form.
  9. `Integrate\`DerivativeDivides[f, x]` — substitution `u(x)`; in the
     cascade the quiet, branch-correct **direct quotient** strategy only.
  10. *(Removed in v0.163.)* `Integrate\`RischNorman` (Bronstein's parallel-Risch
     pmint heuristic) and `Integrate\`RischNormanBlake` (its generalisation to a
     **simple radical extension** `L = K(y)`, `y^m = q(x)`, after S. Blake,
     *Parallel Integration over Simple Radical Extensions*) were removed;
     `Integrate\`ParallelMixedTower` below subsumes both — the transcendental case
     and the simple radical, in a mixed tower (Blake, *Part II*).
  11. `Integrate\`RischTranscendental[f, x]` — the **recursive** transcendental
     Risch algorithm; adds
     closed forms the earlier stages missed.  Correct by construction (no
     differentiation check).  Handles logarithmic polynomials and the
     special-function cases below (Erf, ExpIntegralEi, LogIntegral, PolyLog).
     Its single-extension cases (fractional Rothstein-Trager log-part, Hermite,
     hyperexponential) are decision procedures over the field `C(x)(t)` for the
     kernel `t = Log[u]` or `E^u`, so each first checks — via `rt_is_ratl_in_xt`
     — that the kernelized integrand is genuinely rational in x and t; an
     integrand with a transcendental coefficient of x (`Sin[x]/Log[x]`,
     `Sin[x]/(1+E^x)`, `Gamma[x]/Log[x]`) lies outside that field and is
     **declined**, never mis-certified (this previously produced a wrong `0`).
     Its Risch differential equation is solved by Bronstein's rational one-step
     (SPDE) reduction (polynomial-gcd time, no undetermined-coefficient
     blow-up), closing high-degree `R(x) e^x` forms and — via the
     exponential special-denominator Laurent ansatz and an exact
     tower-variable verification — nested exp/log towers such as
     `∫((e^x−x²+2x)/(x²(x+e^x)²)) e^((x²−1)/x + 1/(x+e^x)) dx = e^(−x+1/(e^x+x)+(x²−1)/x)`
     and `∫(1/(x log(1+e^x)) − …) dx = log(x)/log(1+e^x)`.  Also closes
     polynomial × `Log` with an irreducible-quadratic argument (a new `ArcTan`
     in the answer, e.g. `∫(x⁵−1) log(x²−x+1) dx`), Q-linearly **dependent**
     logarithms (`∫(log(x/(1+x))+log(1+x))/x = log(x)²/2`), mixed **tangent + log**
     towers (`∫12x⁴+5x⁴ log(x¹²cos x)−x⁵ tan x = x⁵ log(x¹²cos x)`), and deeper
     nested exponentials — including a constant-plus-rational inner exponent
     (`∫… = x e^(x e^(1+1/x))`) and a nested exp inside a circular kernel
     (`∫e^(e^x)(1+x e^x) cos(x e^(e^x)) = sin(x e^(e^x))`).  Arithmetic warnings
     from transient internal singular expressions are muted (as in Mathematica).
     When the integrand carries an inverse-trig factor it is exponentialised to
     complex logs; the real reconstruction renders the arctangent part in the
     clean **one-argument** `ArcTan` (not the two-argument `Arg`/`atan2` form
     `ArcTan[1, u]`), so `∫x^(5/2) ArcTan[√x] = (−6x+3x²−2x³+6 log(1+x)+12 x^(7/2) ArcTan[√x])/42`
     and `∫ArcTan[x] = x ArcTan[x] − ½ log(1+x²)`.
  12. `Integrate\`CRCTable[f, x]` — CRC integral table lookup (lazy-loaded
     from `src/internal/CRCMathTablesIntegrals.m` on first call).  Integer-power
     reduction rules cover all six circular and six hyperbolic functions
     (`Sin/Cos/Tan/Cot/Sec/Csc` and `Sinh/…/Csch`, argument `a x`), their
     products (`Sinh^m Cosh^n`), canonical-head reciprocal products
     (`Csc^m Sec^n`, `Csch^m Sech^n`), mixed tangent/secant powers
     (`Tan^m Sec^n`, `Cot^m Csc^n` and hyperbolic analogues — the canonical form
     of `Sin^m/Cos^n` quotients), `E^(a x)` times circular or hyperbolic
     powers (the hyperbolic case guards the `a = n b` resonance), and
     polynomial × hyperbolic (`x^n Sinh/Cosh`, `x Sinh^m`, `x/Sinh^n = x Csch^n`).
  13. `Integrate\`ParallelMixedTower[f, x]` — the parallel
     (Risch-Norman) integrator over a **simple radical in a mixed transcendental
     tower** (S. Blake, *Parallel Integration over Simple Radical Extensions II:
     Mixed Towers*), a Wolfram-language package
     (`src/internal/mixed/ParallelMixed.m`, lazy-loaded on first use).  Integrands
     built from rational operations, one square root `y^m = q`, pure roots that
     flatten (Lemma 3.2), and `Log`/`Exp`/`Tan` of field elements — the radical may
     sit anywhere in the tower (the "mixed" case), and generators above it may
     differentiate through `y`.  A special denominator prime is split into its
     linear factors over the algebraic closure when a single logand does not
     close the integral (SplitSpecials, Theorem 6.1) — e.g. `Sqrt[Tan[x]]`, whose
     `1 + u^4` (with `u = Sqrt[Tan[x]]`) needs distinct coefficients on its four
     factors over `Q(i, Sqrt[2])`.  Every returned antiderivative passes a hard
     verify-or-decline gate — it is re-differentiated and checked against the
     integrand at real points of its domain, so a mis-realised result declines
     rather than escaping as a wrong answer.  As of v0.192 a nested radical whose
     radicand is a genus-0 conic in the tower generators (e.g. `Sqrt[x + Sqrt[x]]`,
     whose radicand is `Sqrt[x]·(1 + Sqrt[x])`) is kept **fused** as the simple
     radical rather than split into a product of independent radicals
     (`Sqrt[a b] = Sqrt[a] Sqrt[b]` holds only for `a, b > 0`), so the
     antiderivative is faithful at negative and complex arguments too, not merely
     on the integrand's real domain; the branch-unsafe split is still used where it
     is genuinely needed to lower the genus of a higher-degree radicand.  A
     non-elementary case, a
     not-yet-ported case, or a run exceeding an internal wall-clock budget all
     decline cleanly (never a wrong answer, never an unbounded hang).
     Cascade-gated to skip pure rational functions (BronsteinRational's job)
     and, as of v0.191, pseudo-elliptic `F/R^p` integrands (a fractional power
     of a degree-≥3 polynomial), which are the following Goursat stage's
     exclusive domain and on which the tower search only grinds before declining.
  14. `Integrate\`GoursatAlgebraic[f, x]` — runs **last** (moved here in v0.191;
     described in full as item 7 above), the deterministic pseudo-elliptic
     specialist that closes the `F/R^p` integrands the general stages leave.
  If every stage gives up the call bubbles back unevaluated.
- `Method -> "<name>"` option (3rd argument) bypasses the cascade and
  dispatches strictly to a single subroutine, with no fallback:
  - `"Automatic"` — default cascade above.
  - `"BronsteinRational"` — `Integrate\`BronsteinRational[f, x]`.
  - `"DerivativeDivides"` — `Integrate\`DerivativeDivides[f, x]` (direct **and**
    the more thorough Eliminate/Solve branch search).  The list form
    `Method -> {"DerivativeDivides", "Substitution" -> u}` **pins** the kernel
    `u(x)`: instead of collecting and trialing every `x`-dependent subexpression,
    only that one substitution is attempted (still both strategies, still strict
    — no fallback to other kernels if `u` does not close the integral). E.g.
    `Integrate[Sqrt[x]/(1 + Sqrt[x]), x, Method -> {"DerivativeDivides", "Substitution" -> Sqrt[x]}]`.
  - `"LinearRadicals"` — `Integrate\`LinearRadicals[f, x]`.
  - `"QuadraticRadicals"` — `Integrate\`QuadraticRadicals[f, x]`.
  - `"LinearRatioRadicals"` — `Integrate\`LinearRatioRadicals[f, x]`.
  - `"ChebychevAlgebraic"` — `Integrate\`ChebychevAlgebraic[f, x]` (Chebychev
    binomial differential `x^p (a x^r + b)^q`).
  - `"GoursatAlgebraic"` — `Integrate\`GoursatAlgebraic[f, x]` (pseudo-elliptic
    `F/R^p`, `p` in `{1/2, 1/3, 2/3, 1/4, 3/4}`, via Mobius eigendescent).
  - `"Weierstrass"` — `Integrate\`Weierstrass[f, x]` (no denominator gate: applies
    to any rational function of the trig/hyperbolic kernels of `x`, including
    polynomial trig).
    - *Automatic cascade only (no pinnable name):* immediately before Weierstrass,
      a **linear-argument substitution** stage (`src/calculus/integrate_linarg.c`)
      reduces a rational trig/hyperbolic integrand whose kernels all share one
      non-trivial linear argument `a·x + b` (`a` a non-zero number, `a ≠ 1` or
      `b ≠ 0`, with a kernel in a denominator) via `u = a·x + b` to the bare
      integral `(1/a)·(Integrate[f(u), u] /. u -> a·x+b)`. Since Weierstrass always
      substitutes `Tan[x/2]`, a scaled argument would otherwise be multiple-angle
      expanded into a high-degree `Tan[x/2]` rational; this keeps `Integrate[Sec[3x]^2,
      x] = Tan[3x]/3` clean. Declines on bare arguments and on polynomial trig (no
      denominator kernel), so other integrands are untouched.
  - `"RischTranscendental"` — `Integrate\`RischTranscendental[f, x]`, the recursive
    transcendental Risch algorithm (`src/calculus/integrate_risch_transcendental.c`).
    A decision procedure over a differential transcendental tower.  (The parallel-Risch
    heuristics `"RischNorman"` and `"RischNormanBlake"` were removed in v0.163 —
    `"ParallelMixedTower"` subsumes both.)  Every case is correct
    by construction — it fires only behind an exact structural certificate, so
    the result is not checked by differentiation.  Cases:
      - rational: delegated to `Integrate\`BronsteinRational`;
      - logarithmic polynomial `P(x, Log[u])`: the recursive primitive-
        polynomial coefficient matching, with a limited-integration oracle that
        folds a would-be new logarithm back into the tower (e.g.
        `Integrate[Log[2 x + 3], x]`, `Integrate[Log[x]/x, x]`);
      - exponential (Laurent) polynomial `sum_i p_i(x) E^(i u)`, `u` polynomial
        in `x`, `i` positive or negative: the powers of `E^u` decouple and each
        `i != 0` term solves the Risch differential equation
        `q_i' + i u' q_i = p_i` by a polynomial ansatz
        (`Integrate[x E^x, x] = (x-1) E^x`, `Integrate[x E^(x^2), x] = E^(x^2)/2`,
        `Integrate[(E^x+E^(-x))/2, x] = Sinh[x]`);
      - Hermite reduction for a repeated pole of `theta = Log[u]` or
        `theta = E^u` (the latter when `D` is coprime to `theta`):
        `Q = H(theta)/Hden(theta) + sum_j c_j Log(g_j)` with
        `Hden = gcd(D, dD/dtheta)`, solved by `SolveAlways` over `theta` and `x`
        (`Integrate[1/(x (1+Log[x])^2), x] = -1/(1+Log[x])`,
        `Integrate[E^x/(1+E^x)^2, x] = -1/(1+E^x)`);
      - a coupled hyperexponential case (a unified ansatz
        `Q = sum_i w_i(x) E^(i u) + H(E^u)/Hden(E^u) + sum_j c_j Log(g_j)` solved
        by `SolveAlways` over `theta` and `x`) that closes mixed
        polynomial-plus-log exponentials such as
        `Integrate[1/(1 + E^x), x] = x - Log[1 + E^x]`, and — with the Hermite
        term `H/Hden` fused in (the `theta`-coprime denominator split into its
        repeated part `Hden = gcd(Dtil, dDtil/dtheta)` and squarefree radical) —
        the repeated / `theta = 0` exponential poles
        `Integrate[1/(1 + E^x)^2, x] = x + 1/(1 + E^x) - Log[1 + E^x]`,
        `Integrate[1/(E^x (1 + E^x)^2), x]`, `Integrate[1/(1 + E^x)^3, x]`;
      - a multi-kernel **sum-of-exponentials** case: an integrand that
        exponentializes to a sum `sum_k p_k(x) E^(W_k)` of NON-commensurate
        exponentials `E^(W_k)` (e.g. the `(1 ± I) x` pair from `E^x Sin[x]`)
        decouples — each term integrates by its own Risch DE
        `q_k' + W_k' q_k = p_k` — closing `Integrate[E^x Sin[x], x]`,
        `Integrate[x E^x Sin[x], x]`, `Integrate[E^(2x) Cos[3x], x]`, ... (via
        the complex exponentials the answer is left in an I-laden `Cosh`/`Sinh`
        form, a `Simplify` opportunity; the diff-back is exactly `0`);
      - nested **logarithmic** and **exponential tower** cases: a rational
        function of a chain of nested logarithms (`Log[x]`, `Log[Log[x]]`, ...) or
        nested exponentials (`E^x`, `E^(E^x)`, ...) is integrated over the tower
        derivation `D = d/dx + sum_i Dt_i d/dt_i` by one unified `SolveAlways`
        ansatz over all tower variables, closing
        `Integrate[1/(x Log[x] Log[Log[x]]), x] = Log[Log[Log[x]]]`,
        `Integrate[Log[Log[x]]/(x Log[x]), x] = Log[Log[x]]^2/2`,
        `Integrate[E^x E^(E^x), x] = E^(E^x)`, ...  The tower cases are
        bounded-ansatz searches and so are diff-back verified (a non-elementary
        integrand declines rather than returning a spurious form); a whole-tower
        rationality gate and a single-kernel nesting gate keep the other cases
        from ever certifying a wrong nested answer;
      - the **genuine one-extension-at-a-time recursion** (Bronstein ch. 5),
        which the flat tower ansatz above cannot express: it builds
        the ordered differential tower (structure-theorem triangularity check) and
        peels one kernel at a time, integrating the polynomial/Laurent part in the
        top kernel *coefficient by coefficient* — each coefficient integral is an
        integration in the LOWER field (the recursion), bottoming out in `C(x)`.
        This closes **mixed exp/log towers** and **rational lower-field
        coefficients** the single-kind flat cases decline, e.g.
        `Integrate[E^x/x + E^x Log[x], x] = E^x Log[x]`,
        `Integrate[Log[1+E^x] + x E^x/(1+E^x), x] = x Log[1+E^x]`,
        `Integrate[1/(x^2 Log[x]) - Log[Log[x]]/x^2, x] = Log[Log[x]]/x`; diff-back
        verified, so non-elementary mixed integrands (`E^x Log[x]`) decline.  A
        **proper rational part** at a logarithmic top level (a genuine `t_n`-pole)
        is integrated by tower **Hermite reduction + Rothstein–Trager**
        (`H(t)/Hden + sum_j c_j Log(g_j)`, `Hden = gcd(den, d den/dt_n)`, constant
        residues), closing `Integrate[1/(x Log[x] (1+Log[Log[x]])^2), x] =
        -1/(1+Log[Log[x]])` and `Integrate[1/(x(1+Log[x])) + E^x, x] = E^x +
        Log[1+Log[x]]` (a non-constant residue declines).  For an **exponential**
        top level the Laurent and log parts couple, so a proper `t_n`-pole is closed
        by a unified coupled-hyperexponential ansatz over the tower derivation, e.g.
        `Integrate[(2 Log[x]/x) E^(Log[x]^2)/(1+E^(Log[x]^2)), x] =
        Log[1+E^(Log[x]^2)]`.  The exponential Laurent step solves a **field Risch
        differential equation** `q_i' + i w_n' q_i = p_i` in the lower field for each
        power; a `q_i` that is rational there — with an arbitrary denominator, by the
        RDE denominator theorem `q_i = h/Denominator[p_i]` with `h` a bounded
        polynomial — is found, closing
        `Integrate[(2/x - 1/(x Log[x]^2)) E^(Log[x]^2), x] = E^(Log[x]^2)/Log[x]`
        (`q = 1/Log[x]`) and
        `Integrate[(2 Log[x]/(x(1+Log[x])) - 1/(x(1+Log[x])^2)) E^(Log[x]^2), x] =
        E^(Log[x]^2)/(1+Log[x])` (`q = 1/(1+Log[x])`, non-monomial denominator).
        A structural pre-pass re-splits an **evaluator-merged exponential monomial**
        `E^(a+b) -> E^a E^b` (the evaluator folds `E^x E^(E^x)` into `E^(x+E^x)`,
        whose exponent carries the foreign kernel `E^x` and breaks tower
        independence), restoring the independent basis `{E^x, E^(E^x)}` and closing
        `Integrate[E^x E^(E^x)/(1+E^(E^x)), x] = Log[1+E^(E^x)]` (the non-elementary
        `E^(E^x)/(1+E^(E^x))` still declines);
      - a trig/hyperbolic front-end (`TrigToExp` -> exponential machinery ->
        `ExpToTrig`) that closes `Sin`, `Cos`, `Sinh`, `Cosh`, `Sin[x]^2`,
        `Sin[x] Cos[x]`, `Tan`, `Tanh`, ...; through the complex substitution
        `Tan`/`Tanh` come out in a correct but I-laden form (e.g.
        `I x - Log[1 + E^(2 I x)] = -Log[Cos[x]]`) that no current simplifier
        reduces to real closed form (a `Simplify` improvement opportunity).
        A structural pre-pass `rt_powers_to_exp` re-exposes a **transcendental
        general power** `b^e` (base `b != e`, non-rational exponent) as a raw
        base-e kernel `E^(e Log b)` — the evaluator collapses `E^(c Log b) -> b^c`,
        so a trig/inverse-trig *of a logarithm* hides its exponential kernels as
        general powers (`Cos[x Log x] -> (x^(I x) + x^(-I x))/2`, kernels
        `E^(±I x Log x)`) and evades every `Exp`/`Power[E, .]`-keyed recognizer.
        Rewriting without going through the evaluator keeps the raw `Power[E, .]`
        spelling the tower collectors/substitution already match, with `Log b` a
        genuine logarithmic sub-kernel; the exponent is split so an
        evaluator-merged polynomial coefficient (`x^2 x^(I x) -> x^(2+I x)`) peels
        back to an algebraic power times the shared primitive kernel `x^(I x)`.
        This closes the elementary mixed log/exp-power towers
        `Integrate[Cos[Log x], x] = (x/2)(Cos[Log x] + Sin[Log x])`,
        `Integrate[Sin[Log x], x] = (x/2)(Sin[Log x] - Cos[Log x])`,
        `Integrate[x^I, x] = (1/2 - I/2) x^(1+I)`,
        `Integrate[x^x (1 + Log x), x] = x^x`, and
        `Integrate[3 x^2 Cos[x Log x] - x^3 (1 + Log x) Sin[x Log x], x] =
        x^3 Cos[x Log x]`; the non-elementary siblings (`Cos[x Log x]`, `x^x`)
        still decline;
      - `K E^(a x^2 + b x + c)` (`a != 0`) → `Erf`/`Erfi` (Gaussian fast-path);
      - `g E^f` with `f = p/q`, `q = s^2` a perfect square → up to two `Erfi`
        error-function terms (Cherry 1989 §3, completing-square), emitted in the
        same solve as the `ei` terms: `Integrate[(1/x+1/x^2)E^(1/x^2),x] =
        -ExpIntegralEi[1/x^2]/2 - Sqrt[Pi] Erfi[1/x]/2`, and (symbolic parameter)
        `Integrate[E^((x^4+a)/x^2),x] = (Sqrt[Pi]/4)(Erfi[(x^2+Sqrt[a])/x]/E^(2Sqrt[a])
        + E^(2Sqrt[a]) Erfi[(x^2-Sqrt[a])/x])`. A non-square `q` gives no erf term;
      - `g E^f` with `g, f` rational in `x` → `ExpIntegralEi` (Cherry 1989, base
        field) — the general rational exponential integral. The `ei` arguments
        `f + alpha_i` are generated by the resultant `Res_x(g1, p + alpha q)`
        (`f = p/q`, `g1` = the `q`-coprime part of `den(g)`) plus the `q`-side
        term `f` itself, and `y` + the constants solved as one linear system over
        the rational matching identity `g = y' + f' y + sum_i c_i f'/(f+alpha_i)`;
        emits `y E^f + sum_i c_i E^(-alpha_i) ExpIntegralEi[f+alpha_i]`. Closes
        e.g. `Integrate[E^(1/x),x] = x E^(1/x) - ExpIntegralEi[1/x]`,
        `Integrate[E^x/(x+1)^2,x] = ExpIntegralEi[x+1]/E - E^x/(x+1)`,
        `Integrate[(x^2+3)E^x/(x^2+3x+2),x] = E^x - 7 ExpIntegralEi[x+2]/E^2 +
        4 ExpIntegralEi[x+1]/E`, and nonlinear exponents
        `Integrate[E^(x^2)/x,x] = ExpIntegralEi[x^2]/2`. A narrow `(M E^(a x + b))/
        (c x + d)` fast-path handles the linear case in O(1). Higher-multiplicity
        poles close (the elementary part `y` carries the pole to multiplicity
        `m-1`), e.g. `Integrate[E^x/(x+1)^3,x] = ExpIntegralEi[x+1]/(2E) -
        E^x(x+2)/(2(x+1)^2)`. **Real algebraic** `ei`-argument constants are
        admitted (Cherry's algebraically-closed constant field): the resultant
        roots may be irrational reals and the coefficients are solved over the
        extension they generate — `Integrate[E^x/(x^2-2),x] =
        (E^Sqrt[2] ExpIntegralEi[x-Sqrt[2]] - E^(-Sqrt[2]) ExpIntegralEi[x+Sqrt[2]])
        /(2 Sqrt[2])` (Cherry p.894), also `E^x/(x^2-x-1)` (golden ratio) etc.
        A single irreducible-quadratic denominator with **complex** conjugate
        roots closes with a conjugate pair of `ExpIntegralEi` over the number field
        the roots generate — `Q(i)` **and** `Q(i sqrt d)`:
        `Integrate[E^x/(x^2+1),x] = (I/2)(ExpIntegralEi[x+I]/E^I - E^I
        ExpIntegralEi[x-I])`, `E^x/(x^2+x+1)` (over `Q(i sqrt3)`), `E^x/(x^2+3)`,
        `E^x/(x^2+2x+5)`, a numerator `(2x+1)E^x/(x^2+x+3)` (over `Q(i sqrt11)`),
        and the polynomial-part case **d12** `Integrate[(x^2+1)E^x/(x^2+x+1),x] =
        E^x + (complex-conjugate ei pair)`. The lone conjugate pair closes over ANY
        `Q(i sqrt d)` — a shifted center (`E^x/(x^2+2x+3)`, `E^x/(x^2+2x+7)`,
        `(3x+1)E^x/(x^2+2x+3)`) that the direct number-field `Solve` cannot crack is
        solved over Q by the `{1, chs}`-basis fallback. Complex poles mixed with a
        `P2`/reciprocal term (`E^(1/x)/(x^2+1)`) and degree-`≥3` constants
        (`E^x/(x^3-2)`) are deferred and decline cleanly. A constant exponent offset
        `E^(c + h(x))` (e.g. `E^(1/x+2) = E^2 E^(1/x)`) is factored out before matching.
        A **multi-term**
        exponential `Sum_i p_i E^(i w)` (rational in a single kernel `E^w`, several
        commensurate Laurent terms) that the single-shape matcher cannot peel is
        integrated term-by-term (Cherry Thm 5.4 case b) —
        `Integrate[(E^x+E^(2x))/(x-1),x] = E ExpIntegralEi[x-1] +
        E^2 ExpIntegralEi[2x-2]`, `E^x(E^x+1)/((x-1)(x-2))`;
      - **error-function integration of transcendental Liouvillian integrands**
        (Knowles 1992/93, the layer above Cherry): an integrand that may itself
        contain `Erf`/`Erfi` is closed in terms of elementary functions plus
        `Erf`/`Erfi` over the `RT_PRIM` Liouvillian-primitive tower. For each
        exponential monomial `e^w` the perfect-square gate (`-w = u^2 -> erf(u)`,
        `+w = u^2 -> erfi(u)`) generates the error-function arguments, and
        `v + Sum_i k_i erf(u_i)` is solved as one linear system over the constants
        (diff-back verified): `Integrate[E^(-x^2-Erf[x]^2),x] = (Pi/4) Erf[Erf[x]]`,
        `Integrate[Erf[x] E^(-x^2-Erf[x]^2),x] = -(Sqrt[Pi]/4) E^(-Erf[x]^2)`,
        `Integrate[2 E^(-x^2) Erf[x] - 3 E^(-1/x^2)/x^2,x] = (Sqrt[Pi]/2) Erf[x]^2 +
        (3 Sqrt[Pi]/2) Erf[1/x]`; `Integrate[x E^(-x^2-Erf[x]^2),x]` correctly
        declines (no erf-elementary antiderivative). The **radical (quasiquadratic)**
        case is also handled (Knowles Part I §6, Ex 8.1): when a top exponential
        hides an algebraic factor `E^(l/2 Log[g]) = g^(l/2)`, it is pulled out first,
        surfacing a half-integer power of a log tower variable; the erf argument is
        then a radical `1/Sqrt[Log g]`, solved in `s = Sqrt[Log g]` — e.g.
        `Integrate[E^(1/2 Log[Log[x]] - 1/Log[x])/(x Log[x]^2),x] =
        -Sqrt[Pi] Erf[1/Sqrt[Log[x]]]` (and its `Erfi` dual). Both **spellings**
        of that integrand reach the engine: the already-reduced
        `Integrate[E^(-1/Log[x])/(x Log[x]^(3/2)),x]` carries the half-integer
        power openly, and the `RischTranscendental` scope gate — which refuses
        every algebraic function of `x`, rightly, since the recursive Risch
        algorithm is not a decision procedure over an algebraic extension —
        admits it as a narrow exception: when *every* algebraic site is a
        fractional power of a transcendental **kernel**, each is rewritten
        `g^(p/q) = g^n E^(r Log g)` (the inverse of the collapse above) and the
        answer is returned only behind a diff-back check, because that rewritten
        tower hides the algebraic relation `(E^(r Log g))^q = g^p`. A genuine
        radical of `x` (`Sqrt[x]`, `Sqrt[1+x^2]`, `Cos[Sqrt[x]]`) is still
        refused and left to the algebraic routes, and the non-elementary
        decision half never runs on the admitted path. `x`-rational
        elementary-part coefficients remain a later increment;
      - `c w^(p-1) w'/Log[w]` → `c LogIntegral[w^p]` (single-li fast path) —
        subsumes `K/Log[x] → K LogIntegral[x]`, a scaled/affine argument
        (`1/Log[2x] → LogIntegral[2x]/2`) and a monomial numerator (`x/Log[x] →
        LogIntegral[x^2]`);
      - general single-log tower `gamma` over `C(x, Log[w])` (w a polynomial) →
        `v(x,Log[w]) + Sum_k d_k LogIntegral[w^k]` (Cherry 1986, MULTI-li via the
        degree-1 Sigma-decomposition over the generator `w`, matched in the tower):
        `Integrate[x^2/Log[x+1],x] = LogIntegral[(1+x)^3] - 2 LogIntegral[(1+x)^2]
        + LogIntegral[1+x]`, `Integrate[x/Log[x]^2,x] = 2 LogIntegral[x^2] -
        x^2/Log[x]`, `Integrate[x^3/Log[x^2-1],x] = LogIntegral[x^2-1]/2 +
        LogIntegral[(x^2-1)^2]/2`. Includes the **transcendental-constant
        rescaling** (a constant root `rho` of the `Log`-denominator →
        `e^rho LogIntegral[e^(-rho) w]`): `Integrate[1/(Log[x]+3),x] =
        LogIntegral[E^3 x]/E^3`, `Integrate[(Log[x]^2+3)/(Log[x]^2+3Log[x]+2),x] =
        x + 4 LogIntegral[E x]/E - 7 LogIntegral[E^2 x]/E^2`. Multi-log towers
        (reducible `w` needing a product decomposition) decline. The tower
        Sigma-decomposition is **Laurent**: a negative-degree (Laurent-polynomial)
        numerator over `Log[w]` admits negative li powers `w^k`, emitted as
        `ExpIntegralEi[k Log[w]]` (`= LogIntegral[w^k]`, the Wolfram form, which
        also avoids `LogIntegral` of an argument `<1`), and higher `Log` poles
        carry an elementary Laurent part `Sum_i b_i(x) Log[w]^i` with `b_i`
        themselves Laurent in `x`: `Integrate[1/(x^4 Log[x]),x] =
        ExpIntegralEi[-3 Log[x]]`, `Integrate[1/(x^4 Log[x]^2),x] =
        -3 ExpIntegralEi[-3 Log[x]] - 1/(x^3 Log[x])`,
        `Integrate[(x-1)^2/(x^4 Log[x]^3),x]` (three `ExpIntegralEi` terms plus a
        Laurent elementary part);
      - `c g'(x)/Log[g]` → `c LogIntegral[g]` for a **mixed exp/log kernel** `g`
        that is not a polynomial in `x` (so `Log[g]` appears as a SUM, not a
        `Log[...]` node, and the polynomial-`w` engines above cannot see it). The
        kernel is discovered from the denominator: for `f = N/D` (`Together`), the
        only elementary `g` with `f ∝ g'/Log[g]` has `Log[g] = D`, i.e. `g =
        Exp[D]`; the certificate `Together[f D / g']` must be a nonzero constant,
        and a `PowerExpand` diff-back reconfirms the branch identity. Closes
        `Integrate[(E^x (1+x))/(x + Log[x]), x] = LogIntegral[x E^x]` (here
        `Log[x E^x] = x + Log[x]`). Gated to a `D` containing a logarithm — a
        `Log`-free `D` is either `Ei`-reducible (`LogIntegral[Exp[u]] =
        ExpIntegralEi[u]`, owned by the `Ei` recognizer) or spurious, so it
        declines fast without normalizing over an `Exp[...]` generator;
      - `K Log[1 + p x]/x` → `PolyLog[2, -p x]` (dilogarithm fast path); and the
        general `R(x) Log[w]` → `Log-Log + PolyLog[2]` form (Cherry degree-2
        Sigma-decomposition, matched in the log tower with dilog arguments the
        interpolants `(x-r_i)/(r_j-r_i)` between the rational roots of the linear
        factors): `Integrate[Log[x]/(1+x),x] = Log[x] Log[1+x] + PolyLog[2,-x]`,
        `Integrate[Log[x]/(1-x),x] = PolyLog[2,1-x]`, `Log[x]/(x^2-1)`,
        `Log[2x+1]/(x+1)`. **Transcendental-constant root spacing** (spacing `!= 1`)
        now closes too — a dilog whose derivative leaves a `Log` of a positive
        constant contributes a real `Log-Log` term:
        `Integrate[Log[2+x]/x,x] = Log[2] Log[x] - PolyLog[2,-x/2]`,
        `Log[2x+3]/(x-1)`. Non-monic linear kernels (`Log[3+2x]/x`) and products
        of logs still decline;
      - **exponential-tower dilogarithm** (Cherry, the exp-tower mirror of the
        `R(x) Log[w]` case above): an integrand that is a rational-in-`E^(c x)`
        combination of the weight-1 logs `{x} ∪ {Log[E^(cx)-ρ_k]}` → `PolyLog[2,
        exp]`. With `θ=E^(cx)` and `x=Log[θ]/c` the "root at 0" tower-log, the
        answer is matched natively in the tower (an ansatz of `PolyLog[2, Möbius(θ)]`,
        `x Log[θ-ρ_k]`, `Log-Log`, and `x^2` terms; coefficients from a linear solve
        over `{θ, x, Log[θ-ρ_k]}`; PowerExpand diff-back verified). Covers both the
        rational-times-`x` forms and the outer-log forms:
        `Integrate[x/(E^x-1),x] = x Log[1-E^-x] - PolyLog[2, E^-x]`,
        `Integrate[x/(1+E^x),x]`, `Integrate[x E^x/(E^x-1),x]`,
        `Integrate[x/(E^(2x)-1),x]` (rate `c=2`), `Integrate[Log[1+E^x],x] =
        -PolyLog[2,-E^x]`, `Integrate[Log[1+E^-x],x]`, and mixtures. Weight-`≥2`
        integrands (`x Log[1+E^x]`, `x^2/(E^x-1)` → `PolyLog[3]`), incommensurate
        exponentials, and algebraic roots decline cleanly. The engine is reached
        automatically from `Integrate` and directly via the debug surface
        `Integrate\`Cherry\`DilogExp[f, x]`;
      - **exponential-tower polylogarithm ladder** (the general-weight,
        algebraic-root generalisation of the exp-tower dilogarithm): an integrand
        `P(x)/Q(E^(c x))` — rational in `θ=E^(c x)`, polynomial in `x` — is
        partial-fractioned over the (rational OR algebraic) roots `ρ` of `Q`, and
        each simple pole `x^n/(θ-ρ)` is closed by the exact Cherry ladder
        `INT x^n/(θ-ρ) = Σ_{k=0}^{n} -(1/ρ)(n!/(n-k)!)/c^{k+1} x^{n-k}
        PolyLog[k+1, ρ/θ]` (from `d/dx PolyLog[k, ρ/θ] = -c PolyLog[k-1, ρ/θ]`,
        `PolyLog[1,z] = -Log[1-z]`), so an `x^n` numerator yields polylogarithms up
        to weight `n+1`: `Integrate[x^2/(E^x-1),x] = x^2 Log[1-E^-x] -
        2 x PolyLog[2,E^-x] - 2 PolyLog[3,E^-x]`,
        `Integrate[x^4/(E^(5x)-1),x] = ... - (24/3125) PolyLog[5, E^-5x]`
        (rate `c=5`), and the **algebraic-root** family
        `Integrate[x/(E^(2x)+E^x-1),x]`, `Integrate[x^2/(E^(2x)+E^x-1),x]` over
        `Q(Sqrt[5])`. Correct by construction (telescoping IBP), PowerExpand-diff-
        back verified over the constant field. Scope: single commensurable kernel,
        `Q` a proper fraction in `θ` with simple nonzero roots; improper fractions,
        `θ=0` (Laurent) poles, and repeated roots decline. Reached automatically
        from `Integrate` (ahead of the weight-2 dilog form) and directly via
        `Integrate\`Cherry\`PolyLogExp[f, x]`;
      - fractional (Rothstein–Trager) log-part: a proper rational function of
        `theta` with squarefree denominator `prod g_i` gives `sum_i c_i Log(g_i)`,
        the constant residues `c_i` solved from `num = sum_i c_i D(g_i)(d/g_i)`
        via `SolveAlways` over `theta` and `x` (`Integrate[1/(x(1+Log[x])),x] =
        Log[1+Log[x]]`, `Integrate[E^x/(1+E^x),x] = Log[1+E^x]`); the residues are
        verified constant (an `x`-dependent `SolveAlways` pseudo-solution is
        rejected so `1/Log[c x+d]` is never mis-closed with a polynomial-coefficient
        logarithm);
      - **pure resultant Lazard–Rioboo–Trager** log-part (when the `SolveAlways`
        path above declines): an irreducible-over-Q denominator factor in `theta`
        with **algebraic (non-rational) residues** is closed by the exact resultant
        `Res_t(a - z D(d), d)` (residues) + Rioboo `LogToReal` (real `Log + ArcTan`
        form), delegated to the internal `Integrate`TranscendentalLogPart` (reuses
        the rational-LRT machinery with the monomial derivation), closing
        `Integrate[1/(x(Log[x]^2+1)),x] = ArcTan[Log[x]]`,
        `Integrate[E^x/(E^(2x)+1),x] = ArcTan[E^x]`, mixed
        `Integrate[(1+Log[x])/(x(Log[x]^2+1)),x] = ArcTan[Log[x]] + Log[1+Log[x]^2]/2`
        and higher-degree denominators; diff-back verified, so non-elementary
        siblings (`1/(Log[x]^2+1)`, `1/(x^2(Log[x]^2+1))`) decline; this same
        resultant LRT is **lifted into the tower recursion** at a logarithmic top,
        so a nested-log proper part with algebraic residues closes too —
        `Integrate[1/(x Log[x] (Log[Log[x]]^2+1)),x] = ArcTan[Log[Log[x]]]` — with
        the residues gated free of *every* lower-field variable (`{x, Log[x], …}`),
        and `1/(Log[x](Log[Log[x]]^2+1))` (x-dependent residues) declining;
    Rational-argument exponents such as `E^(1/x)` are handled (`Integrate[-E^(1/x)/
    x^2, x] = E^(1/x)`) via the `q = h/Denominator[p]` RDE ansatz.  **Multiplicatively
    commensurate merged kernels** are reduced in the tower builder — a collected
    exponential whose exponent is an integer multiple of a class primitive
    (`E^(2 E^x) = (E^(E^x))^2`) becomes a power of the primitive's tower variable
    instead of a spurious extra extension — so towers with `E^(k u)` close and the
    exp-top algebraic-residue LRT is unblocked
    (`Integrate[E^x E^(E^x)/(1+E^(2 E^x)), x] = ArcTan[E^(E^x)]`,
    `Integrate[E^x E^(2 E^x)/(1+E^(E^x)), x] = E^(E^x) - Log[1+E^(E^x)]`).  The class
    primitive is **synthesized** (`p = w_0/lcm(ratio denominators)`) rather than
    required to be a member, so *non-integer* commensurate exponents close too:
    `Integrate[1/(E^(x/2)+E^(x/3)), x]` (primitive `E^(x/6)`) and the nested
    `Integrate[D[E^(E^x/2)/(1+E^(2 E^x/3)), x], x]` (primitive `E^(E^x/6)`); only a
    class whose members are not all rational multiples of one another declines.
    The **RDE-solver degree bounds** are
    Bronstein's exact `RdeBoundDegree` (leading-degree balance,
    `deg_v(q) = deg_v(p) − deg_v(f)` where the exponential dominates), with **no
    arbitrary cap**, so an exponential-Laurent coefficient of any degree closes
    (`Integrate[(6 Log[x]^5 + 2 Log[x]^7)/x E^(Log[x]^2), x] = Log[x]^6 E^(Log[x]^2)`,
    and deg-20 …).  The **flat-tower and proper-part Hermite ansätze are likewise
    cap-free**: exact top-kernel log/exp Laurent bounds and derived inner-exp windows,
    so all top degrees close (`Integrate[Log[Log[x]]^5/(x Log[x]), x] =
    Log[Log[x]]^6/6`, `Integrate[E^x E^(6 E^x)/(1+E^(E^x)), x]`).  The
    **leading-coefficient cancellation / resonance** sub-case of `RdeBoundDegree`
    completes the SPDE degree machinery — the bound is widened monotonically to
    `max(naive, m_res)` at the Bronstein resonance integer `m_res` (detected live for the
    exponential top), correct by the same certification-and-diff-back gate.  Only
    algebraic extensions (`Sqrt`, `RootSum`) remain unimplemented, so integrands needing
    them return unevaluated.
  - `"CRCTable"` — `Integrate\`CRCTable[f, x]`.
  - `"MeijerConvolution"` (definite, (0, Infinity); v0.460): `C x^alpha K1(x) K2(x)` with two kernels from `BesselJ`, `BesselY`, `BesselK`, `E^(-p x)`, `E^(-p x^2)`, `Cos`, `Sin` and `(a^2 + x^2)^-r` at arbitrary scales (`src/internal/MeijerConvolution.m`).
    - Each kernel is written as a Meijer G-function of x^2; exp, cos and sin of x get there through the duplication formula.
    - The x^alpha weight is absorbed into the first G, and the convolution theorem (Prudnikov-Brychkov-Marichev 2.24.1.1) gives a single MeijerG, which Mathilda's MeijerG evaluator reduces where it can (pFq, Bessel, Gamma).
    - Convergence at 0 and Infinity is checked from the kernels' asymptotics (with an equal-frequency rule for oscillatory products). An unproved condition is returned as a `ConditionalExpression`; a refuted one declines. No numeric quadrature is used.
    - Example: `Integrate[BesselJ[nu, a x] E^(-p x), {x, 0, Infinity}, Assumptions -> a > 0 && p > 0 && nu > 0]` closes.
  - `"PolyLogLadder"` — the polylogarithm-ladder stage (`src/internal/PolyLogLadder.m`, v0.458). It covers two classes, both exact constructions accepted only after a differentiation check:
    - `R(x) Log[a x + b]^k`, where R is rational with poles anywhere (algebraic or complex roots of any multiplicity, denominator factors of degree <= 4), k is any power, there is one linear log kernel per term, and parameters may be symbolic. The partial fractions are taken in the kernel variable W = a x + b, so there is no splitting of a logarithm and no branch bookkeeping. The ladder is `INT Log[W]^k/(W - c) = Log[W]^k Log[1 - W/c] + k Sum (-1)^t (k-1)!/(k-1-t)! Log[W]^(k-1-t) PolyLog[t+2, W/c]`.
    - A polynomial in x times a rational function of trig/hyperbolic functions of linear arguments with one base frequency. TrigToExp writes it in the single kernel T = E^(base x), partial fractions are taken in T, and the exponential ladder `INT x^n/(T - rho) = Sum -(1/rho) n!/(n-k)!/base^(k+1) x^(n-k) PolyLog[k+1, rho/T]` applies; repeated roots use rho-derivatives of it.
    - `P(x) Log[G]`, where P is a polynomial and G a rational function of trig/hyperbolic functions of linear arguments (v0.459). With T = E^(base x) and G = K T^e Prod (T - rho)^m, Log[G] = D + e base x + Sum m Log[1 - T/rho]. Here D = Log[K Prod (-rho)^m] + 2 Pi I n is an exact piecewise constant whose integer n is read at a sample point. Each `x^n Log[1 - T/rho]` is the exponential ladder. Examples: `Integrate[Log[Sin[x]], x] = (Log[1/2] + I Pi/2) x - I x^2/2 + I PolyLog[2, E^(I x)] + I PolyLog[2, -E^(I x)]`, and `Integrate[x Log[Sin[x]], {x, 0, Pi}] = -Pi^2 Log[2]/2`. Newton-Leibniz always cross-checks continuity for an antiderivative with PolyLog/Log of an exponential, because a branch jump of D inside the interval is caught there.
    - Examples: `Integrate[Log[x]/(1 + x^2), x]` and `Integrate[x^2 Csc[x], x]` close. `Integrate[Log[x]/(1 + x^2), {x, 0, 1}]` is `-Catalan` and `Integrate[x Csc[x], {x, 0, Pi/2}]` is `2 Catalan`.
    - In the automatic cascade it runs after `CRCTable`, before `LogByParts`.
  - `"ParallelMixedTower"` — `Integrate\`ParallelMixedTower[f, x]` (parallel
    Risch-Norman over a simple radical in a mixed transcendental tower).
  - `"ParallelMixedSpecial"` — `Integrate\`ParallelMixedSpecial[f, x]` (the same
    over a mixed tower with **special-function kernels**, so the answer may carry
    `ExpIntegralEi`, `LogIntegral`, `Erf`/`Erfc`/`Erfi`, the incomplete
    `Gamma[s, z]`, `SinIntegral`/`CosIntegral`, or `EllipticF`/`EllipticE`/`EllipticPi`).  The worker answers with the pair
    `{answer, verified}`; a *verified* pair is unwrapped to the bare
    antiderivative and anything else is handed back as the package said it, so
    `Head[r] === List` remains the caller's test for "no answer".

    **Measured on the 312-case stress corpus of the paper** (120 s cap, against
    the three reference ports, which all score 305 PASS / 7 HONEST):
    Mathilda **247 PASS, 51 HONEST, 11 WEAK, 3 FAIL**, with **no false
    certificate**. Seven groups are at parity. 39 of the 51 HONEST (46 in total, 7 of them HONEST in the reference too) are strict
    mode withholding an answer that is already correct and identical to
    Mathematica's — `Exp[-x^2]` → `-Sqrt[Pi] Erfc[x]/2`, `Exp[-x^3]` →
    `-Gamma[1/3, x^3]/3` — for want of a non-elementarity certificate, and all 11
    WEAK are the partial mode leaving exactly those sub-integrands in the
    remainder for the same reason. The cause is Part II's `q =!= None` guard on
    the holomorphic-remainder certificate (the T2/T10 false-certificate
    protection, absent from the research reference), so it is one root cause for
    50 of the 58 differing cases and it predates this stage. On the 247 cases it does
    close, Mathilda is the fastest of the four ports (0.070 s median against
    Mathematica's 0.110 s).
  - `"Undefined"` — `Integrate\`Undefined[f, x]`.
  - `"Symmetry"` — origin-symmetry reduction for an interval `[-c, c]`
    (`Integrate\`Symmetry[f, {x, -c, c}]`): an odd integrand integrates to `0`,
    an even one to `2 Integrate[f, {x, 0, c}]`. The parity is proved by
    `Simplify`, and a value is claimed only when the half integral converges, so
    a divergent principal value is never reported as `0`. Under Automatic it runs
    after residue and before Newton-Leibniz.
  - `"Beta"` — Euler-Beta reduction on `[0,1]`
    (`Integrate\`Beta[f, {x, 0, 1}]`): `x^(k-1) (1-x)^(l-1) → Beta[k, l]`, with
    `Log[x]^i Log[1-x]^j` weights giving the mixed parameter derivative of
    `Beta`. Gated on `Re[k] > 0 && Re[l] > 0`.
  - `"TrigPower"` — `Sin[x]^m Cos[x]^n` over a canonical trig interval
    (`Integrate\`TrigPower[f, {x, 0, c}]`): over `[0, Pi/2]` it is
    `Beta[(m+1)/2, (n+1)/2]/2`; over `[0, Pi]`/`[0, 2Pi]` the standard parity
    multipliers apply (an odd power integrates to `0`).
    Under Automatic, `Beta`, `TrigPower` and the Mellin/Ramanujan method
    (below) — all deterministic and correct-by-construction, gated by a cheap
    interval/shape test — run **before** Newton-Leibniz (confidence-first
    ordering, after residue/symmetry): a symbolic-exponent integrand they own
    (`x^(s-1) f`, `Sin^a Cos^b`) closes at once, instead of first driving FTC
    into a doomed elementary-antiderivative search. FTC still owns the
    integer-power and genuinely-elementary cases they decline on.
  - `"NewtonLeibniz"` — the real-axis definite-integral mechanism (implicit for
    the `{x, a, b}` form); see **Definite integration** below.
  - `"LineIntegral"` — the complex contour mechanism (implicit for the
    `{x, z0, …, zn}` form); see **Complex line integration** below.
  - `"DiffUnderInt"` (alias `"DifferentiationUnderIntegral"`) — definite
    integration by differentiation under the integral sign (Feynman's trick);
    `Integrate\`DiffUnderInt[f, {x, a, b}]`. Tried last in the definite cascade
    (after residue and Newton-Leibniz). See **Differentiation under the integral
    sign** below.
  - `"SinPowerMonomial"` — `Sin[r x]^k / x^m` on `[0, Infinity)` (the ssp
    family); `Integrate\`SinPowerMonomial[f, {x, 0, Infinity}]`.
  - `"OscillatoryPower"` — Fresnel-type `Cos[b x^n]` / `Sin[b x^n]` on
    `[0, Infinity)`; `Integrate\`OscillatoryPower[f, {x, 0, Infinity}]`.
  - `"RationalLog"` — `R(x) Log[x]^n` on `[0, Infinity)` for a proper rational
    `R` with negative-real-axis poles; `Integrate\`RationalLog[f, {x, 0, Infinity}]`.
  - `"RamanujanMasterTheorem"` (alias `"Mellin"`) — half-line `∫₀^∞ x^{s-1} f(x) dx`
    by the Mellin-transform / Ramanujan Master Theorem method;
    `Integrate\`RamanujanMasterTheorem[f, {x, 0, Infinity}]`. Under Automatic it
    runs before Newton-Leibniz (and before DiffUnderInt). See **Mellin / Ramanujan
    Master Theorem** below.
  - `"IntegralRepresentation"` — a half-line integrand that IS a classical
    integral representation of a special function;
    `Integrate\`IntegralRepresentation[f, {x, 0, Infinity}]`. Families:
    `E^(-p x) BesselJ[nu, q x]` → Laplace-Bessel
    `(Sqrt[q^2+p^2]-p)^nu/(q^nu Sqrt[q^2+p^2])` (`p>0`);
    `BesselJ[0, a x] {Cos,Sin}[b x]` → the Weber–Schafheitlin discontinuous
    integral, a branch chosen by the proved ordering — `Cos`: `1/Sqrt[a^2-b^2]`
    (`a>b>0`) or `0` (`b>a>0`); `Sin`: `0` (`a>b>0`) or `1/Sqrt[b^2-a^2]`
    (`b>a>0`); undecided / `a==b` (divergent) declines;
    `E^(-A Cosh[x]) Cosh[n x]` → `BesselK[n, A]` (`A>0`);
    `x^(nu-1) E^(-A x - B/x)` → `2 (B/A)^(nu/2) BesselK[nu, 2 Sqrt[A B]]`
    (`A,B>0`; `nu` a half-integer gives an elementary form, e.g.
    `E^(-A x - B/x)/Sqrt[x]` → `Sqrt[Pi/A] e^(-2 Sqrt[A B])`);
    `Cos[p x^3 + q x]` → `Pi (3p)^(-1/3) AiryAi[q (3p)^(-1/3)]` (`p>0`, `q` real);
    `x^(s-1) E^(-a x)/(1 - z E^(-c x))` → `c^(-s) Gamma[s] LerchPhi[z, s, a/c]`, and
    `Gamma[s] HurwitzZeta[s, a/c]` when `z = 1` (the Lerch/Hurwitz representation;
    `a,c>0` gated separately so the scaled `a/c>0` discharges, `Re s>1` when `z=1`
    else `Re s>0`) — e.g. `E^(-a x) x^(s-1)/(1-E^(-x))` → `Gamma[s] HurwitzZeta[s, a]`.
    Each gate is proved from Assumptions; correct-by-construction, no NIntegrate.
    Under Automatic it runs after Ramanujan, before Newton-Leibniz.  (The Euler
    finite-interval generalization `x^(a-1)(1-x)^(b-1)(alpha+beta x)^e` on `[0,1]`
    → `alpha^e Beta[a,b] Hypergeometric2F1[-e, a, a+b, -beta/alpha]` runs under the
    `"Beta"` mechanism — e.g. `x^(a-1)(1-x)^(b-1)/(x+c)^(a+b)` on `[0,1]`.)
  - `"Piecewise"` — a piecewise integrand built from `Abs`, `Sign`,
    `HeavisideTheta`, `UnitStep`, `Min`, `Max`, or `Piecewise` whose arguments
    depend on `x` (`Integrate\`Piecewise[f, {x, a, b}]`). It collects the
    breakpoints where a branch switches (root of the `Abs`/`Sign`/step argument,
    `u == v` for `Min`/`Max`, a `Piecewise` condition boundary), keeps the ones
    strictly interior to `(a, b)` — decided assumption-aware with `Refine` — and
    sorts them; on each sub-interval the integrand reduces to its single active
    branch (chosen by a *proved* sign of each governing argument at an interior
    test point), integrates that with the ordinary cascade, and sums. A branch
    is taken only when its governing inequality is proved under the assumptions,
    so an undecidable breakpoint position or branch sign leaves the integral
    unevaluated. Under Automatic it runs early — right after the `DiracDelta`
    sifting, ahead of residue/FTC (none of which model a branch switch) — and
    declines cheaply (no branch-defining head) for every other integrand. E.g.
    `Integrate[Abs[x - a], {x, 0, 2a}, Assumptions -> a > 0] -> a^2`,
    `Integrate[Min[x, a], {x, 0, 2a}, Assumptions -> a > 0] -> 3 a^2/2`,
    `Integrate[HeavisideTheta[a - x], {x, 0, Infinity}, Assumptions -> a > 0] -> a`.
    Two refinements extend its reach: a **perfect-square radical** is first
    normalized to its `Abs` form under a real integration variable
    (`Sqrt[1 - Cos[2x]] -> Sqrt[2] Abs[Sin x]`, so `∫₀^{2π} Sqrt[1-Cos[2x]] dx =
    4 Sqrt[2]`); and a breakpoint `Solve` returns as a **periodic family**
    (`ConditionalExpression[Pi + 2 C[1] Pi, …]`) is enumerated to its concrete
    members inside the interval, so `∫₀^{2π} Abs[Sin x] dx = 4` splits at the
    interior `x = π`. A sub-interval proved to lie in `x <= 0` is integrated under
    the reflection `x -> -x`, keeping a resolved `Abs` branch real
    (`∫_{-1}^{1} Sign[x] Abs[x]^a dx = 0`, not a complex `(-1)^(1+a)` artifact).
  - `"SeriesExpansion"` (alias `"SeriesTermByTerm"`) — a log-weighted rational
    integral `R(x) Log[x]^k` (`k >= 1`) on `[0,1]` with a root-of-unity
    denominator `1 ± x^d` (`Integrate\`SeriesExpansion[f, {x, 0, 1}]`). Expands the
    geometric series `1/(1 ∓ x^d) = Σ (±1)^n x^{dn}`, integrates term-by-term with
    the exact identity `∫₀¹ x^m Log[x]^k dx = (-1)^k k!/(m+1)^(k+1)`, and resums
    each residue class with `Sum[]` (→ `Catalan`, `DirichletBeta`, `Zeta`). The
    general coefficient is produced by exact cyclotomic divisibility, not inferred
    from samples, and the value is cross-checked against `NIntegrate`. It is the
    alternating sibling of the `"Beta"`-family `integrate_ratlog01_try` (which owns
    the `1 - x^m` HurwitzZeta case and runs first); this closes the `1 + x^d`
    alternating case. E.g. `Integrate[Log[x]/(1+x^2), {x,0,1}] = -Catalan`,
    `Integrate[Log[x]^2/(1+x^2), {x,0,1}] = Pi^3/16`,
    `Integrate[Log[x]^2/(1+x), {x,0,1}] = (3/2) Zeta[3]`. Under Automatic it runs in
    the pre-FTC `[0,1]` block (these integrands have non-elementary PolyLog
    antiderivatives).
  - `"King"` (alias `"Reflection"`) — the general King reflection `x -> a+b-x` on a
    finite interval (`Integrate\`King[f, {x, a, b}]`):
    `Integrate[f,{x,a,b}] = (1/2) Integrate[f(x)+f(a+b-x),{x,a,b}]`, and when the
    reflected sum `f(x)+f(a+b-x)` is a constant `c` the value is `c (b-a)/2`. Closes
    e.g. `Integrate[1/(1+Tan[x]^Sqrt[2]), {x,0,Pi/2}] = Pi/4` — and the whole
    `1/(1+Tan[x]^s)` family for any real `s` — which has no elementary
    antiderivative. The constant is read at the midpoint, screened numerically, and
    cross-checked against `NIntegrate` (a mis-detected constancy fails the check and
    declines, so no wrong answer). It is numeric-only: a parametric integrand is out
    of scope. Complements the `x·g(x)`, `g`-invariant reflection already in the
    `"Beta"` file. Under Automatic it runs before Newton-Leibniz.
  - `"FeynmanParameter"` — the parameter / Frullani-logarithm family on `[0,1]`
    (`Integrate\`FeynmanParameter[f, {x, 0, 1}]`), i.e. differentiation under the
    integral sign w.r.t. the exponent:
    `Integrate[(Σ cᵢ x^{eᵢ})/Log[x], {x,0,1}] = Σ cᵢ Log[eᵢ+1]`, valid when
    `Σ cᵢ = 0` (so the `1/Log[x]` singularity at `x=1` is integrable) and each
    `eᵢ+1 > 0` (proved from Assumptions). The inner `x`-power integral is taken
    directly as `1/(e+1)`. E.g. `(x^a-1)/Log[x] → Log[a+1]` (`a>-1`),
    `(x^a-x^b)/Log[x] → Log[(a+1)/(b+1)]`. Numerically cross-checked against
    `NIntegrate`; `1/Log[x]` (coefficient sum `≠ 0`, divergent) declines. Under
    Automatic it runs in the pre-FTC `[0,1]` block.
  The definite mechanisms name themselves only: the actual mechanism is
  chosen from the spec type, so on a definite integral any *other* method name
  is passed through to the inner indefinite integration that produces the
  antiderivative, and either definite-mechanism name on the indefinite
  `Integrate[f, x, …]` form is a no-op (stays unevaluated).
  Unknown method names emit `Integrate::method` and bubble back.
- Universal correctness predicate: `Cancel[Together[D[Integrate[f,x],x] - f]] === 0`.

**Attributes:** `Protected`.

## References

**See also:** [DSolve](../../calculus/DSolve/), [PolynomialQuotientRemainder](../../calculus/PolynomialQuotientRemainder/), [Apart](../../algebra/Apart/), [Log](../../elementary-functions/Log/), [ToRadicals](../../solutions-of-equations/ToRadicals/), [Root](../../solutions-of-equations/Root/), [RootSum](../../solutions-of-equations/RootSum/), [Function](../../functional-programming/Function/)

- Bronstein, "Symbolic Integration I: Transcendental Functions", 2nd ed. (Springer, 2005).
- Geddes, Czapor & Labahn, "Algorithms for Computer Algebra" (Kluwer, 1992), ch. 11–12.
- M. Bronstein, *Symbolic Integration I: Transcendental Functions*, 2nd ed. (Springer, 2005).
- B. M. Trager, "Algebraic Factoring and Rational Function Integration", ACM SYMSAC 1976.
- K. O. Geddes, S. R. Czapor, G. Labahn, *Algorithms for Computer Algebra* (Kluwer, 1992).
- D. J. Jeffrey, A. D. Rich, "The Evaluation of Trigonometric Integrals Avoiding Spurious Discontinuities", ACM TOMS 20(1), 1994, 124–135.
- G. W. Cherry, "Integration in Finite Terms with Special Functions: The Logarithmic Integral", SIAM J. Comput. 1986.
- G. W. Cherry, "An Analysis of the Rational Exponential Integral", SIAM J. Comput. 18(5), 1989.
- S. Blake, "Parallel Integration over Simple Radical Extensions II" (ParallelMixedTower / ParallelMixedSpecial).
- S. Blake, "A Generalisation of Goursat's Algorithm for Integration in Finite Terms" and "A Fourth-Root Generalisation of Goursat's Algorithm".
- K. Roach, "Symbolic Integration", 1992 (undefined-function integrands, §1.7).
- Boulnois, arXiv:2308.09619 — differentiation under the integral sign framed as a first-order ODE in the parameter.
- CRC Standard Mathematical Tables, 31st ed. — reference integral tables.
- Source: [`src/calculus/integrate.c`](https://github.com/stblake/mathilda/blob/main/src/calculus/integrate.c)
- Specification: [`docs/spec/builtins/calculus.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/calculus.md)
- Tests: [`tests/test_cherry_dilog.c`](https://github.com/stblake/mathilda/blob/main/tests/test_cherry_dilog.c)
- Tests: [`tests/test_cherry_dilog_exp.c`](https://github.com/stblake/mathilda/blob/main/tests/test_cherry_dilog_exp.c)
- Tests: [`tests/test_cherry_ei.c`](https://github.com/stblake/mathilda/blob/main/tests/test_cherry_ei.c)
- Tests: [`tests/test_cherry_li.c`](https://github.com/stblake/mathilda/blob/main/tests/test_cherry_li.c)

## Notes & additional examples

### Notes

- **Forms.** `Integrate[f, x]` (indefinite, no constant of integration), `Integrate[f, {x, a, b}]` (definite, including improper bounds), iterated `Integrate[f, {x, a, b}, {y, c, d}]` (innermost spec last), and complex line integrals `Integrate[f, {z, z0, z1, …}]` along a polyline. A `List` integrand threads; `Integrate` itself is not `Listable`, so a range spec is never threaded over. Attributes: `Protected`.
- **Options.** `Method -> Automatic | "name"` (a named method runs that one stage with no fallback), `Assumptions -> cond` (used by the definite stages to prove convergence strips and parameter signs; the result is then `Refine`d under it), `PrincipalValue -> True` (Cauchy principal value at interior poles), and `TimeConstraint -> 3` (seconds for the automatic indefinite cascade; `Infinity` removes the bound; a pinned `Method` is unbounded unless you pass one).
- **Method names.** Indefinite: `"BronsteinRational"`, `"Undefined"`, `"LinearRadicals"`, `"QuadraticRadicals"`, `"LinearRatioRadicals"`, `"ChebychevAlgebraic"`, `"Weierstrass"`, `"DerivativeDivides"` (or `{"DerivativeDivides", "Substitution" -> u}`), `"RischTranscendental"`, `"CRCTable"`, `"PolyLogLadder"`, `"ParallelMixedTower"`, `"GoursatAlgebraic"`, `"ParallelMixedSpecial"`. Definite: `"NewtonLeibniz"`, `"Residue"`, `"Symmetry"`, `"Piecewise"`, `"Beta"`, `"TrigPower"`, `"SeriesExpansion"`, `"FeynmanParameter"`, `"King"`, `"Mellin"` (= `"RamanujanMasterTheorem"`), `"IntegralRepresentation"`, `"MeijerConvolution"`, `"SinPowerMonomial"`, `"OscillatoryPower"`, `"RationalLog"`, `"DiffUnderInt"`, `"LineIntegral"`.
- **Answer quality.** The cascade prefers an elementary answer; special functions (`Erf`, `ExpIntegralEi`, `PolyLog`, `FresnelS`, the elliptic family, …) come only from stages that need them. Antiderivatives are not post-simplified, so the form can differ from Mathematica's by a constant or by an equivalent rearrangement. `Integrate::nonelem` reports a proof of non-elementarity only when no later stage answers.
- **Soundness.** Every definite stage either is correct by construction behind proved gates (residue, Beta, Mellin strips) or is cross-checked numerically; a value still containing the integration variable is discarded, a pole on the path gives `Integrate::idiv` and stays unevaluated, and an undecidable case is left unevaluated rather than guessed.
- **Inexact input.** Machine or arbitrary-precision coefficients are rationalised, integrated exactly, and numericalised back at the input precision. For a purely numeric definite integral, `NIntegrate` is usually the better tool.
- See also: `NIntegrate`, `D`, `Sum`, `Limit`, `Residue`, `DSolve`.
