# DSolve

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`DSolve[eqn, y, x]`**

solves a differential equation for the function y with

**`DSolve[eqn, y[x], x]`**

returns {{y\[x\] -\> expr}} with the solution as an

**`DSolve[{eqn1, ...}, {y1, ...}, x]`**

solves a system; DSolve\[eqn, u, {x, y}\]

<details>
<summary>Notes</summary>

independent variable x, returning {{y -\> Function\[{x}, ...\]}}. expression in x. solves a partial differential equation. Initial/boundary conditions are given as equations at points, e.g. y\[0\] == a, y'\[x0\] == b; DSolve fits the generated constants C\[k\] to them. Every returned branch is verified by back-substitution before it is kept. An automatic DSolve returns only solutions: never a truncated series and never a relation still carrying a derivative of the unknown (those come only from an explicitly pinned series or first-integral Method). DSolve is a cascade polyalgorithm: DSolve\[eqn, y, x, Method -\> "\<Name\>"\] dispatches directly to one method (strict, no fallback), and every method is also callable on its own as DSolve\`\<Name\>\[...\] with its own ?-documentation. Methods (tried roughly in this order): First order:

```text
Quadrature              y^(n)[x] == f(x): integrate n times.
LinearFirstOrder        y'[x] + p y == q: integrating factor.
Separable               y'[x] == g(x)h(y).
Bernoulli               y'[x] == A y[x] + B y^n[x].
Homogeneous             y'[x] == F(y[x]/x).
Exact                   M + N y'[x] == 0, M_y == N_x (+ mu(x)/mu(y) factor search).
Clairaut                y[x] == x y'[x] + f(y'[x]) (+ singular envelope).
Riccati                 y'[x] == q0 + q1 y[x] + q2 y[x]^2: linearize to 2nd-order.
Lagrange                y == x phi(y'[x]) + psi(y'[x]) (parametric solution).
Chini                   y'[x] == f y[x]^n + g y[x] + h (reducible-to-autonomous).
Abel                    y'[x] == f3 y[x]^3 + f2 y[x]^2 + f1 y[x] + f0.
```

FirstOrderSubstitution  y'\[x\] == F(a x + b y\[x\] + c).

```text
Factorable              factor in y'[x], solve each factor.
NthAlgebraic            algebraic in the top derivative.
AlmostLinear            f(x)g(y[x])y'[x] + k(x)l(y[x]) + m(x) == 0.
LinearCoefficients      y'[x] == (a1 x + b1 y[x] + c1)/(a2 x + b2 y[x] + c2).
SeparableReduced        x y'[x]/y == G(x^n y[x]).
```

LieSymmetry (LieGroup)  heuristic Lie point-symmetry backstop. Linear constant-coefficient:

```text
LinearConstantCoefficients   any order, homogeneous + inhomogeneous.
UndeterminedCoefficients     tidy particular for poly/exp/sinusoid forcing.
```

Linear variable-coefficient / 2nd order:

```text
EulerCauchy             a_n x^n y[x]^(n) + ... == g.
ExactODE                higher-order exact linear equation.
NormalForm              y''[x] + P y'[x] + Q y[x] -> z'' == r z (Kovacic prerequisite).
SpecialFunctionForm     Airy/Bessel/Kummer/Gauss recognizers.
Kovacic                 Liouvillian solutions of z''==r z.
FrobeniusSeries         series about a regular singular point (pinned-only).
PowerSeries             series about an ordinary point (pinned-only).
FirstOrderPowerSeries   order-1 power series (pinned-only).
```

OperatorFactor (DFactor)  factor a linear operator into first-order factors. Nonlinear higher-order:

```text
ReductionOfOrder        y''[x] == F(x,y'[x]) (missing y).
AutonomousReduction     y''[x] == f(y[x],y'[x]) (missing x).
Liouville               y''[x] + g(y[x]) y'[x]^2 + h(x) y'[x] == 0.
```

ReducibleFirstIntegral  first integral R(x,y,y') == C\[1\] (pinned-only). Systems (nfun\>1):

```text
DecoupleSystem          each equation involves one function.
TriangularSystem        dependency graph is a DAG.
```

LinearFirstOrderSystem  Y'==A Y+b(x), any constant A (Jordan e^{Ax}). First-order PDE:

```text
PDELinearFirstOrder     a u_{v1}+b u_{v2}+c u==f (characteristics).
```

Options: GeneratedParameters (constant head, default C), Assumptions, Method, IncludeSingularSolutions.  Attributes: Protected.

</details>

## Examples (41)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (15)

```mathematica
In[1]:= DSolve[y''[x] + 4 y'[x] + 5 y[x] == 0, y[x], x]
Out[1]= {{y[x] -> C[1] Cos[x] E^(-2 x) - C[2] Sin[x] E^(-2 x)}}

In[2]:= DSolve[y'''[x] + 4 y'[x] == 5 y[x], y[x], x]
Out[2]= {{y[x] -> C[1] E^x + C[2] E^(-1/2 x) Cos[1/2 Sqrt[19] x] - C[3] E^(-1/2 x) Sin[1/2 Sqrt[19] x]}}
```

BVP

```mathematica
In[3]:= DSolve[{y''[x] + y[x] == 0, y[0] == 0, y[Pi/2] == 1}, y[x], x]
Out[3]= {{y[x] -> Sin[x]}}
```

Kovacic Case 1

```mathematica
In[4]:= DSolve[y''[x] - (x^2 + 3) y[x] == 0, y[x], x]
Out[4]= {{y[x] -> (C[1] (x^2)^(3/4) E^(1/2 x^2))/Sqrt[x] + (C[2] (x^2)^(1/4) HypergeometricPFQ[{1}, {1/2}, x^2] E^(-1/2 x^2))/Sqrt[x]}}
```

No closed form → unevaluated

```mathematica
In[5]:= DSolve[y''[x] + Sin[x] y[x] == 0, y[x], x]
Out[5]= DSolve[Derivative[2][y][x] + Sin[x] y[x] == 0, y[x], x]
```

Euler-Cauchy

```mathematica
In[6]:= DSolve[x^2 y''[x] + 4 x y'[x] + 7 y[x] == 0, y[x], x]
Out[6]= {{y[x] -> (C[1] Cos[1/2 Sqrt[19] Log[x]])/x^(3/2) - (C[2] Sin[1/2 Sqrt[19] Log[x]])/x^(3/2)}}
```

Airy

```mathematica
In[7]:= DSolve[y''[x] - x y[x] == 0, y[x], x]
Out[7]= {{y[x] -> C[1] AiryAi[x] + C[2] AiryBi[x]}}
```

Bessel

```mathematica
In[8]:= DSolve[x^2 y''[x] + x y'[x] + (x^2 - 4) y[x] == 0, y[x], x]
Out[8]= {{y[x] -> C[1] BesselJ[2, x] + C[2] BesselY[2, x]}}
```

First-order substitution

```mathematica
In[9]:= DSolve[y'[x] == (x + y[x])^2, y[x], x]
Out[9]= {{y[x] -> -x - Tan[-C[1] - x]}}
```

2nd-order autonomous, nonlinear

```mathematica
In[10]:= DSolve[y[x] y''[x] == y'[x]^2, y[x], x]
Out[10]= {{y[x] -> C[1] E^(C[2] x)}}
```

System, complex eigenvalues

```mathematica
In[11]:= DSolve[{y'[x] == y[x] - 2 z[x], z'[x] == y[x] - z[x], y[0] == 1, z[0] == 4}, {y[x], z[x]}, x]
Out[11]= {{y[x] -> Cos[x] - 7 Sin[x], z[x] -> Sin[x] + 4 (Cos[x] - Sin[x])}}
```

Defective + singular

```mathematica
In[12]:= DSolve[{y'[t] == 0, x'[t] + y[t] == 0}, {y[t], x[t]}, t]
Out[12]= {{y[t] -> C[1], x[t] -> C[2] - C[1] t}}
```

Defective, x e^{2t}

```mathematica
In[13]:= DSolve[{u'[t] == u[t] - v[t], v'[t] == u[t] + 3 v[t]}, {u[t], v[t]}, t]
Out[13]= {{u[t] -> C[1] E^(2 t) - C[1] t E^(2 t) - C[2] t E^(2 t), v[t] -> C[2] E^(2 t) + C[1] t E^(2 t) + C[2] t E^(2 t)}}
```

Coupled variable-coeff

```mathematica
In[14]:= DSolve[{y'[x] == (2 y[x] + z[x])/x, z'[x] == (y[x] + 2 z[x])/x}, {y, z}, x]
Out[14]= {{y -> Function[{x}, 1/2 x (C[1] - C[2] + C[1] x^2 + C[2] x^2)], z -> Function[{x}, 1/2 x (C[2] - C[1] + C[1] x^2 + C[2] x^2)]}}
```

Inconsistent BVP

```mathematica
In[15]:= DSolve[{y''[x] + y[x] == 0, y[0] == 1, y[Pi] == 1}, y, x]
Out[15]= {}
```

### Scope (14)

Transport

```mathematica
In[16]:= DSolve[D[u[t,x],t] + c D[u[t,x],x] == 0, u, {t,x}]
Out[16]= {{u -> Function[{t, x}, C[1][-c t + x]]}}
```

```mathematica
In[17]:= DSolve[D[u[x,y],x] + 3 D[u[x,y],y] + u[x,y] == 1, u, {x,y}]
Out[17]= {{u -> Function[{x, y}, E^(-x) (E^x + C[1][-3 x + y])]}}
```

Wave/d'Alembert

```mathematica
In[18]:= DSolve[D[u[t,x],{t,2}] == c^2 D[u[t,x],{x,2}], u, {t,x}]
Out[18]= {{u -> Function[{t, x}, C[1][-Sqrt[c^2] t + x] + C[2][Sqrt[c^2] t + x]]}}
```

Laplace/elliptic

```mathematica
In[19]:= DSolve[D[u[x,y],{x,2}] + D[u[x,y],{y,2}] == 0, u, {x,y}]
Out[19]= {{u -> Function[{x, y}, C[1][-I x + y] + C[2][I x + y]]}}
```

```mathematica
In[20]:= DSolve`SeparationOfVariables[D[u[x,t],t] == D[u[x,t],{x,2}], u, {x,t}]
Out[20]= {{u -> Function[{x, t}, E^(-C[3] t) (C[1] E^(-1/2 Sqrt[-4 C[3]] x) + C[2] E^(1/2 Sqrt[-4 C[3]] x))]}}

In[21]:= PDEClassify[D[u[t,x],{t,2}] == c^2 D[u[t,x],{x,2}] /. c -> 2, u, {t,x}]
Out[21]= "Hyperbolic"

In[22]:= {PDEClassify[D[u[x,t],t] == D[u[x,t],{x,2}], u, {x,t}], PDEClassify[D[u[x,y],{x,2}] + D[u[x,y],{y,2}] == 0, u, {x,y}]}
Out[22]= {"Parabolic", "Elliptic"}

In[23]:= DSolve[{D[u[x,t],{t,2}] == c^2 D[u[x,t],{x,2}], u[x,0] == f[x], Derivative[0,1][u][x,0] == g[x]}, u, {x,t}]
Out[23]= {{u -> Function[{x, t}, 1/2 (f[Sqrt[c^2] t + x] + f[-Sqrt[c^2] t + x]) + Integrate[g[K], {K, -Sqrt[c^2] t + x, Sqrt[c^2] t + x}]/(2 Sqrt[c^2])]}}

In[24]:= DSolve[{D[u[x,t],{t,2}] == 4 D[u[x,t],{x,2}], u[x,0] == Sin[x], Derivative[0,1][u][x,0] == 0}, u, {x,t}]
Out[24]= {{u -> Function[{x, t}, 1/2 (Sin[-2 t + x] + Sin[2 t + x])]}}

In[25]:= DSolve[{D[u[x,t],t] == k D[u[x,t],{x,2}], u[x,0] == f[x]}, u, {x,t}]
Out[25]= {{u -> Function[{x, t}, 1/(2 Sqrt[Pi k t]) Integrate[f[K] E^(-(-K + x)^2/(4 k t)), {K, -Infinity, Infinity}]]}}

In[26]:= DSolve[y'[x] + y[x] == a Sin[x], y[x], x]
Out[26]= {{y[x] -> E^(-x) (C[1] + (-1/4 - 1/4*I) a E^((1 + I) x) + (-1/4 + 1/4*I) a E^((1 - I) x))}}

In[27]:= DSolve[{y'[x] == 3 y[x], y[0] == 5}, y[x], x]
Out[27]= {{y[x] -> 5 E^(3 x)}}

In[28]:= DSolve[{y'[x] == -3 y[x]^2, y[1] == 2}, y[x], x]
Out[28]= {{y[x] -> 1/(-5/2 + 3 x)}}

In[29]:= DSolve[y''[x] == 7, y, x]
Out[29]= {{y -> Function[{x}, C[1] + C[2] x + 7/2 x^2]}}
```

### Options (1)

Series only when pinned

```mathematica
In[30]:= DSolve[y''[x] + Sin[x] y[x] == 0, y[x], x, Method -> "PowerSeries"]
Out[30]= {{y[x] -> C[1] + C[2] x - 1/6 C[1] x^3 - 1/12 C[2] x^4 + 1/120 C[1] x^5 + 1/180 (C[1] + C[2]) x^6 + O[x]^7}}
```

### Applications (11)

An initial-value problem: C[1] is fitted and the answer is a pure Function

```mathematica
In[31]:= DSolve[{y'[x] == x + y[x], y[0] == 1}, y, x]
Out[31]= {{y -> Function[{x}, E^x (2 + E^(-x) (-1 - x))]}}
```

Variation of parameters; the particular needs SinIntegral and CosIntegral

```mathematica
In[32]:= DSolve[y''[x] + 4 y[x] == Log[x], y[x], x]
Out[32]= {{y[x] -> C[1] Cos[2 x] - C[2] Sin[2 x] + 1/4 (Log[x] - Cos[2 x] CosIntegral[2 x] - Sin[2 x] SinIntegral[2 x])}}
```

Clairaut: the line family plus the singular envelope -x^2/4

```mathematica
In[33]:= DSolve[y[x] == x y'[x] + y'[x]^2, y[x], x, IncludeSingularSolutions -> True]
Out[33]= {{y[x] -> C[1]^2 + C[1] x}, {y[x] -> -1/4 x^2}}
```

A Bernoulli equation with an irrational exponent, answered as an implicit first integral

```mathematica
In[34]:= DSolve[y'[x] == y[x]^Sqrt[2] + y[x], y[x], x]
Out[34]= {{x + (-Sqrt[2] Log[y[x]] + Log[y[x] + y[x]^Sqrt[2]])/(-1 + Sqrt[2]) == C[1]}}
```

Linear in x(y) once the roles of x and y are swapped

```mathematica
In[35]:= DSolve`InverseOrientation[y'[x] == y[x]/(x + y[x]^3), y[x], x]
Out[35]= {{1/2 (-y[x]^2 + (2 x)/y[x]) == C[1]}}
```

Kovacic: a Liouvillian solution x E^(x^2/2) and its partner

```mathematica
In[36]:= DSolve[y''[x] - (x^2 + 3) y[x] == 0, y[x], x]
Out[36]= {{y[x] -> (C[1] (x^2)^(3/4) E^(1/2 x^2))/Sqrt[x] + (C[2] (x^2)^(1/4) HypergeometricPFQ[{1}, {1/2}, x^2] E^(-1/2 x^2))/Sqrt[x]}}
```

Whittaker's equation, emitted in its confluent 1F1 form

```mathematica
In[37]:= DSolve[x^2 y''[x] + (-x^2/4 + 2 x - 3/16) y[x] == 0, y[x], x]
Out[37]= {{y[x] -> C[1] E^(-1/2 x) x^(1/2 + 1/4 Sqrt[7]) HypergeometricPFQ[{-3/2 + 1/4 Sqrt[7]}, {1 + 1/2 Sqrt[7]}, x] + C[2] E^(-1/2 x) x^(1/2 - 1/4 Sqrt[7]) HypergeometricPFQ[{-3/2 - 1/4 Sqrt[7]}, {1 - 1/2 Sqrt[7]}, x]}}
```

Impulse forcing: the Green's function of the operator

```mathematica
In[38]:= DSolve[{y'[t] + y[t] == DiracDelta[t - 1], y[0] == 0}, y[t], t]
Out[38]= {{y[t] -> HeavisideTheta[-1 + t] E^(1 - t)}}
```

A first-order linear system with initial data

```mathematica
In[39]:= DSolve[{x'[t] == y[t], y'[t] == -x[t], x[0] == 1, y[0] == 0}, {x[t], y[t]}, t]
Out[39]= {{x[t] -> Cos[t], y[t] -> -Sin[t]}}
```

No closed form: left unevaluated rather than truncated

```mathematica
In[40]:= DSolve[y''[x] + Sin[x] y[x] == 0, y[x], x]
Out[40]= DSolve[Derivative[2][y][x] + Sin[x] y[x] == 0, y[x], x]
```

A series answer only on request

```mathematica
In[41]:= DSolve[y''[x] + Sin[x] y[x] == 0, y[x], x, Method -> "PowerSeries"]
Out[41]= {{y[x] -> C[1] + C[2] x - 1/6 C[1] x^3 - 1/12 C[2] x^4 + 1/120 C[1] x^5 + 1/180 (C[1] + C[2]) x^6 + O[x]^7}}
```

## Algorithm

dsolve.c — DSolve dispatcher (cascade polyalgorithm).

Mirrors src/calculus/integrate.c: a Method-option enum selects either the automatic cascade (try each method until one returns a non-NULL result) or a

```text
single pinned method (strict, no fallback).  A per-command fail-memo keyed on
```

eval_toplevel_id() suppresses the fixed-point loop's redundant re-entry on a problem the deterministic cascade already declined, and g_dsolve_depth distinguishes the outermost user call from internal recursions.

The shared problem substrate (parse / verify / fit / assemble) is in dsolve_common.c; each method is one file src/calculus/dsolve_<method>.c.

## Implementation notes

**Algorithm.** `builtin_dsolve` is a cascade polyalgorithm that mirrors
`Integrate`. `dsolve_parse` builds a shared `DSolveProblem`. A `Method` string is
mapped by `ds_method_from_string` to either `DS_AUTOMATIC` or one pinned method,
which runs strictly with no fallback. The dispatcher then branches on the
problem's shape:

- **PDE:** `pde1` (linear first order), then quasilinear Lagrange, Clairaut and
  Charpit, then `pde2` (second order), then the wave (d'Alembert) and heat-kernel
  IVPs.
- **System (`nfun > 1`):** decouple, triangular, higher-order state augmentation
  (`sys_reduce`), constant-coefficient `linsys` and its trig form,
  variable-coefficient, the 2×2 commutative class, and 2-D autonomous via the
  phase plane.
- **Scalar ODE:** the long `DS_AUTOMATIC` chain, where each step is
  `if (!result) result = dsolve_run…(&P, <method>_try)`.

The scalar chain runs roughly in decreasing order of specificity:

1. Factorable and NthAlgebraic split products and powers of the top derivative.
2. The named first-order classes: Quadrature, LinearFirstOrder, Bernoulli,
   Homogeneous, Separable (explicit, then implicit), PolynomialShiftSubstitution,
   Linearizable, Exact (explicit, then implicit), Clairaut, and Lagrange
   (parametric).
3. Linear ODEs: PiecewiseForcing (interval-stitching for `UnitStep`/`Piecewise`
   IVPs), UndeterminedCoefficients, LinearConstantCoefficients, EulerCauchy, and
   ExactODE.
4. Second-order special functions: SpecialFunctionForm (Airy, Bessel, Legendre,
   Kummer/Tricomi, Gauss 2F1, Whittaker via the Liouville normal form
   `DSolve`NormalForm`), PolynomialSolution, Kovacic, ChangeOfVariable, and a
   variation-of-parameters backstop.
5. Higher order: GeneralizedAiry, SymmetricSquare, OperatorFactor,
   ReductionOfOrder, and LowerDerivativeReduction.
6. First-order substitutions: FirstOrderSubstitution, Riccati, Chini, Abel,
   AbelAIR, LinearCoefficients, AlmostLinear, and SeparableReduced.
7. Autonomous, Liouville, and the nonlinear second-order backstops
   ReducibleIntegratingFactor and SecondOrderSymmetry (Lie point symmetry).
8. Implicit fallbacks, then SolvableForY (the "dp" method), HomogeneousG
   (isobaric), and first-order LieSymmetry.
9. The late passes added in v0.41x–v0.45x: nonlinear exactness (the jet peel),
   NthAlgebraicImplicit, BernoulliExpanded (a `yⁿ` term hidden in a product),
   ExactMuU (integrating factor `μ(u)` for `u ∈ {x±y, xy, x²±y²}`), and
   InverseOrientation, which solves `dx/dy = 1/F` for `x(y)`. InverseOrientation
   recurses the whole engine, so it runs last among the deterministic methods.
10. Absolutely last: a VariationOfParameters variant that returns
    non-elementary particular integrals as `Inactive[Integrate]`.

The series methods (`FrobeniusSeries`/`PowerSeries`, `FirstOrderPowerSeries`),
`SolvableForX`, `ReducibleFirstIntegral` and `EigenvalueProblem` are
**pinned-only**. Since v0.361 the automatic answer also has to pass
`ds_automatic_answer_admissible`, which rejects any `SeriesData` result and any
relation that still contains a derivative of the unknown. So an automatic
`DSolve` returns a solution or stays unevaluated. Most methods also register a
qualified head (`DSolve`Kovacic[eqn, y, x]`, `DSolve`InverseOrientation[...]`, …)
that runs that method alone. The `Method` string table does not list every late
pass: InverseOrientation, ExactMuU, BernoulliExpanded and HomogeneousG can be
reached only through their qualified heads.

**Data structures.** The problem substrate (parse, verify, fit, assemble) lives in
`dsolve_common.c`. The dispatcher only sequences the method `*_try` functions,
each of which returns an array of solution branches or `NULL`. These are wrapped
in `dsolve_run`, `dsolve_run_implicit`, `dsolve_run_parametric`,
`dsolve_run_first_integral`, `dsolve_run_system` and `dsolve_run_pde`, which
verify each branch by back-substitution (the implicit-function rule for
`F(x, y) == C`) and fit initial or boundary conditions. A per-command fail-memo
(`ds_fail_tab`, keyed on `eval_toplevel_id()`) stops the fixed-point loop from
re-running a cascade that has already declined. `g_dsolve_depth` separates the
user's call from nested speculative sub-solves; a nested call runs under
`mth_msg_suppress_push` and rolls back its message count, so its probe noise
reaches neither the user nor an enclosing `Check[]`. For the whole cascade,
`Power::infy`/`Infinity::indet` are muted with `arith_warnings_mute_push`.

**Complexity / limits.** The cost is the sum of the declining probes up to the
first method that claims the equation. The order is tuned against Nasser's
DSolve corpus (sections §2.2.37–§2.2.55 added since v0.266), so that each family
is reached before a more general method spins on it. Several methods carry
per-attempt `TimeConstrained` deadlines and decline memos, and the A/B kill
switches (`MATHILDA_NO_INVERSE`, `MATHILDA_NO_HOMOG_G`, …) isolate individual
slots. Options default to `GeneratedParameters -> C`, `Assumptions -> True`,
`Method -> Automatic` and `IncludeSingularSolutions -> False`.

**Attributes:** `Protected`.

## References

**See also:** [Integrate](../../calculus/Integrate/), [HoldAll](../../expression-information/HoldAll/), [SeriesData](../../power-series/SeriesData/), [PolynomialQ](../../algebra/PolynomialQ/), [Solve](../../solutions-of-equations/Solve/), [TimeConstrained](../../time-and-date/TimeConstrained/), [FactorSquareFreeList](../../structural-manipulation/FactorSquareFreeList/), [UnitStep](../../elementary-functions/UnitStep/)

- E. S. Cheb-Terrab and A. D. Roche, *Integrating factors for second-order ODEs* (1999) — the ReducibleIntegratingFactor / ReducibleFirstIntegral classes, cited in the dispatcher.
- J. J. Kovacic, "An algorithm for solving second order linear homogeneous differential equations", J. Symbolic Computation 2 (1986).
- Source: [`src/calculus/dsolve.c`](https://github.com/stblake/mathilda/blob/main/src/calculus/dsolve.c)
- Specification: [`docs/spec/builtins/calculus.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/calculus.md)
- Tests: [`tests/test_dsolve.c`](https://github.com/stblake/mathilda/blob/main/tests/test_dsolve.c)
- Tests: [`tests/test_dsolve_m12_stress.c`](https://github.com/stblake/mathilda/blob/main/tests/test_dsolve_m12_stress.c)
- Tests: [`tests/test_dsolve_m14_stress.c`](https://github.com/stblake/mathilda/blob/main/tests/test_dsolve_m14_stress.c)
- Tests: [`tests/test_dsolve_m17_stress.c`](https://github.com/stblake/mathilda/blob/main/tests/test_dsolve_m17_stress.c)

## Notes & additional examples

### Notes

- `DSolve[eqn, y[x], x]` returns rules for `y[x]`; `DSolve[eqn, y, x]` returns `y -> Function[...]`, which also replaces derivatives, so `eqn /. sol` checks the answer. Constants are `C[1], C[2], ...`; rename them with `GeneratedParameters`.
- Options: `GeneratedParameters -> C`, `Assumptions -> True`, `Method -> Automatic`, `IncludeSingularSolutions -> False`. Attributes: `Protected` only (not `HoldAll`).
- **An automatic answer is always a solution.** Every branch is checked by back-substitution. A truncated series, or a relation that still contains `y'` (a first integral), is never returned from `Method -> Automatic`. Those come only from a pinned method (`"PowerSeries"`/`"FrobeniusSeries"`, `"FirstOrderPowerSeries"`, `"ReducibleFirstIntegral"`). An equation the cascade cannot close stays unevaluated, as in Mathematica.
- Answers may be explicit rules, implicit equations `F(x, y[x]) == C[1]` (separable, exact, Bernoulli, Abel, Lie symmetry), or parametric (Lagrange, SolvableForY). A non-elementary quadrature can stay as `Inactive[Integrate]`.
- `Method -> "<name>"` pins one method with no fallback, e.g. `"Separable"`, `"Bernoulli"`, `"Kovacic"`, `"LieSymmetry"`, `"EulerCauchy"`, `"SpecialFunctionForm"`. The automatic-only late passes, such as `DSolve`InverseOrientation`, `DSolve`ExactMuU`, `DSolve`BernoulliExpanded` and `DSolve`HomogeneousG`, can also be called directly as qualified heads.
- Systems (`{y1, y2, ...}`) and PDEs (`u, {x, t}`) go through their own dispatcher branches: decoupling, triangular, constant- and variable-coefficient linear, 2-D autonomous; PDEs use characteristics, Charpit, d'Alembert and heat-kernel IVPs. `DSolve`EigenvalueProblem` handles Sturm–Liouville problems.
- Known gap: for Laguerre-type `x y'' + (1 - x) y' + n y == 0` with integer `n`, the `HypergeometricU[-n, 1, x]` partner is a polynomial multiple of the first solution, so the returned basis is not independent.
- See also: `NDSolve`, `Integrate`, `Series`, `PDEClassify`.
