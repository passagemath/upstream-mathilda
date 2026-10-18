---
references:
  - "E. S. Cheb-Terrab and A. D. Roche, *Integrating factors for second-order ODEs* (1999) — the ReducibleIntegratingFactor / ReducibleFirstIntegral classes, cited in the dispatcher."
  - "J. J. Kovacic, \"An algorithm for solving second order linear homogeneous differential equations\", J. Symbolic Computation 2 (1986)."
source: src/calculus/dsolve.c
---
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
