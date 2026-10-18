# Minimize

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Minimize[f, x]`**

minimizes f exactly with respect to x.

**`Minimize[f, {x, y, ...}]`**

exact minimum with respect to several variables.

**`Minimize[{f, cons}, vars]`**

exact minimum of f subject to the constraints cons (a single constraint or an And of equalities (==) and inequalities (\<, \<=, \>, \>=); additional list elements {f, c1, c2, ...} are implicitly And-ed).

**`Minimize[..., dom]`**

constrains the variables to the domain dom (Reals, the default).

<details>
<summary>Notes</summary>

Minimize finds a GLOBAL minimum and returns {f\_min, {x -\> x\_min, ...}}.  Exact input yields an exact result (Integer / Rational / radical / Root\[\] values); with inexact (approximate) input it automatically calls NMinimize.  Special return forms: {-Infinity, {x -\> Indeterminate, ...}} when f is unbounded below, {Infinity, {x -\> Indeterminate, ...}} when the constraints are infeasible, and {f\_inf, {x -\> Indeterminate, ...}} with the Minimize::natt message when the infimum f\_inf is not attained. This release handles: univariate polynomial objectives; multivariate polynomial objectives with isolated critical points; and polynomial objectives under polynomial constraints over the Reals (small linear programs included).  Every result is backed by an exact decision procedure (Solve/Reduce over the reals, the exact real-algebraic sign oracle); when a case cannot be decided exactly the expression is left unevaluated rather than guessed.  Minimize is Protected (not HoldAll).

</details>

## Examples (28)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (9)

```mathematica
In[1]:= Minimize[2 x^2 - 3 x + 5, x]
Out[1]= {31/8, {x -> 3/4}}

In[2]:= Minimize[Expand[(x - 1) (x - 2) (x - 4) (x - 7) (x - 9) (x - 10)], x]
Out[2]= {-35721/64, {x -> 11/2}}

In[3]:= Minimize[x^2 + y^2, {x, y}]
Out[3]= {0, {x -> 0, y -> 0}}

In[4]:= Minimize[{x + y, x^2 + y^2 <= 1}, {x, y}]
Out[4]= {-Sqrt[2], {x -> -1/Sqrt[2], y -> -1/Sqrt[2]}}

In[5]:= Minimize[{2 x + 3 y - z, 1 <= x + y + z <= 2 && 1 <= x - y + z <= 2 && x - y - z == 3}, {x, y, z}]
Out[5]= {3, {x -> 2, y -> -1/2, z -> -1/2}}

In[6]:= Minimize[{x + y, x^2 < -1}, {x, y}] Minimize::infeas: The constraints are infeasible; the feasible region is empty.

In[7]:= Maximize[-x^2, x]
Out[7]= {0, {x -> 0}}
```

Inexact input -> NMinimize

```mathematica
In[8]:= Minimize[2.5 x^2 - 3 x, x]
Out[8]= {-0.9, {x -> 0.6}}
```

Flat valley along x y == 3

```mathematica
In[9]:= Minimize[(x y - 3)^2 + 1, {x, y}]
Out[9]= {1, {x -> -1, y -> -3}}
```

### Worked examples (9)

```mathematica
In[10]:= Minimize[{x^2+y^2, (x-2)^2+(y-3)^2==1}, {x,y}]
Out[10]= {1/169 (-39 + 3 Sqrt[13])^2 + 1/169 (-26 + 2 Sqrt[13])^2, {x -> -1/13 (-26 + 2 Sqrt[13]), y -> -1/13 (-39 + 3 Sqrt[13])}}

In[11]:= Power[g, p/q]
Out[11]= g^(p/q)

In[12]:= Minimize[{Sqrt[x+Sqrt[x]]+Sqrt[x-Sqrt[x]], x>=1}, x]
Out[12]= {Sqrt[2], {x -> 1}}

In[13]:= Minimize[Sum[Abs[x-i^2], {i,1,10}], x]
Out[13]= {275, {x -> 25}}

In[14]:= Minimize[x^2+Abs[x-2], x]
Out[14]= {7/4, {x -> 1/2}}

In[15]:= Minimize[{x+y, x^2+y^2==25}, {x,y}, Integers]
Out[15]= {-7, {x -> -4, y -> -3}}

In[16]:= Minimize[{x^2+y^2, x+y>=3 && 0<=x<=5 && 0<=y<=5}, {x,y}, Integers]
Out[16]= {5, {x -> 2, y -> 1}}

In[17]:= Minimize[(x y - 3)^2 + 1, {x, y}]
Out[17]= {1, {x -> -1, y -> -3}}

In[18]:= Minimize[{x + y, x <= y^2}, {x, y}]
Out[18]= {-Infinity, {x -> Indeterminate, y -> Indeterminate}}
```

### Applications (10)

A quartic with two wells: the global minimiser is the least real root of the cubic f' = 0, kept exact as a Root object

```mathematica
In[19]:= Minimize[x^4 - 3 x^2 + x, x]
Out[19]= {Root[1 - 6 #1 + 4 #1^3 &, 1]^4 + Root[1 - 6 #1 + 4 #1^3 &, 1] - 3 Root[1 - 6 #1 + 4 #1^3 &, 1]^2, {x -> Root[1 - 6 #1 + 4 #1^3 &, 1]}}
```

Numerically -3.51391 at x = -1.30084, the deeper well; FindMinimum started at x = 2 stops in the shallower well at -1.07023

```mathematica
In[20]:= N[Minimize[x^4 - 3 x^2 + x, x]]
Out[20]= {-3.51391, {x -> -1.30084}}
```

Mixed equality and inequality constraints: AM-GM on the simplex

```mathematica
In[21]:= Maximize[{x y z, x + y + z == 3 && x >= 0 && y >= 0 && z >= 0}, {x, y, z}]
Out[21]= {1, {x -> 1, y -> 1, z -> 1}}
```

Squared distance from the origin to a circle, (Sqrt[13] - 1)^2; the raw value is an unsimplified radical sum

```mathematica
In[22]:= Simplify[First[Minimize[{x^2 + y^2, (x - 2)^2 + (y - 3)^2 == 1}, {x, y}]]]
Out[22]= 14 - 2 Sqrt[13]
```

Rational objective: the denominator is proved positive on the region, giving 2 Sqrt[5] - 4 at x = Sqrt[5] - 2

```mathematica
In[23]:= Minimize[{(x^2 + 1)/(x + 2), x >= 0}, x]
Out[23]= {-(4 - 2 Sqrt[5]), {x -> -(2 - Sqrt[5])}}
```

Piecewise-linear objective minimised exactly over its sign cells; every x in [25, 36] is optimal

```mathematica
In[24]:= Minimize[Sum[Abs[x - i^2], {i, 1, 10}], x]
Out[24]= {275, {x -> 25}}
```

The Integers domain: lattice points on a circle, found with Solve over the integers

```mathematica
In[25]:= Minimize[{x + y, x^2 + y^2 == 25}, {x, y}, Integers]
Out[25]= {-7, {x -> -4, y -> -3}}
```

Unbounded below: proved with quantifier elimination, returns -Infinity and emits Minimize::natt

```mathematica
In[26]:= Minimize[{x + y, x <= y^2}, {x, y}]
Out[26]= {-Infinity, {x -> Indeterminate, y -> Indeterminate}}
```

An empty region: Infinity with Indeterminate rules and Minimize::infeas

```mathematica
In[27]:= Minimize[{x + y, x^2 < -1}, {x, y}]
Out[27]= {Infinity, {x -> Indeterminate, y -> Indeterminate}}
```

Transcendental objectives are not handled yet, so the call stays unevaluated

```mathematica
In[28]:= Minimize[Exp[x] - x, x]
Out[28]= Minimize[-x + E^x, x]
```

## Algorithm

minimize.c — Minimize / Maximize: exact (symbolic) global optimization.

The exact counterpart to NMinimize/NMaximize (nm_driver.c). Where those run a stochastic numeric global search, Minimize finds a GLOBAL optimum exactly and returns {f_opt, {x -> x_opt, ...}} with exact (Integer / Rational / radical / Root[]) values, delegating to NMinimize only when the input is inexact.

Soundness invariant (every path honours it): a result is emitted only when its value, attainment and feasibility are each backed by an EXACT oracle —

```text
  - Solve[..., Reals] / the zero-dimensional engine for the candidate set,
  - the FLINT qqbar sign oracle (rru_sign_of / rru_sign_compare) for ordering,
  - Reduce[..., Reals] for the global lower-bound certificate and emptiness.
```

Any "don't know" from those oracles (rru_sign == -2, an unclean Solve/Reduce result) is a hard stop: the path is abandoned and the expression is left unevaluated (NULL) rather than guessed. The ONLY licence to return an inexact answer is inexact INPUT, which routes to NMinimize from the top.

Scope this release (see docs/spec/builtins/numerical_calculus.md):

```text
  (a) univariate polynomial objective (unbounded/attainment via the degree
      parity + leading-coefficient tail theorem);
  (b) multivariate polynomial objective with isolated critical points;
  (c) polynomial objective under polynomial constraints over the Reals
      (KKT/active-set enumeration + a Reduce lower-bound certificate; small
      linear programs fall out as the degenerate case);
  (d) unconstrained multivariate polynomial with a POSITIVE-DIMENSIONAL
      (non-isolated) minimizer set — a flat valley the critical-point method
      in (b) cannot resolve (Solve returns Solve::nsdim). The global infimum
      is read directly off real quantifier elimination,
      Reduce[ForAll[{vars}, f >= b], {b}, Reals], and a minimiser is realised
      with FindInstance (mz_qe_infimum). An unbounded objective (Reduce ===
      False) reports -Infinity; an un-attained infimum declines.
```

Deferred (declines, never guesses): transcendental closed forms, parametric Piecewise answers, CONSTRAINED or non-attained positive-dimensional minima, general not-attained-under-constraints via quantifier elimination, exact Integers/ILP optimisation.

Minimize/Maximize are Protected but NOT HoldAll (matching Mathematica): the variables are unbound symbols that evaluate to themselves, so the objective and constraints arrive in symbolic form and the engine proceeds by value-preserving substitution (D, Solve, Reduce).

## Implementation notes

**Algorithm.** `builtin_minimize` and `builtin_maximize` both enter `mz_run` (Maximize minimises `-f` and negates the value with `mz_negate_value`). `mz_run` splits `{f, cons...}` from a bare objective, reads the variable list (bare symbols only) and the domain (`Reals` by default, `Integers`, or `Complexes`, which stays unevaluated for exact input), and rewrites an `Element[{vars}, Integers]` covering every variable into the Integers domain (`mz_strip_domain_element`). Inexact input goes to `NMinimize`/`NMaximize` (`mz_numeric_fallback`), and this is the only path that may return an approximate answer. For exact input over the Reals the order is: (1) radicals and fractional powers are made polynomial by `mz_rad_walk`, which adds a real `u` with `u^q == g && u >= 0` for each `g^(p/q)`; (2) a rational objective `p/q` becomes `Minimize[w, cons && p - w q == 0]` once `q` has a strict sign on the closure (`mz_rational_solve`); (3) without constraints, `mz_unconstrained_min` runs: a univariate `Abs`-piecewise objective is minimised over its sign cells (`mz_univar_piecewise`); a univariate polynomial is decided by degree parity and the sign of the leading coefficient (`mz_univar_poly`, which reports `-Infinity` and `::natt` for an odd degree or a negative leading coefficient); a sum whose terms share no variables is split into blocks with union–find (`mz_separable_min`); then come isolated critical points (`mz_exact_poly`) and finally a flat valley, whose infimum is read off `Reduce[ForAll[vars, f >= b], {b}, Reals]` (`mz_qe_infimum`) and whose minimiser comes from `FindInstance` or a univariate slice; (4) with constraints, the steps are `mz_eliminate_solve` (substitute out any `v == h(others)`), then `mz_exact_poly`, then `mz_qe_infimum` over the region, all three sharing one deadline. `mz_exact_poly` solves the KKT system `∇f − Σ λ_c ∇z_c = 0, z_c = 0` with `Solve[..., Reals]` for every subset of active inequalities (at most `MZ_MAX_INEQ` = 6 inequalities, so at most 64 subsystems). It drops candidates outside the closure and picks the least one by exact comparison. If the closure is provably compact (`mz_bounded_region_cons`: `Resolve[Exists[vars, CL && Σx² > 2^20], Reals] === False`) and the Fritz–John singular points are added (`mz_add_singular_candidates`), the extreme-value theorem decides the minimum. Otherwise the Reduce certificate `closure ⇒ f ≥ f*` must be proved (`mz_prove_lower_bound`; radicals in `f*` are cleared to new variables with their defining relations). Finally the point must satisfy the original strict constraints. An empty region (`Reduce[cons] === False`) gives `{Infinity, {x -> Indeterminate, ...}}` and `::infeas`. For the Integers domain, `mz_integer_opt` tries `Solve[cons, vars, Integers]` (`mz_integer_solve`), then a coercive parametric Diophantine family (`mz_integer_parametric`), then an exact enumeration of a box bounded by the continuous relaxation (`mz_integer_box`).

**Data structures.** Problems are `Expr` trees: a `vars` array of symbols, equality zero-expressions `EQ[]`, and inequalities `MzIneq{z, strict}` meaning `z <= 0` or `z < 0`. Gradient and Jacobian rows are `Expr**` arrays built with `D`. Candidate points are `Expr***` rows of exact values (Integer, Rational, radical, or `Root`) parsed from `Solve` output by `mz_parse_points`, which rejects any parametric or unclean solution set. Lagrange multipliers and auxiliary radical variables are fresh `Minimize$k` symbols that are removed from the reported rules (`mz_strip_aux_rules`). Candidates are ordered by the FLINT qqbar sign oracle `rru_sign_compare` / `rru_sign_of` (`src/solve/reduce_real_util.c`). When that oracle returns "undecided" (`-2`), the call stops and the expression is left unevaluated. The result is `{f_opt, {x -> x_opt, ...}}`.

**Complexity / limits.** The cost is dominated by `Solve[..., Reals]` on the KKT systems (Gröbner/zero-dimensional) and by the CAD inside the `Reduce`/`Resolve` certificates. CAD is doubly exponential in the number of variables, so certificates and flat-valley quantifier elimination are reliable only for about two or three variables and modest degree. The active-set enumeration is `2^nin` and is capped at six inequalities. The integer box is capped at 2·10⁶ points. Every internal probe runs inside `TimeConstrained` with the time left of a per-call budget (`TimeConstraint -> t`, default 30 s, `Infinity` for no limit). The budget is measured with C `clock()`, which is CPU time. A problem that cannot be decided exactly, or that runs out of time, returns unevaluated. Results are never guessed. These cases decline: transcendental objectives, minima approached but not attained on a strict boundary, parametric (`Piecewise`) answers, mixed integer/continuous problems, vector objectives and indexed variables. The minimum value is returned exactly as the objective evaluates at the optimal point, so it can be an unreduced radical or `Root` expression; apply `Simplify` or `RootReduce` to get a canonical form. Minimize is symbolic, so it has no packed-array, NDArray or `Compile[]` surface. Numeric workloads should use `NMinimize`.

**Attributes:** `Protected`.

## References

**See also:** [Maximize](../../calculus/Maximize/), [NMinimize](../../numerical-calculus/NMinimize/), [NMaximize](../../numerical-calculus/NMaximize/), [HoldAll](../../expression-information/HoldAll/), [Solve](../../solutions-of-equations/Solve/), [Reduce](../../solutions-of-equations/Reduce/), [Sqrt](../../arithmetic/Sqrt/), [Abs](../../arithmetic/Abs/)

- H. W. Kuhn and A. W. Tucker, *Nonlinear programming*, Proc. Second Berkeley Symp. on Mathematical Statistics and Probability (1951) 481–492 — the KKT stationarity system enumerated per active set.
- F. John, *Extremum problems with inequalities as subsidiary conditions*, in Studies and Essays Presented to R. Courant (Interscience, 1948) 187–204 — the singular (rank-deficient) candidates that complete the compact-region search.
- G. E. Collins, *Quantifier elimination for real closed fields by cylindrical algebraic decomposition*, LNCS 33 (1975) 134–183 — the CAD behind the Reduce/Resolve certificates.
- Source: [`src/calculus/minimize.c`](https://github.com/stblake/mathilda/blob/main/src/calculus/minimize.c)
- Specification: [`docs/spec/builtins/calculus.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/calculus.md)
- Tests: [`tests/test_minimize.c`](https://github.com/stblake/mathilda/blob/main/tests/test_minimize.c)

## Notes & additional examples

### Notes

- `Minimize[f, x]`, `Minimize[f, {x, y, ...}]`, `Minimize[{f, cons}, vars]` and `Minimize[..., dom]` return `{f_opt, {x -> x_opt, ...}}` with **exact** values (Integer, Rational, radical, or `Root`). Each result is checked with exact oracles: `Solve` over the Reals, exact algebraic sign comparison, and `Reduce`/`Resolve` certificates. If any oracle cannot decide, the call is left unevaluated. Minimize does not guess.
- The minimum value is the objective evaluated at the minimiser. It can come back unreduced (for example a polynomial in a `Root`, or a sum of squared radicals), so wrap it in `Simplify` or `RootReduce` for a canonical form.
- **Inexact input** (`Minimize[2.5 x^2 - 3 x, x]`) is passed straight to `NMinimize`. `NMinimize` is a stochastic global search over machine numbers that always returns *some* point. Minimize either proves its answer or declines. For example, `NMinimize[{x + y, x <= y^2}, {x, y}]` returns a large finite negative number, while Minimize proves the problem is unbounded below (`-Infinity`).
- Special returns: unbounded gives `{-Infinity, {x -> Indeterminate, ...}}` with `Minimize::natt`; infeasible gives `{Infinity, {x -> Indeterminate, ...}}` with `Minimize::infeas`. Both messages go through the message funnel, so `Quiet` and `Check` see them. An infimum that is approached but not attained (for example `Minimize[{x, x > 0}, x]`) currently declines.
- Domains: `Reals` (the default) and `Integers`. The Integers domain can also be written as `Element[{vars}, Integers]` among the constraints, but it must cover every variable; mixed integer/continuous problems decline. `Complexes` with exact input stays unevaluated.
- Options: `TimeConstraint -> t` (seconds of CPU time, default 30; `Infinity` for no limit) bounds all internal `Solve`/`Reduce` probes together.
- Attributes: `{Protected}`. Minimize is not `HoldAll`, as in Mathematica, so a variable that has a value is substituted before Minimize sees it.
- Scaling: the CAD certificates are doubly exponential in the number of variables, and at most six inequality constraints are allowed. Problems in a few variables with moderate degree are the practical range.
- See also: `Maximize`, `NMinimize`, `FindMinimum`, `Reduce`, `Solve`, `FindInstance`.
