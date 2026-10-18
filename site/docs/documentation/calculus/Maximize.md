# Maximize

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Maximize[f, x]`**

maximizes f exactly with respect to x.

**`Maximize[f, {x, y, ...}] / Maximize[{f, cons}, vars] / Maximize[..., dom]`**

as Minimize, for the global maximum.

<details>
<summary>Notes</summary>

Maximize (Protected, not HoldAll) shares Minimize's engine, options, and constraint/domain handling; internally it maximises by minimising -f and negates the objective value in the result.  Returns {f\_max, {x -\> x\_max, ...}}.

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

A double-hump quartic: the two global maxima at x = -1 and x = 1 tie, and one is reported

```mathematica
In[19]:= Maximize[-x^4 + 2 x^2 + 1, x]
Out[19]= {2, {x -> -1}}
```

An equality constraint: largest product for a fixed sum

```mathematica
In[20]:= Maximize[{x y, x + y == 10}, {x, y}]
Out[20]= {25, {x -> 5, y -> 5}}
```

A linear objective on a disk: the maximum is on the boundary in the gradient direction

```mathematica
In[21]:= Maximize[{x + 2 y, x^2 + y^2 <= 5}, {x, y}]
Out[21]= {5, {x -> 1, y -> 2}}
```

A small linear program solved exactly: the optimum is at the vertex (4, 0)

```mathematica
In[22]:= Maximize[{3 x + 2 y, x + y <= 4 && x + 3 y <= 6 && x >= 0 && y >= 0}, {x, y}]
Out[22]= {12, {x -> 4, y -> 0}}
```

A Lagrange-multiplier problem with an irrational optimal point

```mathematica
In[23]:= Maximize[{x^2 y, x^2 + y^2 == 3}, {x, y}]
Out[23]= {2, {x -> -Sqrt[2], y -> 1}}
```

A radical objective, made polynomial with auxiliary variables

```mathematica
In[24]:= Maximize[{Sqrt[x] + Sqrt[1 - x], 0 <= x <= 1}, x]
Out[24]= {Sqrt[2], {x -> 1/2}}
```

The Integers domain: the region is bounded by its continuous relaxation and then enumerated exactly

```mathematica
In[25]:= Maximize[{x y, x^2 + y^2 <= 8 && x >= 0}, {x, y}, Integers]
Out[25]= {4, {x -> 2, y -> 2}}
```

Unbounded above: returns Infinity and emits Maximize::natt

```mathematica
In[26]:= Maximize[x^3 - x, x]
Out[26]= {Infinity, {x -> Indeterminate}}
```

An empty region: -Infinity and Maximize::infeas

```mathematica
In[27]:= Maximize[{x, x^2 + 1 <= 0}, x]
Out[27]= {-Infinity, {x -> Indeterminate}}
```

Transcendental objectives are not handled yet, so the call stays unevaluated

```mathematica
In[28]:= Maximize[Sin[x], x]
Out[28]= Maximize[Sin[x], x]
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

**Algorithm.** `builtin_maximize` calls the shared Minimize driver, `mz_run(res, "Maximize", true)`. The driver forms `fobj = -f`, runs the same exact pipeline as `Minimize` on it, and then negates the reported value with `mz_negate_value`. The pipeline is: make radicals polynomial (`mz_rad_walk`), reformulate a rational objective (`mz_rational_solve`), then dispatch to univariate, piecewise-`Abs`, separable, KKT active-set (`mz_exact_poly`), or quantifier-elimination (`mz_qe_infimum`) solvers, or to the Integers-domain engine (`mz_integer_opt`). Inexact input goes to `NMaximize` through `mz_numeric_fallback`. Messages carry the `Maximize` head and the noun "maximum" (`mz_noun`). An unbounded problem reports `{Infinity, {x -> Indeterminate, ...}}` with `Maximize::natt`, and an infeasible one reports `{-Infinity, ...}` with `Maximize::infeas`. Both are the minimiser's sentinels with the sign flipped.

**Data structures.** Same as `Minimize`: `Expr` trees for the objective and the constraints (equalities `EQ[]`, inequalities `MzIneq{z, strict}`), `Expr***` rows of exact candidate points parsed from `Solve`, and fresh `Minimize$k` symbols for Lagrange multipliers and radical auxiliaries. The fresh symbols keep the `Minimize$` prefix even when called from Maximize, and they are removed from the reported rules. Candidates are ordered by the FLINT qqbar sign oracle `rru_sign_compare` (`src/solve/reduce_real_util.c`). An undecided comparison stops the call and leaves it unevaluated.

**Complexity / limits.** Identical to `Minimize`. The time is spent in `Solve[..., Reals]` on the KKT systems and in the CAD behind the `Reduce`/`Resolve` certificates, which are doubly exponential in the number of variables. There are at most 6 inequalities (64 active sets), and the integer box is capped at 2·10⁶ points. Everything runs under one per-call `TimeConstraint` budget (default 30 s of CPU time). Maximize declines in the same cases as Minimize: transcendental objectives, suprema that are not attained, parametric answers, mixed integer problems, and vector objectives. It is symbolic and has no packed, NDArray or `Compile[]` surface; numeric work belongs in `NMaximize`.

**Attributes:** `Protected`.

## References

**See also:** [Minimize](../../calculus/Minimize/), [NMinimize](../../numerical-calculus/NMinimize/), [NMaximize](../../numerical-calculus/NMaximize/), [HoldAll](../../expression-information/HoldAll/), [Solve](../../solutions-of-equations/Solve/), [Reduce](../../solutions-of-equations/Reduce/), [Sqrt](../../arithmetic/Sqrt/), [Abs](../../arithmetic/Abs/)

- H. W. Kuhn and A. W. Tucker, *Nonlinear programming*, Proc. Second Berkeley Symp. on Mathematical Statistics and Probability (1951) 481–492 — the KKT stationarity system enumerated per active set.
- F. John, *Extremum problems with inequalities as subsidiary conditions*, in Studies and Essays Presented to R. Courant (Interscience, 1948) 187–204 — the singular candidates that complete the compact-region search.
- G. E. Collins, *Quantifier elimination for real closed fields by cylindrical algebraic decomposition*, LNCS 33 (1975) 134–183 — the CAD behind the Reduce/Resolve certificates.
- Source: [`src/calculus/minimize.c`](https://github.com/stblake/mathilda/blob/main/src/calculus/minimize.c)
- Specification: [`docs/spec/builtins/calculus.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/calculus.md)
- Tests: [`tests/test_minimize.c`](https://github.com/stblake/mathilda/blob/main/tests/test_minimize.c)

## Notes & additional examples

### Notes

- `Maximize` is `Minimize` applied to `-f` with the value negated. It has the same argument forms (`f`, `{f, cons}`, variable list, optional domain `Reals`/`Integers`) and returns `{f_max, {x -> x_max, ...}}` with exact values.
- The sentinels are swapped relative to Minimize: unbounded gives `{Infinity, {x -> Indeterminate, ...}}` with `Maximize::natt`; infeasible gives `{-Infinity, ...}` with `Maximize::infeas`. Both go through the message funnel, so `Quiet` and `Check` work.
- Inexact input (`Maximize[1.5 x - x^2, x]`) is passed to `NMaximize`. Compared with `NMaximize`, Maximize proves global optimality exactly (or declines), while `NMaximize` always returns a machine-precision point from a stochastic search.
- `TimeConstraint -> t` (CPU seconds, default 30) bounds the whole call. Attributes: `{Protected}`, not `HoldAll`.
- See also: `Minimize`, `NMaximize`, `FindMaximum`, `Reduce`.
