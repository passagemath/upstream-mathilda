---
references:
  - "G. E. Collins and H. Hong, *Partial Cylindrical Algebraic Decomposition for Quantifier Elimination*, J. Symbolic Computation **12** (1991) 299-328."
  - "S. McCallum, *An Improved Projection Operation for Cylindrical Algebraic Decomposition*, in *Quantifier Elimination and Cylindrical Algebraic Decomposition* (Springer, 1998) 242-268."
  - "S. Basu, R. Pollack and M.-F. Roy, *Algorithms in Real Algebraic Geometry*, 2nd ed. (Springer, 2006), Ch. 5 & 11 (sign determination, CAD)."
  - "A. Schrijver, *Theory of Linear and Integer Programming* (Wiley, 1986), §12.2 (Fourier-Motzkin elimination)."
source:
  - src/solve/reduce.c
  - src/solve/reduce_cad.c
  - src/solve/reduce_sos.c
---
**Algorithm.** `builtin_reduce` wraps `reduce_impl` in an internal-function
message mute (`mth_msg_ifun_suppress`), so the speculative probes the engine
runs internally do not leak diagnostics. `reduce_impl` peels trailing option
`Rule`s, validates the variable spec, and rewrites a `List` in the `expr` slot
into an `And`. In the one-argument form, and when the second slot holds a domain
symbol (`reduce_is_domain_symbol`: `Reals`, `Complexes`, `Integers`,
`Rationals`), `reduce_infer_var_list` takes the variables to be every free
symbol, minus named constants, domain symbols and reserved literals, sorted
alphabetically. So `Reduce[x^2 < 1]` and `Reduce[x^2 < 1, Reals]` reach the
same engine as `Reduce[x^2 < 1, {x}, Reals]`. Two heads short-circuit the pipeline: a top-level `Exists`/`ForAll`
routes to `reduce_qe_dispatch` (quantifier elimination), and a `Modulus -> p`
option routes to `reduce_modular` (residue enumeration over `Z/pZ`). Otherwise the
statement is normalised to a DNF `RForm` by `reduce_form_from_expr` and simplified
(`rform_simplify`); the domain is chosen (default `Complexes`, but `Reals` when an
ordering `<`/`<=` is present with no explicit domain, or when a univariate
real-function / multivariate piecewise/radical statement is force-routed and
pre-rewritten by `reduce_realfn_preprocess` / `reduce_piecewise_preprocess`), and
the per-domain dispatch runs.

The dispatch is: **Integers/Rationals** → `reduce_integers`, which reuses the
`Solve` Diophantine engine and reformats its answer as an `Or` of `And`s carrying
`Element[C[k], dom]` for each free parameter. **Complexes** (equations only) → a
single univariate polynomial equation goes to `reduce_eq_univariate` (which carries
the full leading-coefficient case tree) or, when the residual is transcendental, to
`reduce_eq_transcendental` via `Solve`; a linear system to `reduce_eq_system`; and
anything else (nonlinear, or equations plus `!=`) to the zero-dimensional engine
`reduce_zerodim`. **Reals, one variable** → the exact real-algebraic sign-diagram
engines `reduce_univar` / `reduce_univar_general`, with trig/hyperbolic/exponential
equations pre-empted through `Solve` (the sign diagram is unsound on them) and
periodic-equation-plus-bounds shapes handled by `reduce_periodic_region` /
`reduce_trig_ineq_region`. **Reals, several variables** → `reduce_fm`
(Fourier-Motzkin for linear systems), then `reduce_sos`, then `reduce_cad`
(Cylindrical Algebraic Decomposition, McCallum projection), and finally
`reduce_zerodim` for a zero-dimensional system with irrational fibres. A
rational atom `p/q REL 0` with a variable denominator is cleared first to its
sign-exact polynomial form (`p q < 0`, plus `q != 0` guards), in
`reduce_piecewise_preprocess`.

`reduce_sos` (`reduce_sos.c`, v0.322) is a Positivstellensatz emptiness
certificate. For a region `{g_i >= 0, h_j == 0, q < 0}` it looks for a Putinar
representation `q = σ₀ + Σ σ_i g_i` modulo the equality ideal, with each
`σ` a sum of squares. The Gram matrices come from a numeric SDP solved by
alternating projection, through LAPACK. The interior-zero case uses facial
reduction about a rational zero found by minimisation. The matrices are then
rounded to rationals, and the identity and positive semidefiniteness are
re-verified in exact `mpq` arithmetic by a rational LDLᵀ. Only a verified
certificate returns `False`, and any other shape declines quickly. It runs
ahead of the CAD because the high-degree regions it targets are exactly the
ones on which the CAD grinds.

The n-variable CAD wraps `cad_build` in `cad_build_augmented`, which
quantifier elimination shares. When a fibre polynomial nullifies, that
factor's coefficients are added to the projection level below and the stack is
rebuilt, which repairs McCallum well-orientedness. The rounds are bounded, and
a nullification that persists declines. Every CAD loop head polls `CAD_POLL()`,
one clock read per 64 iterations and only inside a `TimeConstrained` scope.
When the deadline has passed the engine returns `NULL` through its ordinary
decline path, freeing every local, so `TimeConstrained[Reduce[...]]` yields
promptly without leaking.

**Data structures.** The working form is an `RForm`: a DNF of `RConj`
conjunctions of `RAtom`s, each a polynomial (an `Expr`) paired with a relation code
(`R_EQ`, `R_LT`, `R_LE`, `R_ELEM`). The CAD (`reduce_cad.c`) keeps per-level
`PolySet` projection stacks and a `CADRegion` cell tree, and delegates every
numeric step — `Discriminant`, `Resultant`, `FactorList`, `Solve[..,Reals]`
real-root isolation, and the exact `qqbar` sign oracle — to the evaluator and
`reduce_real_util`.

**Complexity / limits.** CAD cost is doubly exponential in the number of
variables (McCallum projection then lifting). The 2-variable case is fully
realised, and the n-variable path is an iterated projection (`cad_project_out` for
each level) with recursive lifting and partial-CAD pruning. The SOS stage needs
LAPACK and declines cleanly without it. The governing invariant is **soundness over completeness**: an
undecidable sign or ordering (`qqbar` returning `-2`), an unwired case, or an
irrational section all return `NULL`, leaving the input unevaluated rather than
emitting a wrong or partial formula. Consequently a returned formula describes the
entire solution set exactly, and a literal `False` is a genuine proof of emptiness.
