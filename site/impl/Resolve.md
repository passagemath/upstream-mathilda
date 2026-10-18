---
references:
  - "G. E. Collins and H. Hong, *Partial Cylindrical Algebraic Decomposition for Quantifier Elimination*, J. Symbolic Computation **12** (1991) 299-328."
  - "S. McCallum, *An Improved Projection Operation for Cylindrical Algebraic Decomposition*, in *Quantifier Elimination and Cylindrical Algebraic Decomposition* (Springer, 1998) 242-268."
source: src/solve/reduce_qe.c
---
**Algorithm.** `Resolve` is the dedicated quantifier-elimination front-end.
`builtin_resolve` accepts one or two arguments, requires the first to be a
quantified head (`Exists` / `ForAll`) — otherwise it declines — and calls
`reduce_qe_dispatch(qexpr, dom)`, the same engine `Reduce` invokes when its input
happens to be top-level quantified. The dispatch normalises the quantifier prefix
(`qe_normalize` peels a maximal same-kind chain, records the bound variables, and
carries any 3-argument restriction as a side condition folded in by
`qe_apply_condition`: `cond && M` for `Exists`, `!cond || M` for `ForAll`). An
alternating prefix is eliminated **inner-block-first** by recursion: the inner
different-kind quantifier is resolved to a quantifier-free `psi`, this block's
restriction is folded onto it, and the outer block is then re-eliminated over the
non-alternating result — composing to arbitrary alternation depth.

Elimination of a single block branches on the number of remaining free variables.
Zero free variables is a decision via `qe_decide`, which re-enters `Reduce` on the
matrix (or its negation, for `ForAll`) and reads emptiness, returning `True` or
`False`. One or more free variables is parametric QE via `qe_parametric` →
`reduce_cad_qe`. Before the CAD sees the matrix, `qe_parametric` runs
`reduce_piecewise_preprocess`, the real-domain preprocessing of the base `Reduce`
engine. It case-splits `Abs`, `Max`/`Min`, `Piecewise`, `Sign` and the other
selectors, rationalises square roots, and (since v0.289) clears variable
denominators to sign-exact polynomial atoms. A chained elimination can therefore carry
an intermediate rational atom such as `m <= 1/del` into the next level, which is what
lets ε–δ limit statements, limits at infinity and infinite limits decide. The
matrix then becomes an `RForm` fed to an iterated McCallum projection with the
bound variables projected out first. Its sign formula is emitted over the
free-variable subspace, with free-variable names sorted alphabetically so the answer
is stable. The CAD is built through `cad_build_augmented` (v0.321), the same
well-orientedness augment-and-retry loop as plain `Reduce`: a positive-dimensional
fibre nullification adds coefficients and rebuilds, instead of declining
immediately. A fully quantified sentence therefore resolves to `True`/`False`; a
partially quantified one resolves to a condition on its free variables.

**Data structures.** `Expr` trees; interned-pointer arrays for the bound and free
variable names; the matrix as an `RForm`. The CAD machinery (`PolySet` projection
stacks, the `CADRegion` cell tree, the `qqbar` sign oracle) lives in
`reduce_cad.c`.

**Complexity / limits.** `Reals` is the default and only supported domain — any
other declines, and a first argument that is not quantified declines too
(`Resolve[x > 0]` stays unevaluated). Cost is that of the CAD, doubly
exponential in the variable count. The CAD's loop-head deadline polls make
`TimeConstrained[Resolve[...], t]` return promptly. A transcendental body
(`Exists[x, Sin[x] == a]`) has no polynomial CAD and declines.
By the engine's hard invariant, an undecidable sign, a non-`Reals` domain, an
unsupported construct, or an inner block whose elimination declines all return
`NULL`, leaving the `Resolve[...]` call unevaluated rather than returning a wrong
formula.
