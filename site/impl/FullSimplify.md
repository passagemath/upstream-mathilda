---
source:
  - src/internal/simp/FullSimplify.m
  - src/internal/simp/transforms/gamma.m
  - src/internal/simp/transforms/logexp.m
  - src/internal/simp/transforms/erf.m
  - src/internal/simp/transforms/polylog.m
  - src/internal/simp/transforms/powerradical.m
---
**Algorithm.** `FullSimplify` is written in the Mathilda language as a
package (`BeginPackage["FullSimplify`"]`, helpers in
``FullSimplify`Private` ``), loaded eagerly at startup. It is a wrapper around
`Simplify`. `Core[expr, opts]` splits trailing `Rule` options from the single
positional assumption. A relevance engine, `Relevant[expr]`, then collects
every symbol that appears in `expr` as a head or atom
(`Cases[expr, _Symbol, {0, Infinity}, Heads -> True]`). For each one it calls
`Ensure[h]`, which consults the manifest `ModuleFile[h]` and `LoadModule`s
that family's identity library the first time the head is seen. It returns
the union of the unary transformation functions the libraries registered
through `RegisterTransforms[h, {f...}]`. The call that follows is
`Simplify[expr, assum, TransformationFunctions -> {Automatic, relevant...,
user...}, ComplexityFunction -> cf]`. Keeping `Automatic` preserves the
whole built-in Simplify pipeline, so the result is never more complex than
`Simplify`'s. A final `Consolidate` sweep applies each relevant identity
once more to the Simplify result and keeps the first rewrite whose score
(`SimplifyCount` or the user `ComplexityFunction`) does not rise. That is
how equal-complexity repackagings such as `x Gamma[x] -> Gamma[1 + x]`
survive the strict-improvement search.

**Identity libraries.** `gamma.m`, triggered by `Gamma`, `LogGamma`,
`PolyGamma`, `Pochhammer`, `Beta` and `Factorial`, holds the `+1`
recurrence, its raising form `z Gamma[z] -> Gamma[z+1]`, and three pair
rules that fire only when both gammas sit in one `Times`. Reflection
`Gamma[z] Gamma[1-z] -> Pi/Sin[Pi z]` is guarded by `!IntegerQ[z]`. The
conjugate pairs `Gamma[1±I b]` and `Gamma[1/2±I b]` are matched through
`PossibleZeroQ[u - Conjugate[v]]`. `Pochhammer`, `Beta` and `Factorial`
are rewritten to Gamma. `erf.m` gives `Erf + Erfc -> 1`; `polylog.m` gives
the dilogarithm duplication `PolyLog[2,z] + PolyLog[2,-z] -> PolyLog[2,z^2]/2`;
`powerradical.m` gives `Surd[x,n]^n -> x`. `logexp.m` (v0.340), triggered
by `Log`, is the real-log antiderivative (Gudermannian) family:
`Log[Sec ± Tan] -> ±ArcTanh[Sin]`, `Log[Csc ± Cot] -> ±ArcTanh[Cos]`,
`Log[Cosh ± Sinh] -> ±x`, `Log[Coth ± Csch] -> ±ArcCoth[Cosh]`, and the
half-angle forms. Each two-term rule comes in a factored version and a
distributed version, so a common coefficient is split off as `Log[a]`, and
each requires `TrueQ[Positive[a]]`.

**Options and limits.** `TimeConstraint -> {tLoc, tTot}` (a bare `t` means
`{t, t}`) wraps each transform in `TimeConstrained[f[e], tLoc, e]`. It also
wraps the whole Simplify call in `TimeConstrained`, whose fallback is the plain
`Simplify` result, so a timeout never returns something worse than
`Simplify`. Identities are Mathilda rewrite rules, not algorithms. The gamma
recurrence matches a literal `+1` shift only, so `Gamma[x+2]/Gamma[x]` stays
as it is. An expression with none of the manifest heads costs exactly one
`Simplify`.
The log rules check only the coefficient, not the sign of the log argument,
so they are valid where that argument is positive (for example
`Sec[x] + Tan[x] > 0` on `-Pi/2 < x < Pi/2`). FullSimplify is purely
symbolic, with no NDArray or `Compile[]` path.
