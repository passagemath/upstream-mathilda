# Refine

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Refine[expr, assum]`**

gives the form of expr that would be obtained if the symbols in it were replaced by explicit values satisfying the assumptions assum.

**`Refine[expr]`**

uses the default assumptions specified by any enclosing Assuming constructs ($Assumptions).

<details>
<summary>Notes</summary>

Options: Assumptions (default $Assumptions) -- default assumptions to append to assum. Given as an option value, it prevents Refine from also using $Assumptions. TimeConstraint (default 30) -- how many seconds to spend on any single condition check before giving up on that transformation. Assumptions can be equations, inequalities, domain specifications such as Element\[x, Integers\], or logical combinations of these. Quantities appearing algebraically in inequalities are assumed to be real. Refine can be applied to expressions, and to equations, inequalities, and domain specifications, which it decides to True or False when they provably follow from (or contradict) the assumptions. Refine is one of the transformations tried by Simplify; use Simplify or FullSimplify for a broader search.

</details>

## Examples (15)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (5)

```mathematica
In[1]:= Refine[Sqrt[x^2], x > 0]
Out[1]= x

In[2]:= Refine[Sqrt[x^2], Element[x, Reals]]
Out[2]= Abs[x]

In[3]:= Refine[Sign[x^2 - x y + y^2 + 1], Element[x | y, Reals]]
Out[3]= 1

In[4]:= Refine[a^2 - b^2 + 1 == 0, a + b == 0]
Out[4]= False

In[5]:= Assuming[x > 0, Refine[Sqrt[x^2 y^2], y < 0]]
Out[5]= -x y
```

### Applications (10)

The square root of a square collapses once the sign is known

```mathematica
In[6]:= Refine[Sqrt[x^2], x > 0]
Out[6]= x
```

Radical reduction over a quotient: the denominator is exposed and |a| pulled out

```mathematica
In[7]:= Refine[Sqrt[1 + b^2/a^2], a > 0]
Out[7]= Sqrt[a^2 + b^2]/a
```

Each squared factor resolved by a CAD sign proof

```mathematica
In[8]:= Refine[Sqrt[(a + b)^2 (a - b)^2], a > b > 0]
Out[8]= (a + b) (a - b)
```

Sign of a compound polynomial on an interval

```mathematica
In[9]:= Refine[Sign[x^3 - x], 0 < x < 1]
Out[9]= -1
```

Integer-assumption structural rules shared with Simplify

```mathematica
In[10]:= Refine[{Floor[n], Mod[2 n, 2], (-1)^(2 n), Sin[n Pi]}, Element[n, Integers]]
Out[10]= {n, 0, 1, 0}
```

The argument of a negative real

```mathematica
In[11]:= Refine[Arg[x], x < 0]
Out[11]= Pi
```

A predicate decided by Reduce/CAD entailment

```mathematica
In[12]:= Refine[x^2 + y^2 >= 2 x y, Element[{x, y}, Reals]]
Out[12]= True
```

An equation refuted by the assumption

```mathematica
In[13]:= Refine[x == y, x > y]
Out[13]= False
```

With no positional assumption Refine reads $Assumptions

```mathematica
In[14]:= Assuming[a > 0 && b > 0, Refine[Sqrt[a^2 b^2]]]
Out[14]= a b
```

Undecided: true for x > 1 only, so the predicate is returned unchanged

```mathematica
In[15]:= Refine[x^2 > 1, x > 0]
Out[15]= x^2 > 1
```

## Algorithm

Mathilda -- Refine implementation.

Refine[expr, assum] gives the form of `expr` that would be obtained if the symbols in it were replaced by explicit values satisfying `assum`. It is a thin orchestrator over machinery that already exists and was deliberately shared for it (simp.h reserves AssumeCtx "so future modules (Refine, ...) can share it"):

```text
  1. Argument / option parsing mirrors builtin_possible_zero_q
     (src/zero_test.c): one positional assumption, an Assumptions -> X option
     that overrides $Assumptions, and a TimeConstraint option. The effective
     assumption is And-combined exactly as PossibleZeroQ/Simplify do.

  2. Predicate expressions are decided to True/False:
       Element[x, dom]              -> element_decide (assumption-aware,
                                       including compound-expression domain
                                       inference via prov_re/prov_int).
       Equal / Unequal              -> zero_test_decide_assuming on lhs - rhs,
                                       with a Reduce fallback.
       Less / ... / logic           -> Reduce/CAD entailment over the reals:
                                       P is True  iff Reduce[A && !P] is False,
                                       P is False iff Reduce[A &&  P] is False.

  3. Every other expression is rewritten by apply_assumption_rules
     (Sqrt[x^2]->x/-x/Abs[x], Log[x^p]->p Log[x], integer-k trig, Abs/Sign/
     Conjugate under sign facts, Floor/Ceiling/Mod under integer/interval
     facts, ...) followed by a deep-positivity post-pass that resolves
     Sign[p]/Abs[p]/Sqrt[p^2] for symbolic polynomials p that the fast
     sign prover cannot settle, using the same Reduce entailment.
```

TimeConstraint (default 30s) is enforced with Simplify's cooperative wall-clock budget (simp_set_time_budget / simp_mono_seconds): Reduce/CAD has no interior abort hook, so the deadline is checked BEFORE each entailment call and a running Reduce is never preempted (best-effort, and never the malloc-lock-prone async TimeConstrained[]).

Refine is a purely symbolic/structural head: it has no element-wise numeric semantics, so it deliberately carries no NDArray/packed kernel and no Compile[] lowering.

Ownership follows the builtin contract: return a new tree or steal from

```text
`res`; never expr_free(res).
```

## Implementation notes

**Algorithm.** `builtin_refine` is a thin orchestrator over the assumption engine
`Simplify` already uses (`AssumeCtx`, `apply_assumption_rules`, declared in `simp.h`,
which reserves that type "so future modules (Refine, …) can share it"). It first splits
trailing options from positional arguments and forms the *effective* assumption with the
identical policy to `PossibleZeroQ`/`Simplify`: an `Assumptions -> X` option replaces
`$Assumptions`, a positional assumption is `And`-conjoined with `$Assumptions`
(`read_dollar_assumptions`), and the result is evaluated so `And[True, p]` canonicalises to
`p`. That expression is parsed into an `AssumeCtx` fact set by `assume_ctx_from_expr`. From
there two paths run. **(1) Predicate head** (`Element`, or a relational/logical head —
`Equal`, `Less`, `And`, …): `refine_decide_predicate` tries to settle it to `True`/`False`.
`Element[x, dom]` goes to the assumption-aware `element_decide`; `Equal`/`Unequal` first
run `zero_test_decide_assuming` on `lhs − rhs` (and an equality-substitution via
`apply_assumption_rules`), then fall back to `Reduce` over the complexes; inequalities and
logic go to a `Reduce`/CAD entailment over the reals, where `P` is `True` iff `Reduce[A &&
!P]` is unsatisfiable (literally `False`) and `False` iff `Reduce[A && P]` is. **(2) Rewrite
path** (everything else): `apply_assumption_rules` rewrites the expression under the facts
(the very rules `Simplify` applies — `Sqrt[x^2]->±x/Abs[x]`, `Log[x^p]->p Log[x]`, integer-`k`
trig, `Abs`/`Sign`/`Conjugate`/`Arg` under sign facts, and the structural rules in
`assume_structural_rewrite`). Those structural rules cover the integer-assumption
reductions (`Floor`/`Mod`/`(-1)^k`/`I^k` parity, `Sin[k Pi]`), the inverse-trig branch
formulas such as `ArcSin[Sin[x]] -> Pi - x` on `(Pi/2, Pi)`, and `Arg` of a canonicalised
exponential. Next, a bottom-up `deep_positivity_walk` resolves `Sign[p]`, `Abs[p]` and
`Sqrt[p^2]` for compound polynomials `p` that the fast sign prover could not settle,
each through a `deep_sign` CAD run. Since v0.379 the walk also reduces radicals over
products and quotients: `sqrt_extract` splits
`Sqrt[Π base_i^(2 k_i) · R] -> (Π |base_i|^k_i) Sqrt[R]`, with each `|base_i|` resolved by
`deep_sign`. A `Plus` radicand is passed through `Together` first to expose its
denominator, and the `^(-1/2)` form is handled too, so `Sqrt[1 + b^2/a^2]` becomes
`Sqrt[a^2 + b^2]/a` under `a > 0`. The rule is sound because every factor pulled out is
non-negative. `Integrate` calls this pass on a definite result whenever it was given
`Assumptions`. The tree is then evaluated. With no usable facts the call is the
identity: it takes the positional expression out of `res` unchanged.

**Data structures.** Everything is `Expr` trees. Assumptions are held as an `AssumeCtx`
(flat `Expr*` fact array) borrowed from `simp`. Entailment is threaded through a
`RefineBudget { double deadline; int calls_left; }`; each `Reduce` query is built as
`Reduce[stmt, varlist, domain]`, where `collect_bare_vars` gathers the distinct bare-symbol
leaves (skipping protected real constants such as `Pi`, and *not* reusing the polynomial
`collect_variables`, which would atomise `x^p`/`Floor[x]` into pseudo-variables `Reduce`
rejects).

**Complexity / limits.** Dominated by the `Reduce`/CAD entailment, which is doubly
exponential in the number of variables; it is bounded two ways — the variable set is capped
at 6 (wider statements decline) and at most `REFINE_MAX_ENTAILMENT_CALLS` (24) `Reduce`
calls run per `Refine`. `TimeConstraint` (default 30 s) is a cooperative wall-clock budget
checked *before* each entailment via `simp_mono_seconds`; a running `Reduce` is never
preempted (CAD has no interior abort hook), and the async, malloc-lock-prone
`TimeConstrained[]` is deliberately avoided. Predicates it cannot decide, and expressions no
rule rewrites, pass through unchanged. `Refine` is purely symbolic/structural: it carries no
NDArray/packed kernel and no `Compile[]` lowering, by design, since it returns symbolic
expressions rather than machine numbers.

- `Protected`. Shares the assumption engine with `Simplify` (`AssumeCtx`,
  `apply_assumption_rules`): every rewrite `Refine` applies is one `Simplify`
  also applies. `Refine` is the rewrite-and-decide pass *without* `Simplify`'s
  complexity-minimising search, so it does not, e.g., factor.
- Assumption-driven rewrites: `Sqrt[x^2] -> x / -x / Abs[x]` (for `x > 0` /
  `x < 0` or `x <= 0` / real), `Abs[x] -> -x` for `x <= 0`,
  `(x^m)^r -> x^(m r)` (for `x >= 0`), `(a^b)^c -> a^(b c)` (for `-1 < b < 1`),
  `a^p b^p -> (a b)^p` (for `a, b > 0`), `Sign[x] -> ±1` and `Arg[x] -> 0 / Pi`
  under sign facts, `Log[x] -> I Pi + Log[-x]` (`x < 0`), `Log[x^p] -> p Log[x]`
  (`x > 0`), `Log[x^2] -> 2 Log[Abs[x]]` and `Log[E^x] -> x` (real `x`),
  `Log[x rest] -> Log[x] + Log[rest]` (positive `x`), `Sin[k Pi] -> 0` and
  `Cos[x + k Pi] -> (-1)^k Cos[x]` (integer `k`), `ArcTan[Tan[x]] -> x` on the
  principal domain, `Re`/`Im`/`Conjugate`/`Arg`/`Abs` of real-symbol expressions
  (e.g. `Abs[a + b I] -> Sqrt[a^2 + b^2]`),
  `Floor`/`Ceiling`/`Round`/`IntegerPart`/`FractionalPart` and `Mod` under
  integer / interval / modular facts. The per-symbol rule synthesis is pruned to
  the symbols the target actually uses, so a large assumption set does not
  overflow the rule buffer.
- Predicate decisions: `Element[x, dom]` via the assumption-aware domain
  prover (including compound expressions, and the sign domains
  `Positive`/`Negative`/`NonNegative`/`NonPositive`); equations via the
  assumption-aware zero test **plus** equality-substitution
  (`Refine[a == b, a - b == 0] -> True`); inequalities and their logical
  combinations via the `Reduce`/CAD entailment (`P` is `True` iff `assum && !P`
  is unsatisfiable over the reals). Quantities appearing algebraically in
  inequalities are assumed real. The zero test is sound under coupling equality
  assumptions: it never reports a genuine identity as non-zero (it downgrades to
  undecided rather than sampling points the assumptions exclude).
- Purely symbolic/structural: no packed/NDArray kernel and no `Compile[]`
  lowering (it returns symbolic expressions, not machine numbers).

**Attributes:** `Protected`.

## References

**See also:** [Assuming](../../simplification/Assuming/), [$Assumptions](../../simplification/$Assumptions/), [Simplify](../../simplification/Simplify/), [PossibleZeroQ](../../expression-information/PossibleZeroQ/), [Reduce](../../solutions-of-equations/Reduce/), [Re](../../arithmetic/Re/), [Im](../../arithmetic/Im/), [Conjugate](../../arithmetic/Conjugate/)

- Source: [`src/refine.c`](https://github.com/stblake/mathilda/blob/main/src/refine.c)
- Specification: [`docs/spec/builtins/simplification.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/simplification.md)
- Tests: [`tests/test_possiblezeroq_assumptions.c`](https://github.com/stblake/mathilda/blob/main/tests/test_possiblezeroq_assumptions.c)
- Tests: [`tests/test_refine.c`](https://github.com/stblake/mathilda/blob/main/tests/test_refine.c)

## Notes & additional examples

### Notes

- `Refine[expr, assum]` gives the form `expr` takes when its symbols are
  given values that satisfy `assum`. It shares the assumption engine of
  `Simplify` (`AssumeCtx`, `apply_assumption_rules`,
  `assume_structural_rewrite`), so every rewrite `Refine` makes is one
  `Simplify` can also make. `Refine` is the rewrite-and-decide pass without
  `Simplify`'s complexity search: it turns `Sqrt[x^2]` into `x` under
  `x > 0`, but it does not factor or expand. Conversely, a correct rewrite
  that `Simplify`'s scorer would reject (`ArcSin[Sin[x]] -> Pi - x` on
  `Pi/2 < x < Pi`) still comes through `Refine`.
- **Predicates.** `Element[...]`, equations, inequalities and logical
  combinations are decided to `True` or `False` when the assumptions prove
  or refute them. Domain membership goes to the assumption-aware prover,
  equations to the assumption-aware zero test plus equality substitution,
  and inequalities to a `Reduce`/CAD entailment over the reals: `P` is
  `True` exactly when `assum && !P` is unsatisfiable. A predicate that is
  neither proved nor refuted comes back unchanged.
- **Radicals.** Under the sign facts, a radical of a product or quotient of
  squares loses each squared factor (as its absolute value, then resolved).
  `Integrate` uses this pass to tidy a definite integral computed with
  `Assumptions`.
- **Defaults and budget.** `Refine[expr]` uses `$Assumptions`, including
  whatever an enclosing `Assuming` adds. `Assumptions -> X` replaces
  `$Assumptions`, and a positional assumption is conjoined with it.
  `TimeConstraint` (default 30 s) caps the time spent on a condition check.
  At most 24 `Reduce` entailments run per call, over at most 6 variables.
- `Protected`. Symbolic only, with no NDArray or `Compile[]` path.
- See also: `Simplify`, `Assuming`, `Reduce`, `PossibleZeroQ`.
