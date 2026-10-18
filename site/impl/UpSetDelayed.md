---
references:
  - "D. Withoff, *Mathematica Internals* (Wolfram Research technical report, 1992), §3.1 (evaluation order: upvalues before downvalues)."
source:
  - src/eval.c
  - src/symtab.c
  - src/parse.c
---
**Algorithm.** `lhs ^:= rhs` parses (`OP_UPSETDELAYED`, precedence 500; the three-character `^:=` token is lexed before bare `^`) to `UpSetDelayed[lhs, rhs]`. `UpSetDelayed` is `HoldAll`, so `rhs` reaches the evaluator's step 6b unevaluated and is stored as written. The path is shared with `UpSet`: `build_assignment_target` evaluates the non-pattern arguments of `lhs`, `collect_level_one_symbols` / `element_dispatch_head` find the distinct level-one tags (argument symbols, argument heads, and the head constraint of a blank such as `a_md`), and `apply_up_assignment` installs the rule on each with `symtab_add_up_value`. One delayed-only step runs first: `move_rhs_condition` rewrites `lhs ^:= body /; test` to the pattern `Condition[lhs, test]` with replacement `body`, exactly as `SetDelayed` does, so a failing test makes the rule not match rather than return a `Condition`. The call returns `Null`. Misuse reports `UpSetDelayed::nosym` (no level-one symbol; the call stays unevaluated) or `UpSetDelayed::write` (a Protected tag, skipped).

**Data structures.** Identical to `UpSet`: a `Rule` on each tag's `SymbolDef.up_values`, ordered by `add_rule`'s specificity sort, counted in `symtab_up_value_count`, with an epoch bump on install. Because the stored replacement is the held `rhs`, `replace_bindings` substitutes the match and the evaluator re-evaluates it on every firing — the difference from `UpSet` is purely when `rhs` is evaluated.

**Complexity / limits.** Installation is linear in the number of level-one arguments plus the list insertion. Firing costs the same scan as any UpValue (see `apply_up_values`): one branch when none exist, otherwise the level-one candidates' rule lists. Mathilda's `Rule` carries no immediate/delayed flag, so `UpValues[]` and `Definition[]` show both kinds with `:>` / `:=`. Compiled in under `UP_VALUES`.
