---
references:
  - "D. Withoff, *Mathematica Internals* (Wolfram Research technical report, 1992), §3.1 (evaluation order: upvalues before downvalues)."
source:
  - src/eval.c
  - src/symtab.c
  - src/parse.c
---
**Algorithm.** `f /: lhs := rhs` parses through the same `OP_TAGSET` led branch as `TagSet`, which inspects the following operator (`=`, `:=` or `=.`) to build `TagSetDelayed[f, lhs, rhs]`. `TagSetDelayed` is `HoldAll` and its `rhs` is never evaluated at assignment time. Step 6c of the evaluator evaluates the non-pattern arguments of `lhs` (`build_assignment_target`) and calls `apply_tag_assignment` with `is_delayed` set: the tag is validated (`TagSetDelayed::sym`, `TagSetDelayed::write`), `classify_tag_position` places the rule as an own-, down- or up-value of the tag, and `move_rhs_condition` first lifts a trailing `body /; test` onto the pattern as `Condition[lhs, test]`, so a failing test makes the rule not apply. `TagSetDelayed::tagnf` (tag absent) and `TagSetDelayed::tagpos` (SubValue position) refuse the install. The call returns `Null`.

**Data structures.** One `Rule` on the tag's `own_values` / `down_values` / `up_values` list, inserted by `add_rule`; the replacement is the held `rhs`, re-evaluated after `replace_bindings` on each firing. For UpValues on `Orderless` heads such as `Plus`/`Times`, `apply_up_values` sorts the call's arguments canonically before matching, and on `Flat` heads `match_flat_leftover` lets a rule for two operands rewrite a pair inside a longer sum and keep the rest.

**Complexity / limits.** Same costs as `TagSet`: a single classification pass and a list insert, then the general UpValue scan at call time. For a non-`Orderless` `Flat` head (`Dot`, `**`) the matched operands must form a leading run. `Definition[f]` shows the rule in `f /: lhs := rhs` form (a lifted condition appears as `lhs /; test`). Compiled in under `UP_VALUES`.
