---
source:
  - src/eval.c
  - src/symtab.c
  - src/parse.c
---
**Algorithm.** `f /: lhs =.` parses (the `OP_TAGSET` branch seeing `=.`) to `TagUnset[f, lhs]`, which is `HoldAll`. The evaluator's step 6d evaluates the arguments of `lhs` with `build_assignment_target` (so the key is built the same way the original assignment built it) and calls `apply_tag_unset`. A non-symbol tag reports `TagUnset::sym`. Otherwise `classify_tag_position` picks the list: OwnValue or DownValue removal goes through `symtab_remove_matching_rule`, UpValue removal through `symtab_remove_matching_up_value`. A SubValue position or a tag not found in `lhs` removes nothing. The call always returns `Null`.

**Data structures.** Matching is by *key*, not by pattern match: both the requested `lhs` and each stored pattern are run through `pattern_canonicalize` and `pattern_alpha_normalize`, then compared with `expr_eq`. Renaming a pattern variable therefore still finds the rule (`f[h[y_]]` removes a rule stored as `f[h[x_]]`), but a stored `Condition` is part of the key, so `lhs /; test` must be repeated to remove a conditioned rule. The first equal node is unlinked and freed, `symtab_up_value_count` is decremented and the evaluation clock bumped.

**Complexity / limits.** Linear in the length of the tag's rule list, with one alpha-normalisation per stored rule. When no rule matches the removal is silent (Mathematica issues `TagUnset::norep`). Plain `lhs =.` (`Unset`) does not reach an UpValue — the tag is required, as in Mathematica. Compiled in under `UP_VALUES`.
