---
source:
  - src/core.c
  - src/symtab.c
  - src/eval.c
---
**Algorithm.** `builtin_up_values` (attributes `HoldAll`, `Protected`) accepts one argument. A symbol reads its list with `symtab_get_up_values`; a string is resolved with the non-creating `symtab_lookup`, so `UpValues /@ Names["s*"]` works, and a string naming no existing symbol issues `UpValues::sym` and stays unevaluated. Any other argument also stays unevaluated. The list is converted by the shared `rules_to_list`, which emits each stored rule as `HoldPattern[lhs] :> rhs`. The assignment `UpValues[s] = {rules}` is handled in the evaluator's `Set` path (`src/eval.c`): after a `Protected` check (`Set::write`), `symtab_set_up_values` frees the old list and installs each `Rule`/`RuleDelayed` element, unwrapping an outer `HoldPattern` and skipping non-rule elements.

**Data structures.** The source is `SymbolDef.up_values`, a singly linked `Rule` list maintained by `add_rule` in specificity order (insertion order breaks ties). Each node is deep-copied into the result, so the returned rules are inert and editing them does not touch the store. A `Condition` lifted from a delayed right side appears inside the `HoldPattern`.

**Complexity / limits.** Linear in the number of stored rules. Immediate (`^=`) and delayed (`^:=`) rules are both returned with `:>` because the rule store does not record which kind created them. `Clear[s]` and `ClearAll[s]` also empty this list. Compiled in under `UP_VALUES`.
