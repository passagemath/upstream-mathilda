---
references:
  - "D. Withoff, *Mathematica Internals* (Wolfram Research technical report, 1992), §3.1 (evaluation order: upvalues before downvalues)."
source:
  - src/eval.c
  - src/symtab.c
  - src/parse.c
---
**Algorithm.** `f /: lhs = rhs` is a ternary operator (`OP_TAGSET`, precedence 500, a dedicated led branch in `src/parse.c`) producing `TagSet[f, lhs, rhs]`. `TagSet` is `HoldAll`; the evaluator's step 6c evaluates the arguments of `lhs` with `build_assignment_target` (as `Set` does), evaluates `rhs` once, and calls `apply_tag_assignment`. That rejects a non-symbol tag (`TagSet::sym`) and a `Protected` tag (`TagSet::write`), then `classify_tag_position` decides where the rule belongs: `lhs` *is* the tag → OwnValue (`symtab_add_own_value`); the tag is the head of `lhs` → DownValue (`symtab_add_down_value`); the tag is the dispatch head of some level-one argument (`element_dispatch_head`: the argument symbol, its head, or a blank's head constraint) → UpValue (`symtab_add_up_value`). A top-level `lhs /; test` is classified on its inner `lhs`. The tag as the head of the head (`f[x][y]`) is a SubValue and is declined with `TagSet::tagpos`; a tag absent from `lhs` gives `TagSet::tagnf`. In every case the call returns `rhs`, matching how `Set` returns its right side on a refused assignment.

**Data structures.** The rule lands on exactly one list of the tag's `SymbolDef` — `own_values`, `down_values` or `up_values` — through the shared `add_rule` insert (specificity-sorted, same-LHS replacement in place). Unlike `UpSet`, which attaches to *every* level-one symbol, `TagSet` attaches only to the named tag, so `a /: fz[a, b] = 1` leaves `b` untouched.

**Complexity / limits.** Classification is one pass over the level-one arguments; installation is a list insert. Firing is the general UpValue scan (`apply_up_values`, before DownValues and builtins). SubValues are not supported. Because the tag's own Protected bit is the only one checked, an UpValue can be attached to a Protected *outer* head (`v /: Length[v] = 42`) without unprotecting it — the intended use. Compiled in under `UP_VALUES`.
