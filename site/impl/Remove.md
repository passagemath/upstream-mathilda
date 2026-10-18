---
source: src/core.c
---
**Algorithm.** `builtin_remove` (`src/core.c`) shares the `core_apply_symbol_action`
walker with `ClearAll`: it visits each argument — a symbol, a string, or a flat
`List` of them — and applies `core_remove_one` to every resolved name, returning
`Null`. `core_remove_one` skips any `Protected` symbol (the guard
that keeps `Remove` from ever deleting a builtin) and otherwise calls
`symtab_remove_symbol(name)`, deleting the symbol's definition from the symbol
table entirely.

This is a stronger erase than `ClearAll`: where `ClearAll` empties a symbol but
leaves the entry in place, `Remove` deletes the entry, so the name no longer
appears in the symbol table until it is next referenced (at which point a fresh,
undefined symbol is created). `Remove` is itself `Protected` so it cannot
be removed.

**Attributes & limits.** `Attributes[Remove]` is `{HoldAll, Protected}`, and
`Protected` alone is what keeps it, and every other builtin, from being removed.
Its arguments arrive unevaluated, and non-symbol/non-string specs are
ignored. A removed symbol's registered default options go with it, whereas
`Clear` and `ClearAll` keep `Options[s]`.
