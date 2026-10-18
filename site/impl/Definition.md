---
source:
  - src/print.c
  - src/core.c
---
**Algorithm.** `Definition` has no evaluator builtin: `Definition[s]` (attributes `HoldAll`, `Protected`) stays as an inert expression, and the printer renders it. In `print_standard`, a `Definition[sym]` node calls `print_definition_body`, which reads `symtab_get_own_values`, `symtab_get_down_values` and `symtab_get_up_values` and prints one line per rule: OwnValues as `lhs = rhs`, DownValues as `lhs := rhs`, UpValues as `s /: lhs := rhs`, in that order and in each list's stored order. A symbol with no values prints `Null`. `render_definitions` in `src/core.c` builds the same text into a string for `?s` / `Information[s]`, which appends it after the usage message.

**Data structures.** It reads the three `Rule` lists on the `SymbolDef` directly and does not copy them. The patterns are printed with the ordinary `print_standard`, so a conditioned rule shows its `Condition` on the left (`lhs /; test := rhs`). DownValues come out in `add_rule`'s specificity order, so a literal `fd[0]` rule appears before the general `fd[x_]` one.

**Complexity / limits.** Linear in the total number of rules. Display only: `FullForm[Definition[s]]` is `Definition[s]`. The rule store has no immediate/delayed flag, so every down- and up-value prints with `:=` and every own-value with `=`. Attributes, Options, messages, NValues and SubValues are not printed (Mathematica includes `Attributes[s] = {...}` lines). Compiled in under `UP_VALUES`.
