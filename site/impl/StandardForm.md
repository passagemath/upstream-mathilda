---
source:
  - src/print.c
  - src/print_latex.c
  - src/repl.c
---
**Semantics.** `StandardForm` is registered to `builtin_display_form`, which always returns `NULL`, so `StandardForm[expr]` evaluates its argument and keeps the wrapper (attribute `Protected`). It is a display directive, not a transformation.

**Rendering.** The plain-text printer `print_standard` (`src/print.c`) prints a `StandardForm[x]` node as `x` itself, the same way it treats `HoldForm`, except under `InputForm`, where the wrapper is printed so the result reads back in. The typeset renderer `src/print_latex.c` keeps a notation flag `LBuf.std`. A `StandardForm[...]` node sets the flag for its subtree and restores it afterwards, giving the language's own names in upright type (`Sin[x]`, `Log[x]`, `x == 1`, `a && b`, ⅇ and ⅈ for `E` and `I`) with two-dimensional fractions, powers and radicals. In the notebook pipe protocol (`src/repl.c`) each result is sent as `slatex` (StandardForm, the default display) and `latex` (TraditionalForm). `slatex` is selectable: each subexpression is tagged `\htmlData{mp=K}`, and its InputForm text is sent in `sparts`. A result with more than 4096 subexpressions or 256 KB of part text is sent untagged.

**Limits.** In the terminal REPL `StandardForm` changes nothing that is visible, because the text printer has a single notation. `ToString[expr, StandardForm]` is an alias for the InputForm string, and `ToExpression[s, StandardForm]` uses the ordinary parser. There is no `MakeBoxes`/box layer behind it.
