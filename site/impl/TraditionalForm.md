---
source:
  - src/print.c
  - src/print_latex.c
  - src/repl.c
---
**Semantics.** `TraditionalForm` shares `builtin_display_form` with `StandardForm`. The builtin always returns `NULL`, so the argument is evaluated and the wrapper is kept (attribute `Protected`). Because the head survives, it stays inside later arithmetic: `TraditionalForm[x^2] + 1` is `Plus[1, TraditionalForm[x^2]]`.

**Rendering.** The text printer (`print_standard`) skips the wrapper and prints its argument, except under `InputForm`, where it is kept. The typeset renderer `src/print_latex.c` (`expr_to_latex_form`) clears the `LBuf.std` flag for the wrapped subtree. That gives textbook notation: `\sin(x)`, `\ln(x)`, `|x|`, `e^{x}`, `x = 1`, `a \land b`, n-th roots as surds, Greek-named symbols as letters, and `C[1]` as c₁. Nesting switches notation for that subtree only, and a `StandardForm[...]` inside switches it back. In the notebook protocol (`src/repl.c`) a result whose outermost head is `TraditionalForm` is sent with `slatex` alone, typeset in that notation. The notebook's *Convert To* menu and its *View › Output Format* setting choose between the forms on the client side, without re-evaluating.

**Limits.** The terminal REPL has no typeset output, so `TraditionalForm` prints like the bare expression there. `TeXForm` uses a separate renderer (`print_tex` in `src/print.c`) and does not read the notation flag. `ToString[expr, TraditionalForm]` is not supported and stays unevaluated. There is no box layer, so `TraditionalForm` output cannot be parsed back.
