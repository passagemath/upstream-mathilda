---
source: src/core.c
---
**Algorithm.** `builtin_tostring` (`src/core.c`) renders its argument to a string through the shared printer in `src/print.c`. It takes 1–2 arguments, and the second must be a form symbol. As in Mathematica, the default form is `OutputForm` (since v0.455; before that it was `InputForm`). A string argument comes back as itself, with no render and no quotes. Anything else goes through `expr_to_string_outputform`, which prints strings anywhere in the tree as their characters (`ToString[{"a", 1}]` is `"{a, 1}"`). The other forms select a printer:

- `InputForm` and `StandardForm` use `expr_to_string`, which keeps string quotes and escapes `"` and `\`, so the result reads back.
- `FullForm` uses `expr_to_string_fullform`.
- `TeXForm` wraps the argument in `TeXForm[...]` and prints that. This sets `C[k]` as `c_k` and infix heads as `\to`, `:\to`, `\land`, …

**Data structures.** The printer writes into a heap buffer, a `tmpfile()` round trip on `_WIN32`. The buffer is copied into a fresh `EXPR_STRING`. Precedence-aware parenthesisation, `Rational`/`Complex` printing and the infix `StringExpression`/`Pattern`/`Repeated` forms (`"a"~~x_`, `x:p`, `p..`) all come from `print_standard`.

**Complexity / limits.** Linear in the size of the printed form. Mathilda's `OutputForm` is one-dimensional (`ToString[1/2]` is `"1/2"`), whereas Mathematica's is a 2-D layout. An unknown form such as `ToString[x, FooForm]` leaves the call unevaluated, so a typo shows at the call site. `Protected`.
