---
source: src/print.c
---
**Mechanism.** `InputForm` is an inert display wrapper. `builtin_inputform` (`src/print.c`) returns `NULL`, so `InputForm[expr]` stays in the tree (`Head[InputForm[x]]` is `InputForm`), and the printer does the work. When `print_standard` meets the head, it increments `g_inputform_depth` and prints the argument in its linear, re-parseable form. While that depth is non-zero several display-only paths are switched off:

- `SeriesData` prints as the literal `SeriesData[x, x0, {...}, nmin, nmax, den]` rather than `1 + x + … + O[x]^4`.
- `TraditionalForm`/`StandardForm` wrappers are kept, not made transparent.
- `NumberForm`, `Row`, `Graphics` and `Graph` are not given their display rendering.

**Strings.** `print_string_literal` quotes strings and escapes `"`, `\`, newline, tab and carriage return. An explicit `InputForm[...]` keeps the quotes even inside `Print`, whose own traversal otherwise prints strings as bare characters (OutputForm). This was fixed in v0.455. `ToString[expr, InputForm]` and `ToString[expr, StandardForm]` go through the same printer (`expr_to_string`).

**Limits.** InputForm round-trips within Mathilda: `ToExpression[ToString[e, InputForm]] === e`. Two spellings differ from Mathematica's. Products print with `*` before `I` (`3*I + Sqrt[2]`), and machine reals in scientific range print as C-style `1.5e+20`, which Mathilda's parser reads back but Mathematica would not (it writes `1.5*^20`). In the notebook, `InputForm` suppresses typesetting and shows the text as code. `Protected`.
