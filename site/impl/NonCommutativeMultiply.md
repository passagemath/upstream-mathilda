---
source:
  - src/parse.c
  - src/attr.c
  - src/print.c
  - src/match.c
---
**Semantics.** `NonCommutativeMultiply` is a purely structural head with no C builtin. Its behaviour comes from its attributes, `Flat`, `OneIdentity` and `Protected`, set in the static table in `src/attr.c`. `Flat` makes the evaluator splice nested `**` calls into one flat call. No other rewriting is done: arguments are not reordered (it is not `Orderless`), numbers are not combined, and identities such as `0 ** a` or `1 ** a` are not applied. A one-argument `NonCommutativeMultiply[a]` stays as it is, because `OneIdentity` only affects pattern matching.

**Parsing and printing.** `src/parse.c` lexes `**` as `OP_NONCOMMULT` at precedence 5900. That is tighter than `.` (`Dot`) and `*` (`Times`) and looser than `^`, so `a b ** c^2` is `Times[a, NonCommutativeMultiply[b, Power[c, 2]]]`. Runs of `**` are flattened when the expression is parsed. `print_standard` prints two or more arguments as infix `a ** b ** c`, and other arities in `FullForm` style.

**Rules and limits.** Since the head is `Flat`, `match_flat_leftover` (`src/match.c`) lets a rule written for two operands rewrite part of a longer product and keep the rest. The head is not `Orderless`, so the matched operands must form a leading run (`op[1] ** op[2] ** x` is rewritten, `x ** op[1] ** op[2]` is not). Its `Protected` bit comes from the static attribute table, which `Unprotect` does not override, so new behaviour is attached through UpValues on the operands (`op /: op[i_] ** op[j_] := ...`) rather than DownValues on `**`. It has no packed-array, NDArray or `Compile[]` path; `**` is symbolic only.
