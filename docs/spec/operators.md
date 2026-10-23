# 3. Operators and Precedence

Precedences are numbered on a **1–10000 ladder** (`;` = 100 at the bottom,
`[[…]]` = 10000 at the top). The levels preserve the relative ordering of
the reference's own scheme but are spread out with generous gaps between adjacent
levels, so new operators can be inserted without renumbering. The values below
are the internal numbers used by the parser (`get_operator` in `src/parse.c`)
and, in lockstep, by the printer's parenthesiser (`get_expr_prec` in
`src/print.c`).

| Operator | FullForm | Precedence | Association |
|----------|----------|------------|-------------|
| `[[...]]`| `Part` | 10000 | Left |
| `f[x]`   | `f[x]` | 9500 | Left |
| `_`, `__`, `___` | `Blank` | 8400 | None |
| `!`, `!!` | `Factorial`, `Factorial2` | 8200 | Left |
| `?`      | `PatternTest` | 7900 | None |
| `'`      | `Derivative` | 7700 | Left |
| `++`, `--` | `Increment`, `Decrement` | 7500 | Left |
| `@*`     | `Composition` | 7200 | Left |
| `@`      | `Prefix` | 7000 | Right |
| `@@`     | `Apply` | 7000 | Right |
| `/@`     | `Map`   | 7000 | Right |
| `<>`     | `StringJoin` | 6700 | Left |
| `^`      | `Power` | 6500 | Right |
| `**`     | `NonCommutativeMultiply` | 5900 | Left |
| `.`      | `Dot` | 5300 | Left |
| `/`      | `Divide`| 5000 | Left |
| `*`      | `Times` | 4500 | Left |
| `+`, `-` | `Plus`  | 3500 | Left |
| `+` (prefix) | `Plus` (one argument: `+x` is `Plus[x]`, which evaluates to `x`) | 3500 | Right |
| `==`, `<`, `<=`, `;;` | `Equal`, `Less`, `LessEqual`, `Span` | 3200 | None |
| `===`    | `SameQ` | 3200 | None |
| `!` (prefix) | `Not` | 3000 | Right |
| `&&`     | `And` | 2800 | Left |
| `\|\|`   | `Or` | 2700 | Left |
| `..`, `...` | `Repeated`, `RepeatedNull` (postfix; a single token that wins over a decimal point, so `0..` is `Repeated[0]`) | 2500 | None |
| `\|`     | `Alternatives` | 2300 | Left |
| `~~`     | `StringExpression` | 2100 | Left |
| `:`      | `Pattern` (`sym:p`) / `Optional` (`p:v`) — see below | 2200 / 2150 | Left |
| `/;`     | `Condition` | 1700 | Left (`a /; b /; c` is `Condition[Condition[a, b], c]`) |
| `->`, `→` | `Rule` | 1500 | Right |
| `:>`     | `RuleDelayed` | 1500 | Right |
| `/.`, `//.` | `ReplaceAll`, `ReplaceRepeated` | 1200 | Left |
| `&`      | `Function`| 1000 | Left |
| `//`     | `Postfix` | 800 | Left |
| `=`      | `Set`   | 500 | Right |
| `:=`     | `SetDelayed` | 500 | Right |
| `=.`     | `Unset` (postfix) | 500 | None |
| `+=`, `-=`, `*=`, `/=` | `AddTo`, `SubtractFrom`, `TimesBy`, `DivideBy` | 500 | Right |
| `>>`, `>>>` | `Put`, `PutAppend` | 300 | Left |
| `;`      | `CompoundExpression` | 100 | Left |

(`_`/`__`/`___` blanks are recognised structurally in `parse_primary`, not via
the operator table; the value above is their documented rank.)

**`:` binds its left operand at atom level.** As in the reference, the
left operand of `:` is the single token right before it: a bare symbol gives
`Pattern` (`a:p`), a blank token gives `Optional` (`x_:v`, `_h:v`). It is never
the expression built so far. The right operand is read at the Pattern (2200) or
Optional (2150) rank, so it takes in everything that binds tighter, `|`, `..`,
`==` and `+` included. One further `:` after `sym:p` supplies a default:

| Input | FullForm |
|-------|----------|
| `n_.*a:_g` | `Times[Optional[Pattern[n, Blank[]]], Pattern[a, Blank[g]]]` |
| `n x_:0` | `Times[n, Optional[Pattern[x, Blank[]], 0]]` |
| `x + a:_g == b` | `Plus[x, Pattern[a, Equal[Blank[g], b]]]` |
| `x?f:_` | `PatternTest[x, Pattern[f, Blank[]]]` |
| `x|a:b` | `Alternatives[x, Pattern[a, b]]` |
| `a:_g:1` | `Optional[Pattern[a, Blank[g]], 1]` |
| `x:a:b:c:d:e` | `Optional[Pattern[x, a], Optional[Pattern[b, c], Pattern[d, e]]]` |
| `a:_g -> 1`, `a:_ /; c` | `Rule[Pattern[a, Blank[g]], 1]`, `Condition[Pattern[a, Blank[]], c]` |

Before v0.475 `:` was an ordinary low-precedence infix operator, so `n_.*a:_g`
parsed as `Optional[n_.*a, _g]`, which is a silently wrong pattern (issue #102).
A left operand that the reference rejects (`(a):_g`, `f[x]:0`) is still accepted, as
`Pattern`/`Optional` of the whole left expression.

The Unicode rightwards arrow `→` (U+2192, the reference `\[Rule]`) is accepted as a
synonym for `->`, so rules and associations written for the reference paste and parse
directly.



## Parenthesised comparisons break chains

A run of comparisons folds into one variadic `Inequality` at parse time, so `a < b <= c` becomes
`Inequality[a, Less, b, LessEqual, c]` rather than a nested pair. **Parentheses stop that**, and are
honoured on both operands:

| input | parse |
|---|---|
| `a >= b == c` | `Inequality[a, GreaterEqual, b, Equal, c]` |
| `(a >= b) == c` | `Equal[GreaterEqual[a, b], c]` |
| `(a > b) == (c > d)` | `Equal[Greater[a, b], Greater[c, d]]` |
| `(a < b < c) == d` | `Equal[Inequality[a, Less, b, Less, c], d]` |

This required tracking whether the left operand came from an **unparenthesised** comparison. Deciding
by inspecting the built subtree — which is what the fold did originally — cannot distinguish the two
cases, because by the time the fold looks, `(a >= b)` and a bare `a >= b` are the same node. The result
was that parentheses were ignored on the left operand only: `(2.0 >= 2.) == (1.0 > 0.)` parsed as
`Inequality[2.0, GreaterEqual, 2.0, Equal, Greater[1.0, 0.0]]` and evaluated to `2.0 == True`. The right
operand was never affected, being parsed as its own subexpression, which is why the asymmetry went
unnoticed: every symmetric-looking test passed.

`Plus`/`Times` flattening in the same part of the parser has the identical shape and is harmless, because
those heads are associative — `(a + b) + c` and `a + b + c` mean the same thing. Chaining changes
*meaning*, so comparisons need provenance tracked rather than guessed.

## Statements in a file (`Get`, `-file`)
A file is read one statement at a time, and each statement is evaluated before
the next is parsed. That way a `BeginPackage`/`Begin` prologue governs the
symbols created by later statements.
- A top-level `;` ends a statement, with one exception. If the next token on
  the same line can only continue an expression, the `;` is the
  `CompoundExpression` operator with an implicit `Null`, and that compound is
  the operator's left operand.
- The tokens that continue a statement are the binary and postfix operators:
  `/;`, `//`, `/.`, `//.`, `->`, `:>`, `:=`, `=`, `==`, `===`, `&&`, `||`, `&`,
  `@@`, `/@`, `@`, `^`, `*`, `<`, `<=`, `>`, `>=`, `<>`, `~~`, `|`, `::`, `-=`
  and `+=`.
- So `f[x_] := g[x];/; x > 0` reads as
  `Condition[CompoundExpression[f[x_] := g[x], Null], x > 0]`. Evaluating it
  defines `f` unconditionally.
- An operator at the start of the next line is a syntax error.
- The same rule holds for `;` inside one expression (`ToExpression`, a
  parenthesised body): `(a; b; -> c)` is
  `Rule[CompoundExpression[a, b, Null], c]`. `;` associates left, so the whole
  chain is built before the operator binds, and a parenthesised compound stays
  nested: `(a; b); c` is `CompoundExpression[CompoundExpression[a, b], c]`.
- On a syntax error, `Get` emits `Syntax::sntx` ("Invalid syntax in or before
  "...` (line n of "file")"), stops reading, and returns `Null`. The message
  goes through the funnel, so `Quiet` and `Check` see it. An expression left
  unfinished at end of file counts as a syntax error.
