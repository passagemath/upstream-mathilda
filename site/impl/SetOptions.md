---
source: src/options_builtin.c
---
**Algorithm.** `builtin_setoptions` (`src/options_builtin.c`) redefines
individual default options of a symbol. The first argument must be a symbol. It
takes a working
copy of the symbol's current option rules (from `symtab_get_options`) as a flat
vector and, for each trailing `name -> value` rule, finds the existing option of
that name (context-insensitive match) and **replaces it in place**, preserving
the option's original position.

`SetOptions` cannot *add* an option. A name not already in the symbol's defaults
raises `SetOptions::optnf` through `mth_message`, leaves the stored options
unchanged, and leaves the call unevaluated. A malformed argument returns `NULL`:
a non-rule, a list of rules, or a string in place of the symbol. On success the new
list is stored with `symtab_set_options`, the evaluation cache is dropped
(`eval_clock_bump`, since option changes can alter results), and a copy of the
updated option list is returned.

**Data structures & limits.** Options live on `SymbolDef.default_options` as
`List[Rule[...]]`; the working vector holds `expr_copy`'d rules so nothing
aliased is mutated before the atomic store. `SetOptions` is `Protected`, but it
works on `Protected` targets: built-ins read their defaults from this store at
call time. For example, the regex string functions (`StringCases`,
`StringPosition`, `StringCount`, `StringSplit`, `StringReplace`, `StringMatchQ`,
`StringContainsQ`, …) seed `IgnoreCase`/`Overlaps` from `Options[head]` before
looking at explicit options. `builtin_setoptions` makes no attribute check on
the target symbol.
