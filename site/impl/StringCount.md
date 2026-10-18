---
source: src/strings/regex/stringcount.c
---
**Algorithm.** `builtin_stringcount` reads its options and builds its rule set exactly as `StringCases` does, seeding `Overlaps -> False` and `IgnoreCase -> False` from `Options[StringCount]` and stripping trailing option rules. Then `sct_scalar` calls `regex_scan` with `want_captures = 0` and returns the number of spans as an `Integer`. No substring is ever created, so each match costs one small span record instead of a string `Expr` and a list slot. A `Rule`/`RuleDelayed` pattern is accepted, but only its left-hand side matters to a count. A list of subjects threads, giving one count each.

**Data structures.** The shared `RegexScan`: a span array and no capture pool. For `Overlaps -> False` the scan uses the `RegexCursor` leftmost-match cursor. Rule sets come from the 64-entry compiled-pattern LRU cache.

**Complexity / limits.** Because the scan is shared, `StringCount[s, p, opts]` always equals `Length[StringCases[s, p, opts]]` and `Length[StringPosition[s, p, opts]]` under the same `Overlaps` setting. Note that `StringPosition` defaults to `Overlaps -> True`. `Overlaps -> All` counts every matching substring at every start, testing `O(len²)` substrings. A positional argument count other than two emits `StringCount::argrx`. Matching is Unicode-aware (PCRE2 UTF + UCP).
