---
source: src/strings/regex/stringposition.c
---
**Algorithm.** `builtin_stringposition` seeds `IgnoreCase` and `Overlaps` from `Options[StringPosition]`, then strips trailing option rules. A single positional argument is the operator form (`src/opform.c`), and an optional third integer `n` keeps the first `n` matches. `sp_scalar` runs the shared `regex_scan_limit`, which `StringCases` and `StringCount` also use. It turns each byte span `[ms, me)` into the 1-based inclusive character pair `{ms + 1, me}`, the form that `StringTake`, `StringDrop` and `StringReplacePart` consume. Only the default differs from its siblings: `Overlaps -> True`, as in the Wolfram Language, so overlapping matches at distinct starts are all listed.

**Data structures.** The `RegexScan` span array, and a `U8View` over the subject that maps byte offsets to character indices. For ASCII the map is the identity, and otherwise it is a binary search over the offset table. Rule sets come from the 64-entry compiled-pattern LRU cache.

**Complexity / limits.** One scan, plus `O(log len)` per reported endpoint on non-ASCII text. Since v0.456 positions count Unicode code points, so they agree with `StringLength`/`StringTake` on accented text. A list of subjects threads, and a wrong argument count emits `StringPosition::argb`.
