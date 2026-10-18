---
source: src/strings/regex/stringcontainsq.c
---
**Algorithm.** The four substring predicates `StringContainsQ`, `StringFreeQ`, `StringStartsQ` and `StringEndsQ` share one core, `sq_dispatch`, and differ only by an `SqKind` tag. It seeds `IgnoreCase` from the registered `Options`, so `SetOptions` takes effect, and strips trailing option rules. With a single positional argument it returns the curried operator form `Function[head[#1, patt, opts...]]`. Otherwise it builds the rule set unanchored and asks `sq_match`: one `regex_match` per rule, stopping at the first hit, since a predicate never needs to know where or how often. `StringStartsQ`/`StringEndsQ` anchor at the pattern level by wrapping it as `StartOfString ~~ p` or `p ~~ EndOfString`. An empty pattern list `{}` matches nothing and is answered directly.

**Data structures.** A `RegexRule` array from the 64-entry compiled-pattern LRU cache. The cache keys on structural equality, so a freshly built anchoring wrapper still hits it. Only the whole-match pair is requested: no `regex_scan` enumeration and no capture pool.

**Complexity / limits.** One PCRE2 probe per rule, with early exit. The predicate is equivalent to `!StringFreeQ[s, p]` and to `StringMatchQ[s, ___ ~~ p ~~ ___]`. A non-string subject, or any non-string element of a subject list, leaves the whole call unevaluated rather than threading a wrong Boolean. Matching is Unicode-aware, including `IgnoreCase` folding of non-ASCII letters.
