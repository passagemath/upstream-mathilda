---
source: src/strings/stringtrim.c
---
**Algorithm.** As in Mathematica, `StringTrim[s, p]` is `StringReplace[s, {StartOfString ~~ p -> "", p ~~ EndOfString -> ""}]`: it removes one match at each end and does not strip to a fixed point. `builtin_stringtrim` takes the subject, an optional pattern (default `Whitespace`, i.e. `\s+`) and `IgnoreCase`. It builds the pattern twice through the shared engine (`regex_rules_build_opts`): once compiled `REGEX_ANCHORED` for the front and once `REGEX_ENDANCHORED` for the back. `st_scalar` runs one anchored match at offset 0 to find the front cut. It then searches for the back match from the end of the front match, which PCRE rejects as soon as a candidate cannot reach the end, so the two cuts never cross (`StringTrim["xxx", "xx"]` is `"x"`).

**Data structures.** Two cached `RegexRule` sets (front and back) and one output buffer for the window `[start, end)`. The full string-pattern language works here, `Condition` and `PatternTest` included.

**Complexity / limits.** Two PCRE2 matches per subject, and the front one never scans. Strip a run with a repeated pattern such as `"x"..`: a single-character pattern removes only one character per end. A list of subjects threads. A non-string element in it emits `StringTrim::strse` and leaves the call unevaluated. A build without PCRE2 warns with `regavail`.
