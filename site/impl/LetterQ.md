---
references:
  - "The Unicode Consortium, *The Unicode Standard, Version 16.0* (2024), §4.5 General Category, and the Unicode Character Database (UnicodeData.txt)."
source: src/strings/stringcase.c, src/strings/unicode_tables.inc, src/strings/regex/string_pattern.c
---
**Algorithm.** `builtin_letterq` calls `charq_builtin(res, "LetterQ", CQ_LETTER)`, which requires one `String` argument (otherwise `LetterQ::string` via `strings_arg_message`, unevaluated), then decodes the UTF-8 code points one at a time with `u8_decode` and returns `False` at the first one for which `char_ok(CQ_LETTER, cp)` fails — so the scan short-circuits — and `True` if none does (vacuously for `""`). `char_ok` answers ASCII code points inline (`A`–`Z`, `a`–`z`) and otherwise binary-searches `k_letter_ranges` with `in_ranges`. Inside a string pattern, `predicate_class` in `string_pattern.c` translates `_?LetterQ` to the PCRE class `[[:alpha:]]` instead of emitting a per-character callout.

**Data structures.** `k_letter_ranges`: sorted, disjoint `[lo, hi]` code-point ranges covering Unicode general categories `L*` (Lu, Ll, Lt, Lm, Lo), generated from Unicode 16.0.0 by `tools/gen_unicode_tables.py` into `src/strings/unicode_tables.inc`.

**Complexity / limits.** `O(L)` for an ASCII string; otherwise `O(L log R)` over the `R` letter ranges. Combining marks (Mn/Mc) and letter-like numbers (Nl) are not letters. Characters are code points, not grapheme clusters. Not `Listable`; no packed-array or `Compile[]` surface (the argument is a `String`).
