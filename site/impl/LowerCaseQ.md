---
references:
  - "The Unicode Consortium, *The Unicode Standard, Version 16.0* (2024), §4.5 General Category, and the Unicode Character Database (UnicodeData.txt)."
source: src/strings/stringcase.c, src/strings/unicode_tables.inc, src/strings/regex/string_pattern.c
---
**Algorithm.** `builtin_lowercaseq` calls `charq_builtin(res, "LowerCaseQ", CQ_LOWER)`, which requires one `String` argument (otherwise `LowerCaseQ::string` via `strings_arg_message`, unevaluated), then decodes the UTF-8 code points one at a time with `u8_decode` and returns `False` at the first one for which `char_ok(CQ_LOWER, cp)` fails — so the scan short-circuits — and `True` if none does (vacuously for `""`). `char_ok` answers ASCII code points inline (`a`–`z`) and otherwise binary-searches `k_lower_ranges` with `in_ranges`. Inside a string pattern, `predicate_class` in `string_pattern.c` translates `_?LowerCaseQ` to the PCRE class `[[:lower:]]` instead of emitting a per-character callout.

**Data structures.** `k_lower_ranges`: sorted, disjoint code-point ranges of Unicode general category Ll, generated from Unicode 16.0.0 by `tools/gen_unicode_tables.py`.

**Complexity / limits.** `O(L)` for an ASCII string; otherwise `O(L log R)`. Any non-letter (digit, space, punctuation) makes the result `False`. Characters are code points, not grapheme clusters. Not `Listable`; no packed-array or `Compile[]` surface (the argument is a `String`).
