---
references:
  - "The Unicode Consortium, *The Unicode Standard, Version 16.0* (2024), §4.5 General Category, and the Unicode Character Database (UnicodeData.txt)."
source: src/strings/stringcase.c, src/strings/unicode_tables.inc, src/strings/regex/string_pattern.c
---
**Algorithm.** `builtin_uppercaseq` calls `charq_builtin(res, "UpperCaseQ", CQ_UPPER)`, which requires one `String` argument (otherwise `UpperCaseQ::string` via `strings_arg_message`, unevaluated), then decodes the UTF-8 code points one at a time with `u8_decode` and returns `False` at the first one for which `char_ok(CQ_UPPER, cp)` fails — so the scan short-circuits — and `True` if none does (vacuously for `""`). `char_ok` answers ASCII code points inline (`A`–`Z`) and otherwise binary-searches `k_upper_ranges` with `in_ranges`. Inside a string pattern, `predicate_class` in `string_pattern.c` translates `_?UpperCaseQ` to the PCRE class `[[:upper:]]` instead of emitting a per-character callout.

**Data structures.** `k_upper_ranges`: sorted, disjoint code-point ranges of Unicode general category Lu, generated from Unicode 16.0.0 by `tools/gen_unicode_tables.py`.

**Complexity / limits.** `O(L)` for an ASCII string; otherwise `O(L log R)`. Titlecase letters (Lt, such as `ǅ`) are not uppercase; any non-letter makes the result `False`. Characters are code points, not grapheme clusters. Not `Listable`; no packed-array or `Compile[]` surface (the argument is a `String`).
