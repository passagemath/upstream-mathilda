---
references:
  - "The Unicode Consortium, *The Unicode Standard, Version 16.0* (2024), §4.5 General Category, and the Unicode Character Database (UnicodeData.txt)."
source: src/strings/stringcase.c, src/strings/unicode_tables.inc, src/strings/regex/string_pattern.c
---
**Algorithm.** `builtin_digitq` calls `charq_builtin(res, "DigitQ", CQ_DIGIT)`, which requires one `String` argument (otherwise `DigitQ::string` via `strings_arg_message`, unevaluated), then decodes the UTF-8 code points one at a time with `u8_decode` and returns `False` at the first one for which `char_ok(CQ_DIGIT, cp)` fails — so the scan short-circuits — and `True` if none does (vacuously for `""`). `char_ok` answers ASCII code points inline (`0`–`9`) and returns 0 for every non-ASCII code point without consulting any table — DigitQ is ASCII-only, as in Mathematica. Inside a string pattern, `predicate_class` in `string_pattern.c` translates `_?DigitQ` to the PCRE class `[[:digit:]]` instead of emitting a per-character callout.

**Data structures.** None beyond the decoder: the test is a pure range check on the code point.

**Complexity / limits.** `O(L)` for an ASCII string; any non-ASCII character fails immediately. Note the asymmetry with the string pattern `DigitCharacter`, which (compiled with PCRE2 UCP) matches every Unicode decimal digit. Characters are code points, not grapheme clusters. Not `Listable`; no packed-array or `Compile[]` surface (the argument is a `String`).
