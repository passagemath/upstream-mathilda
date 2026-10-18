---
references:
  - "The Unicode Consortium, *The Unicode Standard, Version 16.0* (2024), §3.13 Default Case Algorithms, and the Unicode Character Database (UnicodeData.txt, SpecialCasing.txt)."
source: src/strings/stringcase.c, src/strings/unicode_tables.inc, tools/gen_unicode_tables.py
---
**Algorithm.** `builtin_touppercase` calls `case_builtin`, which accepts exactly one `String` argument (the evaluator has already threaded the `Listable` attribute over lists; any other argument returns `NULL`, unevaluated and silent) and hands it to `map_case(s, len, 1, 0, SIZE_MAX)`. `map_case` walks the UTF-8 bytes once: an ASCII byte takes an inline fast path (`a–z`), while a multi-byte sequence is decoded with `u8_decode`, looked up by binary search (`case_find`) in `k_upper_map`, and either re-encoded from the table's 1–3 target code points with `u8_encode` or, when the character has no mapping, copied byte for byte. This is the full (not simple) case mapping, context- and locale-free.

**Data structures.** `k_upper_map` is a sorted array of `U8CaseMap {cp, n, to[3]}` records — every code point whose → mapping differs from itself — generated from Python's `unicodedata` / `str.→()` by `tools/gen_unicode_tables.py` (Unicode 16.0.0) into `src/strings/unicode_tables.inc`. The output buffer starts at input length + 16 bytes and doubles only when fewer than 13 bytes of headroom remain (room for one character expanding to three 4-byte code points), so a same-size mapping never reallocates; the buffer is adopted into the result `String` without a copy.

**Complexity / limits.** `O(L log M)` for a string of `L` characters and `M` table entries (`O(L)` for pure ASCII). No final-sigma or Turkish/Lithuanian tailoring. There is no packed-array or `Compile[]` surface: the argument and result are `String`s, which no machine buffer or compiled register holds.
