---
references:
  - "The Unicode Consortium, *The Unicode Standard, Version 16.0* (2024), §3.13 Default Case Algorithms, and the Unicode Character Database (UnicodeData.txt, SpecialCasing.txt)."
source: src/strings/stringcase.c, src/strings/unicode_tables.inc
---
**Algorithm.** `builtin_decapitalize` calls `capitalize_builtin(res, "Decapitalize", 0)`. A `String` argument goes to the shared case kernel `map_case(s, len, 0, 0, 1)` — the same routine as `ToLowerCase`, restricted by its `from`/`count` window to character 0, so only the first code point is looked up in `k_lower_map` and every later byte is copied verbatim. A flat `List` of strings (checked by `is_string_list`) is mapped element by element into a new `List`; anything else emits `Decapitalize::strse` through `strings_arg_message` (the `mth_message` funnel) and stays unevaluated.

**Data structures.** The generated full-case-mapping table `k_lower_map` (`U8CaseMap` records, Unicode 16.0.0, from `tools/gen_unicode_tables.py`), searched by binary search; ASCII is mapped inline. The first character may expand to up to three code points (`Decapitalize["İ"]` gives `i` plus U+0307`).

**Complexity / limits.** One `O(L)` copy of the string plus one `O(log M)` table lookup. It applies the lowercase mapping, not Unicode titlecasing. Not `Listable`: nested or mixed lists are rejected rather than threaded. No numeric (packed / `Compile[]`) surface.
