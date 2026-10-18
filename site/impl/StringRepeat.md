---
source: src/strings/stringrepeat.c
---
**Algorithm.** `builtin_stringrepeat` takes 2–3 arguments and gives `StringRepeat::argt` otherwise. It works out how many whole copies of the base fit, which is `n`, or fewer when a maximum length `max` caps the result. Each whole copy is one `memcpy`. When the cap cuts the last copy, only that cut needs a character offset (`u8_offset`). The buffer is adopted into the result string (`expr_new_string_adopt`) without a second copy.

**Data structures.** One output buffer of exactly the final byte length. Lengths count Unicode code points (`u8_length`), so `StringRepeat["a\:00e9", 2, 3]` is `"aéa"` and not a string cut in the middle of a character.

**Complexity / limits.** `O(output)`. `n == 0`, an empty base string, or `max == 0` gives `""`. A non-string base gives `StringRepeat::string`, and a negative or non-integer count or `max` gives `StringRepeat::intp`. Both route through `mth_message`, so `Quiet`/`Check` see them, and leave the call unevaluated. A result length that would overflow `size_t` also stays unevaluated.
