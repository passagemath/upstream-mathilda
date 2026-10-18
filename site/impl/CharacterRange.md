---
source: src/strings/charcode.c, src/strings/regex/string_pattern.c
---
**Algorithm.** `builtin_characterrange` requires two arguments and converts each with `range_end`: a machine `Integer` that is a valid code point, or a `String` that decodes (via `u8_decode`) to exactly one code point. If either end is invalid, or the two ends are of different kinds (one string, one integer), it emits `CharacterRange::argtype` through `mth_message` and stays unevaluated. Otherwise `char_list(lo, hi)` encodes each code point from `lo` to `hi` with `u8_encode` into a one-character `String`, giving `{}` when `lo > hi`. As a string *pattern*, `CharacterRange["a", "z"]` is not expanded at all: `string_pattern.c` emits the PCRE character class `[a-z]` (`put_class_char`).

**Data structures.** A `List` of `n = hi - lo + 1` one-character strings built by `strings_char_expr`. A surrogate code point inside an integer range (which `u8_encode` refuses) is rendered as `"?"` rather than producing invalid UTF-8.

**Complexity / limits.** `O(n)` time and memory in the size of the range; there is no lazy form, so a range spanning the whole code space allocates a million-element list. Ranges are by raw code point, with no collation or script awareness.
