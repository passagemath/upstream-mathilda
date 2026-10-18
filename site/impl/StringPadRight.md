---
source: src/strings/stringpad.c
---
**Algorithm.** `StringPadRight` mirrors `StringPadLeft` through the same `pad_dispatch`, called with `left` false. It pads on the right, and when the string is longer than `n` it keeps the first `n` characters. The pad string (default `" "`) cycles one character at a time and is aligned to the result: result character `i` (0-based) is `p[i mod plen]`. So the cycle continues across the original text rather than restarting after it. `StringPadRight["ab", 6, "xyz"]` is `"abzxyz"` and `StringPadRight["x", 6, "ab"]` is `"xbabab"`. Before v0.456 these gave `"xababa"`-style answers. A `List` first argument pads every element to `n`, or, in the one-argument form, to the longest element.

**Data structures.** A `U8View` gives the pad string's code-point offsets. Lengths count Unicode characters. Each string gets one output buffer, and a list input gives a `List`.

**Complexity / limits.** `O(n)` per string. A wrong arity gives `::argb`, a negative or non-integer `n` gives `::intnm`, an empty or non-string pad gives `::stringnz`, and a non-list or non-string element gives `::strlist`. All of these leave the call unevaluated. Not `Listable`.
