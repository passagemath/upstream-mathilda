---
source: src/strings/stringpad.c
---
**Algorithm.** `StringPadLeft` and `StringPadRight` share `pad_dispatch`, with a `left` flag choosing the side. `pad_one` either truncates or pads. Truncation keeps the last `n` characters on the left variant and the first `n` on the right one. Padding cycles the pad string (default `" "`) one character at a time, aligned to the result as Mathematica does: on the left, result character `i` is `p[(i - n) mod plen]`. The cycle is therefore anchored at the result's end, and `StringPadLeft["ab", 7, "xyz"]` is `"zxyzxab"` (v0.456 fixed the alignment). With a `List` first argument every element is padded to the explicit `n`. In the one-argument list form, `n` is the length of the longest element, so all widths come out equal.

**Data structures.** A `U8View` gives the pad string's character offsets, and `u8_length` / `u8_offset` give the subject's. Lengths count Unicode code points, not UTF-8 bytes. Each string gets one output buffer, and a list input gives a `List`.

**Complexity / limits.** `O(n)` per string. A wrong argument count emits `::argb`. Bad arguments give Mathematica's messages and leave the call unevaluated: a negative or non-integer `n` gives `::intnm`, an empty or non-string pad gives `::stringnz`, and a non-list (including a bare string without `n`) or a list holding a non-string gives `::strlist`. `StringPadLeft[{}, n]` is `""` and `StringPadLeft[{}]` is `{}`. The function is not `Listable`, because the one-argument form has to see the whole list.
