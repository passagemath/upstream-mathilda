---
source: src/strings/stringriffle.c
---
**Algorithm.** `builtin_stringriffle` treats the first argument as the data and the rest as separators, one per level, with level 1 outermost. `parse_sep` reads each separator into a `{left, sep, right}` triple. A plain string becomes `{"", sep, ""}`, and a 3-string list is the triple itself. `riffle_build` then recurses over the data. A string leaf is copied verbatim, and any other leaf is rendered in OutputForm (`expr_to_string_outputform`, which is what `ToString` gives). Each level joins its children with its explicit separator. When a level has none, the default scheme applies, chosen by `depth_from_bottom`: one space at the innermost level, and one more newline for each level above it. It is the inverse of `StringSplit`.

**Data structures.** A `SepSpec` array that borrows the argument strings. The whole result is built in one growable `RegexBuf`, so each leaf is copied exactly once whatever the nesting depth, and the buffer is adopted into the result string. This replaced a build-each-level-then-concatenate design in v0.456 (3.6 ms → 1.0 ms on the string benchmark).

**Complexity / limits.** `O(total output)`. No arguments emits `StringRiffle::argm`. A first argument that is neither a list nor a string emits `StringRiffle::list`, and a malformed separator emits `StringRiffle::string`; both leave the call unevaluated. Not `Listable`, so it can see the whole nested structure.
