---
source: src/strings/characters.c
---
**Algorithm.** `builtin_characters` takes one `EXPR_STRING` and returns `NULL` (unevaluated) for anything else. It sizes the result `List` with `u8_length` and then walks the UTF-8 bytes with `u8_char_len`, emitting one string per Unicode code point. A multi-byte character such as `"é"` or `"α"` stays whole. A malformed byte counts as one character of its own, so invalid input is never split mid-sequence and round-trips byte for byte. String literals decode `\[Name]`, `\:hhhh` and the other Mathematica escapes before `Characters` sees them.

**Data structures.** `strings_char_expr` hands out a shared, refcounted string for each ASCII character. The 128 entries are built on first use, and `expr_copy` is a refcount bump. So `Characters`, `StringSplit[s, ""]` and `StringPart` on ASCII text allocate no new string per character, only the list. Before v0.456 this cost 23.7 ms on 200k characters; it now takes about 2 ms.

**Complexity / limits.** `O(len)`. `Length[Characters[s]] == StringLength[s]` for every string, since both count code points. `Listable | Protected`, so a list of strings gives a list of character lists, and `""` gives `{}`.
