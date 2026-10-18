---
source: src/strings/stringrotate.c, src/strings/utf8.c
---
**Algorithm.** `builtin_stringrotateleft` calls the shared kernel `rotate(res, "StringRotateLeft", 1)`. It accepts one or two arguments: the subject must be a `String` (otherwise `StringRotateLeft::string`) and the optional shift a machine `Integer` (otherwise `StringRotateLeft::int`), defaulting to 1; both messages route through `strings_arg_message` / `mth_message` and leave the call unevaluated. The character count `m` comes from `u8_length`; the shift is converted to an equivalent left rotation `k = n mod m`, reduced into `[0, m)`, and its byte offset found with `u8_offset`. The result is the tail from that offset followed by the head — two `memcpy`s into one buffer adopted as the result `String`.

**Data structures.** Only the UTF-8 byte buffer of the subject; positions count Unicode code points, so multi-byte characters are never split.

**Complexity / limits.** `O(B)` in the byte length (one length scan, one offset scan, one copy). Any machine integer shift is accepted (it wraps modulo the length); a bignum shift is rejected with `::int`. A list of strings is not threaded.
