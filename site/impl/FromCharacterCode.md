---
references:
  - "F. Yergeau, *UTF-8, a transformation format of ISO 10646*, RFC 3629 (IETF, 2003)."
source: src/strings/charcode.c, src/strings/utf8.c, src/pack.c, tools/compile_coverage.py
---
**Algorithm.** `builtin_fromcharactercode` takes one argument and dispatches in `fcc_value`: a machine `Integer` becomes a one-character string; a rank-1 or rank-2 `NDT_INT64` `NDArray` (a packed `List` or a visible `NDArray`) is read directly as a code buffer, row by row for rank 2; an `NDArray` of another dtype is unpacked with `pack_unpack` and retried as a `List`; a `List` whose first element is a `List` gives one string per sublist; any other `List` goes to `string_of_list`. Every code is validated by `code_ok` (`0`–`0x10FFFF`, surrogates `0xD800`–`0xDFFF` excluded) and encoded with `u8_encode`. On failure `fcc_bad` emits `::notunicode` (non-negative but out of range, naming the code list and position) or `::intnm` (negative or non-integer) through `mth_message`, and the call stays unevaluated.

**Data structures.** Codes are gathered into a contiguous `int64_t` array (borrowed from the buffer on the packed path), and the output is a single `malloc`'d byte buffer of `4n + 1` bytes — the UTF-8 worst case — adopted into the result `String`. `FromCharacterCode` is on both the `AWARE` and `INT64_OK` lists in `src/pack.c`, so the evaluator's transparency gate passes an `int64` buffer through without boxing; together with `ToCharacterCode`'s packed output, `FromCharacterCode[ToCharacterCode[s]]` never materialises an `Integer` per character.

**Complexity / limits.** `O(n)` in the number of codes. The encoding form `FromCharacterCode[codes, "enc"]` is not implemented. There is no `Compile[]` lowering — the result is a `String`, which no compiled register holds (the documented exemption in `tools/compile_coverage.py`).
