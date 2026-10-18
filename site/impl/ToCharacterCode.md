---
references:
  - "F. Yergeau, *UTF-8, a transformation format of ISO 10646*, RFC 3629 (IETF, 2003)."
source: src/strings/charcode.c, src/strings/utf8.c, src/pack.c
---
**Algorithm.** `builtin_tocharactercode` accepts one argument: a `String` goes to `codes_of`; a flat `List` of strings maps `codes_of` over its elements into a new `List`; anything else emits `ToCharacterCode::strse` (through `strings_arg_message` / `mth_message`) and stays unevaluated. `codes_of` counts the characters with `u8_length`, then asks `ndbuild_open_i64` for a packed `int64` vector of that length; when packing applies it decodes the code points straight into the buffer with `u8_decode`, and otherwise (short strings, packing disabled) builds a `List` of machine `Integer`s.

**Data structures.** Strings are stored as UTF-8, so each code point is decoded (1–4 bytes) on the fly. The packed result is an `EXPR_NDARRAY` of dtype `NDT_INT64` presented as a `List` (`NDA_HEAD_LIST`); `ndbuild_open` returns `NULL` below `pack_min_elements()` (the `MATHILDA_PACK_MIN`-tunable threshold), which selects the boxed path.

**Complexity / limits.** Two linear passes (length, then decode): `O(B)` in the byte length, with no per-character allocation on the packed path. Only Unicode code points are produced; the encoding form `ToCharacterCode[s, "enc"]` is not implemented and stays unevaluated. There is no `Compile[]` lowering because the input is a `String` (recorded as an exemption in the compile-coverage audit).
