# ToLowerCase

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`ToLowerCase["string"]`**

Gives "string" with every letter converted to lowercase.

## Examples (11)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= ToUpperCase["straße"]
Out[1]= "STRASSE"

In[2]:= ToLowerCase[{"ABC", "ΣΟΦΙΑ"}]
Out[2]= {"abc", "σοφια"}

In[3]:= ToUpperCase[x]
Out[3]= ToUpperCase[x]
```

### Applications (8)

```mathematica
In[4]:= ToLowerCase["Hello, World"]
Out[4]= "hello, world"
```

Latin-1 capitals

```mathematica
In[5]:= ToLowerCase["ÀÉÎÕÜ"]
Out[5]= "àéîõü"
```

Greek capitals

```mathematica
In[6]:= ToLowerCase["ΑΘΗΝΑ"]
Out[6]= "αθηνα"
```

Full mapping: capital dotted I becomes i plus a combining dot

```mathematica
In[7]:= ToLowerCase["İ"]
Out[7]= "i̇"
```

The two code points of that result

```mathematica
In[8]:= ToCharacterCode[ToLowerCase["İ"]]
Out[8]= {105, 775}
```

Case mapping is not invertible

```mathematica
In[9]:= ToLowerCase[ToUpperCase["straße"]]
Out[9]= "strasse"
```

Listable

```mathematica
In[10]:= ToLowerCase[{"ABC", "XyZ"}]
Out[10]= {"abc", "xyz"}
```

A non-string stays unevaluated

```mathematica
In[11]:= ToLowerCase[x]
Out[11]= ToLowerCase[x]
```

## Implementation notes

**Algorithm.** `builtin_tolowercase` calls `case_builtin`, which accepts exactly one `String` argument (the evaluator has already threaded the `Listable` attribute over lists; any other argument returns `NULL`, unevaluated and silent) and hands it to `map_case(s, len, 0, 0, SIZE_MAX)`. `map_case` walks the UTF-8 bytes once: an ASCII byte takes an inline fast path (`A–Z`), while a multi-byte sequence is decoded with `u8_decode`, looked up by binary search (`case_find`) in `k_lower_map`, and either re-encoded from the table's 1–3 target code points with `u8_encode` or, when the character has no mapping, copied byte for byte. This is the full (not simple) case mapping, context- and locale-free.

**Data structures.** `k_lower_map` is a sorted array of `U8CaseMap {cp, n, to[3]}` records — every code point whose → mapping differs from itself — generated from Python's `unicodedata` / `str.→()` by `tools/gen_unicode_tables.py` (Unicode 16.0.0) into `src/strings/unicode_tables.inc`. The output buffer starts at input length + 16 bytes and doubles only when fewer than 13 bytes of headroom remain (room for one character expanding to three 4-byte code points), so a same-size mapping never reallocates; the buffer is adopted into the result `String` without a copy.

**Complexity / limits.** `O(L log M)` for a string of `L` characters and `M` table entries (`O(L)` for pure ASCII). No final-sigma or Turkish/Lithuanian tailoring. There is no packed-array or `Compile[]` surface: the argument and result are `String`s, which no machine buffer or compiled register holds.

**Attributes:** `Listable`, `Protected`.

## References

**See also:** [ToUpperCase](../../string-operations/ToUpperCase/)

- The Unicode Consortium, *The Unicode Standard, Version 16.0* (2024), §3.13 Default Case Algorithms, and the Unicode Character Database (UnicodeData.txt, SpecialCasing.txt).
- Source: [`src/strings/stringcase.c`](https://github.com/stblake/mathilda/blob/main/src/strings/stringcase.c)
- Specification: [`docs/spec/builtins/string-operations.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/string-operations.md)
- Tests: [`tests/test_string_review.c`](https://github.com/stblake/mathilda/blob/main/tests/test_string_review.c)

## Notes & additional examples

### Notes

- Full Unicode (16.0) lowercase mapping from the generated `k_lower_map` table; a
  character may map to several code points (`İ → i̇`, length 2).
- No context-sensitive rules: a word-final capital sigma `Σ` always becomes `σ`,
  never final `ς`, and there is no locale tailoring.
- Attributes `{Listable, Protected}`; a non-string is left unevaluated silently.
- See also: `ToUpperCase`, `Decapitalize`, `LowerCaseQ`.
