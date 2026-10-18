# ToUpperCase

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`ToUpperCase["string"]`**

Gives "string" with every letter converted to uppercase (full Unicode case mapping: ToUpperCase\["\\[SZ\]"\] is "SS").

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
In[4]:= ToUpperCase["hello, world"]
Out[4]= "HELLO, WORLD"
```

Accented Latin letters map through the Unicode tables

```mathematica
In[5]:= ToUpperCase["déjà vu"]
Out[5]= "DÉJÀ VU"
```

Full case mapping: one ß becomes the two letters SS

```mathematica
In[6]:= ToUpperCase["straße"]
Out[6]= "STRASSE"
```

So the length can grow

```mathematica
In[7]:= StringLength /@ {"ß", ToUpperCase["ß"]}
Out[7]= {1, 2}
```

Ligatures and other multi-character expansions

```mathematica
In[8]:= ToUpperCase["ﬁ ŉ"]
Out[8]= "FI ʼN"
```

Any cased script, here Greek and Cyrillic

```mathematica
In[9]:= ToUpperCase["αβγ ж"]
Out[9]= "ΑΒΓ Ж"
```

Listable: threads through nested lists

```mathematica
In[10]:= ToUpperCase[{"abc", {"x", "y"}}]
Out[10]= {"ABC", {"X", "Y"}}
```

A non-string element stays unevaluated, with no message

```mathematica
In[11]:= ToUpperCase[{"abc", 42}]
Out[11]= {"ABC", ToUpperCase[42]}
```

## Implementation notes

**Algorithm.** `builtin_touppercase` calls `case_builtin`, which accepts exactly one `String` argument (the evaluator has already threaded the `Listable` attribute over lists; any other argument returns `NULL`, unevaluated and silent) and hands it to `map_case(s, len, 1, 0, SIZE_MAX)`. `map_case` walks the UTF-8 bytes once: an ASCII byte takes an inline fast path (`a–z`), while a multi-byte sequence is decoded with `u8_decode`, looked up by binary search (`case_find`) in `k_upper_map`, and either re-encoded from the table's 1–3 target code points with `u8_encode` or, when the character has no mapping, copied byte for byte. This is the full (not simple) case mapping, context- and locale-free.

**Data structures.** `k_upper_map` is a sorted array of `U8CaseMap {cp, n, to[3]}` records — every code point whose → mapping differs from itself — generated from Python's `unicodedata` / `str.→()` by `tools/gen_unicode_tables.py` (Unicode 16.0.0) into `src/strings/unicode_tables.inc`. The output buffer starts at input length + 16 bytes and doubles only when fewer than 13 bytes of headroom remain (room for one character expanding to three 4-byte code points), so a same-size mapping never reallocates; the buffer is adopted into the result `String` without a copy.

**Complexity / limits.** `O(L log M)` for a string of `L` characters and `M` table entries (`O(L)` for pure ASCII). No final-sigma or Turkish/Lithuanian tailoring. There is no packed-array or `Compile[]` surface: the argument and result are `String`s, which no machine buffer or compiled register holds.

**Attributes:** `Listable`, `Protected`.

## References

**See also:** [ToLowerCase](../../string-operations/ToLowerCase/)

- The Unicode Consortium, *The Unicode Standard, Version 16.0* (2024), §3.13 Default Case Algorithms, and the Unicode Character Database (UnicodeData.txt, SpecialCasing.txt).
- Source: [`src/strings/stringcase.c`](https://github.com/stblake/mathilda/blob/main/src/strings/stringcase.c)
- Specification: [`docs/spec/builtins/string-operations.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/string-operations.md)
- Tests: [`tests/test_string_review.c`](https://github.com/stblake/mathilda/blob/main/tests/test_string_review.c)

## Notes & additional examples

### Notes

- Full Unicode (16.0) uppercase mapping, including the multi-code-point expansions
  (`ß → SS`, `ﬁ → FI`, `ŉ → ʼN`); characters with no uppercase form — digits,
  punctuation, CJK — are copied unchanged.
- The mapping is context-free: there is no locale (Turkish dotless i) or
  language-specific tailoring.
- `ToUpperCase` is not invertible: `ToLowerCase[ToUpperCase["straße"]]` is `"strasse"`.
- Attributes `{Listable, Protected}`. A non-string argument is returned unevaluated
  without a message, as in Mathematica.
- No numeric surface (packed array / `Compile[]`): the result is a `String`.
- See also: `ToLowerCase`, `Capitalize`, `UpperCaseQ`.
