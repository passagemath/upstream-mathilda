# ToCharacterCode

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`ToCharacterCode["string"]`**

Gives the list of the Unicode code points of the characters of "string" (a packed integer array).

**`ToCharacterCode[{s1, s2, ...}]`**

Gives one list of codes per string.

## Examples (12)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= ToCharacterCode["h\:00e9"]
Out[1]= {104, 233}

In[2]:= FromCharacterCode[{72, 105}]
Out[2]= "Hi"

In[3]:= FromCharacterCode[{{72, 105}, {955, 8364}}]
Out[3]= {"Hi", "λ€"}
```

### Applications (9)

```mathematica
In[4]:= ToCharacterCode["Hello"]
Out[4]= {72, 101, 108, 108, 111}
```

Unicode code points, including beyond the BMP

```mathematica
In[5]:= ToCharacterCode["héllo €😀"]
Out[5]= {104, 233, 108, 108, 111, 32, 8364, 128512}
```

A combining accent is its own code point

```mathematica
In[6]:= ToCharacterCode["e\:0301"]
Out[6]= {101, 769}
```

One list per string

```mathematica
In[7]:= ToCharacterCode[{"ab", "xyz"}]
Out[7]= {{97, 98}, {120, 121, 122}}
```

The codes are ordinary integers

```mathematica
In[8]:= ToCharacterCode["abc"] - 96
Out[8]= {1, 2, 3}
```

Long strings give a packed buffer

```mathematica
In[9]:= PackedArrayQ[ToCharacterCode[StringRepeat["a", 100]]]
Out[9]= True
```

ROT13 by arithmetic on codes

```mathematica
In[10]:= FromCharacterCode[Mod[ToCharacterCode["hello"] - 97 + 13, 26] + 97]
Out[10]= "uryyb"
```

Encodings are not supported: unevaluated

```mathematica
In[11]:= ToCharacterCode["abc", "UTF8"]
Out[11]= ToCharacterCode["abc", "UTF8"]
```

ToCharacterCode::strse

```mathematica
In[12]:= ToCharacterCode[5]
Out[12]= ToCharacterCode[5]
```

## Implementation notes

**Algorithm.** `builtin_tocharactercode` accepts one argument: a `String` goes to `codes_of`; a flat `List` of strings maps `codes_of` over its elements into a new `List`; anything else emits `ToCharacterCode::strse` (through `strings_arg_message` / `mth_message`) and stays unevaluated. `codes_of` counts the characters with `u8_length`, then asks `ndbuild_open_i64` for a packed `int64` vector of that length; when packing applies it decodes the code points straight into the buffer with `u8_decode`, and otherwise (short strings, packing disabled) builds a `List` of machine `Integer`s.

**Data structures.** Strings are stored as UTF-8, so each code point is decoded (1–4 bytes) on the fly. The packed result is an `EXPR_NDARRAY` of dtype `NDT_INT64` presented as a `List` (`NDA_HEAD_LIST`); `ndbuild_open` returns `NULL` below `pack_min_elements()` (the `MATHILDA_PACK_MIN`-tunable threshold), which selects the boxed path.

**Complexity / limits.** Two linear passes (length, then decode): `O(B)` in the byte length, with no per-character allocation on the packed path. Only Unicode code points are produced; the encoding form `ToCharacterCode[s, "enc"]` is not implemented and stays unevaluated. There is no `Compile[]` lowering because the input is a `String` (recorded as an exemption in the compile-coverage audit).

**Attributes:** `Protected`.

## References

**See also:** [FromCharacterCode](../../string-operations/FromCharacterCode/), [NDArray](../../linear-algebra/NDArray/), [String](../../other-advanced/String/)

- F. Yergeau, *UTF-8, a transformation format of ISO 10646*, RFC 3629 (IETF, 2003).
- Source: [`src/strings/charcode.c`](https://github.com/stblake/mathilda/blob/main/src/strings/charcode.c)
- Specification: [`docs/spec/builtins/string-operations.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/string-operations.md)
- Tests: [`tests/test_string_review.c`](https://github.com/stblake/mathilda/blob/main/tests/test_string_review.c)

## Notes & additional examples

### Notes

- Returns Unicode code points (decoded from the internal UTF-8), not bytes:
  `ToCharacterCode["é"]` is `{233}`.
- Above the packing threshold the codes are written straight into a packed `int64`
  buffer; short strings give a plain `List`. `FromCharacterCode` reads that buffer
  without boxing, so `FromCharacterCode[ToCharacterCode[s]]` stays on the buffer.
- The two-argument encoding form `ToCharacterCode[s, "enc"]` is not implemented.
- Accepts a string or a flat list of strings; otherwise `::strse`. Not `Listable`.
- No `Compile[]` lowering: the input is a `String`, which no compiled register holds.
- See also: `FromCharacterCode`, `Characters`, `StringLength`.
