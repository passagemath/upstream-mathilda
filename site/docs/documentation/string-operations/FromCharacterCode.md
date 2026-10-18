# FromCharacterCode

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`FromCharacterCode[n]`**

Gives the one-character string with Unicode code point n.

**`FromCharacterCode[{n1, n2, ...}]`**

Gives the string of the characters with those codes.

**`FromCharacterCode[{{n11, ...}, {n21, ...}, ...}]`**

Gives one string per sublist.

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
In[4]:= FromCharacterCode[{72, 105}]
Out[4]= "Hi"
```

A single code gives a one-character string

```mathematica
In[5]:= FromCharacterCode[955]
Out[5]= "λ"
```

A matrix gives one string per row

```mathematica
In[6]:= FromCharacterCode[{{72, 105}, {955, 8364}}]
Out[6]= {"Hi", "λ€"}
```

A code point outside the BMP

```mathematica
In[7]:= FromCharacterCode[128512]
Out[7]= "😀"
```

```mathematica
In[8]:= FromCharacterCode[Reverse[ToCharacterCode["stressed"]]]
Out[8]= "desserts"

In[9]:= FromCharacterCode[Range[97, 102]]
Out[9]= "abcdef"
```

An int64 NDArray is read directly

```mathematica
In[10]:= FromCharacterCode[NDArray[{77, 97, 116}, DataType -> "int64"]]
Out[10]= "Mat"
```

Beyond U+10FFFF: FromCharacterCode::notunicode

```mathematica
In[11]:= FromCharacterCode[{65, 2000000}]
Out[11]= FromCharacterCode[{65, 2000000}]
```

Negative: FromCharacterCode::intnm

```mathematica
In[12]:= FromCharacterCode[-1]
Out[12]= FromCharacterCode[-1]
```

## Implementation notes

**Algorithm.** `builtin_fromcharactercode` takes one argument and dispatches in `fcc_value`: a machine `Integer` becomes a one-character string; a rank-1 or rank-2 `NDT_INT64` `NDArray` (a packed `List` or a visible `NDArray`) is read directly as a code buffer, row by row for rank 2; an `NDArray` of another dtype is unpacked with `pack_unpack` and retried as a `List`; a `List` whose first element is a `List` gives one string per sublist; any other `List` goes to `string_of_list`. Every code is validated by `code_ok` (`0`–`0x10FFFF`, surrogates `0xD800`–`0xDFFF` excluded) and encoded with `u8_encode`. On failure `fcc_bad` emits `::notunicode` (non-negative but out of range, naming the code list and position) or `::intnm` (negative or non-integer) through `mth_message`, and the call stays unevaluated.

**Data structures.** Codes are gathered into a contiguous `int64_t` array (borrowed from the buffer on the packed path), and the output is a single `malloc`'d byte buffer of `4n + 1` bytes — the UTF-8 worst case — adopted into the result `String`. `FromCharacterCode` is on both the `AWARE` and `INT64_OK` lists in `src/pack.c`, so the evaluator's transparency gate passes an `int64` buffer through without boxing; together with `ToCharacterCode`'s packed output, `FromCharacterCode[ToCharacterCode[s]]` never materialises an `Integer` per character.

**Complexity / limits.** `O(n)` in the number of codes. The encoding form `FromCharacterCode[codes, "enc"]` is not implemented. There is no `Compile[]` lowering — the result is a `String`, which no compiled register holds (the documented exemption in `tools/compile_coverage.py`).

**Attributes:** `Protected`.

## References

**See also:** [ToCharacterCode](../../string-operations/ToCharacterCode/), [NDArray](../../linear-algebra/NDArray/), [String](../../other-advanced/String/)

- F. Yergeau, *UTF-8, a transformation format of ISO 10646*, RFC 3629 (IETF, 2003).
- Source: [`src/strings/charcode.c`](https://github.com/stblake/mathilda/blob/main/src/strings/charcode.c)
- Specification: [`docs/spec/builtins/string-operations.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/string-operations.md)
- Tests: [`tests/test_string_review.c`](https://github.com/stblake/mathilda/blob/main/tests/test_string_review.c)

## Notes & additional examples

### Notes

- Codes are Unicode code points `0`–`1114111`, encoded to UTF-8. Surrogate code
  points `U+D800`–`U+DFFF` are rejected as well.
- A rank-1 or rank-2 `int64` buffer (packed `List` or visible `NDArray`) is read
  directly — `FromCharacterCode` is on the packed-aware and `INT64_OK` lists in
  `src/pack.c`. Other `NDArray` dtypes are unpacked first and then fail the integer
  check.
- An invalid code emits `::notunicode` (too large) or `::intnm` (negative or
  non-integer) and leaves the whole call unevaluated; both route through the
  message funnel, so `Check`/`Quiet` work.
- The encoding form `FromCharacterCode[codes, "enc"]` is not implemented.
- See also: `ToCharacterCode`, `CharacterRange`, `StringJoin`.
