# StringPadLeft

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`StringPadLeft["string", n]`**

Makes "string" length n, padding on the left with spaces or truncating (keeping the last n characters) as needed.

**`StringPadLeft["string", n, "padding"]`**

Pads with repeated copies of "padding".

**`StringPadLeft[{s1, s2, ...}]`**

Pads each string on the left with spaces to the length of the longest, making them all the same length.

**`StringPadLeft[{s1, s2, ...}, n, ...]`**

Pads or truncates each string to length n.

## Examples (11)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= StringPadLeft["abcde", 10]
Out[1]= "     abcde"

In[2]:= StringPadRight["abcde", 10, "."]
Out[2]= "abcde....."

In[3]:= StringPadLeft[{"a", "ab", "abc", "abcd", "abcde"}]
Out[3]= {"    a", "   ab", "  abc", " abcd", "abcde"}

In[4]:= StringPadLeft[{"a", "ab", "abc", "abcd", "abcde"}, 3]
Out[4]= {"  a", " ab", "abc", "bcd", "cde"}
```

### Applications (7)

Left-pad to width 5 with zeros

```mathematica
In[5]:= StringPadLeft["42", 5, "0"]
Out[5]= "00042"
```

Too long: truncation keeps the last n characters

```mathematica
In[6]:= StringPadLeft["abcdef", 3]
Out[6]= "def"
```

The one-argument list form pads to the longest element

```mathematica
In[7]:= StringPadLeft[{"a", "bb", "ccc"}]
Out[7]= {"  a", " bb", "ccc"}
```

Right-align a column

```mathematica
In[8]:= StringPadLeft[{"1", "22", "333"}, 4, "."]
Out[8]= {"...1", "..22", ".333"}
```

A multi-character pad is aligned to the end of the result

```mathematica
In[9]:= StringPadLeft["ab", 7, "xyz"]
Out[9]= "zxyzxab"
```

Lengths count characters, not UTF-8 bytes

```mathematica
In[10]:= StringPadLeft["caf\:00e9", 6, "*"]
Out[10]= "**café"
```

A negative length is rejected and stays unevaluated

```mathematica
In[11]:= StringPadLeft["abc", -1]
Out[11]= StringPadLeft["abc", -1]
```

## Algorithm

stringpad.c - StringPadLeft and StringPadRight builtins for Mathilda

```text
StringPadLeft["str", n]        - "str" of length n, padded on the left with
                                 spaces or truncated (keeping the last n chars).
StringPadLeft["str", n, "p"]   - as above, padded by repeating copies of "p".
StringPadLeft["str"]           - unevaluated (StringPadLeft::strlist).
StringPadLeft[{s1, s2, ...}]   - pads each string with spaces to the length of
                                 the longest, so all become the same length.
StringPadLeft[{s1, ...}, n]    - pads or truncates each string to length n.
```

StringPadLeft[{s1, ...}, n, "p"] - as above, using padding string "p".

StringPadRight is identical except padding is added on the right and truncation keeps the first n characters.

Padding is the pad string repeated cyclically and ALIGNED to the result, as in Mathematica: on the right, result character i (0-based) is p[i mod plen]; on the left, it is p[(i - n) mod plen], i.e. the cycle is anchored at the result's end. So StringPadRight["ab", 6, "xyz"] is "abzxyz" and StringPadLeft["ab", 7, "xyz"] is "zxyzxab".

Lengths count characters (Unicode code points); the pad string is cycled a character at a time.

## Implementation notes

**Algorithm.** `StringPadLeft` and `StringPadRight` share `pad_dispatch`, with a `left` flag choosing the side. `pad_one` either truncates or pads. Truncation keeps the last `n` characters on the left variant and the first `n` on the right one. Padding cycles the pad string (default `" "`) one character at a time, aligned to the result as Mathematica does: on the left, result character `i` is `p[(i - n) mod plen]`. The cycle is therefore anchored at the result's end, and `StringPadLeft["ab", 7, "xyz"]` is `"zxyzxab"` (v0.456 fixed the alignment). With a `List` first argument every element is padded to the explicit `n`. In the one-argument list form, `n` is the length of the longest element, so all widths come out equal.

**Data structures.** A `U8View` gives the pad string's character offsets, and `u8_length` / `u8_offset` give the subject's. Lengths count Unicode code points, not UTF-8 bytes. Each string gets one output buffer, and a list input gives a `List`.

**Complexity / limits.** `O(n)` per string. A wrong argument count emits `::argb`. Bad arguments give Mathematica's messages and leave the call unevaluated: a negative or non-integer `n` gives `::intnm`, an empty or non-string pad gives `::stringnz`, and a non-list (including a bare string without `n`) or a list holding a non-string gives `::strlist`. `StringPadLeft[{}, n]` is `""` and `StringPadLeft[{}]` is `{}`. The function is not `Listable`, because the one-argument form has to see the whole list.

**Attributes:** `Protected`.

## References

**See also:** [StringPadRight](../../string-operations/StringPadRight/)

- Source: [`src/strings/stringpad.c`](https://github.com/stblake/mathilda/blob/main/src/strings/stringpad.c)
- Specification: [`docs/spec/builtins/string-operations.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/string-operations.md)
- Tests: [`tests/test_string_review.c`](https://github.com/stblake/mathilda/blob/main/tests/test_string_review.c)
- Tests: [`tests/test_stringpad.c`](https://github.com/stblake/mathilda/blob/main/tests/test_stringpad.c)

## Notes & additional examples

### Notes

- `StringPadLeft[s, n, p]` makes `s` exactly `n` characters long. It pads on the
  left with copies of `p` (default `" "`), or truncates, keeping the last `n`
  characters.
- A multi-character pad is cycled one character at a time and aligned to the
  end of the result, as in Mathematica. `"ab"` padded to 7 with `"xyz"` is
  `"zxyzxab"`.
- On a list of strings every element is padded to `n`. Without `n`, every
  element is padded to the longest, so a ragged column becomes rectangular.
- Since v0.456 lengths count Unicode code points.
- Bad arguments emit Mathematica's messages and leave the call unevaluated:
  `::intnm` for a negative length, `::stringnz` for an empty pad, and
  `::strlist` for a bare string without `n`.
- Not `Listable`, because the one-argument form needs the whole list.
- See also: `StringPadRight`, `StringRepeat`, `IntegerString`.
