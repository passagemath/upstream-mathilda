# StringPadRight

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`StringPadRight["string", n]`**

Makes "string" length n, padding on the right with spaces or truncating (keeping the first n characters) as needed.

**`StringPadRight["string", n, "padding"]`**

Pads with repeated copies of "padding".

**`StringPadRight[{s1, s2, ...}]`**

Pads each string on the right with spaces to the length of the longest, making them all the same length.

**`StringPadRight[{s1, s2, ...}, n, ...]`**

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

Right-pad to width 5 with zeros

```mathematica
In[5]:= StringPadRight["42", 5, "0"]
Out[5]= "42000"
```

Too long: truncation keeps the first n characters

```mathematica
In[6]:= StringPadRight["abcdef", 3]
Out[6]= "abc"
```

The pad cycle is aligned to the result, so it continues past the text

```mathematica
In[7]:= StringPadRight["ab", 6, "xyz"]
Out[7]= "abzxyz"
```

Result character i is the pad's character i mod 2

```mathematica
In[8]:= StringPadRight["x", 6, "ab"]
Out[8]= "xbabab"
```

The one-argument list form pads to the longest element

```mathematica
In[9]:= StringPadRight[{"a", "bb", "ccc"}]
Out[9]= {"a  ", "bb ", "ccc"}
```

Fixed-width fields

```mathematica
In[10]:= StringRiffle[StringPadRight[{"id", "name", "x"}, 6, "."], "|"]
Out[10]= "id....|name..|x....."
```

An empty pad string is rejected

```mathematica
In[11]:= StringPadRight["abc", 5, ""]
Out[11]= StringPadRight["abc", 5, ""]
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

**Algorithm.** `StringPadRight` mirrors `StringPadLeft` through the same `pad_dispatch`, called with `left` false. It pads on the right, and when the string is longer than `n` it keeps the first `n` characters. The pad string (default `" "`) cycles one character at a time and is aligned to the result: result character `i` (0-based) is `p[i mod plen]`. So the cycle continues across the original text rather than restarting after it. `StringPadRight["ab", 6, "xyz"]` is `"abzxyz"` and `StringPadRight["x", 6, "ab"]` is `"xbabab"`. Before v0.456 these gave `"xababa"`-style answers. A `List` first argument pads every element to `n`, or, in the one-argument form, to the longest element.

**Data structures.** A `U8View` gives the pad string's code-point offsets. Lengths count Unicode characters. Each string gets one output buffer, and a list input gives a `List`.

**Complexity / limits.** `O(n)` per string. A wrong arity gives `::argb`, a negative or non-integer `n` gives `::intnm`, an empty or non-string pad gives `::stringnz`, and a non-list or non-string element gives `::strlist`. All of these leave the call unevaluated. Not `Listable`.

**Attributes:** `Protected`.

## References

**See also:** [StringPadLeft](../../string-operations/StringPadLeft/)

- Source: [`src/strings/stringpad.c`](https://github.com/stblake/mathilda/blob/main/src/strings/stringpad.c)
- Specification: [`docs/spec/builtins/string-operations.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/string-operations.md)
- Tests: [`tests/test_string_review.c`](https://github.com/stblake/mathilda/blob/main/tests/test_string_review.c)
- Tests: [`tests/test_stringpad.c`](https://github.com/stblake/mathilda/blob/main/tests/test_stringpad.c)

## Notes & additional examples

### Notes

- `StringPadRight[s, n, p]` mirrors `StringPadLeft`. It pads on the right, and
  when truncating it keeps the first `n` characters.
- The pad is aligned to the result: result character `i` (0-based) is
  `p[[i mod Length[p]]]`. So the cycle does not restart after the original
  text: `StringPadRight["x", 6, "ab"]` is `"xbabab"`. This matches Mathematica.
  Before v0.456 Mathilda gave `"xababa"`.
- Lengths count Unicode characters. Bad arguments emit `::intnm`, `::stringnz`
  or `::strlist` and leave the call unevaluated.
- See also: `StringPadLeft`, `StringRepeat`, `StringRiffle`.
