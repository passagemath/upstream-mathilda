# StringRepeat

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`StringRepeat["string", n]`**

Gives a string with "string" repeated n times.

**`StringRepeat["string", n, max]`**

Gives up to n copies of "string", truncated to a total length of at most max characters.

## Examples (10)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= StringRepeat["a", 50]
Out[1]= "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa"

In[2]:= StringRepeat["abc", 10]
Out[2]= "abcabcabcabcabcabcabcabcabcabc"

In[3]:= StringRepeat["ab", 10, 19]
Out[3]= "abababababababababa"
```

### Applications (7)

Four copies

```mathematica
In[4]:= StringRepeat["ab", 4]
Out[4]= "abababab"
```

A rule of ten dashes

```mathematica
In[5]:= StringRepeat["-", 10]
Out[5]= "----------"
```

Capped at 7 characters, so the last copy is partial

```mathematica
In[6]:= StringRepeat["abc", 5, 7]
Out[6]= "abcabca"
```

A cap larger than the result has no effect

```mathematica
In[7]:= StringRepeat["xy", 3, 100]
Out[7]= "xyxyxy"
```

The cap counts characters, so no character is cut in half

```mathematica
In[8]:= StringRepeat["a\:00e9", 2, 3]
Out[8]= "aéa"
```

Zero copies is the empty string

```mathematica
In[9]:= StringRepeat["abc", 0]
Out[9]= ""
```

A negative count: StringRepeat::intp, left unevaluated

```mathematica
In[10]:= StringRepeat["ab", -1]
Out[10]= StringRepeat["ab", -1]
```

## Algorithm

stringrepeat.c - StringRepeat builtin for Mathilda

```text
StringRepeat["str", n]        - "str" concatenated n times.
StringRepeat["str", n, max]   - up to n copies of "str", truncated so the
                                total length is at most max (a partial final
                                copy is allowed).
```

Lengths count characters (Unicode code points): StringRepeat["a\:00e9", 2, 3] is "a\:00e9a". Invalid arguments emit Mathematica's StringRepeat::string / StringRepeat::intp and leave the call unevaluated.

## Implementation notes

**Algorithm.** `builtin_stringrepeat` takes 2–3 arguments and gives `StringRepeat::argt` otherwise. It works out how many whole copies of the base fit, which is `n`, or fewer when a maximum length `max` caps the result. Each whole copy is one `memcpy`. When the cap cuts the last copy, only that cut needs a character offset (`u8_offset`). The buffer is adopted into the result string (`expr_new_string_adopt`) without a second copy.

**Data structures.** One output buffer of exactly the final byte length. Lengths count Unicode code points (`u8_length`), so `StringRepeat["a\:00e9", 2, 3]` is `"aéa"` and not a string cut in the middle of a character.

**Complexity / limits.** `O(output)`. `n == 0`, an empty base string, or `max == 0` gives `""`. A non-string base gives `StringRepeat::string`, and a negative or non-integer count or `max` gives `StringRepeat::intp`. Both route through `mth_message`, so `Quiet`/`Check` see them, and leave the call unevaluated. A result length that would overflow `size_t` also stays unevaluated.

**Attributes:** `Protected`.

## References

- Source: [`src/strings/stringrepeat.c`](https://github.com/stblake/mathilda/blob/main/src/strings/stringrepeat.c)
- Specification: [`docs/spec/builtins/string-operations.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/string-operations.md)
- Tests: [`tests/test_string_review.c`](https://github.com/stblake/mathilda/blob/main/tests/test_string_review.c)
- Tests: [`tests/test_stringrepeat.c`](https://github.com/stblake/mathilda/blob/main/tests/test_stringrepeat.c)

## Notes & additional examples

### Notes

- `StringRepeat[s, n]` concatenates `n` copies of `s`. `StringRepeat[s, n, max]`
  stops at `max` characters, keeping a partial final copy.
- Since v0.456 lengths count Unicode code points.
- A non-string first argument gives `StringRepeat::string`, and a negative or
  non-integer count or cap gives `StringRepeat::intp`. Both leave the call
  unevaluated.
- Whole copies are `memcpy`'d and the buffer is adopted into the result
  without a second copy.
- See also: `StringJoin`, `StringPadRight`, `Table`.
