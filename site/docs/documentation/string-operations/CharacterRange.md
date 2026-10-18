# CharacterRange

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`CharacterRange["c1", "c2"]`**

Gives the list of characters with codes from that of c1 to that of c2.

**`CharacterRange[n1, n2]`**

Gives the characters with codes n1 through n2.

## Examples (11)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= CharacterRange["a", "e"]
Out[1]= {"a", "b", "c", "d", "e"}

In[2]:= CharacterRange[65, 70]
Out[2]= {"A", "B", "C", "D", "E", "F"}

In[3]:= Length[Alphabet[]]
Out[3]= 26
```

### Applications (8)

```mathematica
In[4]:= CharacterRange["a", "e"]
Out[4]= {"a", "b", "c", "d", "e"}
```

Any Unicode characters, by code point

```mathematica
In[5]:= CharacterRange["α", "ε"]
Out[5]= {"α", "β", "γ", "δ", "ε"}
```

Integer endpoints are code points

```mathematica
In[6]:= CharacterRange[65, 70]
Out[6]= {"A", "B", "C", "D", "E", "F"}
```

Circled digits

```mathematica
In[7]:= CharacterRange[9312, 9316]
Out[7]= {"①", "②", "③", "④", "⑤"}
```

```mathematica
In[8]:= StringJoin[CharacterRange["0", "9"]]
Out[8]= "0123456789"
```

An empty range

```mathematica
In[9]:= CharacterRange["e", "a"]
Out[9]= {}
```

As a string pattern it is a character class

```mathematica
In[10]:= StringCases["x9a3b", CharacterRange["a", "z"]]
Out[10]= {"x", "a", "b"}
```

Mixed endpoint kinds: CharacterRange::argtype

```mathematica
In[11]:= CharacterRange["a", 70]
Out[11]= CharacterRange["a", 70]
```

## Implementation notes

**Algorithm.** `builtin_characterrange` requires two arguments and converts each with `range_end`: a machine `Integer` that is a valid code point, or a `String` that decodes (via `u8_decode`) to exactly one code point. If either end is invalid, or the two ends are of different kinds (one string, one integer), it emits `CharacterRange::argtype` through `mth_message` and stays unevaluated. Otherwise `char_list(lo, hi)` encodes each code point from `lo` to `hi` with `u8_encode` into a one-character `String`, giving `{}` when `lo > hi`. As a string *pattern*, `CharacterRange["a", "z"]` is not expanded at all: `string_pattern.c` emits the PCRE character class `[a-z]` (`put_class_char`).

**Data structures.** A `List` of `n = hi - lo + 1` one-character strings built by `strings_char_expr`. A surrogate code point inside an integer range (which `u8_encode` refuses) is rendered as `"?"` rather than producing invalid UTF-8.

**Complexity / limits.** `O(n)` time and memory in the size of the range; there is no lazy form, so a range spanning the whole code space allocates a million-element list. Ranges are by raw code point, with no collation or script awareness.

**Attributes:** `Protected`.

## References

**See also:** [Alphabet](../../string-operations/Alphabet/)

- Source: [`src/strings/charcode.c`](https://github.com/stblake/mathilda/blob/main/src/strings/charcode.c)
- Specification: [`docs/spec/builtins/string-operations.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/string-operations.md)
- Tests: [`tests/test_string_review.c`](https://github.com/stblake/mathilda/blob/main/tests/test_string_review.c)
- Tests: [`tests/test_stringfns.c`](https://github.com/stblake/mathilda/blob/main/tests/test_stringfns.c)

## Notes & additional examples

### Notes

- The endpoints are both one-character strings or both non-negative integer codes;
  the range is by code point, so `CharacterRange["ά", "ε"]` also picks up the
  accented Greek letters lying between them.
- `{}` when the first code exceeds the second. Anything else emits
  `CharacterRange::argtype` and stays unevaluated.
- Inside a string pattern (`StringCases`, `StringMatchQ`, …) it is translated to a
  PCRE class `[c1-c2]` instead of being expanded.
- Attributes `{Protected}`.
- See also: `Alphabet`, `FromCharacterCode`, `LetterCharacter`.
