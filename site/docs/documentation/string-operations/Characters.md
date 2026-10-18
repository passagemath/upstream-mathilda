# Characters

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Characters["string"]`**

Gives a list of the characters in a string. Each character is given as a length-1 string.

## Examples (13)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (5)

```mathematica
In[1]:= Characters["ABC"]
Out[1]= {"A", "B", "C"}

In[2]:= Characters["A string."]
Out[2]= {"A", " ", "s", "t", "r", "i", "n", "g", "."}

In[3]:= Characters[""]
Out[3]= {}

In[4]:= Characters[{"ABC", "DEF", "XYZ"}]
Out[4]= {{"A", "B", "C"}, {"D", "E", "F"}, {"X", "Y", "Z"}}

In[5]:= Characters[x]
Out[5]= Characters[x]
```

### Applications (8)

A string becomes a list of length-1 strings

```mathematica
In[6]:= Characters["cat"]
Out[6]= {"c", "a", "t"}
```

The inverse is StringJoin

```mathematica
In[7]:= StringJoin[Characters["abcd"]]
Out[7]= "abcd"
```

One element per Unicode character, so é stays whole

```mathematica
In[8]:= Characters["h\:00e9llo"]
Out[8]= {"h", "é", "l", "l", "o"}
```

Named characters decode in string literals

```mathematica
In[9]:= Characters["\[Alpha]\[Beta]"]
Out[9]= {"α", "β"}
```

Both count code points

```mathematica
In[10]:= StringLength["h\:00e9llo"] == Length[Characters["h\:00e9llo"]]
Out[10]= True
```

Listable over a list of strings

```mathematica
In[11]:= Characters[{"ab", "c"}]
Out[11]= {{"a", "b"}, {"c"}}
```

Letter frequencies

```mathematica
In[12]:= Tally[Characters["mississippi"]]
Out[12]= {{"m", 1}, {"i", 4}, {"s", 4}, {"p", 2}}
```

A non-string stays unevaluated

```mathematica
In[13]:= Characters[123]
Out[13]= Characters[123]
```

## Performance

Against other systems, from the benchmark suite (same input, results cross-checked for agreement):

| case | Mathilda | Wolfram | Python |
|---|---:|---:|---:|
| Characters of 200k chars | 4.38 s | 2.54 s | 0.419 s |
| StringSplit on space, 200k chars | 3.12 s | 4.13 s | 0.723 s |
| StringCases regex, 200k chars | 2.36 s | 3.69 s | 3.34 s |
| StringReplace regex, 200k chars | 1.97 s | 4.74 s | 3.35 s |
| StringReplace literal, 200k chars | 0.332 s | 1.17 s | 0.195 s |
| StringCount substring, 200k chars | 0.224 s | 0.364 s | 0.102 s |

## Implementation notes

**Algorithm.** `builtin_characters` takes one `EXPR_STRING` and returns `NULL` (unevaluated) for anything else. It sizes the result `List` with `u8_length` and then walks the UTF-8 bytes with `u8_char_len`, emitting one string per Unicode code point. A multi-byte character such as `"é"` or `"α"` stays whole. A malformed byte counts as one character of its own, so invalid input is never split mid-sequence and round-trips byte for byte. String literals decode `\[Name]`, `\:hhhh` and the other Mathematica escapes before `Characters` sees them.

**Data structures.** `strings_char_expr` hands out a shared, refcounted string for each ASCII character. The 128 entries are built on first use, and `expr_copy` is a refcount bump. So `Characters`, `StringSplit[s, ""]` and `StringPart` on ASCII text allocate no new string per character, only the list. Before v0.456 this cost 23.7 ms on 200k characters; it now takes about 2 ms.

**Complexity / limits.** `O(len)`. `Length[Characters[s]] == StringLength[s]` for every string, since both count code points. `Listable | Protected`, so a list of strings gives a list of character lists, and `""` gives `{}`.

**Attributes:** `Listable`, `Protected`.

## References

**See also:** [List](../../other-advanced/List/)

- Source: [`src/strings/characters.c`](https://github.com/stblake/mathilda/blob/main/src/strings/characters.c)
- Specification: [`docs/spec/builtins/string-operations.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/string-operations.md)
- Tests: [`tests/test_string_review.c`](https://github.com/stblake/mathilda/blob/main/tests/test_string_review.c)
- Tests: [`tests/test_strings.c`](https://github.com/stblake/mathilda/blob/main/tests/test_strings.c)

## Notes & additional examples

### Notes

- `Characters[s]` gives one length-1 string per Unicode code point. Since v0.456
  a multi-byte UTF-8 character stays whole, so
  `Length[Characters[s]] == StringLength[s]` for every string. Earlier versions
  split by byte.
- Each ASCII character is a single shared, refcounted string, so splitting ASCII
  text allocates only the list (200k characters: 23.7 ms → about 2 ms in v0.456).
- `Listable | Protected`, and `Characters[""]` is `{}`.
- See also: `StringJoin`, `StringSplit` (with `""`), `ToCharacterCode`,
  `StringLength`.
