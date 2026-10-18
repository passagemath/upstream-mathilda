# StringTake

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`StringTake["string", n]`**

Gives a string containing the first n characters.

**`StringTake["string", -n]`**

Gives the last n characters.

**`StringTake["string", {n}]`**

Gives the nth character.

**`StringTake["string", {m, n}]`**

Gives characters m through n.

**`StringTake["string", {m, n, s}]`**

Gives characters m through n in steps of s.

**`StringTake["string", UpTo[n]]`**

Gives n characters, or as many as are available.

**`StringTake[{s1, s2, ...}, spec]`**

Gives the list of results for each si.

## Examples (21)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (11)

```mathematica
In[1]:= StringTake["abcdefghijklm", 6]
Out[1]= "abcdef"

In[2]:= StringTake["abcdefghijklm", -4]
Out[2]= "jklm"

In[3]:= StringTake["abcdefghijklm", {5, 10}]
Out[3]= "efghij"

In[4]:= StringTake["abcdefghijklm", {6}]
Out[4]= "f"

In[5]:= StringTake["abcdefghijklm", {1, -1, 2}]
Out[5]= "acegikm"

In[6]:= StringTake[{"abcdef", "stuv", "xyzw"}, -2]
Out[6]= {"ef", "uv", "zw"}

In[7]:= StringTake["abc", UpTo[4]]
Out[7]= "abc"

In[8]:= StringTake["abcdef", {-3, -1}]
Out[8]= "def"

In[9]:= StringTake["abcde", {5, 1, -1}]
Out[9]= "edcba"

In[10]:= StringTake["abc", 0]
Out[10]= ""

In[11]:= StringTake[x, 1]
Out[11]= StringTake[x, 1]
```

### Applications (10)

The first three characters

```mathematica
In[12]:= StringTake["abcdef", 3]
Out[12]= "abc"
```

The last two

```mathematica
In[13]:= StringTake["abcdef", -2]
Out[13]= "ef"
```

Characters 2 through 4

```mathematica
In[14]:= StringTake["abcdef", {2, 4}]
Out[14]= "bcd"
```

Every second character

```mathematica
In[15]:= StringTake["abcdefgh", {1, -1, 2}]
Out[15]= "aceg"
```

A list of specs gives a list of substrings

```mathematica
In[16]:= StringTake["abcdef", {{1, 2}, {4, 5}}]
Out[16]= {"ab", "de"}
```

A list of strings, the same spec for each

```mathematica
In[17]:= StringTake[{"apple", "kiwi"}, 2]
Out[17]= {"ap", "ki"}
```

At most 10 characters, without an error

```mathematica
In[18]:= StringTake["abc", UpTo[10]]
Out[18]= "abc"
```

Positions count characters, not UTF-8 bytes

```mathematica
In[19]:= StringTake["h\:00e9llo", 2]
Out[19]= "hé"
```

The inverse of StringPosition

```mathematica
In[20]:= StringTake["abcdefgh", StringPosition["abcdefgh", "cde"][[1]]]
Out[20]= "cde"
```

Out of range: StringTake::take, left unevaluated

```mathematica
In[21]:= StringTake["abc", 5]
Out[21]= StringTake["abc", 5]
```

## Implementation notes

**Algorithm.** `builtin_stringtake` parses the spec with `st_parse_spec`. The spec is either one sequence spec or a list of them. The one-spec forms are `n`, `-n`, `{n}`, `{m, n}`, `{m, n, s}`, `All`, `None` and `UpTo[n]`, decoded by the shared `strings_seq_parse` into a `StrSeq`. Two Mathematica quirks are kept: `{m, n, 0}` and a list of four or more integers are read as lists of specs, the latter with a `StringTake::ambgsntx` warning. `strings_take_seq` resolves the spec against the character length with `strings_seq_resolve`, which applies the Take range rules exactly (negative positions count from the end). A forward contiguous range is a single `memcpy` of its byte span. Any other step gathers characters through a `U8View` built over just that span. A list of strings maps the scalar kernel over the elements directly, without building a new `StringTake` call for each one.

**Data structures.** `StrChars` / `U8View` from `src/strings/utf8.c`. Positions count Unicode code points. ASCII is detected eight bytes at a time and uses byte arithmetic. Only the prefix a spec needs is scanned (`strings_seq_prefix`), so `StringTake[s, 5]` on a megabyte of non-ASCII text costs microseconds.

**Complexity / limits.** `O(span)` per take. A malformed spec emits `StringTake::mseqs`, an out-of-range one `StringTake::take`, and a non-string `::strse`; each leaves the call unevaluated. In a list of strings or a list of specs, a failing entry becomes an unevaluated `StringTake[si, spec]` inside the result. `Span` (`m ;; n`) is not accepted. `Protected`.

**Attributes:** `Protected`.

## References

- Source: [`src/strings/stringtake.c`](https://github.com/stblake/mathilda/blob/main/src/strings/stringtake.c)
- Specification: [`docs/spec/builtins/string-operations.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/string-operations.md)
- Tests: [`tests/test_graphics.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graphics.c)
- Tests: [`tests/test_string_review.c`](https://github.com/stblake/mathilda/blob/main/tests/test_string_review.c)
- Tests: [`tests/test_stringposition.c`](https://github.com/stblake/mathilda/blob/main/tests/test_stringposition.c)
- Tests: [`tests/test_strings.c`](https://github.com/stblake/mathilda/blob/main/tests/test_strings.c)

## Notes & additional examples

### Notes

- Every Mathematica sequence spec is accepted: `n`, `-n`, `{n}`, `{m, n}`,
  `{m, n, s}`, `All`, `None` and `UpTo[n]`, or a list of them. The Take range
  rules are Mathematica's exactly, so `StringTake["abc", {2, 1}]` is `""`.
- Positions count Unicode code points (since v0.456). Only the prefix a spec
  needs is scanned, so a short take from a very long non-ASCII string is cheap.
- Errors route through `mth_message`, so `Quiet`/`Check` see them:
  `StringTake::take` (out of range), `::mseqs` (not a spec) and `::strse` (not
  a string). Inside a list result, a failing entry stays as an unevaluated
  `StringTake[...]`.
- A list of four or more integers is read as a list of specs, with a
  `StringTake::ambgsntx` warning, as in Mathematica. `Span` (`m ;; n`) is not
  accepted.
- See also: `StringDrop`, `StringTakeDrop`, `StringPart`, `StringPosition`, `Take`.
