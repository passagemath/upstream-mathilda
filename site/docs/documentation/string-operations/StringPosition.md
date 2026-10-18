# StringPosition

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`StringPosition["string", patt]`**

Gives a list of the {start, end} character positions at which substrings matching the string pattern patt occur in "string".

**`StringPosition["string", patt, n]`**

Includes only the first n occurrences.

**`StringPosition["string", {p1, p2, ...}]`**

Gives positions of all the pi.

**`StringPosition[{s1, s2, ...}, patt]`**

Threads over a list of strings. Positions use the form consumed by StringTake / StringReplacePart. Options: Overlaps -\> True (default; overlaps allowed, one substring per start), False (no overlaps), or All (every matching substring); IgnoreCase -\> True treats upper/lowercase as equivalent.

## Examples (13)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= StringPosition["abXYZaaabXYZaaaaXYZXYZ", "XYZ"]
Out[1]= {{3, 5}, {10, 12}, {17, 19}, {20, 22}}

In[2]:= StringPosition["AABBBAABABBCCCBAAA", x_ ~~ x_]
Out[2]= {{1, 2}, {3, 4}, {4, 5}, {6, 7}, {10, 11}, {12, 13}, {13, 14}, {16, 17}, {17, 18}}
```

### Options (2)

```mathematica
In[3]:= StringPosition["AAAAA", "AA", Overlaps -> False]
Out[3]= {{1, 2}, {3, 4}}

In[4]:= StringPosition["abAB", "a", IgnoreCase -> True]
Out[4]= {{1, 1}, {3, 3}}
```

### Applications (9)

{start, end} pairs, 1-based inclusive

```mathematica
In[5]:= StringPosition["abcabc", "bc"]
Out[5]= {{2, 3}, {5, 6}}
```

Overlapping matches are listed by default

```mathematica
In[6]:= StringPosition["aaaa", "aa"]
Out[6]= {{1, 2}, {2, 3}, {3, 4}}
```

Non-overlapping, as StringCases would give

```mathematica
In[7]:= StringPosition["aaaa", "aa", Overlaps -> False]
Out[7]= {{1, 2}, {3, 4}}
```

A character-class pattern

```mathematica
In[8]:= StringPosition["a1b2c3", DigitCharacter]
Out[8]= {{2, 2}, {4, 4}, {6, 6}}
```

Only the first match

```mathematica
In[9]:= StringPosition["abcabc", "bc", 1]
Out[9]= {{2, 3}}
```

Case-insensitive

```mathematica
In[10]:= StringPosition["AbAB", "ab", IgnoreCase -> True]
Out[10]= {{1, 2}, {3, 4}}
```

Positions count characters, so the é is one position

```mathematica
In[11]:= StringPosition["h\:00e9llo", "llo"]
Out[11]= {{3, 5}}
```

Operator form

```mathematica
In[12]:= StringPosition["xyz"]["axyzbxyz"]
Out[12]= {{2, 4}, {6, 8}}
```

Positions feed StringTake directly

```mathematica
In[13]:= StringTake["abcdefgh", StringPosition["abcdefgh", "cde"][[1]]]
Out[13]= "cde"
```

## Algorithm

stringposition.c - StringPosition[subject, pattern, n]

Returns a List of {start, end} character-position pairs at which substrings of

```text
`subject` match the string pattern `pattern`, in the 1-based inclusive form
consumed by StringTake / StringDrop / StringReplacePart.  The pattern may be a
```

literal string, a general string expression (Blank/Pattern/~~/RegularExpression /character classes), or a List of patterns; a List of subjects threads.

Options:

```text
  Overlaps -> True (default) | False | All
    True  - include overlapping substrings, but only the first (natural) match
            starting at each position.
    False - exclude overlapping substrings (greedy left-to-right, global).
    All   - include every matching substring at every start (all lengths).
  IgnoreCase -> True | False (default)
    Treat upper/lowercase as equivalent.
```

A third positional integer argument n keeps only the first n matches.

The match enumeration itself is regex_scan() in regex_common.c, shared with StringCases and StringCount; this file only turns spans into position pairs.

Positions are CHARACTER positions (UTF-8 code points), consistent with StringLength / StringTake: the scan's byte offsets are mapped through the subject's character-offset table (a binary search; identity for ASCII).

## Implementation notes

**Algorithm.** `builtin_stringposition` seeds `IgnoreCase` and `Overlaps` from `Options[StringPosition]`, then strips trailing option rules. A single positional argument is the operator form (`src/opform.c`), and an optional third integer `n` keeps the first `n` matches. `sp_scalar` runs the shared `regex_scan_limit`, which `StringCases` and `StringCount` also use. It turns each byte span `[ms, me)` into the 1-based inclusive character pair `{ms + 1, me}`, the form that `StringTake`, `StringDrop` and `StringReplacePart` consume. Only the default differs from its siblings: `Overlaps -> True`, as in the Wolfram Language, so overlapping matches at distinct starts are all listed.

**Data structures.** The `RegexScan` span array, and a `U8View` over the subject that maps byte offsets to character indices. For ASCII the map is the identity, and otherwise it is a binary search over the offset table. Rule sets come from the 64-entry compiled-pattern LRU cache.

**Complexity / limits.** One scan, plus `O(log len)` per reported endpoint on non-ASCII text. Since v0.456 positions count Unicode code points, so they agree with `StringLength`/`StringTake` on accented text. A list of subjects threads, and a wrong argument count emits `StringPosition::argb`.

**Attributes:** `Protected`.

## References

**See also:** [StringTake](../../string-operations/StringTake/), [StringDrop](../../string-operations/StringDrop/), [StringReplacePart](../../string-operations/StringReplacePart/), [StringCases](../../string-operations/StringCases/), [SetOptions](../../assignment-and-rules/SetOptions/)

- Source: [`src/strings/regex/stringposition.c`](https://github.com/stblake/mathilda/blob/main/src/strings/regex/stringposition.c)
- Specification: [`docs/spec/builtins/string-operations.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/string-operations.md)
- Tests: [`tests/test_string_review.c`](https://github.com/stblake/mathilda/blob/main/tests/test_string_review.c)
- Tests: [`tests/test_stringcontainsq.c`](https://github.com/stblake/mathilda/blob/main/tests/test_stringcontainsq.c)
- Tests: [`tests/test_stringcount.c`](https://github.com/stblake/mathilda/blob/main/tests/test_stringcount.c)
- Tests: [`tests/test_stringfns.c`](https://github.com/stblake/mathilda/blob/main/tests/test_stringfns.c)

## Notes & additional examples

### Notes

- Results are `{start, end}` character positions, 1-based and inclusive: the
  form `StringTake`, `StringDrop` and `StringReplacePart` consume.
- As in the Wolfram Language, the default is `Overlaps -> True`, which lists
  overlapping matches at distinct starts. `StringCases` and `StringCount`
  default to `False`. `Overlaps -> All` lists every matching substring.
- Since v0.456 positions count Unicode code points rather than UTF-8 bytes, so
  they agree with `StringLength` on accented text.
- A third argument `n` keeps the first `n` matches. A one-argument call is the
  operator form. A list of subjects threads.
- See also: `StringCases`, `StringCount`, `StringTake`, `Position`.
