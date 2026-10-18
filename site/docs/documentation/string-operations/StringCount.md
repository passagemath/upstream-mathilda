# StringCount

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`StringCount["string", "sub"]`**

Gives the number of times "sub" appears as a substring of "string".

**`StringCount["string", patt]`**

Gives the number of substrings of "string" matching the string expression patt.

**`StringCount["string", {p1, p2, ...}]`**

Counts the occurrences of any of the pi.

**`StringCount[{s1, s2, ...}, patt]`**

Gives the list of results for each of the si. Equivalent to Length\[StringCases\[...\]\] but does not build the matched substrings. Options: Overlaps -\> False (default; overlapping substrings are not counted as separate), True (overlaps counted separately, one substring per start), or All (every matching substring at every start); IgnoreCase -\> True treats upper/lowercase as equivalent.

## Examples (14)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= StringCount["the cat sat on the mat", "at"]
Out[1]= 3

In[2]:= StringCount["AAAA", "AA"]
Out[2]= 2

In[3]:= StringCount[{"a1", "b22", "ccc"}, DigitCharacter]
Out[3]= {1, 2, 0}

In[4]:= StringCount["x=1.5, y=-2", NumberString]
Out[4]= 2
```

### Options (3)

```mathematica
In[5]:= StringCount["AAAA", "AA", Overlaps -> True]
Out[5]= 3

In[6]:= StringCount["AAAA", x__, Overlaps -> All]
Out[6]= 10

In[7]:= StringCount["aAbB", "a", IgnoreCase -> True]
Out[7]= 2
```

### Applications (7)

Non-overlapping by default

```mathematica
In[8]:= StringCount["mississippi", "ss"]
Out[8]= 2
```

Maximal runs of digits

```mathematica
In[9]:= StringCount["a1 b22 c333", DigitCharacter..]
Out[9]= 3
```

Overlapping occurrences are counted only on request

```mathematica
In[10]:= {StringCount["banana", "ana"], StringCount["banana", "ana", Overlaps -> True]}
Out[10]= {1, 2}
```

Every nonempty substring: 4 + 3 + 2 + 1

```mathematica
In[11]:= StringCount["abcd", __, Overlaps -> All]
Out[11]= 10
```

Case-insensitive

```mathematica
In[12]:= StringCount["Hello HELLO", "hello", IgnoreCase -> True]
Out[12]= 2
```

One count per subject

```mathematica
In[13]:= StringCount[{"aa", "aba", "c"}, "a"]
Out[13]= {2, 2, 0}
```

Non-ASCII patterns

```mathematica
In[14]:= StringCount["h\:00e9h\:00e9", "\:00e9"]
Out[14]= 2
```

## Algorithm

stringcount.c - StringCount[subject, pattern]

```text
Counts the substrings of `subject` that match `pattern`.  This is the
```

counting-only companion of StringCases: it runs the same match enumeration (regex_scan in regex_common.c) but never materialises a substring, so it costs one small span record per match instead of a malloc + Expr + List element.

The pattern may be a literal string, a general string expression (Blank/Pattern/~~/RegularExpression/character classes), a Rule or RuleDelayed (only the LHS matters -- the replacement is irrelevant to a count), or a List

```text
of any of those.  A List of subjects threads, giving one count per subject.
```

Options:

```text
  Overlaps -> False (default) | True | All
    False - overlapping substrings are not counted separately.
    True  - overlapping substrings count separately, but only the first
            matching substring at a given position is counted.
    All   - every matching substring at every position is counted separately.
  IgnoreCase -> True | False (default)
    Treat upper/lowercase as equivalent.
```

Because the scan is shared, StringCount[s, p, opts] is always exactly Length[StringCases[s, p, opts]] and Length[StringPosition[s, p, opts]].

Byte semantics: like the rest of src/strings, offsets are byte offsets (no UTF-8 codepoint decoding), consistent with StringLength / StringPart.

## Implementation notes

**Algorithm.** `builtin_stringcount` reads its options and builds its rule set exactly as `StringCases` does, seeding `Overlaps -> False` and `IgnoreCase -> False` from `Options[StringCount]` and stripping trailing option rules. Then `sct_scalar` calls `regex_scan` with `want_captures = 0` and returns the number of spans as an `Integer`. No substring is ever created, so each match costs one small span record instead of a string `Expr` and a list slot. A `Rule`/`RuleDelayed` pattern is accepted, but only its left-hand side matters to a count. A list of subjects threads, giving one count each.

**Data structures.** The shared `RegexScan`: a span array and no capture pool. For `Overlaps -> False` the scan uses the `RegexCursor` leftmost-match cursor. Rule sets come from the 64-entry compiled-pattern LRU cache.

**Complexity / limits.** Because the scan is shared, `StringCount[s, p, opts]` always equals `Length[StringCases[s, p, opts]]` and `Length[StringPosition[s, p, opts]]` under the same `Overlaps` setting. Note that `StringPosition` defaults to `Overlaps -> True`. `Overlaps -> All` counts every matching substring at every start, testing `O(len²)` substrings. A positional argument count other than two emits `StringCount::argrx`. Matching is Unicode-aware (PCRE2 UTF + UCP).

**Attributes:** `Protected`.

## References

**See also:** [StringCases](../../string-operations/StringCases/), [Rule](../../assignment-and-rules/Rule/), [RuleDelayed](../../assignment-and-rules/RuleDelayed/), [SetOptions](../../assignment-and-rules/SetOptions/), [StringPosition](../../string-operations/StringPosition/)

- Source: [`src/strings/regex/stringcount.c`](https://github.com/stblake/mathilda/blob/main/src/strings/regex/stringcount.c)
- Specification: [`docs/spec/builtins/string-operations.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/string-operations.md)
- Tests: [`tests/test_compile_assoc.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile_assoc.c)
- Tests: [`tests/test_string_review.c`](https://github.com/stblake/mathilda/blob/main/tests/test_string_review.c)
- Tests: [`tests/test_stringcontainsq.c`](https://github.com/stblake/mathilda/blob/main/tests/test_stringcontainsq.c)
- Tests: [`tests/test_stringcount.c`](https://github.com/stblake/mathilda/blob/main/tests/test_stringcount.c)

## Notes & additional examples

### Notes

- `StringCount` is the counting-only companion of `StringCases`. It runs the
  same match scan but never builds a substring.
- Because the scan is shared, `StringCount[s, p, opts]` is always
  `Length[StringCases[s, p, opts]]`. It also equals
  `Length[StringPosition[s, p, opts]]` under the same `Overlaps` setting, but
  note that `StringPosition` defaults to `Overlaps -> True`.
- A rule pattern is accepted, and only its left-hand side matters.
  `Overlaps`/`IgnoreCase` defaults come from `Options[StringCount]`.
- A positional argument count other than two emits `StringCount::argrx`.
- See also: `StringCases`, `StringPosition`, `Count`.
