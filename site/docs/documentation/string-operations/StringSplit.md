# StringSplit

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`StringSplit["string"]`**

Splits "string" at runs of whitespace.

**`StringSplit["string", patt]`**

Splits at delimiters matching the string pattern patt.

**`StringSplit["string", {p1, p2, ...}]`**

Splits at any of the pi.

**`StringSplit["string", patt -> val]`**

Inserts val at the position of each delimiter.

**`StringSplit["string", patt, n]`**

Splits into at most n substrings.

**`StringSplit[{s1, s2, ...}, patt]`**

Gives the list of results for each of the si. Empty substrings between adjacent interior delimiters are kept; those at the start or end are dropped unless All is given as the third argument. "" splits at every character. Option: IgnoreCase -\> True.

## Examples (13)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= StringSplit["a bbb  cccc aa   d"]
Out[1]= {"a", "bbb", "cccc", "aa", "d"}

In[2]:= StringSplit["a-b:c-d:e-f-g", {":", "-"}]
Out[2]= {"a", "b", "c", "d", "e", "f", "g"}

In[3]:= StringSplit["a b::c d::e f g", "::" -> "--"]
Out[3]= {"a b", "--", "c d", "--", "e f g"}

In[4]:= StringSplit["This is a sentence, which goes on.", Except[WordCharacter] ..]
Out[4]= {"This", "is", "a", "sentence", "which", "goes", "on"}
```

### Applications (9)

Default: split at runs of whitespace

```mathematica
In[5]:= StringSplit["a b c"]
Out[5]= {"a", "b", "c"}
```

Interior empty pieces are kept

```mathematica
In[6]:= StringSplit["a,b,,c", ","]
Out[6]= {"a", "b", "", "c"}
```

All keeps the leading and trailing empty pieces too

```mathematica
In[7]:= StringSplit[",a,b,", ",", All]
Out[7]= {"", "a", "b", ""}
```

At most two pieces

```mathematica
In[8]:= StringSplit["a b c d", " ", 2]
Out[8]= {"a", "b c d"}
```

Alternative delimiters

```mathematica
In[9]:= StringSplit["a-b_c", "-" | "_"]
Out[9]= {"a", "b", "c"}
```

A rule inserts a value bound from each delimiter

```mathematica
In[10]:= StringSplit["a1b22c", x:DigitCharacter.. :> FromDigits[x]]
Out[10]= {"a", 1, "b", 22, "c"}
```

Case-insensitive delimiters

```mathematica
In[11]:= StringSplit["aXbxc", "x", IgnoreCase -> True]
Out[11]= {"a", "b", "c"}
```

The empty delimiter splits into Unicode characters

```mathematica
In[12]:= StringSplit["h\:00e9llo", ""]
Out[12]= {"h", "é", "l", "l", "o"}
```

Threads over a list of subjects

```mathematica
In[13]:= StringSplit[{"a b", "c d"}]
Out[13]= {{"a", "b"}, {"c", "d"}}
```

## Algorithm

stringsplit.c - StringSplit[...], the full Wolfram surface.

```text
  StringSplit[s]                 split at runs of whitespace
  StringSplit[s, patt]           split at delimiters matching patt
  StringSplit[s, {p1, p2, ...}]  split at any of the pi
  StringSplit[s, patt -> val]    insert val at each delimiter
  StringSplit[s, patt, n]        at most n substrings
  StringSplit[s, patt, All]      keep leading/trailing empty substrings
  StringSplit[{s1, ...}, patt]   thread over a list of subjects
  IgnoreCase -> True             case-insensitive delimiters
```

The delimiter pattern is translated to PCRE by the shared string-pattern engine (regex_common.c / string_pattern.c), so it accepts literal strings, RegularExpression["re"], the character-class heads (Whitespace, ...), StringExpression (~~), Alternatives (|), Repeated (..), Except, and so on.

Zero-length substrings between two adjacent interior delimiters are kept; empty substrings at the very beginning or end are dropped unless All is given. The empty-string delimiter "" splits at every character.

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

**Algorithm.** `builtin_stringsplit` seeds `IgnoreCase` from `Options[StringSplit]` and strips trailing `IgnoreCase` rules. It takes 1–3 positional arguments and gives `StringSplit::argb` otherwise. The delimiter (default `Whitespace`) is compiled by the shared engine (`regex_rules_build_ex`). `ss_scalar` walks the subject with a `RegexCursor`, taking the leftmost match over all delimiters with ties going to list order. It pushes each inter-delimiter substring, keeping interior empties. For a `patt -> val` or `patt :> val` rule it then pushes the inserted value from `regex_rule_eval_rhs`, the binder shared with `StringCases`/`StringReplace`. So `x:DigitCharacter.. :> FromDigits[x]` inserts integers bound from the match. A third argument `n` stops after `n - 1` cuts, so the remainder is the last piece. `All` keeps the leading and trailing empty substrings, which are dropped by default. The empty delimiter `""` short-circuits to `ss_null`, a split into Unicode characters. A list of subjects threads, and non-string elements are copied through.

**Data structures.** A `PieceVec` of `(Expr, is_sub)` records. The `is_sub` flag marks a substring piece, which is eligible for end-trimming, as distinct from an inserted rule value, which never is. The cursor caches each rule's next match, and rule sets come from the 64-entry compiled-pattern LRU cache.

**Complexity / limits.** One left-to-right pass of `O(len · rules)` probes. A zero-width delimiter advances one code point, so the scan always progresses. The delimiter accepts the full string-pattern language, and PCRE2 runs in UTF + UCP mode. A third argument that is neither an integer nor `All`, or a non-string subject, leaves the call unevaluated.

**Attributes:** `Protected`.

## References

**See also:** [StringReplace](../../string-operations/StringReplace/)

- Source: [`src/strings/stringsplit.c`](https://github.com/stblake/mathilda/blob/main/src/strings/stringsplit.c)
- Specification: [`docs/spec/builtins/string-operations.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/string-operations.md)
- Tests: [`tests/test_string_review.c`](https://github.com/stblake/mathilda/blob/main/tests/test_string_review.c)
- Tests: [`tests/test_stringfns.c`](https://github.com/stblake/mathilda/blob/main/tests/test_stringfns.c)

## Notes & additional examples

### Notes

- `StringSplit` returns the pieces between non-overlapping delimiter matches.
  Interior empty pieces are kept, and leading or trailing ones are dropped
  unless the third argument is `All`. An integer third argument `n` gives at
  most `n` pieces, with the remainder as the last one.
- With a rule `patt -> val` or `patt :> val`, the value is inserted at each
  delimiter. Pattern names are bound to the text they matched (v0.455, issue
  #98), so `x:DigitCharacter.. :> FromDigits[x]` inserts the numbers.
- The delimiter accepts the full string-pattern language, including
  `RegularExpression`. `IgnoreCase` is a registered option.
- `StringSplit[s, ""]` splits into Unicode code points, the same as `Characters`.
- See also: `StringRiffle` (the inverse), `StringCases`, `Characters`,
  `StringPartition`.
