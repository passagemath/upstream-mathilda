# StringReplace

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`StringReplace["string", patt -> rep]`**

Replaces each non-overlapping match of patt in "string" by rep, with the pattern names of patt bound to the text they matched; a non-string rep makes the result a StringExpression. For a RegularExpression pattern a string rep may use $n and $0.

**`StringReplace["string", {patt1 -> rep1, patt2 -> rep2, ...}]`**

Applies a list of replacement rules; at each position the leftmost match wins, ties broken by rule order.

**`StringReplace[{s1, s2, ...}, rules]`**

Gives the list of results for each of the si.

## Examples (16)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (5)

```mathematica
In[1]:= StringReplace["a13b12c1da32efg", RegularExpression["(\\d+)"] -> "[$1]"]
Out[1]= "a[13]b[12]c[1]da[32]efg"

In[2]:= StringReplace["123 45 6 789", RegularExpression["\\b"] :> "X"]
Out[2]= "X123X X45X X6X X789X"

In[3]:= StringReplace["abc123", a__ ~~ "123" :> a <> "!"]
Out[3]= "abc!"

In[4]:= StringReplace["abc123", x:DigitCharacter :> ToExpression[x]]
Out[4]= "abc"~~1~~2~~3

In[5]:= StringReplace["aaaa", "a" -> "b", 2]
Out[5]= "bbaa"
```

### Applications (11)

Every match is replaced

```mathematica
In[6]:= StringReplace["hello world", "o" -> "0"]
Out[6]= "hell0 w0rld"
```

Several rules; the leftmost match wins at each position

```mathematica
In[7]:= StringReplace["the cat", {"the" -> "a", "cat" -> "dog"}]
Out[7]= "a dog"
```

A named piece of the match is bound into the right-hand side

```mathematica
In[8]:= StringReplace["abc123", a__ ~~ "123" :> a <> "!"]
Out[8]= "abc!"
```

The right-hand side is evaluated per match

```mathematica
In[9]:= StringReplace["2+3=5", d:DigitCharacter :> ToString[2 ToExpression[d]]]
Out[9]= "4+6=10"
```

Non-string values turn the result into a StringExpression

```mathematica
In[10]:= StringReplace["abc123", x:DigitCharacter :> ToExpression[x]]
Out[10]= "abc"~~1~~2~~3
```

With RegularExpression a string RHS is a $n template

```mathematica
In[11]:= StringReplace["John Smith", RegularExpression["(\\w+) (\\w+)"] -> "$2 $1"]
Out[11]= "Smith John"
```

Outside RegularExpression, $0 is literal text

```mathematica
In[12]:= StringReplace["abc", "b" -> "$0"]
Out[12]= "a$0c"
```

An occurrence limit: only the first two matches

```mathematica
In[13]:= StringReplace["aaaa", "a" -> "b", 2]
Out[13]= "bbaa"
```

Case-insensitive matching

```mathematica
In[14]:= StringReplace["Hello hello", "hello" -> "bye", IgnoreCase -> True]
Out[14]= "bye bye"
```

Operator form

```mathematica
In[15]:= StringReplace["x" -> "y"]["xox"]
Out[15]= "yoy"
```

A zero-width anchor inserts without consuming text

```mathematica
In[16]:= StringReplace["one\ntwo", StartOfLine -> "> "]
Out[16]= "> one\n> two"
```

## Algorithm

stringreplace.c - StringReplace[subject, rule | {rules...}]

Replaces every non-overlapping match of a rule's pattern by its right-hand

```text
side, scanning left to right; unmatched text is copied verbatim.  Pattern
```

names in the LHS are bound to the text they matched and substituted into the RHS (regex_rule_eval_rhs); a non-string value turns the result into a

```text
StringExpression.  For a RegularExpression pattern a string RHS may use
```

$0/$1... for the whole match and capture groups. With several rules, at each position the leftmost match wins, ties broken by

```text
rule order.  A list of subjects threads.
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

**Algorithm.** `builtin_stringreplace` goes through `sr_dispatch` (shared with `StringDelete` and `StringReplaceList`). `sr_args` seeds `IgnoreCase` from `Options[StringReplace]`, strips trailing `IgnoreCase` rules, and decodes an optional occurrence limit `n` (an integer `>= 0` or `Infinity`); a single positional argument is the operator form, applied by `src/opform.c`. A subject that is not a string or a list of strings emits `StringReplace::strse`. The rules are built unanchored by `regex_rules_build_opts`, and every element must carry an RHS or the call stays unevaluated. `sr_scalar` then walks the subject with a `RegexCursor`. At each position `regex_cursor_next` returns the leftmost match over all rules, ties going to rule order. The literal text before the match is copied, and the match is replaced by `regex_rule_eval_rhs`. That helper is shared with `StringCases`/`StringSplit`. For a `RegularExpression` pattern a string RHS is a `$0`/`$n` template. Otherwise each pattern name is bound to the text it matched (`""` for an untaken alternative), substituted with the matcher's scoping-aware `replace_bindings`, and the RHS evaluated. A zero-width match copies one character (UTF-8 aware, `regex_next_char`) and moves on. Once `n` replacements are done, the rest of the subject is copied unchanged.

**Data structures.** An `SrOut` accumulator: a `RegexBuf` text run, plus a piece array that is only used once a replacement value is not a string. The result is then a `StringExpression` of text runs and values, as in Mathematica: `"abc"~~1~~2~~3`. A `Sequence[...]` value splices in. Compiled rule sets come from a 64-entry refcounted LRU cache keyed by the pattern's structural hash and the compile flags, so mapping one rule over many strings compiles it once. The `RegexCursor` caches each rule's next match, so a scan with several rules costs `O(n · rules)` probes, not `O(n²)`.

**Complexity / limits.** One left-to-right pass. PCRE2 runs in UTF + UCP mode with JIT. Translated patterns are compiled `DOTALL`, so the blanks match newlines, while a bare `RegularExpression` keeps PCRE's own `.`. `Condition` and general `PatternTest` run as PCRE2 callouts that can backtrack. `$n` is literal text outside `RegularExpression` patterns (`"b" -> "$0"` inserts `$0`). A list of subjects threads; a non-string element triggers `strse`.

**Attributes:** `Protected`.

## References

**See also:** [StringCases](../../string-operations/StringCases/), [RegularExpression](../../string-operations/RegularExpression/), [StringExpression](../../other-advanced/StringExpression/)

- Source: [`src/strings/regex/stringreplace.c`](https://github.com/stblake/mathilda/blob/main/src/strings/regex/stringreplace.c)
- Specification: [`docs/spec/builtins/string-operations.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/string-operations.md)
- Tests: [`tests/test_string_review.c`](https://github.com/stblake/mathilda/blob/main/tests/test_string_review.c)
- Tests: [`tests/test_stringcontainsq.c`](https://github.com/stblake/mathilda/blob/main/tests/test_stringcontainsq.c)
- Tests: [`tests/test_stringfns.c`](https://github.com/stblake/mathilda/blob/main/tests/test_stringfns.c)

## Notes & additional examples

### Notes

- `StringReplace` scans left to right and replaces each non-overlapping match
  with its rule's right-hand side, copying the text in between. With several
  rules the leftmost match wins at each position, and ties go to rule order.
- Pattern names on the left (`x:p`, `x_`, `x__`) are bound to the text they
  matched and substituted into the right-hand side before it is evaluated, as in
  Mathematica. This was issue #98, fixed in v0.455; before the fix, names were
  ignored.
- When a replacement value is not a string, the result is a `StringExpression`
  of the pieces (`"abc"~~1~~2~~3`), exactly as in Mathematica. `Sequence[...]`
  values splice in.
- `$0`/`$n` templates work only when the pattern contains a `RegularExpression`.
- `StringReplace[s, rules, n]` replaces at most `n` matches. A one-argument call
  is the operator form. `IgnoreCase` is a registered option, so `SetOptions`
  changes its default.
- A subject that is not a string or a list of strings emits
  `StringReplace::strse` and stays unevaluated. A list of subjects threads.
- Compiled rule sets are cached (a 64-entry LRU), so mapping one replacement
  over many strings compiles the pattern once.
- See also: `StringReplaceList`, `StringDelete`, `StringCases`, `StringSplit`,
  `RegularExpression`.
