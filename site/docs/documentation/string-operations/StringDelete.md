# StringDelete

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`StringDelete["string", patt]`**

Deletes every substring of "string" that matches patt.

**`StringDelete[{s1, s2, ...}, patt]`**

Gives the list of results for each of the si.

**`StringDelete[patt]`**

Represents an operator form that can be applied to a string. Equivalent to StringReplace\["string", patt -\> ""\]. Options: IgnoreCase -\> True.

## Examples (12)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= StringDelete["a1b22c", DigitCharacter..]
Out[1]= "abc"

In[2]:= StringDelete[{"x1", "y22"}, DigitCharacter]
Out[2]= {"x", "y"}
```

### Options (1)

```mathematica
In[3]:= StringDelete["Hello World", "l", IgnoreCase -> True]
Out[3]= "Heo Word"
```

### Applications (9)

```mathematica
In[4]:= StringDelete["a1b22c", DigitCharacter..]
Out[4]= "abc"

In[5]:= StringDelete["The quick brown fox", "o"]
Out[5]= "The quick brwn fx"
```

Deletes l and L

```mathematica
In[6]:= StringDelete["Hello World", "l", IgnoreCase -> True]
Out[6]= "Heo Word"
```

A list of patterns

```mathematica
In[7]:= StringDelete["abcabc", {"a", "c"}]
Out[7]= "bb"
```

Keep only letters, Unicode-aware

```mathematica
In[8]:= StringDelete["naïve café!", Except[LetterCharacter]]
Out[8]= "naïvecafé"
```

```mathematica
In[9]:= StringDelete["foo123bar", RegularExpression["[0-9]+"]]
Out[9]= "foobar"

In[10]:= StringDelete["aXbYc", _?UpperCaseQ]
Out[10]= "abc"
```

One result per string

```mathematica
In[11]:= StringDelete[{"x1", "y22"}, DigitCharacter]
Out[11]= {"x", "y"}
```

Operator form

```mathematica
In[12]:= StringDelete[Whitespace]["a b  c"]
Out[12]= "abc"
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

## Implementation notes

**Algorithm.** `builtin_stringdelete` is `sr_dispatch(res, "StringDelete", SR_DELETE)`, the driver shared with `StringReplace`. `sr_args` strips trailing `IgnoreCase` options (with the default taken from `Options[StringDelete]`); a one-argument call is the operator form `StringDelete[patt][s]`, resolved by `src/opform.c`, and a third positional argument is rejected (no occurrence limit). A subject that is not a string or list of strings emits `StringDelete::strse`. `delete_rules` turns the pattern, or each element of a list of patterns, into `p -> ""`; the rule set is compiled by `regex_rules_build_opts`, which translates each Wolfram string pattern to PCRE2 source in `string_pattern.c` and compiles it (UTF + UCP + DOTALL, plus `REGEX_CASELESS` for `IgnoreCase -> True`). `sr_scalar` then walks the subject with a `RegexCursor`: at each position the leftmost match of any rule wins (ties by rule order), the text before it is kept and the match dropped; a zero-width match advances by one character.

**Data structures.** Compiled rule sets (`RegexRule` arrays) live in a 64-entry LRU cache keyed by the pattern expression and flags, so repeated calls with the same pattern skip translation and compilation. Output accumulates in a growable `RegexBuf`. `Condition` and general `PatternTest` pieces run Mathilda code mid-match through PCRE2 callouts; `_?LetterQ`, `_?UpperCaseQ` etc. become plain character classes.

**Complexity / limits.** One left-to-right PCRE2 scan per subject, linear in the subject for ordinary patterns (backtracking-heavy regexes inherit PCRE2's worst cases). Without PCRE2 in the build (`USE_REGEX` off) the regex layer is a stub and the call stays unevaluated.

**Attributes:** `Protected`.

## References

- Source: [`src/strings/regex/stringreplace.c`](https://github.com/stblake/mathilda/blob/main/src/strings/regex/stringreplace.c)
- Specification: [`docs/spec/builtins/string-operations.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/string-operations.md)
- Tests: [`tests/test_string_review.c`](https://github.com/stblake/mathilda/blob/main/tests/test_string_review.c)

## Notes & additional examples

### Notes

- `StringDelete[s, p]` is `StringReplace[s, p -> ""]`: each pattern is wrapped in a
  rule and the whole set goes through the shared PCRE2-backed `StringReplace`
  scanner (non-overlapping, left to right).
- Any string pattern is accepted — literals, `StringExpression`, `Alternatives`,
  character classes, `RegularExpression`, `PatternTest`, `Condition`.
- Option `IgnoreCase`. Unlike `StringReplace` there is no occurrence-limit argument;
  a third positional argument leaves the call unevaluated.
- A non-string subject emits `StringDelete::strse`. Attributes `{Protected}`.
- See also: `StringReplace`, `StringTrim`, `StringCases`.
