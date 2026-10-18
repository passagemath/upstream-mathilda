# StringMatchQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`StringMatchQ["string", patt]`**

Gives True if the whole "string" matches patt, and False otherwise.

**`StringMatchQ[{s1, s2, ...}, patt]`**

Gives the list of results for each of the si. patt may be any string pattern, RegularExpression\["re"\], a literal string, or a list of alternatives. In a literal string \* matches any run of characters and @ one or more non-uppercase characters (\\\* and \\@ are literal). Options: IgnoreCase -\> True.

## Examples (14)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= StringMatchQ["12345", RegularExpression["\\d+"]]
Out[1]= True

In[2]:= StringMatchQ[{"12", "x"}, RegularExpression["\\d+"]]
Out[2]= {True, False}

In[3]:= StringMatchQ["apple.txt", "*.txt"]
Out[3]= True

In[4]:= StringMatchQ["aXb", x_ ~~ "X" ~~ x_]
Out[4]= False
```

### Applications (10)

The whole string must match; ___ is zero or more characters

```mathematica
In[5]:= StringMatchQ["hello", "h" ~~ ___]
Out[5]= True
```

Letters followed by digits

```mathematica
In[6]:= StringMatchQ["abc123", LetterCharacter .. ~~ DigitCharacter ..]
Out[6]= True
```

Threads over a list of subjects

```mathematica
In[7]:= StringMatchQ[{"cat", "dog"}, "c" ~~ __]
Out[7]= {True, False}
```

The * wildcard inside a literal string

```mathematica
In[8]:= StringMatchQ["apple.txt", "*.txt"]
Out[8]= True
```

The @ wildcard: one or more characters that are not uppercase

```mathematica
In[9]:= StringMatchQ["Hello", "H@"]
Out[9]= True
```

Case-insensitive

```mathematica
In[10]:= StringMatchQ["Hello", "hello", IgnoreCase -> True]
Out[10]= True
```

A repeated name must match the same text

```mathematica
In[11]:= {StringMatchQ["aXa", x_ ~~ "X" ~~ x_], StringMatchQ["aXb", x_ ~~ "X" ~~ x_]}
Out[11]= {True, False}
```

A Condition is evaluated during the match

```mathematica
In[12]:= {StringMatchQ["12", x__ /; ToExpression[x] > 10], StringMatchQ["5", x__ /; ToExpression[x] > 10]}
Out[12]= {True, False}
```

Character classes are Unicode-aware

```mathematica
In[13]:= StringMatchQ["\:00fcber", LetterCharacter..]
Out[13]= True
```

Operator form as a predicate

```mathematica
In[14]:= Select[{"2024", "v2", "17"}, StringMatchQ[DigitCharacter..]]
Out[14]= {"2024", "17"}
```

## Algorithm

stringmatchq.c - StringMatchQ[subject, pattern]

Returns True if the WHOLE subject string matches the pattern, else False. The pattern may be RegularExpression["re"], a literal string, or a List of

```text
alternatives (matches if any one matches).  A list of subjects threads,
```

giving a list of True/False.

## Implementation notes

**Algorithm.** `builtin_stringmatchq` seeds `IgnoreCase` from `Options[StringMatchQ]` and strips a trailing `IgnoreCase` rule. A single positional argument is the operator form, applied by `src/opform.c`. The pattern is built by `regex_rules_build_opts` with `REGEX_RULE_ANCHORED | REGEX_RULE_WILDCARDS`. The anchored flag wraps each rule as `\A(?:...)\z`, with any leading PCRE `(*VERB)` hoisted out of the wrapper, so a match must cover the whole subject. The wildcards flag turns `*` (any run of characters) and `@` (one or more characters that are not uppercase letters) into metacharacters inside literal strings; `\\*` and `\\@` stay literal. `smq_scalar` runs `regex_match` at offset 0 against each rule and stops at the first hit. A `List` pattern is a set of alternatives. A list of subjects threads, and non-string elements are copied through.

**Data structures.** A `RegexRule` array from the shared compiled-pattern cache: a 64-entry refcounted LRU keyed by the pattern's structural hash plus the compile flags. `Map[StringMatchQ[#, p] &, words]` therefore compiles `p` once. Only the whole-match offset pair is requested, so no capture pool and no `regex_scan` enumeration are needed.

**Complexity / limits.** One PCRE2 match per rule per subject, with early exit. PCRE2 runs in UTF + UCP mode with JIT, so `LetterCharacter`, `DigitCharacter` and `IgnoreCase` are Unicode-aware. Pattern names are capture groups, and a repeated name is a backreference (`x_ ~~ "X" ~~ x_`). `Condition` and general `PatternTest` are PCRE2 callouts into the evaluator. A non-string subject leaves the call unevaluated.

**Attributes:** `Protected`.

## References

- Source: [`src/strings/regex/stringmatchq.c`](https://github.com/stblake/mathilda/blob/main/src/strings/regex/stringmatchq.c)
- Specification: [`docs/spec/builtins/string-operations.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/string-operations.md)
- Tests: [`tests/test_mathar_regressions.c`](https://github.com/stblake/mathilda/blob/main/tests/test_mathar_regressions.c)
- Tests: [`tests/test_string_review.c`](https://github.com/stblake/mathilda/blob/main/tests/test_string_review.c)
- Tests: [`tests/test_stringcontainsq.c`](https://github.com/stblake/mathilda/blob/main/tests/test_stringcontainsq.c)
- Tests: [`tests/test_stringfns.c`](https://github.com/stblake/mathilda/blob/main/tests/test_stringfns.c)

## Notes & additional examples

### Notes

- `StringMatchQ` asks whether the *whole* string matches. The pattern is wrapped
  `\A(?:...)\z`. Use `StringContainsQ` for a substring test.
- In a literal-string pattern `*` stands for any run of characters and `@` for one
  or more characters that are not uppercase letters. Write `"\\*"` / `"\\@"` for
  the literal characters.
- The full string-pattern language is accepted: `~~`, `|`, `..`, the character
  classes, `Shortest`/`Longest`, `Except`, named patterns as backreferences,
  `Condition`, and an arbitrary `PatternTest`, which runs as a callout into the
  evaluator.
- A list of patterns is a set of alternatives. A non-string subject leaves the
  call unevaluated, while a non-string element of a subject list is passed
  through unchanged.
- Compiled patterns are cached. Mapping one pattern over 40k words dropped from
  556 ms to about 40 ms in v0.456.
- See also: `StringContainsQ`, `StringStartsQ`, `StringFreeQ`, `MatchQ`.
