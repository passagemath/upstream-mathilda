# RegularExpression

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`RegularExpression["regex"]`**

Represents a class of strings given by the PCRE regular expression "regex", for use in StringMatchQ, StringCases, StringReplace and StringSplit.  It is an inert head: it evaluates to itself. Supported syntax includes . \[c1c2\] \[c1-c2\] \[^...\] p\* p+ p? p{m,n}, non-greedy \*? +? ??, groups (...) and alternation |; the classes \d \D \s \S \w \W and \[\[:name:\]\]; the anchors ^ $ \b \B; and inline options (?i) (?m) (?s).  In a replacement right-hand side $n stands for the n-th captured group and $0 for the whole match.

## Examples (14)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= StringCases["adefgh12c34", RegularExpression["[a-e]+"]]
Out[1]= {"ade", "c"}

In[2]:= StringCases["a23b4222c63333d80", RegularExpression["\\d+"]]
Out[2]= {"23", "4222", "63333", "80"}

In[3]:= StringReplace["John Smith", RegularExpression["(\\w+) (\\w+)"] -> "$2 $1"]
Out[3]= "Smith John"

In[4]:= StringCases["caf\:00e9 na\:00efve", RegularExpression["\\w+"]]
Out[4]= {"café", "naïve"}
```

### Applications (10)

A run of digits

```mathematica
In[5]:= StringCases["x1y22z333", RegularExpression["\\d+"]]
Out[5]= {"1", "22", "333"}
```

$n refers to a captured group

```mathematica
In[6]:= StringReplace["John Smith", RegularExpression["(\\w+) (\\w+)"] -> "$2 $1"]
Out[6]= "Smith John"
```

A backreference inside the regex; $0 is the whole match

```mathematica
In[7]:= StringReplace["hello", RegularExpression["(.)\\1"] -> "<$0>"]
Out[7]= "he<ll>o"
```

Letters then digits, whole string

```mathematica
In[8]:= StringMatchQ["abc123", RegularExpression["[a-z]+\\d+"]]
Out[8]= True
```

The digits are the delimiters

```mathematica
In[9]:= StringSplit["a1b2c3", RegularExpression["\\d"]]
Out[9]= {"a", "b", "c"}
```

Word boundaries and counted repeats

```mathematica
In[10]:= StringCases["The year 1999 and 2024", RegularExpression["\\b\\d{4}\\b"]]
Out[10]= {"1999", "2024"}
```

\w is Unicode-aware

```mathematica
In[11]:= StringCases["caf\:00e9 na\:00efve", RegularExpression["\\w+"]]
Out[11]= {"café", "naïve"}
```

A regex dot skips newlines; the blank _ does not

```mathematica
In[12]:= {StringCases["one\ntwo", RegularExpression["e.t"]], StringCases["one\ntwo", "e" ~~ _ ~~ "t"]}
Out[12]= {{}, {"e\nt"}}
```

A regex mixes with symbolic string patterns

```mathematica
In[13]:= StringCases["ab12cd", RegularExpression["[a-z]+"] ~~ DigitCharacter]
Out[13]= {"ab1"}
```

A bad pattern emits RegularExpression::regex

```mathematica
In[14]:= Check[RegularExpression["a(b"], "invalid"]
Out[14]= "invalid"
```

## Algorithm

regularexpression.c - the RegularExpression[...] head.

RegularExpression["re"] is an inert data head: it carries a PCRE pattern string for use by StringMatchQ / StringCases / StringReplace / StringSplit and does not evaluate to anything else (exactly like Wolfram Language). The builtin therefore validates its argument and otherwise returns NULL so

```text
the expression survives evaluation unchanged.  When the pattern does not
```

compile a RegularExpression::regex diagnostic is emitted (still inert), and when Mathilda was built without PCRE2 a RegularExpression::regavail note is printed the first time an invalid-but-present call is seen.

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

**Algorithm.** `RegularExpression["re"]` is an inert data head. `builtin_regularexpression` always returns `NULL`, so the expression survives evaluation unchanged and carries raw PCRE source to the string functions that consume it (`StringMatchQ`, `StringCases`, `StringReplace`, `StringSplit`, `StringPosition`, `StringCount`, …). For a single string argument it does a syntax check through the shared rule builder (`regex_rules_build_opts`). A bad pattern emits `RegularExpression::regex` through `mth_message`, which `Check`/`Quiet` see, and the head stays inert. A valid pattern is compiled once and left in the compiled-pattern cache under exactly the key the string functions look up. So `StringCases[s, RegularExpression[re]]` does not compile `re` a second time.

**Data structures.** The translator (`string_pattern.c`) passes the source through verbatim. When the whole pattern is the regex, it is flagged `raw_regex` and compiled without `DOTALL`, so `.` keeps PCRE's meaning and does not match a newline. Inside a larger string expression the regex is wrapped `(?-s:...)`. `regex_compile_ex` (`regex_engine.c`) calls `pcre2_compile` with `PCRE2_UTF | PCRE2_UCP | PCRE2_MATCH_INVALID_UTF`, so `\w`, `\d` and caseless matching are Unicode-aware. It then runs a best-effort `pcre2_jit_compile`. Leading `(*VERB)`s are hoisted out of any anchoring wrapper, and `IgnoreCase` is a real `PCRE2_CASELESS` flag, not a `(?i)` prefix. A string RHS of a rule whose pattern contains a `RegularExpression` is a template: `$0` is the whole match, `$n` the n-th group, and `$$` a literal `$` (`regex_expand_template`).

**Complexity / limits.** Matching cost is PCRE2's. Up to 64 offset pairs are exposed (`REGEX_MAX_PAIRS`). Backslashes are doubled in a Mathilda string literal (`"\\d"` is the regex `\d`). `$n` templates apply only to patterns that contain a `RegularExpression`. Elsewhere `"$1"` is literal text, and named string-pattern captures (`x:p :> f[x]`) are the way to reuse a match. Without PCRE2 (`USE_REGEX` unset) the engine is stubbed out and every string-pattern builtin warns `regavail` and stays unevaluated.

**Attributes:** `Protected`.

## References

**See also:** [StringMatchQ](../../string-operations/StringMatchQ/), [StringCases](../../string-operations/StringCases/), [StringReplace](../../string-operations/StringReplace/), [StringSplit](../../string-operations/StringSplit/)

- P. Hazel, *PCRE2 — Perl-Compatible Regular Expressions*, 10.x, https://www.pcre.org/.
- Source: [`src/strings/regex/regularexpression.c`](https://github.com/stblake/mathilda/blob/main/src/strings/regex/regularexpression.c)
- Specification: [`docs/spec/builtins/string-operations.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/string-operations.md)
- Tests: [`tests/test_names.c`](https://github.com/stblake/mathilda/blob/main/tests/test_names.c)
- Tests: [`tests/test_regex.c`](https://github.com/stblake/mathilda/blob/main/tests/test_regex.c)
- Tests: [`tests/test_stringcontainsq.c`](https://github.com/stblake/mathilda/blob/main/tests/test_stringcontainsq.c)
- Tests: [`tests/test_stringcount.c`](https://github.com/stblake/mathilda/blob/main/tests/test_stringcount.c)

## Notes & additional examples

### Notes

- `RegularExpression["re"]` is an inert head carrying PCRE2 source. Evaluating
  it compiles the pattern once, which is a syntax check that also warms the
  compiled-pattern cache the string functions use.
- In the right-hand side of a rule whose pattern contains a `RegularExpression`,
  `$0` is the whole match, `$n` the n-th group, and `$$` a literal `$`. With
  purely symbolic patterns, use named pieces (`x:p :> ...`) instead. There
  `"$1"` is literal text (v0.455).
- PCRE2 runs in UTF + UCP mode, so `\w`, `\d`, `[[:alpha:]]` and `(?i)` are
  Unicode-aware. A bare `RegularExpression` keeps PCRE's own `.`, which does not
  match a newline. The symbolic blanks `_`/`__` do match newlines, as in
  Mathematica. Use `(?s)` to change the regex behaviour.
- Backslashes are doubled in a Mathilda string literal (`"\\d"` is the regex
  `\d`). An invalid pattern emits `RegularExpression::regex` through the message
  funnel, so `Check` and `Quiet` see it, and stays inert.
- See also: `StringCases`, `StringReplace`, `StringMatchQ`, `StringSplit`,
  `StringExpression`.
