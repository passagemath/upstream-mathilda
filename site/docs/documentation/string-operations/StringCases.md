# StringCases

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`StringCases["string", patt]`**

Gives the list of non-overlapping substrings of "string" that match patt, from left to right.

**`StringCases["string", patt -> rhs]`**

Gives rhs for each match, with the pattern names of patt bound to the text they matched (x\_\_ ~~ "1" :\> x). For a RegularExpression pattern a string rhs may use $n (n-th group) and $0 (whole match).

**`StringCases["string", {p1, p2, ...}]`**

Gives the matches of any of the pi.

**`StringCases[{s1, s2, ...}, patt]`**

Gives the list of results for each of the si. Options: Overlaps -\> False (default; overlapping substrings are not treated as separate), True (overlaps separate, one substring per start), or All (every matching substring at every start); IgnoreCase -\> True treats upper/lowercase as equivalent.

## Examples (18)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= StringCases["a13b12c17a32", RegularExpression["[^a1]"]]
Out[1]= {"3", "b", "2", "c", "7", "3", "2"}

In[2]:= StringCases["AaBBccDDeefG", RegularExpression["[[:upper:]]+"]]
Out[2]= {"A", "BB", "DD", "G"}

In[3]:= StringCases["abc123", (x__ ~~ "123") :> x]
Out[3]= {"abc"}

In[4]:= StringCases["a1b22c333", x:DigitCharacter.. /; StringLength[x] > 1 :> ToExpression[x]]
Out[4]= {22, 333}
```

### Options (3)

```mathematica
In[5]:= StringCases["AAAA", "AA", Overlaps -> True]
Out[5]= {"AA", "AA", "AA"}

In[6]:= StringCases["AAA", x__, Overlaps -> All]
Out[6]= {"AAA", "AA", "A", "AA", "A", "A"}

In[7]:= StringCases["aAbB", "a", IgnoreCase -> True]
Out[7]= {"a", "A"}
```

### Applications (11)

Runs of letters, left to right

```mathematica
In[8]:= StringCases["the cat sat", LetterCharacter ..]
Out[8]= {"the", "cat", "sat"}
```

A named piece of each match, bound into the rule

```mathematica
In[9]:= StringCases["abc123", (x__ ~~ "123") :> x]
Out[9]= {"abc"}
```

The right-hand side is evaluated, so these are integers

```mathematica
In[10]:= StringCases["a1b22c333", x:DigitCharacter.. :> ToExpression[x]]
Out[10]= {1, 22, 333}
```

Parse key=value pairs

```mathematica
In[11]:= StringCases["x=1, y=22", (k:LetterCharacter ~~ "=" ~~ v:DigitCharacter..) :> k -> ToExpression[v]]
Out[11]= {"x" -> 1, "y" -> 22}
```

With RegularExpression a string RHS may use $n

```mathematica
In[12]:= StringCases["a1b2c3", RegularExpression["[a-z](\\d)"] -> "$1"]
Out[12]= {"1", "2", "3"}
```

An occurrence limit

```mathematica
In[13]:= StringCases["one two three four", WordCharacter.., 2]
Out[13]= {"one", "two"}
```

Overlapping matches, one per start

```mathematica
In[14]:= StringCases["AAAA", "AA", Overlaps -> True]
Out[14]= {"AA", "AA", "AA"}
```

Every matching substring

```mathematica
In[15]:= StringCases["abcd", __, Overlaps -> All]
Out[15]= {"abcd", "abc", "ab", "a", "bcd", "bc", "b", "cd", "c", "d"}
```

Greedy by default; Shortest is lazy

```mathematica
In[16]:= {StringCases["aab abb", "a" ~~ __ ~~ "b"], StringCases["aab abb", Shortest["a" ~~ __ ~~ "b"]]}
Out[16]= {{"aab abb"}, {"aab", "abb"}}
```

Letter classes are Unicode-aware

```mathematica
In[17]:= StringCases["h\:00e9llo w\:00f6rld", LetterCharacter..]
Out[17]= {"héllo", "wörld"}
```

Operator form

```mathematica
In[18]:= StringCases[LetterCharacter..]["hi there"]
Out[18]= {"hi", "there"}
```

## Algorithm

stringcases.c - StringCases[subject, pattern]

Returns a List of the substrings of `subject` that match `pattern`, left to

```text
right.  With a rule pattern (patt -> rhs / patt :> rhs) each match is replaced
```

by the rhs with the pattern's names bound to the pieces they matched (StringCases["abc123", x__ ~~ "123" :> x] = {"abc"}); for a RegularExpression

```text
pattern a string rhs may use $0/$1... for the whole match and capture groups.  The
```

pattern may also be a List of alternatives/rules; at each position the

```text
leftmost match wins, ties broken by rule order.  A list of subjects threads.
```

Options:

```text
  Overlaps -> False (default) | True | All
    False - non-overlapping, greedy left-to-right.
    True  - overlapping substrings count separately, but only the first match
            starting at each position.
    All   - every matching substring at every start (all lengths).
  IgnoreCase -> True | False (default)
    Treat upper/lowercase as equivalent.
```

The match enumeration itself is regex_scan() in regex_common.c, shared with StringCount and StringPosition, so StringCount[s, p, opts] always equals Length[StringCases[s, p, opts]].

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

**Algorithm.** `builtin_stringcases` seeds `Overlaps`/`IgnoreCase` from `Options[StringCases]`, so `SetOptions` takes effect, then strips trailing option rules. A single positional argument is the operator form, applied by `src/opform.c`. An optional third argument `n` (an integer `>= 0` or `Infinity`) keeps the first `n` matches, and `n = 0` gives `{}` at once. The rule set is built by `regex_rules_build_ex`, anchored only in `Overlaps -> All` mode, where whole candidate substrings are tested. `sc_scalar` calls `regex_scan_limit` once per subject. A bare pattern gives the matched substring. A `patt -> rhs` / `patt :> rhs` rule goes through `regex_rule_eval_rhs`. For a `RegularExpression` pattern a string RHS is a `$0`/`$n` template. Otherwise each pattern name is bound to the text its capture group matched (`""` for an untaken alternative) through `replace_bindings`, and the RHS is evaluated. So `(x__ ~~ "123") :> x` gives `"abc"`, and `x:DigitCharacter.. :> ToExpression[x]` gives integers. This is the issue #98 fix in v0.455: before it, names were dropped and the whole match came back.

**Data structures.** `regex_scan` fills a `RegexScan`: a span array, plus a flat capture pool that spans index by offset. `want_captures` is set only when some rule has an RHS, so plain extraction allocates no pool. For `Overlaps -> False` the scan uses the shared `RegexCursor`, which caches each rule's next match, and an occurrence limit stops the scan early. Rule sets come from the 64-entry compiled-pattern LRU cache.

**Complexity / limits.** `False`/`True` modes are one pass. `All` mode probes `O(len² · rules)` anchored substrings. `StringCount` and `StringPosition` share `regex_scan`, so the three agree on every `Overlaps` policy. PCRE2 runs in UTF + UCP mode, and offsets are byte offsets on code-point boundaries. A non-string subject is left unevaluated; inside a subject list, non-string elements are copied through.

**Attributes:** `Protected`.

## References

**See also:** [RegularExpression](../../string-operations/RegularExpression/), [SetOptions](../../assignment-and-rules/SetOptions/)

- Source: [`src/strings/regex/stringcases.c`](https://github.com/stblake/mathilda/blob/main/src/strings/regex/stringcases.c)
- Specification: [`docs/spec/builtins/string-operations.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/string-operations.md)
- Tests: [`tests/test_string_review.c`](https://github.com/stblake/mathilda/blob/main/tests/test_string_review.c)
- Tests: [`tests/test_stringcontainsq.c`](https://github.com/stblake/mathilda/blob/main/tests/test_stringcontainsq.c)
- Tests: [`tests/test_stringcount.c`](https://github.com/stblake/mathilda/blob/main/tests/test_stringcount.c)
- Tests: [`tests/test_stringfns.c`](https://github.com/stblake/mathilda/blob/main/tests/test_stringfns.c)

## Notes & additional examples

### Notes

- With a rule `patt :> rhs`, every pattern name in `patt` is bound to the text
  it matched, and `rhs` is evaluated. A name in an untaken alternative is bound
  to `""`. Until v0.455 the names were ignored and the whole match was returned
  (issue #98).
- `$0`/`$n` templates apply only to `RegularExpression` patterns. Elsewhere a
  string RHS is literal.
- `Overlaps -> False` (the default) is non-overlapping and greedy. `True` gives
  overlapping matches, the first one at each start. `All` gives every matching
  substring of every length.
- `StringCases[s, p, n]` keeps the first `n` matches, and the scan stops once it
  has them.
- `IgnoreCase` and `Overlaps` defaults come from `Options[StringCases]`, so
  `SetOptions` changes them.
- `StringCount` and `StringPosition` use the same scanner, so
  `StringCount[s, p, opts] == Length[StringCases[s, p, opts]]` always holds.
- A list of subjects threads, and non-string elements pass through unchanged.
- See also: `StringCount`, `StringPosition`, `StringReplace`, `Cases`.
