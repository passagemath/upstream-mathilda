# StringReplaceList

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`StringReplaceList["string", rules]`**

Gives the list of strings obtained by replacing ONE match of the rules at a time, for every match start of every rule, in order of position.

**`StringReplaceList["string", rules, n]`**

Gives at most n results.

**`StringReplaceList[{s1, s2, ...}, rules]`**

Gives the list of results for each of the si. Options: IgnoreCase -\> True.

## Examples (13)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= StringReplaceList["aaa", "a" -> "b"]
Out[1]= {"baa", "aba", "aab"}

In[2]:= StringReplaceList["aaa", "aa" -> "X"]
Out[2]= {"Xa", "aX"}

In[3]:= StringReplaceList["abc", {"a" -> "1", "b" -> "2"}]
Out[3]= {"1bc", "a2c"}
```

### Applications (10)

Each single replacement separately

```mathematica
In[4]:= StringReplaceList["aaa", "a" -> "b"]
Out[4]= {"baa", "aba", "aab"}
```

Overlapping matches each give a result

```mathematica
In[5]:= StringReplaceList["aaa", "aa" -> "X"]
Out[5]= {"Xa", "aX"}
```

Several rules, ordered by position

```mathematica
In[6]:= StringReplaceList["abc", {"a" -> "1", "b" -> "2"}]
Out[6]= {"1bc", "a2c"}
```

At most two results

```mathematica
In[7]:= StringReplaceList["aaaa", "a" -> "b", 2]
Out[7]= {"baaa", "abaa"}
```

A named pattern used on the right-hand side

```mathematica
In[8]:= StringReplaceList["a1b2", d:DigitCharacter :> "<" <> d <> ">"]
Out[8]= {"a<1>b2", "a1b<2>"}
```

```mathematica
In[9]:= StringReplaceList["aAa", "a" -> "-", IgnoreCase -> True]
Out[9]= {"-Aa", "a-a", "aA-"}

In[10]:= StringReplaceList["αβα", "α" -> "a"]
Out[10]= {"aβα", "αβa"}
```

One list per string

```mathematica
In[11]:= StringReplaceList[{"ab", "ba"}, "a" -> "X"]
Out[11]= {{"Xb"}, {"bX"}}
```

No match gives an empty list

```mathematica
In[12]:= StringReplaceList["abc", "z" -> "y"]
Out[12]= {}
```

A non-string value gives a StringExpression

```mathematica
In[13]:= StringReplaceList["ab", "b" -> f[1]]
Out[13]= {"a"~~f[1]}
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

**Algorithm.** `builtin_stringreplacelist` is `sr_dispatch(res, "StringReplaceList", SR_LIST)`. `sr_args` strips `IgnoreCase` options and decodes the optional limit `n` (integer ≥ 0 or `Infinity`; anything else leaves the call unevaluated). The subject must be a string or list of strings (`::strse` otherwise), and every element of the rule set must carry a right-hand side — a bare pattern returns unevaluated. The rules are compiled by `regex_rules_build_opts` (Wolfram string pattern → PCRE2, UTF/UCP). For each subject, `sr_list_scalar` calls `regex_scan_limit` in `REGEX_OV_TRUE` (overlapping) mode with captures: one pass per rule finds a match at every start position, the per-rule runs are merged by a stable sort (`scan_stable_sort`) into position-then-rule order, and the list is cut to `n`. Each span then yields a fresh string — the prefix, the rule's right-hand side instantiated with the match's named captures by `regex_rule_eval_rhs`, and the suffix. A non-string value turns that result into a `StringExpression` (`sr_finish`), and `Sequence[...]` values splice.

**Data structures.** A `RegexScan` of `RegexSpan {ms, me, rule, caps_off, npairs}` records with a shared capture-offset array, and a per-result `SrOut` builder. Compiled rule sets come from the same 64-entry LRU cache as the rest of the string-pattern family.

**Complexity / limits.** The overlapping scan does not stop early: the limit `n` is applied after sorting, so a small `n` on a long subject still pays the full scan, and the output is `O(k · L)` for `k` matches in a string of length `L`. `n = 0` short-circuits to `{}`.

**Attributes:** `Protected`.

## References

- Source: [`src/strings/regex/stringreplace.c`](https://github.com/stblake/mathilda/blob/main/src/strings/regex/stringreplace.c)
- Specification: [`docs/spec/builtins/string-operations.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/string-operations.md)
- Tests: [`tests/test_string_review.c`](https://github.com/stblake/mathilda/blob/main/tests/test_string_review.c)

## Notes & additional examples

### Notes

- Enumerates every match start of every rule with overlapping allowed
  (`regex_scan_limit` in `Overlaps -> True` mode), stably sorted by position and then
  rule order, and builds one string per match with only that match replaced.
- The limit `n` (integer ≥ 0 or `Infinity`) truncates after sorting, so it costs a
  full scan; `n = 0` returns `{}` immediately.
- The second argument must be rules: a bare pattern leaves the call unevaluated.
- Option `IgnoreCase`; a non-string subject emits `::strse`. Attributes `{Protected}`.
- See also: `StringReplace`, `ReplaceList`, `StringCases`.
