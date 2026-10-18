# StringContainsQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`StringContainsQ["string", patt]`**

Gives True if any substring of "string" matches the string expression patt, and False otherwise.

**`StringContainsQ["string", {p1, p2, ...}]`**

Gives True if any substring matches any of the pi.

**`StringContainsQ[{s1, s2, ...}, patt]`**

Gives the list of results for each of the si.

**`StringContainsQ[patt]`**

Represents an operator form that can be applied to a string. Equivalent to !StringFreeQ\["string", patt\], and to StringMatchQ\["string", \_\_\_ ~~ patt ~~ \_\_\_\]. Options: IgnoreCase -\> True treats upper/lowercase as equivalent.

## Examples (14)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= StringContainsQ["bcde", "b" ~~ __ ~~ "e"]
Out[1]= True

In[2]:= StringContainsQ[{"a", "b", "ab", "abcd", "bcde"}, "a"]
Out[2]= {True, False, True, True, False}

In[3]:= StringContainsQ["bac 123", RegularExpression["a.*"] ~~ DigitCharacter ..]
Out[3]= True

In[4]:= Select[{"abc", "xyz", "bat"}, StringContainsQ["a"]]
Out[4]= {"abc", "bat"}
```

### Options (1)

```mathematica
In[5]:= StringContainsQ["abcd", "BC", IgnoreCase -> True]
Out[5]= True
```

### Worked examples (1)

```mathematica
In[6]:= Options[StringContainsQ]
Out[6]= {IgnoreCase -> False}
```

### Applications (8)

Some substring matches

```mathematica
In[7]:= StringContainsQ["abcdef", "cd"]
Out[7]= True
```

Any string pattern

```mathematica
In[8]:= StringContainsQ["abc123", DigitCharacter]
Out[8]= True
```

Alternatives

```mathematica
In[9]:= StringContainsQ["abc", "x" | "c"]
Out[9]= True
```

Threads over a list of subjects

```mathematica
In[10]:= StringContainsQ[{"apple", "banana", "cherry"}, "an"]
Out[10]= {False, True, False}
```

Case folded

```mathematica
In[11]:= StringContainsQ["Mathilda", "math", IgnoreCase -> True]
Out[11]= True
```

Case folding is Unicode-aware

```mathematica
In[12]:= StringContainsQ["caf\:00e9", "\:00c9", IgnoreCase -> True]
Out[12]= True
```

A regex anchor still means the start of the string

```mathematica
In[13]:= StringContainsQ["abc", RegularExpression["^b"]]
Out[13]= False
```

Operator form as a predicate

```mathematica
In[14]:= Select[{"cat", "dog", "cow"}, StringContainsQ["o"]]
Out[14]= {"dog", "cow"}
```

## Implementation notes

**Algorithm.** The four substring predicates `StringContainsQ`, `StringFreeQ`, `StringStartsQ` and `StringEndsQ` share one core, `sq_dispatch`, and differ only by an `SqKind` tag. It seeds `IgnoreCase` from the registered `Options`, so `SetOptions` takes effect, and strips trailing option rules. With a single positional argument it returns the curried operator form `Function[head[#1, patt, opts...]]`. Otherwise it builds the rule set unanchored and asks `sq_match`: one `regex_match` per rule, stopping at the first hit, since a predicate never needs to know where or how often. `StringStartsQ`/`StringEndsQ` anchor at the pattern level by wrapping it as `StartOfString ~~ p` or `p ~~ EndOfString`. An empty pattern list `{}` matches nothing and is answered directly.

**Data structures.** A `RegexRule` array from the 64-entry compiled-pattern LRU cache. The cache keys on structural equality, so a freshly built anchoring wrapper still hits it. Only the whole-match pair is requested: no `regex_scan` enumeration and no capture pool.

**Complexity / limits.** One PCRE2 probe per rule, with early exit. The predicate is equivalent to `!StringFreeQ[s, p]` and to `StringMatchQ[s, ___ ~~ p ~~ ___]`. A non-string subject, or any non-string element of a subject list, leaves the whole call unevaluated rather than threading a wrong Boolean. Matching is Unicode-aware, including `IgnoreCase` folding of non-ASCII letters.

**Attributes:** `Protected`.

## References

**See also:** [StringMatchQ](../../string-operations/StringMatchQ/), [SetOptions](../../assignment-and-rules/SetOptions/)

- Source: [`src/strings/regex/stringcontainsq.c`](https://github.com/stblake/mathilda/blob/main/src/strings/regex/stringcontainsq.c)
- Specification: [`docs/spec/builtins/string-operations.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/string-operations.md)
- Tests: [`tests/test_parallelmixedspecial.c`](https://github.com/stblake/mathilda/blob/main/tests/test_parallelmixedspecial.c)
- Tests: [`tests/test_string_review.c`](https://github.com/stblake/mathilda/blob/main/tests/test_string_review.c)
- Tests: [`tests/test_stringcontainsq.c`](https://github.com/stblake/mathilda/blob/main/tests/test_stringcontainsq.c)
- Tests: [`tests/test_stringfns.c`](https://github.com/stblake/mathilda/blob/main/tests/test_stringfns.c)

## Notes & additional examples

### Notes

- `StringContainsQ` searches inside the string, where `StringMatchQ` must match
  the whole of it. It is equivalent to `StringMatchQ[s, ___ ~~ p ~~ ___]` and to
  `!StringFreeQ[s, p]`.
- The one-argument operator form `StringContainsQ[p]` (with any
  `IgnoreCase -> ...` it was given) curries into a predicate, which suits
  `Select`.
- It stops at the first match, with one PCRE2 probe per rule. It shares its core
  with `StringFreeQ`, `StringStartsQ` and `StringEndsQ`.
- A non-string subject, or a list holding one, leaves the whole call unevaluated
  rather than threading a wrong Boolean. An empty pattern list matches nothing.
- See also: `StringMatchQ`, `StringFreeQ`, `StringStartsQ`, `StringEndsQ`.
