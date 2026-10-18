# StringTrim

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`StringTrim["string"]`**

Trims whitespace from the beginning and end of "string".

**`StringTrim["string", patt]`**

Removes one substring matching the string pattern patt from the beginning and one from the end (use patt.. to remove runs).

**`StringTrim[{s1, s2, ...}, ...]`**

Gives the list of results for each of the si. Whitespace covers runs of spaces, tabs, and newlines. Options: IgnoreCase -\> True.

## Examples (12)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= StringTrim["   aaa bbb ccc   "]
Out[1]= "aaa bbb ccc"

In[2]:= StringTrim["++++aaa bbb ccc----", ("+" | "-") ..]
Out[2]= "aaa bbb ccc"

In[3]:= StringTrim["   aaa bbb ccc   ", RegularExpression["^ *"]]
Out[3]= "aaa bbb ccc   "

In[4]:= StringTrim["007bond007", DigitCharacter ..]
Out[4]= "bond"
```

### Applications (8)

Default: whitespace from both ends

```mathematica
In[5]:= StringTrim["   hello   "]
Out[5]= "hello"
```

Tabs and newlines are whitespace too

```mathematica
In[6]:= StringTrim["\n\t text \n"]
Out[6]= "text"
```

One match is removed at each end, as in Mathematica

```mathematica
In[7]:= StringTrim["xxabxx", "x"]
Out[7]= "xabx"
```

Use a repeated pattern to strip a whole run

```mathematica
In[8]:= StringTrim["xxabxx", "x"..]
Out[8]= "ab"
```

Any string pattern works

```mathematica
In[9]:= StringTrim["123abc456", DigitCharacter..]
Out[9]= "abc"
```

The two cuts never overlap

```mathematica
In[10]:= StringTrim["xxx", "xx"]
Out[10]= "x"
```

Case-insensitive

```mathematica
In[11]:= StringTrim["XXhelloxx", "x", IgnoreCase -> True]
Out[11]= "Xhellox"
```

Threads over a list of strings

```mathematica
In[12]:= StringTrim[{"  a ", " b  "}]
Out[12]= {"a", "b"}
```

## Algorithm

stringtrim.c - StringTrim[...], trims matching substrings from both ends.

```text
  StringTrim["string"]           trims whitespace runs from start and end
  StringTrim["string", patt]     trims a substring matching patt from each end
  StringTrim[{s1, s2, ...}, ...] threads over a list of subject strings
  IgnoreCase -> True             case-insensitive pattern
```

As in Mathematica, ONE match is removed at each end -- StringTrim is StringReplace[s, {StartOfString ~~ patt -> "", patt ~~ EndOfString -> ""}]: StringTrim["xxabcxx", "x"] is "xabcx" (use "x".. to strip runs). The back match must start at or after the end of the front match, so StringTrim["xxx", "xx"] is "x".

Both ends go through the shared string-pattern engine (regex_common.c), so the full pattern language works, Condition and PatternTest included:

```text
  - front: the pattern compiled REGEX_ANCHORED, one match at offset 0 -- no
    scan of the subject at all;
  - back: the pattern compiled REGEX_ENDANCHORED, searched from the end of the
    front match; PCRE rejects a candidate as soon as it cannot reach the end.
```

## Implementation notes

**Algorithm.** As in Mathematica, `StringTrim[s, p]` is `StringReplace[s, {StartOfString ~~ p -> "", p ~~ EndOfString -> ""}]`: it removes one match at each end and does not strip to a fixed point. `builtin_stringtrim` takes the subject, an optional pattern (default `Whitespace`, i.e. `\s+`) and `IgnoreCase`. It builds the pattern twice through the shared engine (`regex_rules_build_opts`): once compiled `REGEX_ANCHORED` for the front and once `REGEX_ENDANCHORED` for the back. `st_scalar` runs one anchored match at offset 0 to find the front cut. It then searches for the back match from the end of the front match, which PCRE rejects as soon as a candidate cannot reach the end, so the two cuts never cross (`StringTrim["xxx", "xx"]` is `"x"`).

**Data structures.** Two cached `RegexRule` sets (front and back) and one output buffer for the window `[start, end)`. The full string-pattern language works here, `Condition` and `PatternTest` included.

**Complexity / limits.** Two PCRE2 matches per subject, and the front one never scans. Strip a run with a repeated pattern such as `"x"..`: a single-character pattern removes only one character per end. A list of subjects threads. A non-string element in it emits `StringTrim::strse` and leaves the call unevaluated. A build without PCRE2 warns with `regavail`.

**Attributes:** `Protected`.

## References

- Source: [`src/strings/stringtrim.c`](https://github.com/stblake/mathilda/blob/main/src/strings/stringtrim.c)
- Specification: [`docs/spec/builtins/string-operations.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/string-operations.md)
- Tests: [`tests/test_string_review.c`](https://github.com/stblake/mathilda/blob/main/tests/test_string_review.c)
- Tests: [`tests/test_stringfns.c`](https://github.com/stblake/mathilda/blob/main/tests/test_stringfns.c)

## Notes & additional examples

### Notes

- `StringTrim[s, p]` is `StringReplace[s, {StartOfString ~~ p -> "", p ~~ EndOfString -> ""}]`.
  It removes **one** match at each end and never touches the middle. Before
  v0.455 it stripped each end to a fixed point, so `StringTrim["xxabxx", "x"]`
  wrongly gave `"ab"`.
- The default pattern `Whitespace` is itself a run (`\s+`), so the
  one-argument form still removes all leading and trailing whitespace.
- The front match is compiled anchored and the back match end-anchored, so
  neither end needs a scan of the subject.
- A list holding a non-string emits `StringTrim::strse` and stays unevaluated.
- See also: `StringReplace`, `StringDelete`, `StringPadLeft`.
