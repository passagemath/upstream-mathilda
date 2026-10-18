# LowerCaseQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`LowerCaseQ["string"]`**

Gives True if every character of "string" is a lowercase letter.

## Examples (12)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= LetterQ["αβγ"]
Out[1]= True

In[2]:= DigitQ["\:0662"]
Out[2]= False

In[3]:= UpperCaseQ["ABC"]
Out[3]= True

In[4]:= LowerCaseQ["abc1"]
Out[4]= False
```

### Applications (8)

```mathematica
In[5]:= LowerCaseQ["hello"]
Out[5]= True

In[6]:= LowerCaseQ["Hello"]
Out[6]= False
```

A digit is not a lowercase letter

```mathematica
In[7]:= LowerCaseQ["abc1"]
Out[7]= False
```

The sharp s is lowercase, though it has no single uppercase form

```mathematica
In[8]:= LowerCaseQ["straße"]
Out[8]= True
```

```mathematica
In[9]:= LowerCaseQ["σοφία"]
Out[9]= True

In[10]:= Select[{"abc", "ABC", "Abc", "x y"}, LowerCaseQ]
Out[10]= {"abc"}
```

Runs of lowercase letters

```mathematica
In[11]:= StringCases["iPhone eBay", _?LowerCaseQ ..]
Out[11]= {"i", "hone", "e", "ay"}
```

Not a string: LowerCaseQ::string

```mathematica
In[12]:= LowerCaseQ[{"a"}]
Out[12]= LowerCaseQ[{"a"}]
```

## Implementation notes

**Algorithm.** `builtin_lowercaseq` calls `charq_builtin(res, "LowerCaseQ", CQ_LOWER)`, which requires one `String` argument (otherwise `LowerCaseQ::string` via `strings_arg_message`, unevaluated), then decodes the UTF-8 code points one at a time with `u8_decode` and returns `False` at the first one for which `char_ok(CQ_LOWER, cp)` fails — so the scan short-circuits — and `True` if none does (vacuously for `""`). `char_ok` answers ASCII code points inline (`a`–`z`) and otherwise binary-searches `k_lower_ranges` with `in_ranges`. Inside a string pattern, `predicate_class` in `string_pattern.c` translates `_?LowerCaseQ` to the PCRE class `[[:lower:]]` instead of emitting a per-character callout.

**Data structures.** `k_lower_ranges`: sorted, disjoint code-point ranges of Unicode general category Ll, generated from Unicode 16.0.0 by `tools/gen_unicode_tables.py`.

**Complexity / limits.** `O(L)` for an ASCII string; otherwise `O(L log R)`. Any non-letter (digit, space, punctuation) makes the result `False`. Characters are code points, not grapheme clusters. Not `Listable`; no packed-array or `Compile[]` surface (the argument is a `String`).

**Attributes:** `Protected`.

## References

**See also:** [LetterQ](../../string-operations/LetterQ/), [DigitQ](../../string-operations/DigitQ/), [UpperCaseQ](../../string-operations/UpperCaseQ/)

- The Unicode Consortium, *The Unicode Standard, Version 16.0* (2024), §4.5 General Category, and the Unicode Character Database (UnicodeData.txt).
- Source: [`src/strings/stringcase.c`](https://github.com/stblake/mathilda/blob/main/src/strings/stringcase.c)
- Specification: [`docs/spec/builtins/string-operations.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/string-operations.md)
- Tests: [`tests/test_string_review.c`](https://github.com/stblake/mathilda/blob/main/tests/test_string_review.c)

## Notes & additional examples

### Notes

- `True` when every character is in Unicode category Ll (lowercase letter), using
  Unicode 16.0 ranges; `True` for `""`.
- A non-string emits `LowerCaseQ::string` and stays unevaluated; not `Listable`.
- See also: `UpperCaseQ`, `ToLowerCase`, `LetterQ`.
