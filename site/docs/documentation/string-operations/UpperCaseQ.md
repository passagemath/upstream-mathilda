# UpperCaseQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`UpperCaseQ["string"]`**

Gives True if every character of "string" is an uppercase letter.

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
In[5]:= UpperCaseQ["ABC"]
Out[5]= True

In[6]:= UpperCaseQ["AbC"]
Out[6]= False
```

A space is not an uppercase letter

```mathematica
In[7]:= UpperCaseQ["A B"]
Out[7]= False
```

Any script with case

```mathematica
In[8]:= UpperCaseQ["ΔΣÉ"]
Out[8]= True
```

Titlecase digraphs are category Lt, not Lu

```mathematica
In[9]:= UpperCaseQ["ǅ"]
Out[9]= False
```

The 30 Latin-1 capitals, U+00C0-U+00DE minus ×

```mathematica
In[10]:= Length[Select[CharacterRange[192, 255], UpperCaseQ]]
Out[10]= 30
```

Capitalised words

```mathematica
In[11]:= StringCases["Hello World", _?UpperCaseQ ~~ (_?LowerCaseQ)..]
Out[11]= {"Hello", "World"}
```

```mathematica
In[12]:= UpperCaseQ[""]
Out[12]= True
```

## Implementation notes

**Algorithm.** `builtin_uppercaseq` calls `charq_builtin(res, "UpperCaseQ", CQ_UPPER)`, which requires one `String` argument (otherwise `UpperCaseQ::string` via `strings_arg_message`, unevaluated), then decodes the UTF-8 code points one at a time with `u8_decode` and returns `False` at the first one for which `char_ok(CQ_UPPER, cp)` fails — so the scan short-circuits — and `True` if none does (vacuously for `""`). `char_ok` answers ASCII code points inline (`A`–`Z`) and otherwise binary-searches `k_upper_ranges` with `in_ranges`. Inside a string pattern, `predicate_class` in `string_pattern.c` translates `_?UpperCaseQ` to the PCRE class `[[:upper:]]` instead of emitting a per-character callout.

**Data structures.** `k_upper_ranges`: sorted, disjoint code-point ranges of Unicode general category Lu, generated from Unicode 16.0.0 by `tools/gen_unicode_tables.py`.

**Complexity / limits.** `O(L)` for an ASCII string; otherwise `O(L log R)`. Titlecase letters (Lt, such as `ǅ`) are not uppercase; any non-letter makes the result `False`. Characters are code points, not grapheme clusters. Not `Listable`; no packed-array or `Compile[]` surface (the argument is a `String`).

**Attributes:** `Protected`.

## References

**See also:** [LetterQ](../../string-operations/LetterQ/), [DigitQ](../../string-operations/DigitQ/), [LowerCaseQ](../../string-operations/LowerCaseQ/)

- The Unicode Consortium, *The Unicode Standard, Version 16.0* (2024), §4.5 General Category, and the Unicode Character Database (UnicodeData.txt).
- Source: [`src/strings/stringcase.c`](https://github.com/stblake/mathilda/blob/main/src/strings/stringcase.c)
- Specification: [`docs/spec/builtins/string-operations.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/string-operations.md)
- Tests: [`tests/test_string_review.c`](https://github.com/stblake/mathilda/blob/main/tests/test_string_review.c)

## Notes & additional examples

### Notes

- `True` when every character is in Unicode category Lu (uppercase letter); digits,
  spaces and punctuation make it `False`, so it is not "has no lowercase letters".
- Titlecase letters (Lt, e.g. `ǅ`) are neither upper- nor lowercase.
- `True` for `""`; a non-string emits `UpperCaseQ::string`. Not `Listable`.
- See also: `LowerCaseQ`, `ToUpperCase`, `LetterQ`.
