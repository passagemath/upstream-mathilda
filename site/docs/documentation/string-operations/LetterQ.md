# LetterQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`LetterQ["string"]`**

Gives True if every character of "string" is a letter (any script), and False otherwise. LetterQ\[""\] is True.

## Examples (14)

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

### Applications (10)

```mathematica
In[5]:= LetterQ["abc"]
Out[5]= True
```

One non-letter makes it False

```mathematica
In[6]:= LetterQ["abc1"]
Out[6]= False
```

Letters of any script count

```mathematica
In[7]:= LetterQ["αβγ日本語"]
Out[7]= True
```

Unicode category Lo: the feminine ordinal is a letter

```mathematica
In[8]:= LetterQ["ª"]
Out[8]= True
```

A Roman-numeral letterlike symbol is Nl, not a letter

```mathematica
In[9]:= LetterQ["Ⅻ"]
Out[9]= False
```

A combining accent is a separate, non-letter character

```mathematica
In[10]:= LetterQ["e\:0301"]
Out[10]= False
```

Vacuously True

```mathematica
In[11]:= LetterQ[""]
Out[11]= True
```

```mathematica
In[12]:= Select[{"cat", "r2d2", "naïve", "hello world"}, LetterQ]
Out[12]= {"cat", "naïve"}
```

Inside a string pattern it becomes a character class

```mathematica
In[13]:= StringCases["a1B2ç3", _?LetterQ]
Out[13]= {"a", "B", "ç"}
```

Not Listable: a list is not a string

```mathematica
In[14]:= LetterQ[{"a", "b"}]
Out[14]= LetterQ[{"a", "b"}]
```

## Implementation notes

**Algorithm.** `builtin_letterq` calls `charq_builtin(res, "LetterQ", CQ_LETTER)`, which requires one `String` argument (otherwise `LetterQ::string` via `strings_arg_message`, unevaluated), then decodes the UTF-8 code points one at a time with `u8_decode` and returns `False` at the first one for which `char_ok(CQ_LETTER, cp)` fails — so the scan short-circuits — and `True` if none does (vacuously for `""`). `char_ok` answers ASCII code points inline (`A`–`Z`, `a`–`z`) and otherwise binary-searches `k_letter_ranges` with `in_ranges`. Inside a string pattern, `predicate_class` in `string_pattern.c` translates `_?LetterQ` to the PCRE class `[[:alpha:]]` instead of emitting a per-character callout.

**Data structures.** `k_letter_ranges`: sorted, disjoint `[lo, hi]` code-point ranges covering Unicode general categories `L*` (Lu, Ll, Lt, Lm, Lo), generated from Unicode 16.0.0 by `tools/gen_unicode_tables.py` into `src/strings/unicode_tables.inc`.

**Complexity / limits.** `O(L)` for an ASCII string; otherwise `O(L log R)` over the `R` letter ranges. Combining marks (Mn/Mc) and letter-like numbers (Nl) are not letters. Characters are code points, not grapheme clusters. Not `Listable`; no packed-array or `Compile[]` surface (the argument is a `String`).

**Attributes:** `Protected`.

## References

**See also:** [DigitQ](../../string-operations/DigitQ/), [UpperCaseQ](../../string-operations/UpperCaseQ/), [LowerCaseQ](../../string-operations/LowerCaseQ/)

- The Unicode Consortium, *The Unicode Standard, Version 16.0* (2024), §4.5 General Category, and the Unicode Character Database (UnicodeData.txt).
- Source: [`src/strings/stringcase.c`](https://github.com/stblake/mathilda/blob/main/src/strings/stringcase.c)
- Specification: [`docs/spec/builtins/string-operations.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/string-operations.md)
- Tests: [`tests/test_string_review.c`](https://github.com/stblake/mathilda/blob/main/tests/test_string_review.c)

## Notes & additional examples

### Notes

- `True` when every character (code point) belongs to a Unicode letter category
  `L*` (Lu, Ll, Lt, Lm, Lo), using Unicode 16.0 ranges; ASCII is tested inline.
- Characters are code points, not grapheme clusters, so a decomposed `e` + U+0301
  is not all letters.
- A non-string argument emits `LetterQ::string` and stays unevaluated; it is not
  `Listable` (use `Map` or `AllTrue`).
- In a string pattern `_?LetterQ` is translated directly to a PCRE character class
  rather than a per-character callout.
- See also: `LetterCharacter`, `DigitQ`, `UpperCaseQ`, `LowerCaseQ`.
