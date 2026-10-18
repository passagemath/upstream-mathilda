# DigitQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`DigitQ["string"]`**

Gives True if every character of "string" is a digit 0-9, and False otherwise. DigitQ\[""\] is True.

## Examples (13)

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

### Applications (9)

```mathematica
In[5]:= DigitQ["0123456789"]
Out[5]= True

In[6]:= DigitQ["12a"]
Out[6]= False
```

A sign is not a digit

```mathematica
In[7]:= DigitQ["-12"]
Out[7]= False
```

Nor is a decimal point

```mathematica
In[8]:= DigitQ["1.5"]
Out[8]= False
```

Arabic-Indic digits: DigitQ is ASCII-only

```mathematica
In[9]:= DigitQ["١٢٣"]
Out[9]= False
```

While the pattern DigitCharacter accepts them

```mathematica
In[10]:= StringMatchQ["١٢٣", DigitCharacter..]
Out[10]= True
```

Vacuously True

```mathematica
In[11]:= DigitQ[""]
Out[11]= True
```

```mathematica
In[12]:= Select[{"2024", "v2", "007", ""}, DigitQ]
Out[12]= {"2024", "007", ""}
```

A number is not a string: DigitQ::string

```mathematica
In[13]:= DigitQ[123]
Out[13]= DigitQ[123]
```

## Implementation notes

**Algorithm.** `builtin_digitq` calls `charq_builtin(res, "DigitQ", CQ_DIGIT)`, which requires one `String` argument (otherwise `DigitQ::string` via `strings_arg_message`, unevaluated), then decodes the UTF-8 code points one at a time with `u8_decode` and returns `False` at the first one for which `char_ok(CQ_DIGIT, cp)` fails — so the scan short-circuits — and `True` if none does (vacuously for `""`). `char_ok` answers ASCII code points inline (`0`–`9`) and returns 0 for every non-ASCII code point without consulting any table — DigitQ is ASCII-only, as in Mathematica. Inside a string pattern, `predicate_class` in `string_pattern.c` translates `_?DigitQ` to the PCRE class `[[:digit:]]` instead of emitting a per-character callout.

**Data structures.** None beyond the decoder: the test is a pure range check on the code point.

**Complexity / limits.** `O(L)` for an ASCII string; any non-ASCII character fails immediately. Note the asymmetry with the string pattern `DigitCharacter`, which (compiled with PCRE2 UCP) matches every Unicode decimal digit. Characters are code points, not grapheme clusters. Not `Listable`; no packed-array or `Compile[]` surface (the argument is a `String`).

**Attributes:** `Protected`.

## References

**See also:** [LetterQ](../../string-operations/LetterQ/), [UpperCaseQ](../../string-operations/UpperCaseQ/), [LowerCaseQ](../../string-operations/LowerCaseQ/)

- The Unicode Consortium, *The Unicode Standard, Version 16.0* (2024), §4.5 General Category, and the Unicode Character Database (UnicodeData.txt).
- Source: [`src/strings/stringcase.c`](https://github.com/stblake/mathilda/blob/main/src/strings/stringcase.c)
- Specification: [`docs/spec/builtins/string-operations.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/string-operations.md)
- Tests: [`tests/test_string_review.c`](https://github.com/stblake/mathilda/blob/main/tests/test_string_review.c)

## Notes & additional examples

### Notes

- `True` exactly when every character is an ASCII digit `0`–`9`, as in Mathematica —
  deliberately narrower than the string pattern `DigitCharacter`, which matches any
  Unicode decimal digit.
- `True` for `""`. A non-string emits `DigitQ::string` and stays unevaluated; not
  `Listable`.
- `DigitQ` is a test on text, not a number test: use `IntegerQ` / `NumberQ` for
  expressions.
- See also: `DigitCharacter`, `LetterQ`, `ToExpression`.
