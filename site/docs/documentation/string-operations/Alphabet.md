# Alphabet

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Alphabet[]`**

Gives the list of the lowercase letters a through z.

## Examples (9)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= CharacterRange["a", "e"]
Out[1]= {"a", "b", "c", "d", "e"}

In[2]:= CharacterRange[65, 70]
Out[2]= {"A", "B", "C", "D", "E", "F"}

In[3]:= Length[Alphabet[]]
Out[3]= 26
```

### Applications (6)

```mathematica
In[4]:= Alphabet[]
Out[4]= {"a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o", "p", "q", "r", "s", "t", "u", "v", "w", "x", "y", "z"}

In[5]:= Length[Alphabet[]]
Out[5]= 26

In[6]:= StringJoin[ToUpperCase[Alphabet[]]]
Out[6]= "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
```

The alphabet index of a letter

```mathematica
In[7]:= Position[Alphabet[], "m"]
Out[7]= {{13}}
```

```mathematica
In[8]:= Alphabet[] === CharacterRange["a", "z"]
Out[8]= True
```

Language arguments are not implemented: unevaluated

```mathematica
In[9]:= Alphabet["Greek"]
Out[9]= Alphabet["Greek"]
```

## Implementation notes

**Algorithm.** `builtin_alphabet` accepts only the zero-argument call and returns `char_list('a', 'z')` — the same helper as `CharacterRange`, which encodes each code point `0x61`–`0x7A` as a one-character `String`. Any argument (a language name, a type spec) makes it return `NULL`, so the call stays unevaluated.

**Data structures.** A fresh 26-element `List` of one-character strings on every call.

**Complexity / limits.** Constant time. Mathematica's language-specific forms (`Alphabet["Greek"]`, `Alphabet[lang, type]`) and the language-data they rely on are not implemented.

**Attributes:** `Protected`.

## References

**See also:** [CharacterRange](../../string-operations/CharacterRange/)

- Source: [`src/strings/charcode.c`](https://github.com/stblake/mathilda/blob/main/src/strings/charcode.c)
- Specification: [`docs/spec/builtins/string-operations.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/string-operations.md)
- Tests: [`tests/test_string_review.c`](https://github.com/stblake/mathilda/blob/main/tests/test_string_review.c)

## Notes & additional examples

### Notes

- Only the zero-argument form is implemented: the 26 lowercase English letters,
  built by the same code path as `CharacterRange["a", "z"]`.
- Mathematica's `Alphabet["language"]` and `Alphabet[spec, type]` forms are not
  supported and stay unevaluated.
- Attributes `{Protected}`.
- See also: `CharacterRange`, `LetterQ`.
