### Worked examples

```mathematica
In[1]:= LetterQ["abc"]
In[2]:= LetterQ["abc1"]  (* one non-letter makes it False *)
In[3]:= LetterQ["αβγ日本語"]  (* letters of any script count *)
In[4]:= LetterQ["ª"]  (* Unicode category Lo: the feminine ordinal is a letter *)
In[5]:= LetterQ["Ⅻ"]  (* a Roman-numeral letterlike symbol is Nl, not a letter *)
In[6]:= LetterQ["e\:0301"]  (* a combining accent is a separate, non-letter character *)
In[7]:= LetterQ[""]  (* vacuously True *)
In[8]:= Select[{"cat", "r2d2", "naïve", "hello world"}, LetterQ]
In[9]:= StringCases["a1B2ç3", _?LetterQ]  (* inside a string pattern it becomes a character class *)
In[10]:= LetterQ[{"a", "b"}]  (* not Listable: a list is not a string *)
```

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
