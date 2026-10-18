### Worked examples

```mathematica
In[1]:= CharacterRange["a", "e"]
In[2]:= CharacterRange["α", "ε"]  (* any Unicode characters, by code point *)
In[3]:= CharacterRange[65, 70]  (* integer endpoints are code points *)
In[4]:= CharacterRange[9312, 9316]  (* circled digits *)
In[5]:= StringJoin[CharacterRange["0", "9"]]
In[6]:= CharacterRange["e", "a"]  (* an empty range *)
In[7]:= StringCases["x9a3b", CharacterRange["a", "z"]]  (* as a string pattern it is a character class *)
In[8]:= CharacterRange["a", 70]  (* mixed endpoint kinds: CharacterRange::argtype *)
```

### Notes

- The endpoints are both one-character strings or both non-negative integer codes;
  the range is by code point, so `CharacterRange["ά", "ε"]` also picks up the
  accented Greek letters lying between them.
- `{}` when the first code exceeds the second. Anything else emits
  `CharacterRange::argtype` and stays unevaluated.
- Inside a string pattern (`StringCases`, `StringMatchQ`, …) it is translated to a
  PCRE class `[c1-c2]` instead of being expanded.
- Attributes `{Protected}`.
- See also: `Alphabet`, `FromCharacterCode`, `LetterCharacter`.
