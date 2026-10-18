### Worked examples

```mathematica
In[1]:= UpperCaseQ["ABC"]
In[2]:= UpperCaseQ["AbC"]
In[3]:= UpperCaseQ["A B"]  (* a space is not an uppercase letter *)
In[4]:= UpperCaseQ["ΔΣÉ"]  (* any script with case *)
In[5]:= UpperCaseQ["ǅ"]  (* titlecase digraphs are category Lt, not Lu *)
In[6]:= Length[Select[CharacterRange[192, 255], UpperCaseQ]]  (* the 30 Latin-1 capitals, U+00C0-U+00DE minus × *)
In[7]:= StringCases["Hello World", _?UpperCaseQ ~~ (_?LowerCaseQ)..]  (* capitalised words *)
In[8]:= UpperCaseQ[""]
```

### Notes

- `True` when every character is in Unicode category Lu (uppercase letter); digits,
  spaces and punctuation make it `False`, so it is not "has no lowercase letters".
- Titlecase letters (Lt, e.g. `ǅ`) are neither upper- nor lowercase.
- `True` for `""`; a non-string emits `UpperCaseQ::string`. Not `Listable`.
- See also: `LowerCaseQ`, `ToUpperCase`, `LetterQ`.
