### Worked examples

```mathematica
In[1]:= Alphabet[]
In[2]:= Length[Alphabet[]]
In[3]:= StringJoin[ToUpperCase[Alphabet[]]]
In[4]:= Position[Alphabet[], "m"]  (* the alphabet index of a letter *)
In[5]:= Alphabet[] === CharacterRange["a", "z"]
In[6]:= Alphabet["Greek"]  (* language arguments are not implemented: unevaluated *)
```

### Notes

- Only the zero-argument form is implemented: the 26 lowercase English letters,
  built by the same code path as `CharacterRange["a", "z"]`.
- Mathematica's `Alphabet["language"]` and `Alphabet[spec, type]` forms are not
  supported and stay unevaluated.
- Attributes `{Protected}`.
- See also: `CharacterRange`, `LetterQ`.
