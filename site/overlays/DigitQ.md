### Worked examples

```mathematica
In[1]:= DigitQ["0123456789"]
In[2]:= DigitQ["12a"]
In[3]:= DigitQ["-12"]  (* a sign is not a digit *)
In[4]:= DigitQ["1.5"]  (* nor is a decimal point *)
In[5]:= DigitQ["١٢٣"]  (* Arabic-Indic digits: DigitQ is ASCII-only *)
In[6]:= StringMatchQ["١٢٣", DigitCharacter..]  (* while the pattern DigitCharacter accepts them *)
In[7]:= DigitQ[""]  (* vacuously True *)
In[8]:= Select[{"2024", "v2", "007", ""}, DigitQ]
In[9]:= DigitQ[123]  (* a number is not a string: DigitQ::string *)
```

### Notes

- `True` exactly when every character is an ASCII digit `0`–`9`, as in Mathematica —
  deliberately narrower than the string pattern `DigitCharacter`, which matches any
  Unicode decimal digit.
- `True` for `""`. A non-string emits `DigitQ::string` and stays unevaluated; not
  `Listable`.
- `DigitQ` is a test on text, not a number test: use `IntegerQ` / `NumberQ` for
  expressions.
- See also: `DigitCharacter`, `LetterQ`, `ToExpression`.
