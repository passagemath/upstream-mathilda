### Worked examples

```mathematica
In[1]:= StringDelete["a1b22c", DigitCharacter..]
In[2]:= StringDelete["The quick brown fox", "o"]
In[3]:= StringDelete["Hello World", "l", IgnoreCase -> True]  (* deletes l and L *)
In[4]:= StringDelete["abcabc", {"a", "c"}]  (* a list of patterns *)
In[5]:= StringDelete["naïve café!", Except[LetterCharacter]]  (* keep only letters, Unicode-aware *)
In[6]:= StringDelete["foo123bar", RegularExpression["[0-9]+"]]
In[7]:= StringDelete["aXbYc", _?UpperCaseQ]
In[8]:= StringDelete[{"x1", "y22"}, DigitCharacter]  (* one result per string *)
In[9]:= StringDelete[Whitespace]["a b  c"]  (* operator form *)
```

### Notes

- `StringDelete[s, p]` is `StringReplace[s, p -> ""]`: each pattern is wrapped in a
  rule and the whole set goes through the shared PCRE2-backed `StringReplace`
  scanner (non-overlapping, left to right).
- Any string pattern is accepted — literals, `StringExpression`, `Alternatives`,
  character classes, `RegularExpression`, `PatternTest`, `Condition`.
- Option `IgnoreCase`. Unlike `StringReplace` there is no occurrence-limit argument;
  a third positional argument leaves the call unevaluated.
- A non-string subject emits `StringDelete::strse`. Attributes `{Protected}`.
- See also: `StringReplace`, `StringTrim`, `StringCases`.
