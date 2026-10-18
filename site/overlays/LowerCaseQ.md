### Worked examples

```mathematica
In[1]:= LowerCaseQ["hello"]
In[2]:= LowerCaseQ["Hello"]
In[3]:= LowerCaseQ["abc1"]  (* a digit is not a lowercase letter *)
In[4]:= LowerCaseQ["straße"]  (* the sharp s is lowercase, though it has no single uppercase form *)
In[5]:= LowerCaseQ["σοφία"]
In[6]:= Select[{"abc", "ABC", "Abc", "x y"}, LowerCaseQ]
In[7]:= StringCases["iPhone eBay", _?LowerCaseQ ..]  (* runs of lowercase letters *)
In[8]:= LowerCaseQ[{"a"}]  (* not a string: LowerCaseQ::string *)
```

### Notes

- `True` when every character is in Unicode category Ll (lowercase letter), using
  Unicode 16.0 ranges; `True` for `""`.
- A non-string emits `LowerCaseQ::string` and stays unevaluated; not `Listable`.
- See also: `UpperCaseQ`, `ToLowerCase`, `LetterQ`.
