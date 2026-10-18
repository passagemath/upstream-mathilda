### Worked examples

```mathematica
In[1]:= Capitalize["hello world"]  (* only the first character changes *)
In[2]:= Capitalize["élan vital"]
In[3]:= Capitalize["ßx"]  (* the first character uses the full mapping, so it can expand *)
In[4]:= Capitalize[{"alpha", "beta", "gamma"}]  (* a list of strings maps *)
In[5]:= Capitalize[""]
In[6]:= StringRiffle[Capitalize /@ StringSplit["the quick brown fox"]]  (* title case, word by word *)
In[7]:= Capitalize[5]  (* not a string: Capitalize::strse, unevaluated *)
```

### Notes

- Maps only the first character, with the same full Unicode uppercase table as
  `ToUpperCase`; the rest of the string is copied byte for byte.
- It uses the *uppercase*, not the Unicode titlecase, mapping.
- Accepts a string or a flat list of strings; it is not `Listable`, so a mixed
  or nested list emits `Capitalize::strse` and stays unevaluated.
- Attributes `{Protected}`.
- See also: `Decapitalize`, `ToUpperCase`.
