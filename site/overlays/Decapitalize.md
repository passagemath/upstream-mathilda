### Worked examples

```mathematica
In[1]:= Decapitalize["Hello World"]
In[2]:= Decapitalize["ABC"]  (* only the first character is lowered *)
In[3]:= Decapitalize["ÉCOLE"]
In[4]:= Decapitalize["Ωmega"]
In[5]:= Decapitalize[{"Alpha", "Beta"}]  (* a list of strings maps *)
In[6]:= Decapitalize["123abc"]  (* a non-letter first character is unchanged *)
In[7]:= Decapitalize[{"Alpha", 1}]  (* a list with a non-string: Decapitalize::strse *)
```

### Notes

- Lowercases the first character with the full Unicode table used by `ToLowerCase`;
  everything after it is copied unchanged.
- Accepts a string or a flat list of strings (not `Listable`); anything else
  emits `Decapitalize::strse` and stays unevaluated.
- Attributes `{Protected}`.
- See also: `Capitalize`, `ToLowerCase`.
