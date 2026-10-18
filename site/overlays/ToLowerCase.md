### Worked examples

```mathematica
In[1]:= ToLowerCase["Hello, World"]
In[2]:= ToLowerCase["ÀÉÎÕÜ"]  (* Latin-1 capitals *)
In[3]:= ToLowerCase["ΑΘΗΝΑ"]  (* Greek capitals *)
In[4]:= ToLowerCase["İ"]  (* full mapping: capital dotted I becomes i plus a combining dot *)
In[5]:= ToCharacterCode[ToLowerCase["İ"]]  (* the two code points of that result *)
In[6]:= ToLowerCase[ToUpperCase["straße"]]  (* case mapping is not invertible *)
In[7]:= ToLowerCase[{"ABC", "XyZ"}]  (* Listable *)
In[8]:= ToLowerCase[x]  (* a non-string stays unevaluated *)
```

### Notes

- Full Unicode (16.0) lowercase mapping from the generated `k_lower_map` table; a
  character may map to several code points (`İ → i̇`, length 2).
- No context-sensitive rules: a word-final capital sigma `Σ` always becomes `σ`,
  never final `ς`, and there is no locale tailoring.
- Attributes `{Listable, Protected}`; a non-string is left unevaluated silently.
- See also: `ToUpperCase`, `Decapitalize`, `LowerCaseQ`.
