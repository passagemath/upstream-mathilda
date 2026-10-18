### Worked examples

```mathematica
In[1]:= ToUpperCase["hello, world"]
In[2]:= ToUpperCase["déjà vu"]  (* accented Latin letters map through the Unicode tables *)
In[3]:= ToUpperCase["straße"]  (* full case mapping: one ß becomes the two letters SS *)
In[4]:= StringLength /@ {"ß", ToUpperCase["ß"]}  (* so the length can grow *)
In[5]:= ToUpperCase["ﬁ ŉ"]  (* ligatures and other multi-character expansions *)
In[6]:= ToUpperCase["αβγ ж"]  (* any cased script, here Greek and Cyrillic *)
In[7]:= ToUpperCase[{"abc", {"x", "y"}}]  (* Listable: threads through nested lists *)
In[8]:= ToUpperCase[{"abc", 42}]  (* a non-string element stays unevaluated, with no message *)
```

### Notes

- Full Unicode (16.0) uppercase mapping, including the multi-code-point expansions
  (`ß → SS`, `ﬁ → FI`, `ŉ → ʼN`); characters with no uppercase form — digits,
  punctuation, CJK — are copied unchanged.
- The mapping is context-free: there is no locale (Turkish dotless i) or
  language-specific tailoring.
- `ToUpperCase` is not invertible: `ToLowerCase[ToUpperCase["straße"]]` is `"strasse"`.
- Attributes `{Listable, Protected}`. A non-string argument is returned unevaluated
  without a message, as in Mathematica.
- No numeric surface (packed array / `Compile[]`): the result is a `String`.
- See also: `ToLowerCase`, `Capitalize`, `UpperCaseQ`.
