### Worked examples

```mathematica
In[1]:= FromCharacterCode[{72, 105}]
In[2]:= FromCharacterCode[955]  (* a single code gives a one-character string *)
In[3]:= FromCharacterCode[{{72, 105}, {955, 8364}}]  (* a matrix gives one string per row *)
In[4]:= FromCharacterCode[128512]  (* a code point outside the BMP *)
In[5]:= FromCharacterCode[Reverse[ToCharacterCode["stressed"]]]
In[6]:= FromCharacterCode[Range[97, 102]]
In[7]:= FromCharacterCode[NDArray[{77, 97, 116}, DataType -> "int64"]]  (* an int64 NDArray is read directly *)
In[8]:= FromCharacterCode[{65, 2000000}]  (* beyond U+10FFFF: FromCharacterCode::notunicode *)
In[9]:= FromCharacterCode[-1]  (* negative: FromCharacterCode::intnm *)
```

### Notes

- Codes are Unicode code points `0`–`1114111`, encoded to UTF-8. Surrogate code
  points `U+D800`–`U+DFFF` are rejected as well.
- A rank-1 or rank-2 `int64` buffer (packed `List` or visible `NDArray`) is read
  directly — `FromCharacterCode` is on the packed-aware and `INT64_OK` lists in
  `src/pack.c`. Other `NDArray` dtypes are unpacked first and then fail the integer
  check.
- An invalid code emits `::notunicode` (too large) or `::intnm` (negative or
  non-integer) and leaves the whole call unevaluated; both route through the
  message funnel, so `Check`/`Quiet` work.
- The encoding form `FromCharacterCode[codes, "enc"]` is not implemented.
- See also: `ToCharacterCode`, `CharacterRange`, `StringJoin`.
