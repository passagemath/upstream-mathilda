### Worked examples

```mathematica
In[1]:= ToCharacterCode["Hello"]
In[2]:= ToCharacterCode["héllo €😀"]  (* Unicode code points, including beyond the BMP *)
In[3]:= ToCharacterCode["e\:0301"]  (* a combining accent is its own code point *)
In[4]:= ToCharacterCode[{"ab", "xyz"}]  (* one list per string *)
In[5]:= ToCharacterCode["abc"] - 96  (* the codes are ordinary integers *)
In[6]:= PackedArrayQ[ToCharacterCode[StringRepeat["a", 100]]]  (* long strings give a packed buffer *)
In[7]:= FromCharacterCode[Mod[ToCharacterCode["hello"] - 97 + 13, 26] + 97]  (* ROT13 by arithmetic on codes *)
In[8]:= ToCharacterCode["abc", "UTF8"]  (* encodings are not supported: unevaluated *)
In[9]:= ToCharacterCode[5]  (* ToCharacterCode::strse *)
```

### Notes

- Returns Unicode code points (decoded from the internal UTF-8), not bytes:
  `ToCharacterCode["é"]` is `{233}`.
- Above the packing threshold the codes are written straight into a packed `int64`
  buffer; short strings give a plain `List`. `FromCharacterCode` reads that buffer
  without boxing, so `FromCharacterCode[ToCharacterCode[s]]` stays on the buffer.
- The two-argument encoding form `ToCharacterCode[s, "enc"]` is not implemented.
- Accepts a string or a flat list of strings; otherwise `::strse`. Not `Listable`.
- No `Compile[]` lowering: the input is a `String`, which no compiled register holds.
- See also: `FromCharacterCode`, `Characters`, `StringLength`.
