### Worked examples

```mathematica
In[1]:= StringRotateLeft["abcdef", 2]
In[2]:= StringRotateLeft["abcdef"]  (* one place by default *)
In[3]:= StringRotateLeft["abcdef", -1]  (* a negative shift rotates right *)
In[4]:= StringRotateLeft["abcdef", 13]  (* shifts wrap modulo the length *)
In[5]:= StringRotateLeft["αβγδ", 1]  (* positions count characters, not bytes *)
In[6]:= Table[StringRotateLeft["abcd", k], {k, 0, 4}]
In[7]:= StringRotateLeft["", 3]
In[8]:= StringRotateLeft["abc", x]  (* a symbolic shift: StringRotateLeft::int *)
```

### Notes

- The shift is reduced modulo the character count, the cut point located with
  `u8_offset`, and the result assembled with two `memcpy`s of the UTF-8 bytes.
- Takes a single string; a list of strings is not threaded (emits `::string`),
  so use `Map`. A non-integer shift emits `::int`.
- Attributes `{Protected}`.
- See also: `StringRotateRight`, `RotateLeft`, `StringTake`.
