### Worked examples

```mathematica
In[1]:= StringTakeDrop["abcdef", 2]
In[2]:= StringTakeDrop["abcdef", -2]  (* the last two characters *)
In[3]:= StringTakeDrop["abcdef", {2, 4}]
In[4]:= StringTakeDrop["abcdefgh", {1, -1, 2}]  (* a stepped span splits odd and even positions *)
In[5]:= StringTakeDrop["abcdef", UpTo[10]]  (* clamping with UpTo instead of failing *)
In[6]:= StringTakeDrop["héllo", 2]  (* positions are characters *)
In[7]:= StringTakeDrop[{"abcd", "wxyz"}, 2]  (* one pair per string *)
In[8]:= StringTakeDrop["abc", 5]  (* out of range: StringTakeDrop::take *)
```

### Notes

- Equivalent to `{StringTake[s, spec], StringDrop[s, spec]}` and built on the same
  scalar kernels (`strings_take_seq` / `strings_drop_seq`), so the two halves always
  partition the characters of `s`.
- Accepts a single sequence spec: `n`, `-n`, `{n}`, `{m, n}`, `{m, n, s}`, `UpTo[n]`,
  `All`, `None`. Positions are validated against every string first, so an
  out-of-range spec reports `StringTakeDrop::take` (not `StringTake::take`).
- Messages: `::strse` (subject), `::seqs` (spec), `::take` (range), `::argr`
  (argument count). Attributes `{Protected}`.
- See also: `StringTake`, `StringDrop`, `TakeDrop`.
