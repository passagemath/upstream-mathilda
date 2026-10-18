### Worked examples

```mathematica
In[1]:= StringTake["abcdef", 3]  (* the first three characters *)
```

```mathematica
In[1]:= StringTake["abcdef", -2]  (* the last two *)
```

```mathematica
In[1]:= StringTake["abcdef", {2, 4}]  (* characters 2 through 4 *)
```

```mathematica
In[1]:= StringTake["abcdefgh", {1, -1, 2}]  (* every second character *)
```

```mathematica
In[1]:= StringTake["abcdef", {{1, 2}, {4, 5}}]  (* a list of specs gives a list of substrings *)
```

```mathematica
In[1]:= StringTake[{"apple", "kiwi"}, 2]  (* a list of strings, the same spec for each *)
```

```mathematica
In[1]:= StringTake["abc", UpTo[10]]  (* at most 10 characters, without an error *)
```

```mathematica
In[1]:= StringTake["h\:00e9llo", 2]  (* positions count characters, not UTF-8 bytes *)
```

```mathematica
In[1]:= StringTake["abcdefgh", StringPosition["abcdefgh", "cde"][[1]]]  (* the inverse of StringPosition *)
```

```mathematica
In[1]:= StringTake["abc", 5]  (* out of range: StringTake::take, left unevaluated *)
```

### Notes

- Every Mathematica sequence spec is accepted: `n`, `-n`, `{n}`, `{m, n}`,
  `{m, n, s}`, `All`, `None` and `UpTo[n]`, or a list of them. The Take range
  rules are Mathematica's exactly, so `StringTake["abc", {2, 1}]` is `""`.
- Positions count Unicode code points (since v0.456). Only the prefix a spec
  needs is scanned, so a short take from a very long non-ASCII string is cheap.
- Errors route through `mth_message`, so `Quiet`/`Check` see them:
  `StringTake::take` (out of range), `::mseqs` (not a spec) and `::strse` (not
  a string). Inside a list result, a failing entry stays as an unevaluated
  `StringTake[...]`.
- A list of four or more integers is read as a list of specs, with a
  `StringTake::ambgsntx` warning, as in Mathematica. `Span` (`m ;; n`) is not
  accepted.
- See also: `StringDrop`, `StringTakeDrop`, `StringPart`, `StringPosition`, `Take`.
