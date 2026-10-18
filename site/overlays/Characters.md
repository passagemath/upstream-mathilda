### Worked examples

```mathematica
In[1]:= Characters["cat"]  (* a string becomes a list of length-1 strings *)
```

```mathematica
In[1]:= StringJoin[Characters["abcd"]]  (* the inverse is StringJoin *)
```

```mathematica
In[1]:= Characters["h\:00e9llo"]  (* one element per Unicode character, so é stays whole *)
```

```mathematica
In[1]:= Characters["\[Alpha]\[Beta]"]  (* named characters decode in string literals *)
```

```mathematica
In[1]:= StringLength["h\:00e9llo"] == Length[Characters["h\:00e9llo"]]  (* both count code points *)
```

```mathematica
In[1]:= Characters[{"ab", "c"}]  (* listable over a list of strings *)
```

```mathematica
In[1]:= Tally[Characters["mississippi"]]  (* letter frequencies *)
```

```mathematica
In[1]:= Characters[123]  (* a non-string stays unevaluated *)
```

### Notes

- `Characters[s]` gives one length-1 string per Unicode code point. Since v0.456
  a multi-byte UTF-8 character stays whole, so
  `Length[Characters[s]] == StringLength[s]` for every string. Earlier versions
  split by byte.
- Each ASCII character is a single shared, refcounted string, so splitting ASCII
  text allocates only the list (200k characters: 23.7 ms → about 2 ms in v0.456).
- `Listable | Protected`, and `Characters[""]` is `{}`.
- See also: `StringJoin`, `StringSplit` (with `""`), `ToCharacterCode`,
  `StringLength`.
