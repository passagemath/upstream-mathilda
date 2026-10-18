### Worked examples

```mathematica
In[1]:= StringPadLeft["42", 5, "0"]  (* left-pad to width 5 with zeros *)
```

```mathematica
In[1]:= StringPadLeft["abcdef", 3]  (* too long: truncation keeps the last n characters *)
```

```mathematica
In[1]:= StringPadLeft[{"a", "bb", "ccc"}]  (* the one-argument list form pads to the longest element *)
```

```mathematica
In[1]:= StringPadLeft[{"1", "22", "333"}, 4, "."]  (* right-align a column *)
```

```mathematica
In[1]:= StringPadLeft["ab", 7, "xyz"]  (* a multi-character pad is aligned to the end of the result *)
```

```mathematica
In[1]:= StringPadLeft["caf\:00e9", 6, "*"]  (* lengths count characters, not UTF-8 bytes *)
```

```mathematica
In[1]:= StringPadLeft["abc", -1]  (* a negative length is rejected and stays unevaluated *)
```

### Notes

- `StringPadLeft[s, n, p]` makes `s` exactly `n` characters long. It pads on the
  left with copies of `p` (default `" "`), or truncates, keeping the last `n`
  characters.
- A multi-character pad is cycled one character at a time and aligned to the
  end of the result, as in Mathematica. `"ab"` padded to 7 with `"xyz"` is
  `"zxyzxab"`.
- On a list of strings every element is padded to `n`. Without `n`, every
  element is padded to the longest, so a ragged column becomes rectangular.
- Since v0.456 lengths count Unicode code points.
- Bad arguments emit Mathematica's messages and leave the call unevaluated:
  `::intnm` for a negative length, `::stringnz` for an empty pad, and
  `::strlist` for a bare string without `n`.
- Not `Listable`, because the one-argument form needs the whole list.
- See also: `StringPadRight`, `StringRepeat`, `IntegerString`.
