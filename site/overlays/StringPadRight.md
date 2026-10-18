### Worked examples

```mathematica
In[1]:= StringPadRight["42", 5, "0"]  (* right-pad to width 5 with zeros *)
```

```mathematica
In[1]:= StringPadRight["abcdef", 3]  (* too long: truncation keeps the first n characters *)
```

```mathematica
In[1]:= StringPadRight["ab", 6, "xyz"]  (* the pad cycle is aligned to the result, so it continues past the text *)
```

```mathematica
In[1]:= StringPadRight["x", 6, "ab"]  (* result character i is the pad's character i mod 2 *)
```

```mathematica
In[1]:= StringPadRight[{"a", "bb", "ccc"}]  (* the one-argument list form pads to the longest element *)
```

```mathematica
In[1]:= StringRiffle[StringPadRight[{"id", "name", "x"}, 6, "."], "|"]  (* fixed-width fields *)
```

```mathematica
In[1]:= StringPadRight["abc", 5, ""]  (* an empty pad string is rejected *)
```

### Notes

- `StringPadRight[s, n, p]` mirrors `StringPadLeft`. It pads on the right, and
  when truncating it keeps the first `n` characters.
- The pad is aligned to the result: result character `i` (0-based) is
  `p[[i mod Length[p]]]`. So the cycle does not restart after the original
  text: `StringPadRight["x", 6, "ab"]` is `"xbabab"`. This matches Mathematica.
  Before v0.456 Mathilda gave `"xababa"`.
- Lengths count Unicode characters. Bad arguments emit `::intnm`, `::stringnz`
  or `::strlist` and leave the call unevaluated.
- See also: `StringPadLeft`, `StringRepeat`, `StringRiffle`.
