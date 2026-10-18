### Worked examples

```mathematica
In[1]:= StringPosition["abcabc", "bc"]  (* {start, end} pairs, 1-based inclusive *)
```

```mathematica
In[1]:= StringPosition["aaaa", "aa"]  (* overlapping matches are listed by default *)
```

```mathematica
In[1]:= StringPosition["aaaa", "aa", Overlaps -> False]  (* non-overlapping, as StringCases would give *)
```

```mathematica
In[1]:= StringPosition["a1b2c3", DigitCharacter]  (* a character-class pattern *)
```

```mathematica
In[1]:= StringPosition["abcabc", "bc", 1]  (* only the first match *)
```

```mathematica
In[1]:= StringPosition["AbAB", "ab", IgnoreCase -> True]  (* case-insensitive *)
```

```mathematica
In[1]:= StringPosition["h\:00e9llo", "llo"]  (* positions count characters, so the é is one position *)
```

```mathematica
In[1]:= StringPosition["xyz"]["axyzbxyz"]  (* operator form *)
```

```mathematica
In[1]:= StringTake["abcdefgh", StringPosition["abcdefgh", "cde"][[1]]]  (* positions feed StringTake directly *)
```

### Notes

- Results are `{start, end}` character positions, 1-based and inclusive: the
  form `StringTake`, `StringDrop` and `StringReplacePart` consume.
- As in the Wolfram Language, the default is `Overlaps -> True`, which lists
  overlapping matches at distinct starts. `StringCases` and `StringCount`
  default to `False`. `Overlaps -> All` lists every matching substring.
- Since v0.456 positions count Unicode code points rather than UTF-8 bytes, so
  they agree with `StringLength` on accented text.
- A third argument `n` keeps the first `n` matches. A one-argument call is the
  operator form. A list of subjects threads.
- See also: `StringCases`, `StringCount`, `StringTake`, `Position`.
