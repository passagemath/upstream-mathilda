### Worked examples

```mathematica
In[1]:= StringRepeat["ab", 4]  (* four copies *)
```

```mathematica
In[1]:= StringRepeat["-", 10]  (* a rule of ten dashes *)
```

```mathematica
In[1]:= StringRepeat["abc", 5, 7]  (* capped at 7 characters, so the last copy is partial *)
```

```mathematica
In[1]:= StringRepeat["xy", 3, 100]  (* a cap larger than the result has no effect *)
```

```mathematica
In[1]:= StringRepeat["a\:00e9", 2, 3]  (* the cap counts characters, so no character is cut in half *)
```

```mathematica
In[1]:= StringRepeat["abc", 0]  (* zero copies is the empty string *)
```

```mathematica
In[1]:= StringRepeat["ab", -1]  (* a negative count: StringRepeat::intp, left unevaluated *)
```

### Notes

- `StringRepeat[s, n]` concatenates `n` copies of `s`. `StringRepeat[s, n, max]`
  stops at `max` characters, keeping a partial final copy.
- Since v0.456 lengths count Unicode code points.
- A non-string first argument gives `StringRepeat::string`, and a negative or
  non-integer count or cap gives `StringRepeat::intp`. Both leave the call
  unevaluated.
- Whole copies are `memcpy`'d and the buffer is adopted into the result
  without a second copy.
- See also: `StringJoin`, `StringPadRight`, `Table`.
