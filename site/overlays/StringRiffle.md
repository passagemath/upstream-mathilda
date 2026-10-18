### Worked examples

```mathematica
In[1]:= StringRiffle[{"a", "b", "c"}]  (* default scheme: a single space *)
```

```mathematica
In[1]:= StringRiffle[{"2024", "01", "02"}, "-"]  (* a custom separator *)
```

```mathematica
In[1]:= StringRiffle[{"x", "y", "z"}, {"(", ", ", ")"}]  (* a {left, sep, right} triple *)
```

```mathematica
In[1]:= StringRiffle[{1, 2, 3}, ", "]  (* non-string elements are converted as ToString would *)
```

```mathematica
In[1]:= StringRiffle[{{"a", "b"}, {"c", "d"}}]  (* nested: rows by newline, cells by space *)
```

```mathematica
In[1]:= StringRiffle[{{1, 2}, {3, 4}}, "; ", ","]  (* one separator per level, outermost first *)
```

```mathematica
In[1]:= StringRiffle[StringSplit["a  b   c"], "_"]  (* the inverse of StringSplit *)
```

```mathematica
In[1]:= StringRiffle[5]  (* the first argument must be a list: StringRiffle::list *)
```

### Notes

- `StringRiffle[list, sep1, sep2, ...]` joins a possibly nested list. `sep_i`
  goes between elements at level `i`, counting from the outside. A separator is
  a string or a `{left, sep, right}` triple.
- Levels without an explicit separator use the default scheme: one space at the
  innermost level and one more newline per level above it.
- String leaves are used verbatim, and any other leaf appears in its OutputForm.
- The result is built in one buffer, so each leaf is copied once whatever the
  depth (v0.456: 3.6 ms → 1.0 ms on the string benchmark).
- See also: `StringJoin`, `StringSplit`, `Riffle`.
