### Worked examples

```mathematica
In[1]:= StringCount["mississippi", "ss"]  (* non-overlapping by default *)
```

```mathematica
In[1]:= StringCount["a1 b22 c333", DigitCharacter..]  (* maximal runs of digits *)
```

```mathematica
In[1]:= {StringCount["banana", "ana"], StringCount["banana", "ana", Overlaps -> True]}  (* overlapping occurrences are counted only on request *)
```

```mathematica
In[1]:= StringCount["abcd", __, Overlaps -> All]  (* every nonempty substring: 4 + 3 + 2 + 1 *)
```

```mathematica
In[1]:= StringCount["Hello HELLO", "hello", IgnoreCase -> True]  (* case-insensitive *)
```

```mathematica
In[1]:= StringCount[{"aa", "aba", "c"}, "a"]  (* one count per subject *)
```

```mathematica
In[1]:= StringCount["h\:00e9h\:00e9", "\:00e9"]  (* non-ASCII patterns *)
```

### Notes

- `StringCount` is the counting-only companion of `StringCases`. It runs the
  same match scan but never builds a substring.
- Because the scan is shared, `StringCount[s, p, opts]` is always
  `Length[StringCases[s, p, opts]]`. It also equals
  `Length[StringPosition[s, p, opts]]` under the same `Overlaps` setting, but
  note that `StringPosition` defaults to `Overlaps -> True`.
- A rule pattern is accepted, and only its left-hand side matters.
  `Overlaps`/`IgnoreCase` defaults come from `Options[StringCount]`.
- A positional argument count other than two emits `StringCount::argrx`.
- See also: `StringCases`, `StringPosition`, `Count`.
