### Worked examples

```mathematica
In[1]:= StringPartition["abcdefgh", 3]  (* non-overlapping blocks; a short tail is dropped *)
```

```mathematica
In[1]:= StringPartition["abcdef", 2, 1]  (* length 2, offset 1: overlapping blocks *)
```

```mathematica
In[1]:= StringPartition["abcdef", 2, 3]  (* an offset larger than the block skips characters *)
```

```mathematica
In[1]:= StringPartition["abcdefg", UpTo[3]]  (* allow a shorter final block *)
```

```mathematica
In[1]:= StringPartition["abcdefgh", UpTo[3], 2]  (* with an offset, UpTo keeps the clipped tail *)
```

```mathematica
In[1]:= StringPartition["h\:00e9llo", 2]  (* blocks count characters, not bytes *)
```

```mathematica
In[1]:= StringPartition[{"ab", "cd"}, 1]  (* as in Mathematica, a list of strings is not threaded *)
```

### Notes

- `StringPartition[s, n, d]` gives the length-`n` substrings that start every `d`
  characters, with `d` defaulting to `n`. A tail too short for a whole block is
  dropped, unless the length is given as `UpTo[n]`.
- Since v0.456 lengths count Unicode code points, and a list first argument is
  rejected with `StringPartition::string` instead of being threaded, matching
  Mathematica. Map it over the list if needed.
- Other messages are `::intpm` (block length not a positive integer), `::intnm`
  (once an offset is given, `n = 0` is allowed) and `::intp` (bad offset). All
  leave the call unevaluated.
- See also: `Partition`, `StringTake`, `Characters`.
