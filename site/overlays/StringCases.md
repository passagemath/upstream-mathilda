### Worked examples

```mathematica
In[1]:= StringCases["the cat sat", LetterCharacter ..]  (* runs of letters, left to right *)
```

```mathematica
In[1]:= StringCases["abc123", (x__ ~~ "123") :> x]  (* a named piece of each match, bound into the rule *)
```

```mathematica
In[1]:= StringCases["a1b22c333", x:DigitCharacter.. :> ToExpression[x]]  (* the right-hand side is evaluated, so these are integers *)
```

```mathematica
In[1]:= StringCases["x=1, y=22", (k:LetterCharacter ~~ "=" ~~ v:DigitCharacter..) :> k -> ToExpression[v]]  (* parse key=value pairs *)
```

```mathematica
In[1]:= StringCases["a1b2c3", RegularExpression["[a-z](\\d)"] -> "$1"]  (* with RegularExpression a string RHS may use $n *)
```

```mathematica
In[1]:= StringCases["one two three four", WordCharacter.., 2]  (* an occurrence limit *)
```

```mathematica
In[1]:= StringCases["AAAA", "AA", Overlaps -> True]  (* overlapping matches, one per start *)
```

```mathematica
In[1]:= StringCases["abcd", __, Overlaps -> All]  (* every matching substring *)
```

```mathematica
In[1]:= {StringCases["aab abb", "a" ~~ __ ~~ "b"], StringCases["aab abb", Shortest["a" ~~ __ ~~ "b"]]}  (* greedy by default; Shortest is lazy *)
```

```mathematica
In[1]:= StringCases["h\:00e9llo w\:00f6rld", LetterCharacter..]  (* letter classes are Unicode-aware *)
```

```mathematica
In[1]:= StringCases[LetterCharacter..]["hi there"]  (* operator form *)
```

### Notes

- With a rule `patt :> rhs`, every pattern name in `patt` is bound to the text
  it matched, and `rhs` is evaluated. A name in an untaken alternative is bound
  to `""`. Until v0.455 the names were ignored and the whole match was returned
  (issue #98).
- `$0`/`$n` templates apply only to `RegularExpression` patterns. Elsewhere a
  string RHS is literal.
- `Overlaps -> False` (the default) is non-overlapping and greedy. `True` gives
  overlapping matches, the first one at each start. `All` gives every matching
  substring of every length.
- `StringCases[s, p, n]` keeps the first `n` matches, and the scan stops once it
  has them.
- `IgnoreCase` and `Overlaps` defaults come from `Options[StringCases]`, so
  `SetOptions` changes them.
- `StringCount` and `StringPosition` use the same scanner, so
  `StringCount[s, p, opts] == Length[StringCases[s, p, opts]]` always holds.
- A list of subjects threads, and non-string elements pass through unchanged.
- See also: `StringCount`, `StringPosition`, `StringReplace`, `Cases`.
