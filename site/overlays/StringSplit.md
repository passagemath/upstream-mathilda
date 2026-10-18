### Worked examples

```mathematica
In[1]:= StringSplit["a b c"]  (* default: split at runs of whitespace *)
```

```mathematica
In[1]:= StringSplit["a,b,,c", ","]  (* interior empty pieces are kept *)
```

```mathematica
In[1]:= StringSplit[",a,b,", ",", All]  (* All keeps the leading and trailing empty pieces too *)
```

```mathematica
In[1]:= StringSplit["a b c d", " ", 2]  (* at most two pieces *)
```

```mathematica
In[1]:= StringSplit["a-b_c", "-" | "_"]  (* alternative delimiters *)
```

```mathematica
In[1]:= StringSplit["a1b22c", x:DigitCharacter.. :> FromDigits[x]]  (* a rule inserts a value bound from each delimiter *)
```

```mathematica
In[1]:= StringSplit["aXbxc", "x", IgnoreCase -> True]  (* case-insensitive delimiters *)
```

```mathematica
In[1]:= StringSplit["h\:00e9llo", ""]  (* the empty delimiter splits into Unicode characters *)
```

```mathematica
In[1]:= StringSplit[{"a b", "c d"}]  (* threads over a list of subjects *)
```

### Notes

- `StringSplit` returns the pieces between non-overlapping delimiter matches.
  Interior empty pieces are kept, and leading or trailing ones are dropped
  unless the third argument is `All`. An integer third argument `n` gives at
  most `n` pieces, with the remainder as the last one.
- With a rule `patt -> val` or `patt :> val`, the value is inserted at each
  delimiter. Pattern names are bound to the text they matched (v0.455, issue
  #98), so `x:DigitCharacter.. :> FromDigits[x]` inserts the numbers.
- The delimiter accepts the full string-pattern language, including
  `RegularExpression`. `IgnoreCase` is a registered option.
- `StringSplit[s, ""]` splits into Unicode code points, the same as `Characters`.
- See also: `StringRiffle` (the inverse), `StringCases`, `Characters`,
  `StringPartition`.
