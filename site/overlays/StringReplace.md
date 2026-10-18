### Worked examples

```mathematica
In[1]:= StringReplace["hello world", "o" -> "0"]  (* every match is replaced *)
```

```mathematica
In[1]:= StringReplace["the cat", {"the" -> "a", "cat" -> "dog"}]  (* several rules; the leftmost match wins at each position *)
```

```mathematica
In[1]:= StringReplace["abc123", a__ ~~ "123" :> a <> "!"]  (* a named piece of the match is bound into the right-hand side *)
```

```mathematica
In[1]:= StringReplace["2+3=5", d:DigitCharacter :> ToString[2 ToExpression[d]]]  (* the right-hand side is evaluated per match *)
```

```mathematica
In[1]:= StringReplace["abc123", x:DigitCharacter :> ToExpression[x]]  (* non-string values turn the result into a StringExpression *)
```

```mathematica
In[1]:= StringReplace["John Smith", RegularExpression["(\\w+) (\\w+)"] -> "$2 $1"]  (* with RegularExpression a string RHS is a $n template *)
```

```mathematica
In[1]:= StringReplace["abc", "b" -> "$0"]  (* outside RegularExpression, $0 is literal text *)
```

```mathematica
In[1]:= StringReplace["aaaa", "a" -> "b", 2]  (* an occurrence limit: only the first two matches *)
```

```mathematica
In[1]:= StringReplace["Hello hello", "hello" -> "bye", IgnoreCase -> True]  (* case-insensitive matching *)
```

```mathematica
In[1]:= StringReplace["x" -> "y"]["xox"]  (* operator form *)
```

```mathematica
In[1]:= StringReplace["one\ntwo", StartOfLine -> "> "]  (* a zero-width anchor inserts without consuming text *)
```

### Notes

- `StringReplace` scans left to right and replaces each non-overlapping match
  with its rule's right-hand side, copying the text in between. With several
  rules the leftmost match wins at each position, and ties go to rule order.
- Pattern names on the left (`x:p`, `x_`, `x__`) are bound to the text they
  matched and substituted into the right-hand side before it is evaluated, as in
  Mathematica. This was issue #98, fixed in v0.455; before the fix, names were
  ignored.
- When a replacement value is not a string, the result is a `StringExpression`
  of the pieces (`"abc"~~1~~2~~3`), exactly as in Mathematica. `Sequence[...]`
  values splice in.
- `$0`/`$n` templates work only when the pattern contains a `RegularExpression`.
- `StringReplace[s, rules, n]` replaces at most `n` matches. A one-argument call
  is the operator form. `IgnoreCase` is a registered option, so `SetOptions`
  changes its default.
- A subject that is not a string or a list of strings emits
  `StringReplace::strse` and stays unevaluated. A list of subjects threads.
- Compiled rule sets are cached (a 64-entry LRU), so mapping one replacement
  over many strings compiles the pattern once.
- See also: `StringReplaceList`, `StringDelete`, `StringCases`, `StringSplit`,
  `RegularExpression`.
