### Worked examples

```mathematica
In[1]:= StringContainsQ["abcdef", "cd"]  (* some substring matches *)
```

```mathematica
In[1]:= StringContainsQ["abc123", DigitCharacter]  (* any string pattern *)
```

```mathematica
In[1]:= StringContainsQ["abc", "x" | "c"]  (* alternatives *)
```

```mathematica
In[1]:= StringContainsQ[{"apple", "banana", "cherry"}, "an"]  (* threads over a list of subjects *)
```

```mathematica
In[1]:= StringContainsQ["Mathilda", "math", IgnoreCase -> True]  (* case folded *)
```

```mathematica
In[1]:= StringContainsQ["caf\:00e9", "\:00c9", IgnoreCase -> True]  (* case folding is Unicode-aware *)
```

```mathematica
In[1]:= StringContainsQ["abc", RegularExpression["^b"]]  (* a regex anchor still means the start of the string *)
```

```mathematica
In[1]:= Select[{"cat", "dog", "cow"}, StringContainsQ["o"]]  (* operator form as a predicate *)
```

### Notes

- `StringContainsQ` searches inside the string, where `StringMatchQ` must match
  the whole of it. It is equivalent to `StringMatchQ[s, ___ ~~ p ~~ ___]` and to
  `!StringFreeQ[s, p]`.
- The one-argument operator form `StringContainsQ[p]` (with any
  `IgnoreCase -> ...` it was given) curries into a predicate, which suits
  `Select`.
- It stops at the first match, with one PCRE2 probe per rule. It shares its core
  with `StringFreeQ`, `StringStartsQ` and `StringEndsQ`.
- A non-string subject, or a list holding one, leaves the whole call unevaluated
  rather than threading a wrong Boolean. An empty pattern list matches nothing.
- See also: `StringMatchQ`, `StringFreeQ`, `StringStartsQ`, `StringEndsQ`.
