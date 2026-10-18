### Worked examples

```mathematica
In[1]:= StringQ["AbC"]  (* a string *)
```

```mathematica
In[1]:= StringQ[""]  (* the empty string still counts *)
```

```mathematica
In[1]:= StringQ[123]  (* a number is not a string *)
```

```mathematica
In[1]:= StringQ[{"a", "b"}]  (* not Listable: a list is not a string *)
```

```mathematica
In[1]:= StringQ /@ {"a", 1, x, "b"}  (* map it for an element-wise test *)
```

```mathematica
In[1]:= Select[{1, "a", x, "bc"}, StringQ]  (* as a predicate *)
```

```mathematica
In[1]:= StringQ[ToString[42]]  (* ToString always returns a String *)
```

### Notes

- `StringQ[expr]` is `True` exactly when `expr` is a string, including `""`, and
  `False` for everything else. The test is a single type-tag check.
- It is deliberately **not** `Listable`, so a list is tested as a single object.
  Map it for an element-wise answer.
- Any arity other than one emits `StringQ::argx` through the message funnel and
  leaves the call unevaluated. `Protected`.
- See also: `StringMatchQ`, `NumberQ`, `MatchQ` with `_String`.
