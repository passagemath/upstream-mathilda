### Worked examples

```mathematica
In[1]:= StringReplaceList["aaa", "a" -> "b"]  (* each single replacement separately *)
In[2]:= StringReplaceList["aaa", "aa" -> "X"]  (* overlapping matches each give a result *)
In[3]:= StringReplaceList["abc", {"a" -> "1", "b" -> "2"}]  (* several rules, ordered by position *)
In[4]:= StringReplaceList["aaaa", "a" -> "b", 2]  (* at most two results *)
In[5]:= StringReplaceList["a1b2", d:DigitCharacter :> "<" <> d <> ">"]  (* a named pattern used on the right-hand side *)
In[6]:= StringReplaceList["aAa", "a" -> "-", IgnoreCase -> True]
In[7]:= StringReplaceList["αβα", "α" -> "a"]
In[8]:= StringReplaceList[{"ab", "ba"}, "a" -> "X"]  (* one list per string *)
In[9]:= StringReplaceList["abc", "z" -> "y"]  (* no match gives an empty list *)
In[10]:= StringReplaceList["ab", "b" -> f[1]]  (* a non-string value gives a StringExpression *)
```

### Notes

- Enumerates every match start of every rule with overlapping allowed
  (`regex_scan_limit` in `Overlaps -> True` mode), stably sorted by position and then
  rule order, and builds one string per match with only that match replaced.
- The limit `n` (integer ≥ 0 or `Infinity`) truncates after sorting, so it costs a
  full scan; `n = 0` returns `{}` immediately.
- The second argument must be rules: a bare pattern leaves the call unevaluated.
- Option `IgnoreCase`; a non-string subject emits `::strse`. Attributes `{Protected}`.
- See also: `StringReplace`, `ReplaceList`, `StringCases`.
