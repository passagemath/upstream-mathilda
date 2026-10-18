### Worked examples

```mathematica
In[1]:= mod /: mod[a_, p_] + mod[b_, p_] := mod[Mod[a + b, p], p]; mod[2, 5] + mod[4, 5]   (* modular addition, attached to mod and not to the Protected Plus *)
```

```mathematica
In[1]:= mod /: i_Integer mod[a_, p_] := mod[Mod[i a, p], p]; 3 mod[4, 5]   (* scalar multiplication through an UpValue for Times; 12 mod 5 = 2 *)
```

```mathematica
In[1]:= mod[2, 5] + 3 mod[3, 5] + x   (* both rules work together, and the unrelated term x is kept *)
```

```mathematica
In[1]:= cnt = 0; tk /: hit[tk] := ++cnt; {hit[tk], hit[tk]}   (* delayed: the right side runs on every use *)
```

```mathematica
In[1]:= pt /: norm[pt[x_, y_]] := Sqrt[x^2 + y^2] /; NumericQ[x]; {norm[pt[3, 4]], norm[pt[a, b]]}   (* a trailing condition guards the rule *)
```

```mathematica
In[1]:= Definition[pt]   (* the stored rule, with the condition moved onto the left side *)
```

```mathematica
In[1]:= Check[Protect[pq]; pq /: q[pq] := 1, failed]   (* a Protected tag refuses the assignment: TagSetDelayed::write *)
```

### Notes

- `f /: lhs := rhs` stores `rhs` unevaluated on the tag `f`, as an own-, down- or up-value depending on where `f` occurs in `lhs`, and returns `Null`.
- For UpValues on `Orderless` heads (`Plus`, `Times`) the arguments are put in canonical order before matching, so `3 mod[4, 5]` and `mod[4, 5] 3` behave the same. On `Flat` heads a two-operand rule rewrites a matching pair and keeps the other operands.
- `lhs := body /; test` is stored as `lhs /; test`, which is how `Definition` shows it.
- Messages: `TagSetDelayed::sym`, `::write`, `::tagnf`, `::tagpos`, the same set as `TagSet`.
- Attributes `{HoldAll, Protected, SequenceHold}`.
- See also: `TagSet`, `TagUnset`, `UpSetDelayed`, `Definition`.
