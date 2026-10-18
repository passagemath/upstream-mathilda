### Worked examples

```mathematica
In[1]:= x /: f2[x] = 1; {f2[x], f2[y]}   (* the rule is associated with the tag x *)
```

```mathematica
In[1]:= UpValues[x]   (* x is an argument of the left side, so the rule is an UpValue *)
```

```mathematica
In[1]:= tg /: tg[0] = zero; {DownValues[tg], UpValues[tg]}   (* when the tag is the head of the left side, the rule is an ordinary DownValue *)
```

```mathematica
In[1]:= fdv[x_] := down; gg /: fdv[gg] = up; {fdv[gg], fdv[hh]}   (* the UpValue on gg is tried before the DownValue on fdv *)
```

```mathematica
In[1]:= vv /: Length[vv] = 42; {Length[vv], Length[{1, 2}]}   (* a rule for the Protected head Length, stored on vv *)
```

```mathematica
In[1]:= xx = 5; yy /: gg2[xx, yy] = 9; UpValues[yy]   (* other arguments of the left side are evaluated before the rule is stored *)
```

```mathematica
In[1]:= Check[hh /: k[z] = 3, failed]   (* the tag must occur in the left side: TagSet::tagnf *)
```

```mathematica
In[1]:= Check[uu /: uu[a][b] = 1, failed]   (* a tag in head-of-head position would make a SubValue, which is not supported: TagSet::tagpos *)
```

### Notes

- `f /: lhs = rhs` evaluates `rhs` once, stores the rule on the **single** tag `f`, and returns `rhs`. Where it goes depends on where `f` occurs: as `lhs` itself it is an OwnValue, as the head of `lhs` a DownValue, as an argument or the head of an argument an UpValue.
- `UpSet` attaches a rule to every level-one symbol; `TagSet` attaches it to the named tag only. With `a /: fz[a, b] = 1`, `b` gets nothing.
- Messages: `TagSet::sym` (the tag is not a symbol), `TagSet::write` (the tag is Protected), `TagSet::tagnf` (the tag is not in `lhs`), `TagSet::tagpos` (SubValues are not supported). In each case nothing is stored and `rhs` is still returned.
- Only the tag's Protected bit is checked. Attaching a rule to a Protected *outer* head through an unprotected tag is the intended use.
- Attributes `{HoldAll, Protected, SequenceHold}`.
- See also: `TagSetDelayed`, `TagUnset`, `UpSet`, `UpValues`.
