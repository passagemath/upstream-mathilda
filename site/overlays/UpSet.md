### Worked examples

```mathematica
In[1]:= area[square] ^= 4; area[square]   (* the value is attached to square, the argument, not to area *)
```

```mathematica
In[1]:= {UpValues[square], DownValues[area]}   (* stored as an UpValue of square; area has no DownValues *)
```

```mathematica
In[1]:= prop[a, b[c]] ^= val; {UpValues[a], UpValues[b]}   (* every level-one symbol gets the rule: a as an argument, b as the head of an argument *)
```

```mathematica
In[1]:= Length[vec] ^= 3; {Length[vec], Length[{1, 2}]}   (* a rule for the Protected head Length, attached to vec, so Length need not be unprotected *)
```

```mathematica
In[1]:= ctr = 0; stamp[t] ^= ctr++; {stamp[t], stamp[t], ctr}   (* immediate: the right side was evaluated once, when the assignment was made *)
```

```mathematica
In[1]:= fdv[x_] := down; fdv[gup] ^= up; {fdv[gup], fdv[other]}   (* an UpValue is tried before the DownValues of the enclosing head *)
```

```mathematica
In[1]:= Check[fx[x_] ^= 1, failed]   (* no symbol at level one to attach to: UpSet::nosym, nothing is stored *)
```

```mathematica
In[1]:= Check[Sin[Pi] ^= 2, failed]   (* Pi is Protected: UpSet::write *)
```

### Notes

- `lhs ^= rhs` evaluates `rhs` once (`UpSet` is `HoldFirst`) and returns it. The rule is attached to **every** distinct symbol at level one of `lhs`: each argument that is a symbol, the head of each argument, and the head a blank names (`a_h`). The head of `lhs` itself gets nothing.
- The arguments of `lhs` are evaluated before storing, as with `Set`; pattern-bearing arguments stay as written.
- UpValues fire before the enclosing head's DownValues and built-in code, and not inside `HoldAllComplete`. This is how to give a Protected built-in (`Length`, `Plus`, `N`) a meaning for your own objects.
- Messages: `UpSet::nosym` (no level-one symbol; the call stays unevaluated) and `UpSet::write` (a Protected tag is skipped, and any other tags still get the rule). Both go through the message system, so `Quiet` and `Check` work.
- Attributes `{HoldFirst, Protected, SequenceHold}`. It is part of the `UP_VALUES` subsystem, which is on by default.
- See also: `UpSetDelayed`, `TagSet`, `UpValues`, `TagUnset`, `Set`.
