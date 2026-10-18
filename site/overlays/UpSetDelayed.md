### Worked examples

```mathematica
In[1]:= area[sq[s_]] ^:= s^2; area[sq[3]]   (* the rule lives on sq and fires on area[sq[...]] *)
```

```mathematica
In[1]:= UpValues[sq]   (* stored as a delayed HoldPattern rule *)
```

```mathematica
In[1]:= perim[sq[s_]] ^:= 4 s; {perim[sq[a]], perim[circle[a]]}   (* only calls with an sq argument are touched *)
```

```mathematica
In[1]:= ctr = 0; tick[t] ^:= ++ctr; {tick[t], tick[t], tick[t]}   (* delayed: the right side is re-evaluated each time the rule fires *)
```

```mathematica
In[1]:= vol[cube[s_]] ^:= s^3 /; s > 0; {vol[cube[2]], vol[cube[-2]]}   (* a trailing condition is moved onto the pattern, so a failed test leaves the call unevaluated *)
```

```mathematica
In[1]:= md[a_] + md[b_] ^:= md[a + b]; md[1] + md[2] + md[3] + z   (* on the Flat, Orderless head Plus, a two-term rule merges pairs inside a longer sum and keeps the rest *)
```

```mathematica
In[1]:= Check[Length[x_] ^:= 0, failed]   (* a bare blank names no symbol: UpSetDelayed::nosym *)
```

### Notes

- `lhs ^:= rhs` keeps `rhs` unevaluated (`UpSetDelayed` is `HoldAll`) and returns `Null`. It attaches the rule to the same level-one symbols as `UpSet`.
- `lhs ^:= body /; test` is stored as the pattern `lhs /; test`, as with `SetDelayed`.
- `UpValues` and `Definition` show immediate and delayed rules the same way (`:>` / `:=`), because the rule store does not record which kind made them.
- Limitation: on a `Flat` head that is not `Orderless` (`Dot`, `**`), the matched operands must form a leading run.
- Attributes `{HoldAll, Protected, SequenceHold}`.
- See also: `UpSet`, `TagSetDelayed`, `UpValues`, `SetDelayed`.
