### Worked examples

```mathematica
In[1]:= v = 10
In[2]:= v =.  (* the postfix form of Unset *)
In[3]:= v
```

```mathematica
In[1]:= g[1] = a; g[2] = b; g[1] =.; DownValues[g]  (* removes only the g[1] rule *)
```

```mathematica
In[1]:= gg[x_] := x + 1; gg[y_] =.; gg[2]  (* pattern names may differ: the rule is found up to renaming *)
```

```mathematica
In[1]:= fact[1] = 1; fact[n_] := n fact[n - 1]; fact[1] =.; DownValues[fact]  (* the general rule survives *)
```

```mathematica
In[1]:= mm[x_] /; x > 0 := 1; mm[x_] /; x > 0 =.; DownValues[mm]  (* a conditional rule *)
```

```mathematica
In[1]:= Check[Unset[Pi], "refused"]  (* a Protected symbol: Unset::wrsym *)
```

```mathematica
In[1]:= undefinedthing =.  (* no matching rule: still Null *)
```

### Notes

- `lhs =.` (`Unset`) removes the one definition whose left-hand side is `lhs`,
  rather than every rule on a symbol. A bare symbol loses its OwnValue, and
  `f[args] =.` loses exactly that DownValue.
- Left-hand sides are compared after canonicalization and alpha-renaming of
  pattern variables, so `f[y_] =.` removes `f[x_] := ...`. A
  `Condition`-wrapped left-hand side is matched on its inner head.
- `HoldFirst | Protected`: the target is not evaluated before the rule is
  located. The result is `Null` whether or not a rule was found.
- A `Protected` owner is refused with `Unset::wrsym`, routed through the message
  funnel so `Check` sees it.
- See also: `Clear`, `ClearAll`, `Set`, `TagUnset`, `DownValues`.
