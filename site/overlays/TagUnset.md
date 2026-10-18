### Worked examples

```mathematica
In[1]:= h /: f[h[x_]] := fh[x]; f[h[5]]   (* an UpValue on h *)
```

```mathematica
In[1]:= f[h[x_]] =.; f[h[5]]   (* plain Unset does not reach an UpValue, so the rule still fires *)
```

```mathematica
In[1]:= h /: f[h[y_]] =.; f[h[5]]   (* with the tag the rule is removed; the left side is compared up to renaming of pattern variables *)
```

```mathematica
In[1]:= UpValues[h]   (* nothing left on h *)
```

```mathematica
In[1]:= tt /: tt[1] = a; tt /: tt[1] =.; tt[1]   (* when the tag is the head, a DownValue is removed *)
```

```mathematica
In[1]:= sq /: area[sq[s_]] := s^2 /; s > 0; sq /: area[sq[s_]] =.; Length[UpValues[sq]]   (* a stored condition is part of the key, so this removes nothing *)
```

```mathematica
In[1]:= sq /: area[sq[s_]] /; s > 0 =.; UpValues[sq]   (* repeat the condition to remove the rule *)
```

```mathematica
In[1]:= FullForm[Hold[h /: f[h] =.]]   (* the operator form f /: lhs =. is TagUnset *)
```

### Notes

- `f /: lhs =.` removes the own-, down- or up-value of `f` whose left side equals `lhs`, after canonicalising and renaming pattern variables, and returns `Null`.
- The arguments of `lhs` are evaluated the same way the original assignment evaluated them, so the key matches what was stored.
- When no rule matches, nothing happens and no message is given (Mathematica issues `TagUnset::norep`). A non-symbol tag gives `TagUnset::sym`.
- Attributes `{HoldAll, Protected, SequenceHold}`.
- See also: `TagSet`, `TagSetDelayed`, `Unset`, `UpValues`, `Clear`.
