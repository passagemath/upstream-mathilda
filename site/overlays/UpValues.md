### Worked examples

```mathematica
In[1]:= g /: f[g[x_]] := h[x]; UpValues[g]   (* one rule, stored on g *)
```

```mathematica
In[1]:= UpValues[neverdefined]   (* an empty list for a symbol with no UpValues *)
```

```mathematica
In[1]:= UpValues["g"]   (* a string naming the symbol also works *)
```

```mathematica
In[1]:= UpValues[k] = {HoldPattern[w[k]] :> done}; w[k]   (* assigning a list replaces all UpValues of k *)
```

```mathematica
In[1]:= UpValues[k] = {}; w[k]   (* and an empty list removes them *)
```

```mathematica
In[1]:= area[sq[s_]] ^:= s^2 /; s > 0; UpValues[sq]   (* a moved condition shows up inside the HoldPattern *)
```

```mathematica
In[1]:= Clear[g]; UpValues[g]   (* Clear removes UpValues as well *)
```

```mathematica
In[1]:= Attributes[UpValues]   (* the argument is held, so a symbol that has a value is not replaced by it *)
```

### Notes

- Returns `{HoldPattern[lhs] :> rhs, ...}` for the rules made by `^=`, `^:=`, `f /: ... =` and `f /: ... :=`, in matching order, and `{}` when there are none. Immediate rules are shown with `:>` too.
- Accepts a symbol, or a string naming an existing symbol (so `UpValues /@ Names["s*"]` works). A string naming no symbol gives `UpValues::sym` and stays unevaluated.
- `UpValues[s] = {rules}` replaces the whole list; `HoldPattern` wrappers on the left sides are optional. A Protected `s` is refused with `Set::write`.
- Attributes `{HoldAll, Protected}`.
- See also: `DownValues`, `OwnValues`, `Definition`, `UpSet`, `TagSet`.
