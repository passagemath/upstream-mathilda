### Worked examples

```mathematica
In[1]:= SetAttributes[g, Orderless]; g = 5  (* a value and an attribute *)
In[2]:= ClearAll[g]
In[3]:= {g, Attributes[g]}  (* both are gone *)
```

```mathematica
In[1]:= h[x_] := x^2; SetAttributes[h, Listable]; ClearAll[h]; {DownValues[h], Attributes[h]}  (* DownValues and attributes together *)
```

```mathematica
In[1]:= w = 3; SetAttributes[w, Flat]; Clear[w]; {w, Attributes[w]}  (* compare Clear, which keeps the attributes *)
```

```mathematica
In[1]:= p = 1; q = 2; ClearAll[{p, "q"}]; {p, q}  (* a list of specs, mixing symbols and strings *)
```

```mathematica
In[1]:= Options[o1] = {opt -> 1}; ClearAll[o1]; Options[o1]  (* default options survive; Remove deletes them *)
```

```mathematica
In[1]:= k = 1; Protect[k]; ClearAll[k]; k  (* a Protected symbol is skipped *)
```

```mathematica
In[1]:= ClearAll[Sin]; Sin[0]  (* so built-ins are never wiped *)
```

### Notes

- `ClearAll[s]` removes `s`'s OwnValues, DownValues, attributes and usage
  message, returning the symbol to its undefined state. Compare `Clear[s]`,
  which drops only the values. Registered default options are kept. `Remove`
  deletes them along with the symbol.
- Arguments may be symbols, strings naming symbols, or a flat list of them. The
  result is `Null`.
- `HoldAll | Protected`, so the symbols themselves are cleared, not their values.
- A `Protected` symbol is skipped silently, which is what keeps `ClearAll` from
  wiping a built-in. To clear one of your own protected symbols, `Unprotect` it
  first. `Protected` is the only guard, and every protection can be lifted.
- See also: `Clear`, `Remove`, `Unset`, `ClearAttributes`.
