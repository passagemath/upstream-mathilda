### Worked examples

```mathematica
In[1]:= tmpvar = 7  (* a user symbol *)
In[2]:= Remove[tmpvar]
In[3]:= tmpvar  (* a fresh, undefined symbol *)
```

```mathematica
In[1]:= r = 7; Remove[r]; Names["r"]  (* the entry is gone from the symbol table *)
```

```mathematica
In[1]:= Options[o2] = {opt -> 1}; Remove[o2]; Options[o2]  (* default options go with it *)
```

```mathematica
In[1]:= Remove[Plus]; 1 + 1  (* built-ins are Protected, so they are skipped *)
```

```mathematica
In[1]:= Attributes[Remove]  (* protected like every built-in *)
```

### Notes

- `Remove[s]` deletes the symbol `s` from the symbol table entirely. This is
  stronger than `ClearAll`, which empties a symbol but keeps its entry and its
  default options.
- After removal the name no longer exists, and the next reference to it creates
  a fresh, undefined symbol. The result is `Null`.
- Arguments may be symbols, strings, or a flat list of them.
- A `Protected` symbol is skipped, so `Remove` never deletes a built-in.
  `Attributes[Remove]` is `{HoldAll, Protected}`.
- See also: `ClearAll`, `Clear`, `Names`.
