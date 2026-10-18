### Worked examples

```mathematica
In[1]:= Protect[foo, bar]  (* returns the names whose Protected attribute was newly set *)
In[2]:= MemberQ[Attributes[foo], Protected]
```

```mathematica
In[1]:= Protect[foo]  (* already protected: nothing changed *)
```

```mathematica
In[1]:= Check[foo = 6, "refused"]  (* assignment to a protected symbol emits Set::wrsym *)
```

```mathematica
In[1]:= foo  (* the value was not assigned *)
```

```mathematica
In[1]:= Protect[{"k1", k2}]  (* a list of specs; strings name symbols *)
```

```mathematica
In[1]:= Attributes[Protect]
```

### Notes

- `Protect[s1, s2, ...]` sets the `Protected` attribute on each symbol. That is
  what `Set`, `SetDelayed`, `Unset`, `ClearAll`, `Remove` and `SetAttributes`
  consult before changing a symbol.
- It returns the names, as strings, whose state actually changed, so protecting
  an already-protected symbol gives `{}`.
- Arguments may be symbols, strings, or a flat list of them. `HoldAll | Protected`.
- Every built-in carries `Protected`, and `Unprotect` always lifts it. No
  attribute seals a symbol permanently (since v0.301).
- See also: `Unprotect`, `Attributes`, `SetAttributes`.
