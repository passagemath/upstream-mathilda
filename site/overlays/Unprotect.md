### Worked examples

```mathematica
In[1]:= Protect[baz]
In[2]:= Unprotect[baz]  (* returns the names whose Protected attribute was actually cleared *)
In[3]:= MemberQ[Attributes[baz], Protected]
```

```mathematica
In[1]:= Unprotect[baz]  (* not protected any more: nothing changed *)
```

```mathematica
In[1]:= baz = 6; baz  (* assignment works again *)
```

```mathematica
In[1]:= Unprotect[Sin]  (* built-ins can be unprotected too *)
```

```mathematica
In[1]:= Protect[Sin]  (* and re-protected after extending them *)
```


### Notes

- `Unprotect[s1, s2, ...]` is the inverse of `Protect`. It clears the
  `Protected` attribute and returns the names, as strings, whose state changed.
  Unprotecting a symbol that was not protected gives `{}`.
- Once a symbol is unprotected, `Set`/`SetDelayed` and the clearing heads
  accept it again. The usual way to extend a built-in is
  `Unprotect[f]; new rule; Protect[f]`.
- Nothing is sealed beyond `Unprotect`'s reach: `Protected` is the only
  sealing attribute (since v0.301).
- Arguments may be symbols, strings, or a flat list. `HoldAll | Protected`.
- See also: `Protect`, `Attributes`, `ClearAttributes`.
