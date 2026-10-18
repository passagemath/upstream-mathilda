### Worked examples

```mathematica
In[1]:= c0 = 3; Definition[c0]   (* an OwnValue *)
```

```mathematica
In[1]:= fd[x_] := x^2; fd[0] = 1; Definition[fd]   (* DownValues, the special case first; the immediate fd[0] = 1 is shown with := *)
```

```mathematica
In[1]:= g /: f[g[x_]] := h[x]; Definition[g]   (* an UpValue, shown in tag form *)
```

```mathematica
In[1]:= gq = 7; gq[0] = z; gq /: f[gq[x_]] := h[x]; Definition[gq]   (* all three kinds: own, then down, then up *)
```

```mathematica
In[1]:= Definition[neverdefined]   (* nothing defined *)
```

```mathematica
In[1]:= FullForm[Definition[fd]]   (* the head is only displayed this way; the expression is still Definition[fd] *)
```

```mathematica
In[1]:= SetAttributes[fl, Listable]; fl[x_] := x + 1; Definition[fl]   (* attributes are not listed, unlike Mathematica *)
```

### Notes

- `Definition[s]` stays unevaluated and is displayed as one line per stored rule: OwnValues as `s = v`, DownValues as `lhs := rhs`, UpValues as `s /: lhs := rhs`. With no rules it prints `Null`. `?s` and `Information[s]` add the same lines after the usage message.
- The rule store does not record immediate versus delayed, so every down- and up-value is shown with `:=`.
- Unlike Mathematica, it does not print `Attributes[s] = {...}`, Options, messages, NValues or SubValues.
- Attributes `{HoldAll, Protected}`.
- See also: `OwnValues`, `DownValues`, `UpValues`, `Information`, `Attributes`.
