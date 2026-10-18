### Worked examples

```mathematica
In[1]:= InputForm[1/2 + x^2]  (* linear syntax *)
```

```mathematica
In[1]:= InputForm[{1, 1/2, "a", x}]  (* strings keep their quotes *)
```

```mathematica
In[1]:= InputForm["say \"hi\""]  (* embedded quotes are escaped *)
```

```mathematica
In[1]:= InputForm[Sqrt[2] + 3 I]  (* explicit * before I *)
```

```mathematica
In[1]:= InputForm[Series[Exp[x], {x, 0, 3}]]  (* a series shows its SeriesData form *)
```

```mathematica
In[1]:= InputForm[TraditionalForm[Sin[x]]]  (* display wrappers are kept, so the text re-reads *)
```

```mathematica
In[1]:= InputForm["a" ~~ x_ ~~ "b"]  (* string patterns print infix *)
```

```mathematica
In[1]:= Head[InputForm[x]]  (* an inert wrapper that stays in the expression *)
```

```mathematica
In[1]:= ToExpression[ToString[(a + b)^2/c, InputForm]] === (a + b)^2/c  (* round trip through text *)
```

### Notes

- `InputForm[expr]` prints `expr` in a linear form the parser can read back:
  operator syntax (`a + b`, not `Plus[a, b]`) with explicit string quotes. Use
  `FullForm` to see the internal tree.
- It is an inert wrapper: the head stays in the expression, and only the printer
  acts on it. Under it, `SeriesData`, `NumberForm`, `Row` and the
  `TraditionalForm`/`StandardForm` wrappers lose their display-only rendering.
- An explicit `InputForm[...]` inside `Print` keeps string quotes (v0.455).
  `Print` otherwise shows strings as bare characters. `ToString[e, InputForm]`
  gives the same text as a `String`.
- Two spellings differ from Mathematica's. Machine reals in scientific range
  print as `1.5e+20`, which Mathilda's parser reads back but Mathematica writes
  as `1.5*^20`. Products put an explicit `*` before `I`.
- In the notebook, `InputForm` suppresses typesetting and shows the text as code.
- See also: `FullForm`, `ToString`, `ToExpression`, `StandardForm`.
