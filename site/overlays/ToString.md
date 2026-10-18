### Worked examples

```mathematica
In[1]:= ToString[x^2 + y^3]  (* the default is OutputForm *)
```

```mathematica
In[1]:= ToString["hello"]  (* a string is its own OutputForm: no added quotes *)
```

```mathematica
In[1]:= ToString[{"a", 1}]  (* strings inside an expression print as their characters *)
```

```mathematica
In[1]:= ToString["a\"b", InputForm]  (* under InputForm, quotes and backslashes are escaped *)
```

```mathematica
In[1]:= ToString[{1, "a"}, FullForm]  (* the internal tree *)
```

```mathematica
In[1]:= ToString[x^2 + y^3, TeXForm]  (* TeX for typesetting *)
```

```mathematica
In[1]:= ToString[a :> b, TeXForm]  (* rule-delayed stays distinct from rule *)
```

```mathematica
In[1]:= ToString[C[10], TeXForm]  (* a generated constant becomes a subscript *)
```

```mathematica
In[1]:= StringLength[ToString[2^100]]  (* the decimal digits of 2^100 *)
```

```mathematica
In[1]:= ToExpression[ToString[(a + b)^2/c, InputForm]] === (a + b)^2/c  (* the InputForm text reads back *)
```

```mathematica
In[1]:= ToString[x, FooForm]  (* an unknown form stays unevaluated *)
```

### Notes

- Since v0.455 `ToString[expr]` uses `OutputForm` by default, as Mathematica
  does, so strings come out without quotes: `ToString["hello"]` is `"hello"`.
  It used to be `"\"hello\""`.
- The supported forms are `OutputForm` (the default), `InputForm`, `FullForm`
  and `TeXForm`. `StandardForm` is accepted as an alias of `InputForm`.
  `InputForm` keeps string quotes and escapes `"` and `\`, so the result can be
  read back with `ToExpression`.
- Mathilda's `OutputForm` is one-dimensional (`ToString[1/2]` is `"1/2"`),
  whereas Mathematica's is a 2-D layout.
- Under `TeXForm` the generated constants `C[k]` render as `c_k` (single
  characters bare, longer subscripts braced), and `Rule`/`RuleDelayed` render
  as `\to` / `:\to`.
- An unknown form leaves the call unevaluated, so a typo shows at the call site.
  `Protected`.
- See also: `ToExpression`, `InputForm`, `FullForm`, `TeXForm`, `StringJoin`.
