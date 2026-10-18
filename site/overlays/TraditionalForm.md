### Worked examples

```mathematica
In[1]:= TraditionalForm[Sin[x]^2 + Log[x]]   (* the terminal prints the expression; the notebook typesets it as ln(x) + sin²(x) *)
```

```mathematica
In[1]:= TraditionalForm[Integrate[x^2, x]]   (* the argument is evaluated before display *)
```

```mathematica
In[1]:= FullForm[TraditionalForm[x^2] + 1]   (* the wrapper stays in the expression and in any later arithmetic *)
```

```mathematica
In[1]:= InputForm[TraditionalForm[Sin[x]]]   (* InputForm keeps it *)
```

```mathematica
In[1]:= {Sin[x], TraditionalForm[Cos[x]]}   (* nested, it changes only that part (cos(x) in the notebook) *)
```

```mathematica
In[1]:= ToString[x^2, TraditionalForm]   (* not supported as a ToString form *)
```

```mathematica
In[1]:= Attributes[TraditionalForm]   (* only Protected *)
```

### Notes

- A display wrapper for textbook notation: `sin(x)`, `ln(x)`, `|x|`, `e^x`, `x = 1`, `a ∧ b`, surds for n-th roots, Greek letters for Greek-named symbols, `c₁` for `C[1]`. In the notebook, an outermost `TraditionalForm[...]` sets how the whole result is shown, and a nested one sets it for that part.
- The terminal REPL has no typeset output, so it prints the wrapped expression as plain text. Only `InputForm` and `FullForm` show the head.
- The head is kept, so wrap a result only for display. `TraditionalForm[x^2] + 1` is `Plus[1, TraditionalForm[x^2]]`.
- There is no box layer: TraditionalForm output cannot be parsed back, and `ToString[expr, TraditionalForm]` stays unevaluated. Use `TeXForm` for LaTeX text.
- Attributes `{Protected}`.
- See also: `StandardForm`, `TeXForm`, `InputForm`, `FullForm`.
