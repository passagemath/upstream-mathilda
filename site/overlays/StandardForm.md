### Worked examples

```mathematica
In[1]:= StandardForm[Sin[x]^2 + Log[x]]   (* in the terminal the wrapper is not shown; the notebook typesets it in StandardForm *)
```

```mathematica
In[1]:= StandardForm[1 + 1]   (* the argument is evaluated *)
```

```mathematica
In[1]:= FullForm[StandardForm[Sin[x]]]   (* the wrapper stays in the expression *)
```

```mathematica
In[1]:= InputForm[StandardForm[x^2]]   (* InputForm keeps the wrapper so the text reads back in *)
```

```mathematica
In[1]:= First[StandardForm[a + b]]   (* it is an ordinary head *)
```

```mathematica
In[1]:= ToString[x^2, StandardForm]   (* as a ToString form it is the same as InputForm *)
```

```mathematica
In[1]:= TraditionalForm[{Sin[x], StandardForm[Cos[x]]}]   (* nested inside TraditionalForm, it switches its part back; the terminal shows plain text *)
```

```mathematica
In[1]:= Attributes[StandardForm]   (* only Protected *)
```

### Notes

- A display wrapper. In the notebook it is the default output form: the language's own names (`Sin[x]`, `Log[x]`, `x == 1`) set upright in Arial, with two-dimensional fractions, powers and radicals, and ⅇ / ⅈ for `E` / `I`. Parts of a StandardForm result can be selected and copied as InputForm.
- The terminal REPL has only text output, so `StandardForm[expr]` prints as `expr`. Only `InputForm` and `FullForm` show the wrapper.
- There is no box layer: `MakeBoxes`, `ToBoxes` and `ToExpression` of box structures are not supported. `ToExpression[s, StandardForm]` parses `s` as ordinary input.
- Attributes `{Protected}`.
- See also: `TraditionalForm`, `InputForm`, `FullForm`, `TeXForm`, `HoldForm`.
