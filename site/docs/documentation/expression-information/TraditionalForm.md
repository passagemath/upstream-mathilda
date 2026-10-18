# TraditionalForm

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`TraditionalForm[expr]`**

displays expr in textbook notation: sin(x), log(x), |x|, e^x, x = 1.

<details>
<summary>Notes</summary>

TraditionalForm is a display wrapper: the notebook typesets its argument in that notation (in the math font); as text it prints as expr itself. TraditionalForm can be nested to switch notation for one subexpression.

</details>

## Examples (10)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

The wrapper survives arithmetic

```mathematica
In[1]:= FullForm[TraditionalForm[x^2] + 1]
Out[1]= Plus[1, TraditionalForm[Power[x, 2]]]
```

```mathematica
In[2]:= InputForm[TraditionalForm[Sin[x]]]
Out[2]= TraditionalForm[Sin[x]]
```

Not a ToString form

```mathematica
In[3]:= ToString[x^2, TraditionalForm]
Out[3]= ToString[x^2, TraditionalForm]
```

### Applications (7)

The terminal prints the expression; the notebook typesets it as ln(x) + sin²(x)

```mathematica
In[4]:= TraditionalForm[Sin[x]^2 + Log[x]]
Out[4]= Log[x] + Sin[x]^2
```

The argument is evaluated before display

```mathematica
In[5]:= TraditionalForm[Integrate[x^2, x]]
Out[5]= 1/3 x^3
```

The wrapper stays in the expression and in any later arithmetic

```mathematica
In[6]:= FullForm[TraditionalForm[x^2] + 1]
Out[6]= Plus[1, TraditionalForm[Power[x, 2]]]
```

InputForm keeps it

```mathematica
In[7]:= InputForm[TraditionalForm[Sin[x]]]
Out[7]= TraditionalForm[Sin[x]]
```

Nested, it changes only that part (cos(x) in the notebook)

```mathematica
In[8]:= {Sin[x], TraditionalForm[Cos[x]]}
Out[8]= {Sin[x], Cos[x]}
```

Not supported as a ToString form

```mathematica
In[9]:= ToString[x^2, TraditionalForm]
Out[9]= ToString[x^2, TraditionalForm]
```

Only Protected

```mathematica
In[10]:= Attributes[TraditionalForm]
Out[10]= {Protected}
```

## Implementation notes

**Semantics.** `TraditionalForm` shares `builtin_display_form` with `StandardForm`. The builtin always returns `NULL`, so the argument is evaluated and the wrapper is kept (attribute `Protected`). Because the head survives, it stays inside later arithmetic: `TraditionalForm[x^2] + 1` is `Plus[1, TraditionalForm[x^2]]`.

**Rendering.** The text printer (`print_standard`) skips the wrapper and prints its argument, except under `InputForm`, where it is kept. The typeset renderer `src/print_latex.c` (`expr_to_latex_form`) clears the `LBuf.std` flag for the wrapped subtree. That gives textbook notation: `\sin(x)`, `\ln(x)`, `|x|`, `e^{x}`, `x = 1`, `a \land b`, n-th roots as surds, Greek-named symbols as letters, and `C[1]` as c₁. Nesting switches notation for that subtree only, and a `StandardForm[...]` inside switches it back. In the notebook protocol (`src/repl.c`) a result whose outermost head is `TraditionalForm` is sent with `slatex` alone, typeset in that notation. The notebook's *Convert To* menu and its *View › Output Format* setting choose between the forms on the client side, without re-evaluating.

**Limits.** The terminal REPL has no typeset output, so `TraditionalForm` prints like the bare expression there. `TeXForm` uses a separate renderer (`print_tex` in `src/print.c`) and does not read the notation flag. `ToString[expr, TraditionalForm]` is not supported and stays unevaluated. There is no box layer, so `TraditionalForm` output cannot be parsed back.

**Attributes:** `Protected`.

## References

**See also:** [StandardForm](../../expression-information/StandardForm/), [InputForm](../../expression-information/InputForm/)

- Source: [`src/print.c`](https://github.com/stblake/mathilda/blob/main/src/print.c)
- Specification: [`docs/spec/builtins/expression-information.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/expression-information.md)
- Tests: [`tests/test_print.c`](https://github.com/stblake/mathilda/blob/main/tests/test_print.c)

## Notes & additional examples

### Notes

- A display wrapper for textbook notation: `sin(x)`, `ln(x)`, `|x|`, `e^x`, `x = 1`, `a ∧ b`, surds for n-th roots, Greek letters for Greek-named symbols, `c₁` for `C[1]`. In the notebook, an outermost `TraditionalForm[...]` sets how the whole result is shown, and a nested one sets it for that part.
- The terminal REPL has no typeset output, so it prints the wrapped expression as plain text. Only `InputForm` and `FullForm` show the head.
- The head is kept, so wrap a result only for display. `TraditionalForm[x^2] + 1` is `Plus[1, TraditionalForm[x^2]]`.
- There is no box layer: TraditionalForm output cannot be parsed back, and `ToString[expr, TraditionalForm]` stays unevaluated. Use `TeXForm` for LaTeX text.
- Attributes `{Protected}`.
- See also: `StandardForm`, `TeXForm`, `InputForm`, `FullForm`.
