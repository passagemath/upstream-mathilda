# StandardForm

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`StandardForm[expr]`**

displays expr in Mathilda's own typeset notation: Sin\[x\], Log\[x\], Abs\[x\], with 2D fractions, powers and radicals.

<details>
<summary>Notes</summary>

StandardForm is the notebook's default output form, set in Arial. As a display wrapper it switches the notation back inside TraditionalForm\[...\]; as text it prints as expr itself.

</details>

## Examples (11)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= StandardForm[1 + 1]
Out[1]= 2

In[2]:= FullForm[StandardForm[Sin[x]]]
Out[2]= StandardForm[Sin[x]]

In[3]:= InputForm[StandardForm[x^2]]
Out[3]= StandardForm[x^2]
```

### Applications (8)

In the terminal the wrapper is not shown; the notebook typesets it in StandardForm

```mathematica
In[4]:= StandardForm[Sin[x]^2 + Log[x]]
Out[4]= Log[x] + Sin[x]^2
```

The argument is evaluated

```mathematica
In[5]:= StandardForm[1 + 1]
Out[5]= 2
```

The wrapper stays in the expression

```mathematica
In[6]:= FullForm[StandardForm[Sin[x]]]
Out[6]= StandardForm[Sin[x]]
```

InputForm keeps the wrapper so the text reads back in

```mathematica
In[7]:= InputForm[StandardForm[x^2]]
Out[7]= StandardForm[x^2]
```

It is an ordinary head

```mathematica
In[8]:= First[StandardForm[a + b]]
Out[8]= a + b
```

As a ToString form it is the same as InputForm

```mathematica
In[9]:= ToString[x^2, StandardForm]
Out[9]= "x^2"
```

Nested inside TraditionalForm, it switches its part back; the terminal shows plain text

```mathematica
In[10]:= TraditionalForm[{Sin[x], StandardForm[Cos[x]]}]
Out[10]= {Sin[x], Cos[x]}
```

Only Protected

```mathematica
In[11]:= Attributes[StandardForm]
Out[11]= {Protected}
```

## Implementation notes

**Semantics.** `StandardForm` is registered to `builtin_display_form`, which always returns `NULL`, so `StandardForm[expr]` evaluates its argument and keeps the wrapper (attribute `Protected`). It is a display directive, not a transformation.

**Rendering.** The plain-text printer `print_standard` (`src/print.c`) prints a `StandardForm[x]` node as `x` itself, the same way it treats `HoldForm`, except under `InputForm`, where the wrapper is printed so the result reads back in. The typeset renderer `src/print_latex.c` keeps a notation flag `LBuf.std`. A `StandardForm[...]` node sets the flag for its subtree and restores it afterwards, giving the language's own names in upright type (`Sin[x]`, `Log[x]`, `x == 1`, `a && b`, ⅇ and ⅈ for `E` and `I`) with two-dimensional fractions, powers and radicals. In the notebook pipe protocol (`src/repl.c`) each result is sent as `slatex` (StandardForm, the default display) and `latex` (TraditionalForm). `slatex` is selectable: each subexpression is tagged `\htmlData{mp=K}`, and its InputForm text is sent in `sparts`. A result with more than 4096 subexpressions or 256 KB of part text is sent untagged.

**Limits.** In the terminal REPL `StandardForm` changes nothing that is visible, because the text printer has a single notation. `ToString[expr, StandardForm]` is an alias for the InputForm string, and `ToExpression[s, StandardForm]` uses the ordinary parser. There is no `MakeBoxes`/box layer behind it.

**Attributes:** `Protected`.

## References

**See also:** [Sqrt](../../arithmetic/Sqrt/), [Pi](../../mathematical-constants/Pi/), [E](../../mathematical-constants/E/), [I](../../mathematical-constants/I/), [EulerGamma](../../mathematical-constants/EulerGamma/), [InputForm](../../expression-information/InputForm/)

- Source: [`src/print.c`](https://github.com/stblake/mathilda/blob/main/src/print.c)
- Specification: [`docs/spec/builtins/expression-information.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/expression-information.md)
- Tests: [`tests/test_print.c`](https://github.com/stblake/mathilda/blob/main/tests/test_print.c)

## Notes & additional examples

### Notes

- A display wrapper. In the notebook it is the default output form: the language's own names (`Sin[x]`, `Log[x]`, `x == 1`) set upright in Arial, with two-dimensional fractions, powers and radicals, and ⅇ / ⅈ for `E` / `I`. Parts of a StandardForm result can be selected and copied as InputForm.
- The terminal REPL has only text output, so `StandardForm[expr]` prints as `expr`. Only `InputForm` and `FullForm` show the wrapper.
- There is no box layer: `MakeBoxes`, `ToBoxes` and `ToExpression` of box structures are not supported. `ToExpression[s, StandardForm]` parses `s` as ordinary input.
- Attributes `{Protected}`.
- See also: `TraditionalForm`, `InputForm`, `FullForm`, `TeXForm`, `HoldForm`.
