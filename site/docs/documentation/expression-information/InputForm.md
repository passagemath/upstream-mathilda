# InputForm

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`InputForm[expr]`**

prints expr in a form suitable to be re-read by the parser, using operator syntax (a + b, not Plus\[a, b\]) and explicit string quotes.

<details>
<summary>Notes</summary>

Like FullForm, InputForm is a printer wrapper: it is consumed during output and does not appear in the printed result.

</details>

## Examples (12)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= InputForm[{1, 1/2, "a", x}]
Out[1]= {1, 1/2, "a", x}

In[2]:= InputForm["say \"hi\""]
Out[2]= "say \"hi\""

In[3]:= InputForm[Series[Exp[x], {x, 0, 3}]]
Out[3]= SeriesData[x, 0, {1, 1, 1/2, 1/6}, 0, 4, 1]
```

### Applications (9)

Linear syntax

```mathematica
In[4]:= InputForm[1/2 + x^2]
Out[4]= 1/2 + x^2
```

Strings keep their quotes

```mathematica
In[5]:= InputForm[{1, 1/2, "a", x}]
Out[5]= {1, 1/2, "a", x}
```

Embedded quotes are escaped

```mathematica
In[6]:= InputForm["say \"hi\""]
Out[6]= "say \"hi\""
```

Explicit * before I

```mathematica
In[7]:= InputForm[Sqrt[2] + 3 I]
Out[7]= 3*I + Sqrt[2]
```

A series shows its SeriesData form

```mathematica
In[8]:= InputForm[Series[Exp[x], {x, 0, 3}]]
Out[8]= SeriesData[x, 0, {1, 1, 1/2, 1/6}, 0, 4, 1]
```

Display wrappers are kept, so the text re-reads

```mathematica
In[9]:= InputForm[TraditionalForm[Sin[x]]]
Out[9]= TraditionalForm[Sin[x]]
```

String patterns print infix

```mathematica
In[10]:= InputForm["a" ~~ x_ ~~ "b"]
Out[10]= "a"~~x_~~"b"
```

An inert wrapper that stays in the expression

```mathematica
In[11]:= Head[InputForm[x]]
Out[11]= InputForm
```

Round trip through text

```mathematica
In[12]:= ToExpression[ToString[(a + b)^2/c, InputForm]] === (a + b)^2/c
Out[12]= True
```

## Implementation notes

**Mechanism.** `InputForm` is an inert display wrapper. `builtin_inputform` (`src/print.c`) returns `NULL`, so `InputForm[expr]` stays in the tree (`Head[InputForm[x]]` is `InputForm`), and the printer does the work. When `print_standard` meets the head, it increments `g_inputform_depth` and prints the argument in its linear, re-parseable form. While that depth is non-zero several display-only paths are switched off:

- `SeriesData` prints as the literal `SeriesData[x, x0, {...}, nmin, nmax, den]` rather than `1 + x + … + O[x]^4`.
- `TraditionalForm`/`StandardForm` wrappers are kept, not made transparent.
- `NumberForm`, `Row`, `Graphics` and `Graph` are not given their display rendering.

**Strings.** `print_string_literal` quotes strings and escapes `"`, `\`, newline, tab and carriage return. An explicit `InputForm[...]` keeps the quotes even inside `Print`, whose own traversal otherwise prints strings as bare characters (OutputForm). This was fixed in v0.455. `ToString[expr, InputForm]` and `ToString[expr, StandardForm]` go through the same printer (`expr_to_string`).

**Limits.** InputForm round-trips within Mathilda: `ToExpression[ToString[e, InputForm]] === e`. Two spellings differ from Mathematica's. Products print with `*` before `I` (`3*I + Sqrt[2]`), and machine reals in scientific range print as C-style `1.5e+20`, which Mathilda's parser reads back but Mathematica would not (it writes `1.5*^20`). In the notebook, `InputForm` suppresses typesetting and shows the text as code. `Protected`.

**Attributes:** `Protected`.

## References

**See also:** [FullForm](../../expression-information/FullForm/), [TeXForm](../../expression-information/TeXForm/), [NumberForm](../../expression-information/NumberForm/), [StandardForm](../../expression-information/StandardForm/), [TraditionalForm](../../expression-information/TraditionalForm/), [HoldForm](../../expression-information/HoldForm/), [Print](../../expression-information/Print/), [SeriesData](../../power-series/SeriesData/)

- Source: [`src/print.c`](https://github.com/stblake/mathilda/blob/main/src/print.c)
- Specification: [`docs/spec/builtins/expression-information.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/expression-information.md)
- Tests: [`tests/test_graph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph.c)
- Tests: [`tests/test_graph_ops.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_ops.c)
- Tests: [`tests/test_hypergraph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_hypergraph.c)
- Tests: [`tests/test_numeric_largearg.c`](https://github.com/stblake/mathilda/blob/main/tests/test_numeric_largearg.c)

## Notes & additional examples

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
