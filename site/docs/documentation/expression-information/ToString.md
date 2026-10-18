# ToString

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`ToString[expr]`**

gives the printed form of expr (as OutputForm) as a String; strings inside expr appear without quotes, ToString\["s"\] is "s".

**`ToString[expr, form]`**

uses the specified output form.

<details>
<summary>Notes</summary>

Supported forms: OutputForm (default, linear), InputForm, FullForm, TeXForm.

</details>

## Examples (16)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (5)

```mathematica
In[1]:= ToString[x^2 + y^3]
Out[1]= "x^2 + y^3"

In[2]:= ToString[x^2 + y^3, FullForm]
Out[2]= "Plus[Power[x, 2], Power[y, 3]]"

In[3]:= ToString[x^2 + y^3, TeXForm]
Out[3]= "x^{2}+y^{3}"

In[4]:= ToString[{"a", 1}]
Out[4]= "{a, 1}"

In[5]:= ToString["a\"b", InputForm]
Out[5]= "\"a\\\"b\""
```

### Applications (11)

The default is OutputForm

```mathematica
In[6]:= ToString[x^2 + y^3]
Out[6]= "x^2 + y^3"
```

A string is its own OutputForm: no added quotes

```mathematica
In[7]:= ToString["hello"]
Out[7]= "hello"
```

Strings inside an expression print as their characters

```mathematica
In[8]:= ToString[{"a", 1}]
Out[8]= "{a, 1}"
```

Under InputForm, quotes and backslashes are escaped

```mathematica
In[9]:= ToString["a\"b", InputForm]
Out[9]= "\"a\\\"b\""
```

The internal tree

```mathematica
In[10]:= ToString[{1, "a"}, FullForm]
Out[10]= "List[1, \"a\"]"
```

TeX for typesetting

```mathematica
In[11]:= ToString[x^2 + y^3, TeXForm]
Out[11]= "x^{2}+y^{3}"
```

Rule-delayed stays distinct from rule

```mathematica
In[12]:= ToString[a :> b, TeXForm]
Out[12]= "a:\\to b"
```

A generated constant becomes a subscript

```mathematica
In[13]:= ToString[C[10], TeXForm]
Out[13]= "c_{10}"
```

The decimal digits of 2^100

```mathematica
In[14]:= StringLength[ToString[2^100]]
Out[14]= 31
```

The InputForm text reads back

```mathematica
In[15]:= ToExpression[ToString[(a + b)^2/c, InputForm]] === (a + b)^2/c
Out[15]= True
```

An unknown form stays unevaluated

```mathematica
In[16]:= ToString[x, FooForm]
Out[16]= ToString[x, FooForm]
```

## Implementation notes

**Algorithm.** `builtin_tostring` (`src/core.c`) renders its argument to a string through the shared printer in `src/print.c`. It takes 1–2 arguments, and the second must be a form symbol. As in Mathematica, the default form is `OutputForm` (since v0.455; before that it was `InputForm`). A string argument comes back as itself, with no render and no quotes. Anything else goes through `expr_to_string_outputform`, which prints strings anywhere in the tree as their characters (`ToString[{"a", 1}]` is `"{a, 1}"`). The other forms select a printer:

- `InputForm` and `StandardForm` use `expr_to_string`, which keeps string quotes and escapes `"` and `\`, so the result reads back.
- `FullForm` uses `expr_to_string_fullform`.
- `TeXForm` wraps the argument in `TeXForm[...]` and prints that. This sets `C[k]` as `c_k` and infix heads as `\to`, `:\to`, `\land`, …

**Data structures.** The printer writes into a heap buffer, a `tmpfile()` round trip on `_WIN32`. The buffer is copied into a fresh `EXPR_STRING`. Precedence-aware parenthesisation, `Rational`/`Complex` printing and the infix `StringExpression`/`Pattern`/`Repeated` forms (`"a"~~x_`, `x:p`, `p..`) all come from `print_standard`.

**Complexity / limits.** Linear in the size of the printed form. Mathilda's `OutputForm` is one-dimensional (`ToString[1/2]` is `"1/2"`), whereas Mathematica's is a 2-D layout. An unknown form such as `ToString[x, FooForm]` leaves the call unevaluated, so a typo shows at the call site. `Protected`.

- `Protected`.
- An unsupported form leaves the call unevaluated (e.g. `ToString[x, FooForm]` returns `ToString[x, FooForm]`), so a typo is visible at the call site rather than silently downgraded.
- `TeXForm` renders the generated constant `C[k]` (from `DSolve`, `Reduce`, `Integrate`) as the subscripted `c_k`, matching Mathematica: `ToString[C[1], TeXForm]` is `"c_1"`, and `ToString[C[10], TeXForm]` is `"c_{10}"` (single-character subscripts bare, longer ones braced). The same holds for the notebook LaTeX renderer.
- `TeXForm` sets the infix heads infix: `Rule` as `\to`, `RuleDelayed` as `:\to`, `Set` as `=`, `SetDelayed` as `:=`, `And`/`Or` as `\land`/`\lor`, the relations as `=`, `\neq`, `<`, `>`, `\leq`, `\geq`, `SameQ`/`UnsameQ` as `\equiv`/`\not\equiv`, `Not` as `\neg`, and a chained `Inequality` as the chain (`1<x<2`). `RuleDelayed` is deliberately distinguished from `Rule` — one glyph for both would make `a -> b` and `a :> b` indistinguishable once typeset, the same reason `SetDelayed` is `:=`. The notebook LaTeX renderer agrees, and additionally brackets a nested operator that binds looser than the one holding it (`(a || b) && c`, but not `a && b || c`).

**Attributes:** `Protected`.

## References

**See also:** [String](../../other-advanced/String/), [InputForm](../../expression-information/InputForm/), [FullForm](../../expression-information/FullForm/), [TeXForm](../../expression-information/TeXForm/), [StandardForm](../../expression-information/StandardForm/), [DSolve](../../calculus/DSolve/), [Reduce](../../solutions-of-equations/Reduce/), [Integrate](../../calculus/Integrate/)

- Source: [`src/core.c`](https://github.com/stblake/mathilda/blob/main/src/core.c)
- Specification: [`docs/spec/builtins/expression-information.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/expression-information.md)
- Tests: [`tests/test_assoc_forms.c`](https://github.com/stblake/mathilda/blob/main/tests/test_assoc_forms.c)
- Tests: [`tests/test_assoc_ops.c`](https://github.com/stblake/mathilda/blob/main/tests/test_assoc_ops.c)
- Tests: [`tests/test_blas.c`](https://github.com/stblake/mathilda/blob/main/tests/test_blas.c)
- Tests: [`tests/test_compile_assoc.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile_assoc.c)

## Notes & additional examples

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
