# Definition

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Definition[s]`**

prints the own-, down- and up-values associated with s.

## Examples (10)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= g /: f[g[x_]] := h[x]; Definition[g]
Out[1]= g /: f[g[x_]] := h[x]

In[2]:= c0 = 3; Definition[c0]
Out[2]= c0 = 3
```

Display only

```mathematica
In[3]:= fd[x_] := x^2; FullForm[Definition[fd]]
Out[3]= Definition[fd]
```

### Applications (7)

An OwnValue

```mathematica
In[4]:= c0 = 3; Definition[c0]
Out[4]= c0 = 3
```

DownValues, the special case first; the immediate fd[0] = 1 is shown with :=

```mathematica
In[5]:= fd[x_] := x^2; fd[0] = 1; Definition[fd]
Out[5]= fd[0] := 1
fd[x_] := x^2
```

An UpValue, shown in tag form

```mathematica
In[6]:= g /: f[g[x_]] := h[x]; Definition[g]
Out[6]= g /: f[g[x_]] := h[x]
```

All three kinds: own, then down, then up

```mathematica
In[7]:= gq = 7; gq[0] = z; gq /: f[gq[x_]] := h[x]; Definition[gq]
Out[7]= gq = 7
gq[0] := z
gq /: f[gq[x_]] := h[x]
```

Nothing defined

```mathematica
In[8]:= Definition[neverdefined]
```

The head is only displayed this way; the expression is still Definition[fd]

```mathematica
In[9]:= FullForm[Definition[fd]]
Out[9]= Definition[fd]
```

Attributes are not listed, unlike Mathematica

```mathematica
In[10]:= SetAttributes[fl, Listable]; fl[x_] := x + 1; Definition[fl]
Out[10]= fl[x_] := x + 1
```

## Implementation notes

**Algorithm.** `Definition` has no evaluator builtin: `Definition[s]` (attributes `HoldAll`, `Protected`) stays as an inert expression, and the printer renders it. In `print_standard`, a `Definition[sym]` node calls `print_definition_body`, which reads `symtab_get_own_values`, `symtab_get_down_values` and `symtab_get_up_values` and prints one line per rule: OwnValues as `lhs = rhs`, DownValues as `lhs := rhs`, UpValues as `s /: lhs := rhs`, in that order and in each list's stored order. A symbol with no values prints `Null`. `render_definitions` in `src/core.c` builds the same text into a string for `?s` / `Information[s]`, which appends it after the usage message.

**Data structures.** It reads the three `Rule` lists on the `SymbolDef` directly and does not copy them. The patterns are printed with the ordinary `print_standard`, so a conditioned rule shows its `Condition` on the left (`lhs /; test := rhs`). DownValues come out in `add_rule`'s specificity order, so a literal `fd[0]` rule appears before the general `fd[x_]` one.

**Complexity / limits.** Linear in the total number of rules. Display only: `FullForm[Definition[s]]` is `Definition[s]`. The rule store has no immediate/delayed flag, so every down- and up-value prints with `:=` and every own-value with `=`. Attributes, Options, messages, NValues and SubValues are not printed (Mathematica includes `Attributes[s] = {...}` lines). Compiled in under `UP_VALUES`.

**Attributes:** `HoldAll`, `Protected`.

## References

**See also:** [FullForm](../../expression-information/FullForm/)

- Source: [`src/print.c`](https://github.com/stblake/mathilda/blob/main/src/print.c)
- Specification: [`docs/spec/builtins/expression-information.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/expression-information.md)
- Tests: [`tests/test_upvalues.c`](https://github.com/stblake/mathilda/blob/main/tests/test_upvalues.c)

## Notes & additional examples

### Notes

- `Definition[s]` stays unevaluated and is displayed as one line per stored rule: OwnValues as `s = v`, DownValues as `lhs := rhs`, UpValues as `s /: lhs := rhs`. With no rules it prints `Null`. `?s` and `Information[s]` add the same lines after the usage message.
- The rule store does not record immediate versus delayed, so every down- and up-value is shown with `:=`.
- Unlike Mathematica, it does not print `Attributes[s] = {...}`, Options, messages, NValues or SubValues.
- Attributes `{HoldAll, Protected}`.
- See also: `OwnValues`, `DownValues`, `UpValues`, `Information`, `Attributes`.
