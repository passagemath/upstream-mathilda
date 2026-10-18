# Remove

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Remove[s1, s2, ...]`**

removes the named symbols completely, deleting their definitions from the symbol table. Remove\[{s1, s2, ...}\] accepts a list of specs.

<details>
<summary>Notes</summary>

Remove has attribute HoldAll; symbols with attribute Protected are not affected.

</details>

## Examples (10)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= f[x_] := x^2; SetAttributes[f, Listable]; Attributes[f]
Out[1]= {Listable}

In[2]:= ClearAll[f]; {Attributes[f], DownValues[f]}
Out[2]= {{}, {}}

In[3]:= x = 2; Remove[x]; x
Out[3]= x
```

### Applications (7)

A user symbol

```mathematica
In[4]:= tmpvar = 7
Out[4]= 7
```

```mathematica
In[5]:= Remove[tmpvar]
```

A fresh, undefined symbol

```mathematica
In[6]:= tmpvar
Out[6]= tmpvar
```

The entry is gone from the symbol table

```mathematica
In[7]:= r = 7; Remove[r]; Names["r"]
Out[7]= {}
```

Default options go with it

```mathematica
In[8]:= Options[o2] = {opt -> 1}; Remove[o2]; Options[o2]
Out[8]= {}
```

Built-ins are Protected, so they are skipped

```mathematica
In[9]:= Remove[Plus]; 1 + 1
Out[9]= 2
```

Protected like every built-in

```mathematica
In[10]:= Attributes[Remove]
Out[10]= {HoldAll, Protected}
```

## Implementation notes

**Algorithm.** `builtin_remove` (`src/core.c`) shares the `core_apply_symbol_action`
walker with `ClearAll`: it visits each argument — a symbol, a string, or a flat
`List` of them — and applies `core_remove_one` to every resolved name, returning
`Null`. `core_remove_one` skips any `Protected` symbol (the guard
that keeps `Remove` from ever deleting a builtin) and otherwise calls
`symtab_remove_symbol(name)`, deleting the symbol's definition from the symbol
table entirely.

This is a stronger erase than `ClearAll`: where `ClearAll` empties a symbol but
leaves the entry in place, `Remove` deletes the entry, so the name no longer
appears in the symbol table until it is next referenced (at which point a fresh,
undefined symbol is created). `Remove` is itself `Protected` so it cannot
be removed.

**Attributes & limits.** `Attributes[Remove]` is `{HoldAll, Protected}`, and
`Protected` alone is what keeps it, and every other builtin, from being removed.
Its arguments arrive unevaluated, and non-symbol/non-string specs are
ignored. A removed symbol's registered default options go with it, whereas
`Clear` and `ClearAll` keep `Options[s]`.

- `ClearAll` and `Remove` both have attributes `{HoldAll, Protected}`. Both hold
  their arguments, so they operate on the symbol, not its current value.
- Neither affects symbols with the attribute `Protected`. This is what prevents
  `Remove`/`ClearAll` from ever deleting or wiping a built-in.
- `ClearAll`, unlike `Clear`, also removes attributes and the usage message.
- Both return `Null`.

**Attributes:** `HoldAll`, `Protected`.

## References

**See also:** [ClearAll](../../assignment-and-rules/ClearAll/), [Clear](../../assignment-and-rules/Clear/)

- Source: [`src/core.c`](https://github.com/stblake/mathilda/blob/main/src/core.c)
- Specification: [`docs/spec/builtins/assignment-and-rules.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/assignment-and-rules.md)
- Tests: [`tests/test_clearall_remove_protect.c`](https://github.com/stblake/mathilda/blob/main/tests/test_clearall_remove_protect.c)

## Notes & additional examples

### Notes

- `Remove[s]` deletes the symbol `s` from the symbol table entirely. This is
  stronger than `ClearAll`, which empties a symbol but keeps its entry and its
  default options.
- After removal the name no longer exists, and the next reference to it creates
  a fresh, undefined symbol. The result is `Null`.
- Arguments may be symbols, strings, or a flat list of them.
- A `Protected` symbol is skipped, so `Remove` never deletes a built-in.
  `Attributes[Remove]` is `{HoldAll, Protected}`.
- See also: `ClearAll`, `Clear`, `Names`.
