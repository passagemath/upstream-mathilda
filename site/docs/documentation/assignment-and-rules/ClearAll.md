# ClearAll

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`ClearAll[s1, s2, ...]`**

clears all values, definitions, attributes and messages for the named symbols. ClearAll\[{s1, s2, ...}\] accepts a list of specs.

<details>
<summary>Notes</summary>

ClearAll has attribute HoldAll; symbols with attribute Protected are not affected.

</details>

## Examples (12)

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

### Applications (9)

A value and an attribute

```mathematica
In[4]:= SetAttributes[g, Orderless]; g = 5
Out[4]= 5
```

```mathematica
In[5]:= ClearAll[g]
```

Both are gone

```mathematica
In[6]:= {g, Attributes[g]}
Out[6]= {g, {}}
```

DownValues and attributes together

```mathematica
In[7]:= h[x_] := x^2; SetAttributes[h, Listable]; ClearAll[h]; {DownValues[h], Attributes[h]}
Out[7]= {{}, {}}
```

Compare Clear, which keeps the attributes

```mathematica
In[8]:= w = 3; SetAttributes[w, Flat]; Clear[w]; {w, Attributes[w]}
Out[8]= {w, {Flat}}
```

A list of specs, mixing symbols and strings

```mathematica
In[9]:= p = 1; q = 2; ClearAll[{p, "q"}]; {p, q}
Out[9]= {p, q}
```

Default options survive; Remove deletes them

```mathematica
In[10]:= Options[o1] = {opt -> 1}; ClearAll[o1]; Options[o1]
Out[10]= {opt -> 1}
```

A Protected symbol is skipped

```mathematica
In[11]:= k = 1; Protect[k]; ClearAll[k]; k
Out[11]= 1
```

So built-ins are never wiped

```mathematica
In[12]:= ClearAll[Sin]; Sin[0]
Out[12]= 0
```

## Implementation notes

**Algorithm.** `builtin_clear_all` (`src/core.c`) is one of the four
symbol-management heads driven by `core_apply_symbol_action`, which walks the
argument list and applies a per-symbol action to each spec — a bare symbol, a
string naming a symbol, or a flat `List` of such specs (so `ClearAll[{a, b}]`
works), with the name read out by `core_symbol_name_of`. The action here is
`core_clear_all_one`.

For each name, `core_clear_all_one` first skips any `Protected`
symbol (which is what shields every builtin), then does the full erase that
distinguishes `ClearAll` from `Clear`: `symtab_clear_symbol` drops the
OwnValues/DownValues, the attribute word is zeroed (bumping the rule epoch so
the evaluation cache is invalidated), and the docstring (usage message) is
freed. The C builtin function pointer, if any, is left intact — the `Protected`
guard already keeps `ClearAll` away from builtins. The head returns `Null`.

**Attributes & limits.** `ClearAll` carries `HoldAll | Protected`, so the
symbols arrive unevaluated rather than being replaced by their current values.
Non-symbol and non-string specs are ignored without a message. A `Protected`
user symbol is also skipped silently, and keeps its values until it is
`Unprotect`ed. `Protected` is the only guard (since v0.301), and every protection
can be lifted with `Unprotect`.

- `ClearAll` and `Remove` both have attributes `{HoldAll, Protected}`. Both hold
  their arguments, so they operate on the symbol, not its current value.
- Neither affects symbols with the attribute `Protected`. This is what prevents
  `Remove`/`ClearAll` from ever deleting or wiping a built-in.
- `ClearAll`, unlike `Clear`, also removes attributes and the usage message.
- Both return `Null`.

**Attributes:** `HoldAll`, `Protected`.

## References

**See also:** [Remove](../../assignment-and-rules/Remove/), [Clear](../../assignment-and-rules/Clear/)

- Source: [`src/core.c`](https://github.com/stblake/mathilda/blob/main/src/core.c)
- Specification: [`docs/spec/builtins/assignment-and-rules.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/assignment-and-rules.md)
- Tests: [`tests/test_clearall_remove_protect.c`](https://github.com/stblake/mathilda/blob/main/tests/test_clearall_remove_protect.c)
- Tests: [`tests/test_condition_downvalue.c`](https://github.com/stblake/mathilda/blob/main/tests/test_condition_downvalue.c)
- Tests: [`tests/test_core.c`](https://github.com/stblake/mathilda/blob/main/tests/test_core.c)
- Tests: [`tests/test_evaluate.c`](https://github.com/stblake/mathilda/blob/main/tests/test_evaluate.c)

## Notes & additional examples

### Notes

- `ClearAll[s]` removes `s`'s OwnValues, DownValues, attributes and usage
  message, returning the symbol to its undefined state. Compare `Clear[s]`,
  which drops only the values. Registered default options are kept. `Remove`
  deletes them along with the symbol.
- Arguments may be symbols, strings naming symbols, or a flat list of them. The
  result is `Null`.
- `HoldAll | Protected`, so the symbols themselves are cleared, not their values.
- A `Protected` symbol is skipped silently, which is what keeps `ClearAll` from
  wiping a built-in. To clear one of your own protected symbols, `Unprotect` it
  first. `Protected` is the only guard, and every protection can be lifted.
- See also: `Clear`, `Remove`, `Unset`, `ClearAttributes`.
