# UpValues

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`UpValues[s]`**

gives a list of up-value rules for s (values for g\[..,s,..\] or g\[..,s\[..\],..\]).

## Examples (12)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= g /: f[g[x_]] := h[x]; UpValues[g]
Out[1]= {HoldPattern[f[g[x_]]] :> h[x]}
```

Set the whole list

```mathematica
In[2]:= UpValues[k] = {HoldPattern[w[k]] :> done}; w[k]
Out[2]= done
```

```mathematica
In[3]:= UpValues[neverdefined]
Out[3]= {}
```

The moved condition

```mathematica
In[4]:= area[sq[s_]] ^:= s^2 /; s > 0; UpValues[sq]
Out[4]= {HoldPattern[area[sq[s_]] /; s > 0] :> s^2}
```

### Applications (8)

One rule, stored on g

```mathematica
In[5]:= g /: f[g[x_]] := h[x]; UpValues[g]
Out[5]= {HoldPattern[f[g[x_]]] :> h[x]}
```

An empty list for a symbol with no UpValues

```mathematica
In[6]:= UpValues[neverdefined]
Out[6]= {}
```

A string naming the symbol also works

```mathematica
In[7]:= UpValues["g"]
Out[7]= {HoldPattern[f[g[x_]]] :> h[x]}
```

Assigning a list replaces all UpValues of k

```mathematica
In[8]:= UpValues[k] = {HoldPattern[w[k]] :> done}; w[k]
Out[8]= done
```

And an empty list removes them

```mathematica
In[9]:= UpValues[k] = {}; w[k]
Out[9]= w[k]
```

A moved condition shows up inside the HoldPattern

```mathematica
In[10]:= area[sq[s_]] ^:= s^2 /; s > 0; UpValues[sq]
Out[10]= {HoldPattern[area[sq[s_]] /; s > 0] :> s^2}
```

Clear removes UpValues as well

```mathematica
In[11]:= Clear[g]; UpValues[g]
Out[11]= {}
```

The argument is held, so a symbol that has a value is not replaced by it

```mathematica
In[12]:= Attributes[UpValues]
Out[12]= {HoldAll, Protected}
```

## Implementation notes

**Algorithm.** `builtin_up_values` (attributes `HoldAll`, `Protected`) accepts one argument. A symbol reads its list with `symtab_get_up_values`; a string is resolved with the non-creating `symtab_lookup`, so `UpValues /@ Names["s*"]` works, and a string naming no existing symbol issues `UpValues::sym` and stays unevaluated. Any other argument also stays unevaluated. The list is converted by the shared `rules_to_list`, which emits each stored rule as `HoldPattern[lhs] :> rhs`. The assignment `UpValues[s] = {rules}` is handled in the evaluator's `Set` path (`src/eval.c`): after a `Protected` check (`Set::write`), `symtab_set_up_values` frees the old list and installs each `Rule`/`RuleDelayed` element, unwrapping an outer `HoldPattern` and skipping non-rule elements.

**Data structures.** The source is `SymbolDef.up_values`, a singly linked `Rule` list maintained by `add_rule` in specificity order (insertion order breaks ties). Each node is deep-copied into the result, so the returned rules are inert and editing them does not touch the store. A `Condition` lifted from a delayed right side appears inside the `HoldPattern`.

**Complexity / limits.** Linear in the number of stored rules. Immediate (`^=`) and delayed (`^:=`) rules are both returned with `:>` because the rule store does not record which kind created them. `Clear[s]` and `ClearAll[s]` also empty this list. Compiled in under `UP_VALUES`.

**Attributes:** `HoldAll`, `Protected`.

## References

**See also:** [UpSet](../../assignment-and-rules/UpSet/), [UpSetDelayed](../../assignment-and-rules/UpSetDelayed/), [TagSet](../../assignment-and-rules/TagSet/), [TagSetDelayed](../../assignment-and-rules/TagSetDelayed/), [DownValues](../../assignment-and-rules/DownValues/), [Rule](../../assignment-and-rules/Rule/), [RuleDelayed](../../assignment-and-rules/RuleDelayed/)

- Source: [`src/core.c`](https://github.com/stblake/mathilda/blob/main/src/core.c)
- Specification: [`docs/spec/builtins/expression-information.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/expression-information.md)
- Tests: [`tests/test_upvalues.c`](https://github.com/stblake/mathilda/blob/main/tests/test_upvalues.c)

## Notes & additional examples

### Notes

- Returns `{HoldPattern[lhs] :> rhs, ...}` for the rules made by `^=`, `^:=`, `f /: ... =` and `f /: ... :=`, in matching order, and `{}` when there are none. Immediate rules are shown with `:>` too.
- Accepts a symbol, or a string naming an existing symbol (so `UpValues /@ Names["s*"]` works). A string naming no symbol gives `UpValues::sym` and stays unevaluated.
- `UpValues[s] = {rules}` replaces the whole list; `HoldPattern` wrappers on the left sides are optional. A Protected `s` is refused with `Set::write`.
- Attributes `{HoldAll, Protected}`.
- See also: `DownValues`, `OwnValues`, `Definition`, `UpSet`, `TagSet`.
