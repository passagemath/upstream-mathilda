# Unset

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Unset[lhs]`**

or lhs =. removes any rule whose left-hand side is lhs, up to renaming of pattern variables. A bare symbol clears its value; a function form clears the matching definition on the head symbol.

<details>
<summary>Notes</summary>

Unset has attribute HoldFirst; Protected symbols are not affected.

</details>

## Examples (12)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= x = 5; x =.; x
Out[1]= x

In[2]:= f[x_] := x^2; f[x_] =.; f[3]
Out[2]= f[3]

In[3]:= fact[1] = 1; fact[n_] := n fact[n - 1]; fact[1] =.; DownValues[fact]
Out[3]= {HoldPattern[fact[n_]] :> n fact[n - 1]}
```

### Applications (9)

```mathematica
In[4]:= v = 10
Out[4]= 10
```

The postfix form of Unset

```mathematica
In[5]:= v =.
```

```mathematica
In[6]:= v
Out[6]= v
```

Removes only the g[1] rule

```mathematica
In[7]:= g[1] = a; g[2] = b; g[1] =.; DownValues[g]
Out[7]= {HoldPattern[g[2]] :> b}
```

Pattern names may differ: the rule is found up to renaming

```mathematica
In[8]:= gg[x_] := x + 1; gg[y_] =.; gg[2]
Out[8]= gg[2]
```

The general rule survives

```mathematica
In[9]:= fact[1] = 1; fact[n_] := n fact[n - 1]; fact[1] =.; DownValues[fact]
Out[9]= {HoldPattern[fact[n_]] :> n fact[n - 1]}
```

A conditional rule

```mathematica
In[10]:= mm[x_] /; x > 0 := 1; mm[x_] /; x > 0 =.; DownValues[mm]
Out[10]= {}
```

A Protected symbol: Unset::wrsym

```mathematica
In[11]:= Check[Unset[Pi], "refused"]
Out[11]= "refused"
```

No matching rule: still Null

```mathematica
In[12]:= undefinedthing =.
```

## Implementation notes

**Algorithm.** `builtin_unset` (`src/core.c`) implements `lhs =.`, removing the
single definition whose left-hand side is `lhs` rather than every rule on a
symbol. It first works out which symbol owns the rule and whether it is an
OwnValue or a DownValue: a bare symbol is an OwnValue on itself; `f[...]` is a
DownValue keyed by the head `f`; and because `f[x_] /; cond =.` parses to
`Unset[Condition[f[x_], cond]]`, a `Condition` wrapper is unwrapped to find the
inner head. A non-assignable left-hand side (e.g. `Unset[5]`) returns `NULL`.

A `Protected` owner is refused with `Unset::wrsym` (mirroring `Set`),
returning `Null`. Otherwise `symtab_remove_matching_rule(name, lhs, own_value)`
deletes exactly one rule, leaving the symbol's other definitions, attributes and
remaining rules intact. To find it, the function canonicalizes the query the way
`add_rule` canonicalized the stored patterns (`pattern_canonicalize`), then
alpha-normalizes both sides (`pattern_alpha_normalize`). So `f[y_] =.` removes
a rule defined as `f[x_] := ...`. A removal bumps the evaluation clock so that
cached results which relied on the rule are dropped. The head returns `Null`
whether or not a rule was found.

**Attributes & limits.** `Unset` carries `HoldFirst | Protected`, so the target
is not evaluated to its current value before the rule is located. It takes one
argument.

- `=.` is a low-precedence postfix operator (precedence 40, like `Set`), so it
  captures the whole preceding expression: `a b =.` parses as `Unset[a b]`. The
  guard against a trailing digit keeps `k =.5` parsing as `Set[k, 0.5]`.
- `Unset` has attributes `{HoldFirst, Protected}`; it holds `lhs`, so the symbol
  (not its value) is operated on. `Protected` symbols are not affected.
- Always returns `Null`, whether or not a matching rule was found.

**Attributes:** `HoldFirst`, `Protected`.

## References

**See also:** [Set](../../assignment-and-rules/Set/)

- Source: [`src/core.c`](https://github.com/stblake/mathilda/blob/main/src/core.c)
- Specification: [`docs/spec/builtins/assignment-and-rules.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/assignment-and-rules.md)
- Tests: [`tests/test_unset.c`](https://github.com/stblake/mathilda/blob/main/tests/test_unset.c)

## Notes & additional examples

### Notes

- `lhs =.` (`Unset`) removes the one definition whose left-hand side is `lhs`,
  rather than every rule on a symbol. A bare symbol loses its OwnValue, and
  `f[args] =.` loses exactly that DownValue.
- Left-hand sides are compared after canonicalization and alpha-renaming of
  pattern variables, so `f[y_] =.` removes `f[x_] := ...`. A
  `Condition`-wrapped left-hand side is matched on its inner head.
- `HoldFirst | Protected`: the target is not evaluated before the rule is
  located. The result is `Null` whether or not a rule was found.
- A `Protected` owner is refused with `Unset::wrsym`, routed through the message
  funnel so `Check` sees it.
- See also: `Clear`, `ClearAll`, `Set`, `TagUnset`, `DownValues`.
