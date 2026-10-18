# Protect

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Protect[s1, s2, ...]`**

sets the attribute Protected for the named symbols and returns the list of their names. Protect\[{s1, s2, ...}\] accepts a list of specs.

<details>
<summary>Notes</summary>

Protect has attribute HoldAll.

</details>

## Examples (12)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (5)

```mathematica
In[1]:= f[x_] := x^2; Protect[f]
Out[1]= {"f"}

In[2]:= f[x_, y_] := x + y

In[3]:= Unprotect[f]
Out[3]= {"f"}

In[4]:= Protect[f]; Protect[f]
Out[4]= {}

In[5]:= Check[f[x_] := 0, "refused"]
Out[5]= "refused"
```

### Applications (7)

Returns the names whose Protected attribute was newly set

```mathematica
In[6]:= Protect[foo, bar]
Out[6]= {"foo", "bar"}
```

```mathematica
In[7]:= MemberQ[Attributes[foo], Protected]
Out[7]= True
```

Already protected: nothing changed

```mathematica
In[8]:= Protect[foo]
Out[8]= {}
```

Assignment to a protected symbol emits Set::wrsym

```mathematica
In[9]:= Check[foo = 6, "refused"]
Out[9]= "refused"
```

The value was not assigned

```mathematica
In[10]:= foo
Out[10]= foo
```

A list of specs; strings name symbols

```mathematica
In[11]:= Protect[{"k1", k2}]
Out[11]= {"k1", "k2"}
```

```mathematica
In[12]:= Attributes[Protect]
Out[12]= {HoldAll, Protected}
```

## Implementation notes

**Algorithm.** `builtin_protect` (`src/core.c`) calls the shared driver
`core_protect_unprotect(res, protecting = true)`. That driver walks the argument
list — each spec being a symbol, a string, or a flat `List` of them — and applies
`core_protect_one` to every name. `core_protect_one` does nothing if the symbol
is already `Protected`; otherwise it
sets the `ATTR_PROTECTED` bit and bumps the rule epoch (invalidating the
evaluation cache), returning `true` only when the bit was *newly* set.

The driver collects the names whose state actually changed into a growable
`Expr**` buffer and returns them as a `List` of strings — matching the Wolfram
convention where `Protect` reports exactly what it altered, so re-protecting an
already-protected symbol yields `{}`. The `Protected` attribute is what the
evaluator and `Set` consult to refuse redefinition of a symbol.

**Attributes & limits.** `Protect` carries `HoldAll | Protected`, so its
arguments reach the handler as unevaluated symbol names. Assigning to a protected
symbol emits `Set::wrsym` / `SetDelayed::wrsym` through the message funnel, so
`Check` sees the refusal. `Protected` is the only sealing attribute (since
v0.301), and `Unprotect` always lifts it.

- Both have attributes `{HoldAll, Protected}` and hold their arguments.
- The typical sequence for adding rules to an existing symbol is
  `Unprotect[f]; definition; Protect[f]`.

**Attributes:** `HoldAll`, `Protected`.

## References

**See also:** [Unprotect](../../assignment-and-rules/Unprotect/)

- Source: [`src/core.c`](https://github.com/stblake/mathilda/blob/main/src/core.c)
- Specification: [`docs/spec/builtins/assignment-and-rules.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/assignment-and-rules.md)
- Tests: [`tests/test_clearall_remove_protect.c`](https://github.com/stblake/mathilda/blob/main/tests/test_clearall_remove_protect.c)

## Notes & additional examples

### Notes

- `Protect[s1, s2, ...]` sets the `Protected` attribute on each symbol. That is
  what `Set`, `SetDelayed`, `Unset`, `ClearAll`, `Remove` and `SetAttributes`
  consult before changing a symbol.
- It returns the names, as strings, whose state actually changed, so protecting
  an already-protected symbol gives `{}`.
- Arguments may be symbols, strings, or a flat list of them. `HoldAll | Protected`.
- Every built-in carries `Protected`, and `Unprotect` always lifts it. No
  attribute seals a symbol permanently (since v0.301).
- See also: `Unprotect`, `Attributes`, `SetAttributes`.
