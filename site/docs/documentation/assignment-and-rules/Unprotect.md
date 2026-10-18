# Unprotect

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Unprotect[s1, s2, ...]`**

removes the attribute Protected from the named symbols and returns the list of their names. Unprotect\[{s1, ...}\] accepts a list of specs.

<details>
<summary>Notes</summary>

Unprotect has attribute HoldAll.

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

```mathematica
In[6]:= Protect[baz]
Out[6]= {"baz"}
```

Returns the names whose Protected attribute was actually cleared

```mathematica
In[7]:= Unprotect[baz]
Out[7]= {"baz"}
```

```mathematica
In[8]:= MemberQ[Attributes[baz], Protected]
Out[8]= False
```

Not protected any more: nothing changed

```mathematica
In[9]:= Unprotect[baz]
Out[9]= {}
```

Assignment works again

```mathematica
In[10]:= baz = 6; baz
Out[10]= 6
```

Built-ins can be unprotected too

```mathematica
In[11]:= Unprotect[Sin]
Out[11]= {"Sin"}
```

And re-protected after extending them

```mathematica
In[12]:= Protect[Sin]
Out[12]= {"Sin"}
```

## Implementation notes

**Algorithm.** `builtin_unprotect` (`src/core.c`) is the mirror of `Protect`: it
calls the shared driver `core_protect_unprotect(res, protecting = false)`, which
walks the argument list (symbols, strings, or a flat `List` of them) and applies
`core_unprotect_one` to each name. `core_unprotect_one` does nothing if the
symbol is not currently `Protected`; otherwise it
clears the `ATTR_PROTECTED` bit and bumps the rule epoch, returning `true` only
when the bit was actually cleared.

The driver returns a `List` of the names (as strings) whose protection state
changed, so `Unprotect` of a symbol that was never protected yields `{}`. Once
the bit is cleared, `Set`/`SetDelayed` and the clearing heads will again accept
the symbol, which is the usual prelude to redefining or extending a built-in's
behaviour.

**Attributes & limits.** `Unprotect` carries `HoldAll | Protected`, so its
arguments arrive as unevaluated names. It works on built-ins too
(`Unprotect[Sin]` gives `{"Sin"}`). Since v0.301 `Protected` is the only
sealing attribute, so nothing is beyond `Unprotect`'s reach.

- Both have attributes `{HoldAll, Protected}` and hold their arguments.
- The typical sequence for adding rules to an existing symbol is
  `Unprotect[f]; definition; Protect[f]`.

**Attributes:** `HoldAll`, `Protected`.

## References

**See also:** [Protect](../../assignment-and-rules/Protect/)

- Source: [`src/core.c`](https://github.com/stblake/mathilda/blob/main/src/core.c)
- Specification: [`docs/spec/builtins/assignment-and-rules.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/assignment-and-rules.md)
- Tests: [`tests/test_clearall_remove_protect.c`](https://github.com/stblake/mathilda/blob/main/tests/test_clearall_remove_protect.c)

## Notes & additional examples

### Notes

- `Unprotect[s1, s2, ...]` is the inverse of `Protect`. It clears the
  `Protected` attribute and returns the names, as strings, whose state changed.
  Unprotecting a symbol that was not protected gives `{}`.
- Once a symbol is unprotected, `Set`/`SetDelayed` and the clearing heads
  accept it again. The usual way to extend a built-in is
  `Unprotect[f]; new rule; Protect[f]`.
- Nothing is sealed beyond `Unprotect`'s reach: `Protected` is the only
  sealing attribute (since v0.301).
- Arguments may be symbols, strings, or a flat list. `HoldAll | Protected`.
- See also: `Protect`, `Attributes`, `ClearAttributes`.
