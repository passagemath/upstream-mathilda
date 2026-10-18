# TagSet

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

f /: lhs = rhs assigns rhs and associates the assignment with the symbol f (as an up-, down- or own-value as appropriate).

## Examples (11)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= mod /: mod[a_, p_] + mod[b_, p_] := mod[Mod[a + b, p], p]; mod /: i_Integer mod[a_, p_] := mod[Mod[i a, p], p]; mod[2, 5] + 3 mod[3, 5] - mod[1, 5]
Out[1]= mod[0, 5]

In[2]:= h /: f[h[x_]] := fh[x]; f[h[5]]
Out[2]= fh[5]

In[3]:= h /: f[h[x_]] =.; f[h[5]]
Out[3]= f[h[5]]
```

### Applications (8)

The rule is associated with the tag x

```mathematica
In[4]:= x /: f2[x] = 1; {f2[x], f2[y]}
Out[4]= {1, f2[y]}
```

X is an argument of the left side, so the rule is an UpValue

```mathematica
In[5]:= UpValues[x]
Out[5]= {HoldPattern[f2[x]] :> 1}
```

When the tag is the head of the left side, the rule is an ordinary DownValue

```mathematica
In[6]:= tg /: tg[0] = zero; {DownValues[tg], UpValues[tg]}
Out[6]= {{HoldPattern[tg[0]] :> zero}, {}}
```

The UpValue on gg is tried before the DownValue on fdv

```mathematica
In[7]:= fdv[x_] := down; gg /: fdv[gg] = up; {fdv[gg], fdv[hh]}
Out[7]= {up, down}
```

A rule for the Protected head Length, stored on vv

```mathematica
In[8]:= vv /: Length[vv] = 42; {Length[vv], Length[{1, 2}]}
Out[8]= {42, 2}
```

Other arguments of the left side are evaluated before the rule is stored

```mathematica
In[9]:= xx = 5; yy /: gg2[xx, yy] = 9; UpValues[yy]
Out[9]= {HoldPattern[gg2[5, yy]] :> 9}
```

The tag must occur in the left side: TagSet::tagnf

```mathematica
In[10]:= Check[hh /: k[z] = 3, failed]
Out[10]= failed
```

A tag in head-of-head position would make a SubValue, which is not supported: TagSet::tagpos

```mathematica
In[11]:= Check[uu /: uu[a][b] = 1, failed]
Out[11]= failed
```

## Implementation notes

**Algorithm.** `f /: lhs = rhs` is a ternary operator (`OP_TAGSET`, precedence 500, a dedicated led branch in `src/parse.c`) producing `TagSet[f, lhs, rhs]`. `TagSet` is `HoldAll`; the evaluator's step 6c evaluates the arguments of `lhs` with `build_assignment_target` (as `Set` does), evaluates `rhs` once, and calls `apply_tag_assignment`. That rejects a non-symbol tag (`TagSet::sym`) and a `Protected` tag (`TagSet::write`), then `classify_tag_position` decides where the rule belongs: `lhs` *is* the tag → OwnValue (`symtab_add_own_value`); the tag is the head of `lhs` → DownValue (`symtab_add_down_value`); the tag is the dispatch head of some level-one argument (`element_dispatch_head`: the argument symbol, its head, or a blank's head constraint) → UpValue (`symtab_add_up_value`). A top-level `lhs /; test` is classified on its inner `lhs`. The tag as the head of the head (`f[x][y]`) is a SubValue and is declined with `TagSet::tagpos`; a tag absent from `lhs` gives `TagSet::tagnf`. In every case the call returns `rhs`, matching how `Set` returns its right side on a refused assignment.

**Data structures.** The rule lands on exactly one list of the tag's `SymbolDef` — `own_values`, `down_values` or `up_values` — through the shared `add_rule` insert (specificity-sorted, same-LHS replacement in place). Unlike `UpSet`, which attaches to *every* level-one symbol, `TagSet` attaches only to the named tag, so `a /: fz[a, b] = 1` leaves `b` untouched.

**Complexity / limits.** Classification is one pass over the level-one arguments; installation is a list insert. Firing is the general UpValue scan (`apply_up_values`, before DownValues and builtins). SubValues are not supported. Because the tag's own Protected bit is the only one checked, an UpValue can be attached to a Protected *outer* head (`v /: Length[v] = 42`) without unprotecting it — the intended use. Compiled in under `UP_VALUES`.

- `TagSet`, `TagSetDelayed` and `TagUnset` are `HoldAll` and `Protected`; for the
  immediate `TagSet` the `rhs` is evaluated before the rule is stored.
- If the tag does not occur in `lhs`, a `tagnf` message is issued and nothing is
  installed. `SubValues` (`f[x][y]`, the tag as the head of the head) are not
  supported and are declined with a message.

**Attributes:** `HoldAll`, `Protected`, `SequenceHold`.

## References

**See also:** [TagSetDelayed](../../assignment-and-rules/TagSetDelayed/), [TagUnset](../../assignment-and-rules/TagUnset/), [HoldAll](../../expression-information/HoldAll/)

- D. Withoff, *Mathematica Internals* (Wolfram Research technical report, 1992), §3.1 (evaluation order: upvalues before downvalues).
- Source: [`src/eval.c`](https://github.com/stblake/mathilda/blob/main/src/eval.c)
- Specification: [`docs/spec/builtins/assignment-and-rules.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/assignment-and-rules.md)

## Notes & additional examples

### Notes

- `f /: lhs = rhs` evaluates `rhs` once, stores the rule on the **single** tag `f`, and returns `rhs`. Where it goes depends on where `f` occurs: as `lhs` itself it is an OwnValue, as the head of `lhs` a DownValue, as an argument or the head of an argument an UpValue.
- `UpSet` attaches a rule to every level-one symbol; `TagSet` attaches it to the named tag only. With `a /: fz[a, b] = 1`, `b` gets nothing.
- Messages: `TagSet::sym` (the tag is not a symbol), `TagSet::write` (the tag is Protected), `TagSet::tagnf` (the tag is not in `lhs`), `TagSet::tagpos` (SubValues are not supported). In each case nothing is stored and `rhs` is still returned.
- Only the tag's Protected bit is checked. Attaching a rule to a Protected *outer* head through an unprotected tag is the intended use.
- Attributes `{HoldAll, Protected, SequenceHold}`.
- See also: `TagSetDelayed`, `TagUnset`, `UpSet`, `UpValues`.
