# TagUnset

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

f /: lhs =. removes the rule for lhs associated with the symbol f.

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

An UpValue on h

```mathematica
In[4]:= h /: f[h[x_]] := fh[x]; f[h[5]]
Out[4]= fh[5]
```

Plain Unset does not reach an UpValue, so the rule still fires

```mathematica
In[5]:= f[h[x_]] =.; f[h[5]]
Out[5]= fh[5]
```

With the tag the rule is removed; the left side is compared up to renaming of pattern variables

```mathematica
In[6]:= h /: f[h[y_]] =.; f[h[5]]
Out[6]= f[h[5]]
```

Nothing left on h

```mathematica
In[7]:= UpValues[h]
Out[7]= {}
```

When the tag is the head, a DownValue is removed

```mathematica
In[8]:= tt /: tt[1] = a; tt /: tt[1] =.; tt[1]
Out[8]= tt[1]
```

A stored condition is part of the key, so this removes nothing

```mathematica
In[9]:= sq /: area[sq[s_]] := s^2 /; s > 0; sq /: area[sq[s_]] =.; Length[UpValues[sq]]
Out[9]= 1
```

Repeat the condition to remove the rule

```mathematica
In[10]:= sq /: area[sq[s_]] /; s > 0 =.; UpValues[sq]
Out[10]= {}
```

The operator form f /: lhs =. is TagUnset

```mathematica
In[11]:= FullForm[Hold[h /: f[h] =.]]
Out[11]= Hold[TagUnset[h, f[h]]]
```

## Implementation notes

**Algorithm.** `f /: lhs =.` parses (the `OP_TAGSET` branch seeing `=.`) to `TagUnset[f, lhs]`, which is `HoldAll`. The evaluator's step 6d evaluates the arguments of `lhs` with `build_assignment_target` (so the key is built the same way the original assignment built it) and calls `apply_tag_unset`. A non-symbol tag reports `TagUnset::sym`. Otherwise `classify_tag_position` picks the list: OwnValue or DownValue removal goes through `symtab_remove_matching_rule`, UpValue removal through `symtab_remove_matching_up_value`. A SubValue position or a tag not found in `lhs` removes nothing. The call always returns `Null`.

**Data structures.** Matching is by *key*, not by pattern match: both the requested `lhs` and each stored pattern are run through `pattern_canonicalize` and `pattern_alpha_normalize`, then compared with `expr_eq`. Renaming a pattern variable therefore still finds the rule (`f[h[y_]]` removes a rule stored as `f[h[x_]]`), but a stored `Condition` is part of the key, so `lhs /; test` must be repeated to remove a conditioned rule. The first equal node is unlinked and freed, `symtab_up_value_count` is decremented and the evaluation clock bumped.

**Complexity / limits.** Linear in the length of the tag's rule list, with one alpha-normalisation per stored rule. When no rule matches the removal is silent (Mathematica issues `TagUnset::norep`). Plain `lhs =.` (`Unset`) does not reach an UpValue — the tag is required, as in Mathematica. Compiled in under `UP_VALUES`.

- `TagSet`, `TagSetDelayed` and `TagUnset` are `HoldAll` and `Protected`; for the
  immediate `TagSet` the `rhs` is evaluated before the rule is stored.
- If the tag does not occur in `lhs`, a `tagnf` message is issued and nothing is
  installed. `SubValues` (`f[x][y]`, the tag as the head of the head) are not
  supported and are declined with a message.

**Attributes:** `HoldAll`, `Protected`, `SequenceHold`.

## References

**See also:** [TagSet](../../assignment-and-rules/TagSet/), [TagSetDelayed](../../assignment-and-rules/TagSetDelayed/), [HoldAll](../../expression-information/HoldAll/)

- Source: [`src/eval.c`](https://github.com/stblake/mathilda/blob/main/src/eval.c)
- Specification: [`docs/spec/builtins/assignment-and-rules.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/assignment-and-rules.md)

## Notes & additional examples

### Notes

- `f /: lhs =.` removes the own-, down- or up-value of `f` whose left side equals `lhs`, after canonicalising and renaming pattern variables, and returns `Null`.
- The arguments of `lhs` are evaluated the same way the original assignment evaluated them, so the key matches what was stored.
- When no rule matches, nothing happens and no message is given (Mathematica issues `TagUnset::norep`). A non-symbol tag gives `TagUnset::sym`.
- Attributes `{HoldAll, Protected, SequenceHold}`.
- See also: `TagSet`, `TagSetDelayed`, `Unset`, `UpValues`, `Clear`.
