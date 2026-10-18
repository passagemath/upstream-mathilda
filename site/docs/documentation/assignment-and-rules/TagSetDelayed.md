# TagSetDelayed

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

f /: lhs := rhs assigns rhs as a delayed value associated with the symbol f.

## Examples (10)

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

### Applications (7)

Modular addition, attached to mod and not to the Protected Plus

```mathematica
In[4]:= mod /: mod[a_, p_] + mod[b_, p_] := mod[Mod[a + b, p], p]; mod[2, 5] + mod[4, 5]
Out[4]= mod[1, 5]
```

Scalar multiplication through an UpValue for Times; 12 mod 5 = 2

```mathematica
In[5]:= mod /: i_Integer mod[a_, p_] := mod[Mod[i a, p], p]; 3 mod[4, 5]
Out[5]= mod[2, 5]
```

Both rules work together, and the unrelated term x is kept

```mathematica
In[6]:= mod[2, 5] + 3 mod[3, 5] + x
Out[6]= mod[1, 5] + x
```

Delayed: the right side runs on every use

```mathematica
In[7]:= cnt = 0; tk /: hit[tk] := ++cnt; {hit[tk], hit[tk]}
Out[7]= {1, 2}
```

A trailing condition guards the rule

```mathematica
In[8]:= pt /: norm[pt[x_, y_]] := Sqrt[x^2 + y^2] /; NumericQ[x]; {norm[pt[3, 4]], norm[pt[a, b]]}
Out[8]= {5, norm[pt[a, b]]}
```

The stored rule, with the condition moved onto the left side

```mathematica
In[9]:= Definition[pt]
Out[9]= pt /: norm[pt[x_, y_]] /; NumericQ[x] := Sqrt[x^2 + y^2]
```

A Protected tag refuses the assignment: TagSetDelayed::write

```mathematica
In[10]:= Check[Protect[pq]; pq /: q[pq] := 1, failed]
Out[10]= failed
```

## Implementation notes

**Algorithm.** `f /: lhs := rhs` parses through the same `OP_TAGSET` led branch as `TagSet`, which inspects the following operator (`=`, `:=` or `=.`) to build `TagSetDelayed[f, lhs, rhs]`. `TagSetDelayed` is `HoldAll` and its `rhs` is never evaluated at assignment time. Step 6c of the evaluator evaluates the non-pattern arguments of `lhs` (`build_assignment_target`) and calls `apply_tag_assignment` with `is_delayed` set: the tag is validated (`TagSetDelayed::sym`, `TagSetDelayed::write`), `classify_tag_position` places the rule as an own-, down- or up-value of the tag, and `move_rhs_condition` first lifts a trailing `body /; test` onto the pattern as `Condition[lhs, test]`, so a failing test makes the rule not apply. `TagSetDelayed::tagnf` (tag absent) and `TagSetDelayed::tagpos` (SubValue position) refuse the install. The call returns `Null`.

**Data structures.** One `Rule` on the tag's `own_values` / `down_values` / `up_values` list, inserted by `add_rule`; the replacement is the held `rhs`, re-evaluated after `replace_bindings` on each firing. For UpValues on `Orderless` heads such as `Plus`/`Times`, `apply_up_values` sorts the call's arguments canonically before matching, and on `Flat` heads `match_flat_leftover` lets a rule for two operands rewrite a pair inside a longer sum and keep the rest.

**Complexity / limits.** Same costs as `TagSet`: a single classification pass and a list insert, then the general UpValue scan at call time. For a non-`Orderless` `Flat` head (`Dot`, `**`) the matched operands must form a leading run. `Definition[f]` shows the rule in `f /: lhs := rhs` form (a lifted condition appears as `lhs /; test`). Compiled in under `UP_VALUES`.

- `TagSet`, `TagSetDelayed` and `TagUnset` are `HoldAll` and `Protected`; for the
  immediate `TagSet` the `rhs` is evaluated before the rule is stored.
- If the tag does not occur in `lhs`, a `tagnf` message is issued and nothing is
  installed. `SubValues` (`f[x][y]`, the tag as the head of the head) are not
  supported and are declined with a message.

**Attributes:** `HoldAll`, `Protected`, `SequenceHold`.

## References

**See also:** [TagSet](../../assignment-and-rules/TagSet/), [TagUnset](../../assignment-and-rules/TagUnset/), [HoldAll](../../expression-information/HoldAll/)

- D. Withoff, *Mathematica Internals* (Wolfram Research technical report, 1992), §3.1 (evaluation order: upvalues before downvalues).
- Source: [`src/eval.c`](https://github.com/stblake/mathilda/blob/main/src/eval.c)
- Specification: [`docs/spec/builtins/assignment-and-rules.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/assignment-and-rules.md)

## Notes & additional examples

### Notes

- `f /: lhs := rhs` stores `rhs` unevaluated on the tag `f`, as an own-, down- or up-value depending on where `f` occurs in `lhs`, and returns `Null`.
- For UpValues on `Orderless` heads (`Plus`, `Times`) the arguments are put in canonical order before matching, so `3 mod[4, 5]` and `mod[4, 5] 3` behave the same. On `Flat` heads a two-operand rule rewrites a matching pair and keeps the other operands.
- `lhs := body /; test` is stored as `lhs /; test`, which is how `Definition` shows it.
- Messages: `TagSetDelayed::sym`, `::write`, `::tagnf`, `::tagpos`, the same set as `TagSet`.
- Attributes `{HoldAll, Protected, SequenceHold}`.
- See also: `TagSet`, `TagUnset`, `UpSetDelayed`, `Definition`.
