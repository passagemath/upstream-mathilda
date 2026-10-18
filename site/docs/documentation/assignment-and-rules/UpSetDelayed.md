# UpSetDelayed

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

lhs ^:= rhs assigns rhs as a delayed value and associates it with the symbols at level one of lhs (UpValues).

## Examples (10)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= area[square] ^= s^2; UpValues[square]
Out[1]= {HoldPattern[area[square]] :> s^2}

In[2]:= prop[a, b[c]] ^= value; {UpValues[a], UpValues[b]}
Out[2]= {{HoldPattern[prop[a, b[c]]] :> value}, {HoldPattern[prop[a, b[c]]] :> value}}

In[3]:= f[g[x_]] ^:= h[x]; {f[g[2]], f[k[2]]}
Out[3]= {h[2], f[k[2]]}
```

### Applications (7)

The rule lives on sq and fires on area[sq[...]]

```mathematica
In[4]:= area[sq[s_]] ^:= s^2; area[sq[3]]
Out[4]= 9
```

Stored as a delayed HoldPattern rule

```mathematica
In[5]:= UpValues[sq]
Out[5]= {HoldPattern[area[sq[s_]]] :> s^2}
```

Only calls with an sq argument are touched

```mathematica
In[6]:= perim[sq[s_]] ^:= 4 s; {perim[sq[a]], perim[circle[a]]}
Out[6]= {4 a, perim[circle[a]]}
```

Delayed: the right side is re-evaluated each time the rule fires

```mathematica
In[7]:= ctr = 0; tick[t] ^:= ++ctr; {tick[t], tick[t], tick[t]}
Out[7]= {1, 2, 3}
```

A trailing condition is moved onto the pattern, so a failed test leaves the call unevaluated

```mathematica
In[8]:= vol[cube[s_]] ^:= s^3 /; s > 0; {vol[cube[2]], vol[cube[-2]]}
Out[8]= {8, vol[cube[-2]]}
```

On the Flat, Orderless head Plus, a two-term rule merges pairs inside a longer sum and keeps the rest

```mathematica
In[9]:= md[a_] + md[b_] ^:= md[a + b]; md[1] + md[2] + md[3] + z
Out[9]= md[6] + z
```

A bare blank names no symbol: UpSetDelayed::nosym

```mathematica
In[10]:= Check[Length[x_] ^:= 0, failed]
Out[10]= failed
```

## Implementation notes

**Algorithm.** `lhs ^:= rhs` parses (`OP_UPSETDELAYED`, precedence 500; the three-character `^:=` token is lexed before bare `^`) to `UpSetDelayed[lhs, rhs]`. `UpSetDelayed` is `HoldAll`, so `rhs` reaches the evaluator's step 6b unevaluated and is stored as written. The path is shared with `UpSet`: `build_assignment_target` evaluates the non-pattern arguments of `lhs`, `collect_level_one_symbols` / `element_dispatch_head` find the distinct level-one tags (argument symbols, argument heads, and the head constraint of a blank such as `a_md`), and `apply_up_assignment` installs the rule on each with `symtab_add_up_value`. One delayed-only step runs first: `move_rhs_condition` rewrites `lhs ^:= body /; test` to the pattern `Condition[lhs, test]` with replacement `body`, exactly as `SetDelayed` does, so a failing test makes the rule not match rather than return a `Condition`. The call returns `Null`. Misuse reports `UpSetDelayed::nosym` (no level-one symbol; the call stays unevaluated) or `UpSetDelayed::write` (a Protected tag, skipped).

**Data structures.** Identical to `UpSet`: a `Rule` on each tag's `SymbolDef.up_values`, ordered by `add_rule`'s specificity sort, counted in `symtab_up_value_count`, with an epoch bump on install. Because the stored replacement is the held `rhs`, `replace_bindings` substitutes the match and the evaluator re-evaluates it on every firing — the difference from `UpSet` is purely when `rhs` is evaluated.

**Complexity / limits.** Installation is linear in the number of level-one arguments plus the list insertion. Firing costs the same scan as any UpValue (see `apply_up_values`): one branch when none exist, otherwise the level-one candidates' rule lists. Mathilda's `Rule` carries no immediate/delayed flag, so `UpValues[]` and `Definition[]` show both kinds with `:>` / `:=`. Compiled in under `UP_VALUES`.

- `UpSet` is `HoldFirst` (its `rhs` is evaluated once), `UpSetDelayed` is
  `HoldAll` (its `rhs` is held and re-evaluated on each use); both are `Protected`.
- The head of a pattern-constrained argument also counts: `a_mod + b_mod ^:= …`
  keys on `mod` (the head of the `Blank`), so the rule fires for `mod[…] + mod[…]`.
- An UpValue on a `Protected` symbol is refused with a message.

**Attributes:** `HoldAll`, `Protected`, `SequenceHold`.

## References

**See also:** [UpSet](../../assignment-and-rules/UpSet/), [HoldFirst](../../other-advanced/HoldFirst/), [HoldAll](../../expression-information/HoldAll/), [Blank](../../pattern-matching/Blank/)

- D. Withoff, *Mathematica Internals* (Wolfram Research technical report, 1992), §3.1 (evaluation order: upvalues before downvalues).
- Source: [`src/eval.c`](https://github.com/stblake/mathilda/blob/main/src/eval.c)
- Specification: [`docs/spec/builtins/assignment-and-rules.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/assignment-and-rules.md)

## Notes & additional examples

### Notes

- `lhs ^:= rhs` keeps `rhs` unevaluated (`UpSetDelayed` is `HoldAll`) and returns `Null`. It attaches the rule to the same level-one symbols as `UpSet`.
- `lhs ^:= body /; test` is stored as the pattern `lhs /; test`, as with `SetDelayed`.
- `UpValues` and `Definition` show immediate and delayed rules the same way (`:>` / `:=`), because the rule store does not record which kind made them.
- Limitation: on a `Flat` head that is not `Orderless` (`Dot`, `**`), the matched operands must form a leading run.
- Attributes `{HoldAll, Protected, SequenceHold}`.
- See also: `UpSet`, `TagSetDelayed`, `UpValues`, `SetDelayed`.
