# UpSet

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

lhs ^= rhs assigns rhs and associates the assignment with the symbols at level one of lhs (UpValues).

## Examples (11)

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

### Applications (8)

The value is attached to square, the argument, not to area

```mathematica
In[4]:= area[square] ^= 4; area[square]
Out[4]= 4
```

Stored as an UpValue of square; area has no DownValues

```mathematica
In[5]:= {UpValues[square], DownValues[area]}
Out[5]= {{HoldPattern[area[square]] :> 4}, {}}
```

Every level-one symbol gets the rule: a as an argument, b as the head of an argument

```mathematica
In[6]:= prop[a, b[c]] ^= val; {UpValues[a], UpValues[b]}
Out[6]= {{HoldPattern[prop[a, b[c]]] :> val}, {HoldPattern[prop[a, b[c]]] :> val}}
```

A rule for the Protected head Length, attached to vec, so Length need not be unprotected

```mathematica
In[7]:= Length[vec] ^= 3; {Length[vec], Length[{1, 2}]}
Out[7]= {3, 2}
```

Immediate: the right side was evaluated once, when the assignment was made

```mathematica
In[8]:= ctr = 0; stamp[t] ^= ctr++; {stamp[t], stamp[t], ctr}
Out[8]= {0, 0, 1}
```

An UpValue is tried before the DownValues of the enclosing head

```mathematica
In[9]:= fdv[x_] := down; fdv[gup] ^= up; {fdv[gup], fdv[other]}
Out[9]= {up, down}
```

No symbol at level one to attach to: UpSet::nosym, nothing is stored

```mathematica
In[10]:= Check[fx[x_] ^= 1, failed]
Out[10]= failed
```

Pi is Protected: UpSet::write

```mathematica
In[11]:= Check[Sin[Pi] ^= 2, failed]
Out[11]= failed
```

## Implementation notes

**Algorithm.** `lhs ^= rhs` parses (`OP_UPSET`, precedence 500, right-associative, `src/parse.c`) to `UpSet[lhs, rhs]`. `UpSet` is `HoldFirst`, so the evaluator's assignment branch (step 6b in `evaluate_step`) receives an already-evaluated `rhs`. `build_assignment_target` evaluates the *arguments* of `lhs` exactly as `Set` does (pattern-bearing and Hold-attributed positions stay held), then `apply_up_assignment` calls `collect_level_one_symbols`: a top-level `Condition` is unwrapped, and for each argument `element_dispatch_head` yields the symbol an occurrence would consult — the argument itself when it is a symbol, the head of a call `g[...]`, the head named by a constrained blank (`a_h`, `_h`, and through `PatternTest`/`Condition`/`Optional` wrappers), and nothing for a bare `_` or an `Alternatives`. The distinct symbols (up to 64) each receive the rule via `symtab_add_up_value`; a `Protected` one is skipped with `UpSet::write`, and an `lhs` with no such symbol issues `UpSet::nosym` and leaves the `UpSet` unevaluated. The call returns `rhs`.

**Data structures.** An UpValue is an ordinary `Rule` node (pattern, replacement, dispatch arity and first-argument-head keys) on the tag's `SymbolDef.up_values` list, inserted by the same `add_rule` specificity sort that DownValues use; a same-LHS rule is replaced in place. A global `symtab_up_value_count` tracks the total, and each install advances the rule epoch (`eval_rule_epoch_mark`) so memoised GROUND nodes re-check. At call time `apply_up_values` (step 3b, before the head's DownValues and builtin) collects the distinct UpValue-bearing symbols at level one of the call (≤ 32 candidates), canonically sorts the arguments first when the head is `Orderless`, and tries each candidate's list with `match`, falling back to `match_flat_leftover` when the enclosing head is `Flat`.

**Complexity / limits.** With no UpValues anywhere, every call pays one load-and-branch (`symtab_up_value_count == 0`). Otherwise the cost is a scan of the level-one arguments plus the rule lists of the candidates found, filtered by arity and first-argument head. UpValues are skipped under `HoldAllComplete` and inside the assignment primitives themselves, so an `_[g[x_]]` rule cannot hijack a held assignment LHS. For a non-`Orderless` `Flat` head the matched operands must form a leading run. The whole subsystem is compiled in under the `UP_VALUES` build flag (default on).

- `UpSet` is `HoldFirst` (its `rhs` is evaluated once), `UpSetDelayed` is
  `HoldAll` (its `rhs` is held and re-evaluated on each use); both are `Protected`.
- The head of a pattern-constrained argument also counts: `a_mod + b_mod ^:= …`
  keys on `mod` (the head of the `Blank`), so the rule fires for `mod[…] + mod[…]`.
- An UpValue on a `Protected` symbol is refused with a message.

**Attributes:** `HoldFirst`, `Protected`, `SequenceHold`.

## References

**See also:** [UpSetDelayed](../../assignment-and-rules/UpSetDelayed/), [HoldFirst](../../other-advanced/HoldFirst/), [HoldAll](../../expression-information/HoldAll/), [Blank](../../pattern-matching/Blank/)

- D. Withoff, *Mathematica Internals* (Wolfram Research technical report, 1992), §3.1 (evaluation order: upvalues before downvalues).
- Source: [`src/eval.c`](https://github.com/stblake/mathilda/blob/main/src/eval.c)
- Specification: [`docs/spec/builtins/assignment-and-rules.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/assignment-and-rules.md)

## Notes & additional examples

### Notes

- `lhs ^= rhs` evaluates `rhs` once (`UpSet` is `HoldFirst`) and returns it. The rule is attached to **every** distinct symbol at level one of `lhs`: each argument that is a symbol, the head of each argument, and the head a blank names (`a_h`). The head of `lhs` itself gets nothing.
- The arguments of `lhs` are evaluated before storing, as with `Set`; pattern-bearing arguments stay as written.
- UpValues fire before the enclosing head's DownValues and built-in code, and not inside `HoldAllComplete`. This is how to give a Protected built-in (`Length`, `Plus`, `N`) a meaning for your own objects.
- Messages: `UpSet::nosym` (no level-one symbol; the call stays unevaluated) and `UpSet::write` (a Protected tag is skipped, and any other tags still get the rule). Both go through the message system, so `Quiet` and `Check` work.
- Attributes `{HoldFirst, Protected, SequenceHold}`. It is part of the `UP_VALUES` subsystem, which is on by default.
- See also: `UpSetDelayed`, `TagSet`, `UpValues`, `TagUnset`, `Set`.
