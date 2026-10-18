# SetOptions

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`SetOptions[s, name -> value, ...]`**

sets default options for the symbol s and returns the new Options\[s\].  It can change Protected symbols, and only changes existing options -- an unknown name raises SetOptions::optnf.  Use AppendTo\[Options\[s\], ...\] to add one.

## Examples (16)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (6)

```mathematica
In[1]:= Options[LinearSolve]
Out[1]= {Method -> Automatic, Modulus -> 0, ZeroTest -> Automatic}

In[2]:= Options[f] = {a -> 1, b -> 2}; f[OptionsPattern[]] := {OptionValue[a], OptionValue[b]}

In[3]:= {f[], f[a -> 17], f[b -> 18], f[a -> 17, b -> 18]}
Out[3]= {{1, 2}, {17, 2}, {1, 18}, {17, 18}}

In[4]:= SetOptions[f, b -> 20]; f[]
Out[4]= {1, 20}

In[5]:= SetOptions[f, c -> 3]
Out[5]= SetOptions[f, c -> 3]

In[6]:= AppendTo[Options[f], c -> 3]
Out[6]= {a -> 1, b -> 20, c -> 3}
```

### Applications (10)

```mathematica
In[7]:= Options[FactorInteger]
Out[7]= {GaussianIntegers -> False, Method -> Automatic}
```

Returns the updated list

```mathematica
In[8]:= SetOptions[FactorInteger, GaussianIntegers -> True]
Out[8]= {GaussianIntegers -> True, Method -> Automatic}
```

```mathematica
In[9]:= Options[FactorInteger]
Out[9]= {GaussianIntegers -> True, Method -> Automatic}
```

Restore the default

```mathematica
In[10]:= SetOptions[FactorInteger, GaussianIntegers -> False]
Out[10]= {GaussianIntegers -> False, Method -> Automatic}
```

The option keeps its position

```mathematica
In[11]:= Options[mf] = {alpha -> 1, beta -> 2}; SetOptions[mf, beta -> 20]
Out[11]= {alpha -> 1, beta -> 20}
```

Option-reading code sees the new default

```mathematica
In[12]:= mf[OptionsPattern[]] := {OptionValue[alpha], OptionValue[beta]}; mf[]
Out[12]= {1, 20}
```

An unknown option: SetOptions::optnf, left unevaluated

```mathematica
In[13]:= SetOptions[mf, gamma -> 3]
Out[13]= SetOptions[mf, gamma -> 3]
```

Adding an option goes through Options[...] instead

```mathematica
In[14]:= AppendTo[Options[mf], gamma -> 3]
Out[14]= {alpha -> 1, beta -> 20, gamma -> 3}
```

Built-ins read their defaults at call time

```mathematica
In[15]:= SetOptions[StringCases, IgnoreCase -> True]; StringCases["aAa", "a"]
Out[15]= {"a", "A", "a"}
```

Restore the default

```mathematica
In[16]:= SetOptions[StringCases, IgnoreCase -> False]
Out[16]= {IgnoreCase -> False, Overlaps -> False}
```

## Algorithm

options_builtin.c — Options, SetOptions, OptionValue and the registry of default option settings for option-accepting builtins.

Mathilda stores a symbol's default options as a List[Rule[name, val], ...] on SymbolDef.default_options (the DefaultValues-equivalent), reached through symtab_set_options / symtab_get_options. This file implements:

```text
  Options[...]      query a symbol's defaults or an expression's explicit
                    options, optionally selected by name.
  SetOptions[...]   redefine individual default options of a symbol.
  OptionValue[...]  resolve a single option value from explicit options plus
                    defaults; the bare/2-arg forms are resolved inside a rule
                    by optionvalue_inject_context (see apply_down_values).
  options_register_defaults  the comprehensive table wired from core_init.
```

Memory: every result is freshly built. Sub-expressions taken from `res` or from the stored options are duplicated with expr_copy (a refcount bump), so nothing aliased is mutated in place and the evaluator remains free to release

```text
`res` after a builtin returns.
```

## Implementation notes

**Algorithm.** `builtin_setoptions` (`src/options_builtin.c`) redefines
individual default options of a symbol. The first argument must be a symbol. It
takes a working
copy of the symbol's current option rules (from `symtab_get_options`) as a flat
vector and, for each trailing `name -> value` rule, finds the existing option of
that name (context-insensitive match) and **replaces it in place**, preserving
the option's original position.

`SetOptions` cannot *add* an option. A name not already in the symbol's defaults
raises `SetOptions::optnf` through `mth_message`, leaves the stored options
unchanged, and leaves the call unevaluated. A malformed argument returns `NULL`:
a non-rule, a list of rules, or a string in place of the symbol. On success the new
list is stored with `symtab_set_options`, the evaluation cache is dropped
(`eval_clock_bump`, since option changes can alter results), and a copy of the
updated option list is returned.

**Data structures & limits.** Options live on `SymbolDef.default_options` as
`List[Rule[...]]`; the working vector holds `expr_copy`'d rules so nothing
aliased is mutated before the atomic store. `SetOptions` is `Protected`, but it
works on `Protected` targets: built-ins read their defaults from this store at
call time. For example, the regex string functions (`StringCases`,
`StringPosition`, `StringCount`, `StringSplit`, `StringReplace`, `StringMatchQ`,
`StringContainsQ`, …) seed `IgnoreCase`/`Overlaps` from `Options[head]` before
looking at explicit options. `builtin_setoptions` makes no attribute check on
the target symbol.

- `Options`, `SetOptions`, and `OptionValue` all have attribute `{Protected}`.
- Default options survive `Clear[f]` (only rules are cleared) and are removed
  with the symbol by `Remove[f]`.
- The registered defaults mirror the options each builtin's evaluator actually
  honors, with the value it falls back to when the option is absent. The sweep
  covers, among others: `Integrate` (`Method`), `Limit` (`Direction`,
  `Assumptions`), `Series`/`PowerExpand` (`Assumptions`), `D`/`Dt`
  (`NonConstants`), `Sum`/`Product` (`Method`), `Simplify`/`FullSimplify`
  (`Assumptions`, `ComplexityFunction`, `TransformationFunctions`),
  `GroebnerBasis` (`MonomialOrder -> Lexicographic`, `CoefficientDomain ->
  Rationals`, `Method`, `Sort`, `Modulus`), `Factor`/`Together`/`Cancel`/`Apart`
  and the `Polynomial*` family (`Extension`, `Modulus`), `IrreduciblePolynomialQ`
  /`SquareFreeQ` (`GaussianIntegers`, ...), the eigen/linear-algebra heads
  (`Eigenvalues`/`Eigenvectors` → `Cubics`, `Quartics`; `LinearSolve`,
  `LeastSquares`, `NullSpace`, `MatrixRank`, `PseudoInverse`,
  `SingularValueDecomposition`), `PrimeQ`/`CoprimeQ`/`FactorInteger`
  (`GaussianIntegers`), the `Random*` heads (`WorkingPrecision`), the numerical
  calculus/root-finding heads (`NIntegrate`, `NSum`, `NProduct`, `NLimit`, `ND`,
  `NSeries`, `NResidue`, `FindRoot`, `FindMinimum`, `FindMaximum`, `NRoots`,
  `NSolve`, `Solve`, `Fit`), `Plot`, and the structural family that reads
  `Heads` (`Cases`/`Count`/`DeleteCases`/`MemberQ`/`Map`/`Apply`/`MapAll`/`Level`
  /`Depth`/`LeafCount` default `Heads -> False`; `Position`/`FreeQ` default
  `Heads -> True`).

**Attributes:** `Protected`.

## References

**See also:** [Options](../../assignment-and-rules/Options/), [OptionValue](../../assignment-and-rules/OptionValue/), [Set](../../assignment-and-rules/Set/), [Hold](../../expression-information/Hold/), [Integrate](../../calculus/Integrate/), [Limit](../../calculus/Limit/), [Series](../../power-series/Series/), [PowerExpand](../../algebra/PowerExpand/)

- Source: [`src/options_builtin.c`](https://github.com/stblake/mathilda/blob/main/src/options_builtin.c)
- Specification: [`docs/spec/builtins/assignment-and-rules.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/assignment-and-rules.md)
- Tests: [`tests/test_compiledfunction.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compiledfunction.c)
- Tests: [`tests/test_image.c`](https://github.com/stblake/mathilda/blob/main/tests/test_image.c)
- Tests: [`tests/test_numberform.c`](https://github.com/stblake/mathilda/blob/main/tests/test_numberform.c)
- Tests: [`tests/test_options.c`](https://github.com/stblake/mathilda/blob/main/tests/test_options.c)

## Notes & additional examples

### Notes

- `SetOptions[s, name -> value, ...]` changes a symbol's default option settings
  and returns the full updated option list. Each changed option keeps its
  original position. Later `Options[s]` calls, `OptionValue`, and the built-ins
  that read their own defaults all see the new value.
- It works on `Protected` symbols. Built-ins such as `StringCases`,
  `StringPosition`, `StringSplit` and `StringMatchQ` seed `IgnoreCase` /
  `Overlaps` from `Options[head]`.
- It can only change options a symbol already has. An unknown name raises
  `SetOptions::optnf` and leaves both the settings and the call unchanged. Use
  `AppendTo[Options[s], ...]` or `Options[s] = {...}` to add options.
- The first argument must be a symbol, and each option a separate rule argument.
  A list of rules is not accepted.
- The evaluation cache is dropped on every change.
- See also: `Options`, `OptionValue`, `OptionsPattern`.
