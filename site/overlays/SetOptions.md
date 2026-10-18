### Worked examples

```mathematica
In[1]:= Options[FactorInteger]
In[2]:= SetOptions[FactorInteger, GaussianIntegers -> True]  (* returns the updated list *)
In[3]:= Options[FactorInteger]
```

```mathematica
In[1]:= SetOptions[FactorInteger, GaussianIntegers -> False]  (* restore the default *)
```

```mathematica
In[1]:= Options[mf] = {alpha -> 1, beta -> 2}; SetOptions[mf, beta -> 20]  (* the option keeps its position *)
```

```mathematica
In[1]:= mf[OptionsPattern[]] := {OptionValue[alpha], OptionValue[beta]}; mf[]  (* option-reading code sees the new default *)
```

```mathematica
In[1]:= SetOptions[mf, gamma -> 3]  (* an unknown option: SetOptions::optnf, left unevaluated *)
```

```mathematica
In[1]:= AppendTo[Options[mf], gamma -> 3]  (* adding an option goes through Options[...] instead *)
```

```mathematica
In[1]:= SetOptions[StringCases, IgnoreCase -> True]; StringCases["aAa", "a"]  (* built-ins read their defaults at call time *)
```

```mathematica
In[1]:= SetOptions[StringCases, IgnoreCase -> False]  (* restore the default *)
```

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
