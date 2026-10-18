### Worked examples

```mathematica
In[1]:= StringCases["x1y22z333", RegularExpression["\\d+"]]  (* a run of digits *)
```

```mathematica
In[1]:= StringReplace["John Smith", RegularExpression["(\\w+) (\\w+)"] -> "$2 $1"]  (* $n refers to a captured group *)
```

```mathematica
In[1]:= StringReplace["hello", RegularExpression["(.)\\1"] -> "<$0>"]  (* a backreference inside the regex; $0 is the whole match *)
```

```mathematica
In[1]:= StringMatchQ["abc123", RegularExpression["[a-z]+\\d+"]]  (* letters then digits, whole string *)
```

```mathematica
In[1]:= StringSplit["a1b2c3", RegularExpression["\\d"]]  (* the digits are the delimiters *)
```

```mathematica
In[1]:= StringCases["The year 1999 and 2024", RegularExpression["\\b\\d{4}\\b"]]  (* word boundaries and counted repeats *)
```

```mathematica
In[1]:= StringCases["caf\:00e9 na\:00efve", RegularExpression["\\w+"]]  (* \w is Unicode-aware *)
```

```mathematica
In[1]:= {StringCases["one\ntwo", RegularExpression["e.t"]], StringCases["one\ntwo", "e" ~~ _ ~~ "t"]}  (* a regex dot skips newlines; the blank _ does not *)
```

```mathematica
In[1]:= StringCases["ab12cd", RegularExpression["[a-z]+"] ~~ DigitCharacter]  (* a regex mixes with symbolic string patterns *)
```

```mathematica
In[1]:= Check[RegularExpression["a(b"], "invalid"]  (* a bad pattern emits RegularExpression::regex *)
```

### Notes

- `RegularExpression["re"]` is an inert head carrying PCRE2 source. Evaluating
  it compiles the pattern once, which is a syntax check that also warms the
  compiled-pattern cache the string functions use.
- In the right-hand side of a rule whose pattern contains a `RegularExpression`,
  `$0` is the whole match, `$n` the n-th group, and `$$` a literal `$`. With
  purely symbolic patterns, use named pieces (`x:p :> ...`) instead. There
  `"$1"` is literal text (v0.455).
- PCRE2 runs in UTF + UCP mode, so `\w`, `\d`, `[[:alpha:]]` and `(?i)` are
  Unicode-aware. A bare `RegularExpression` keeps PCRE's own `.`, which does not
  match a newline. The symbolic blanks `_`/`__` do match newlines, as in
  Mathematica. Use `(?s)` to change the regex behaviour.
- Backslashes are doubled in a Mathilda string literal (`"\\d"` is the regex
  `\d`). An invalid pattern emits `RegularExpression::regex` through the message
  funnel, so `Check` and `Quiet` see it, and stays inert.
- See also: `StringCases`, `StringReplace`, `StringMatchQ`, `StringSplit`,
  `StringExpression`.
