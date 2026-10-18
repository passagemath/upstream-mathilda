### Worked examples

```mathematica
In[1]:= StringMatchQ["hello", "h" ~~ ___]  (* the whole string must match; ___ is zero or more characters *)
```

```mathematica
In[1]:= StringMatchQ["abc123", LetterCharacter .. ~~ DigitCharacter ..]  (* letters followed by digits *)
```

```mathematica
In[1]:= StringMatchQ[{"cat", "dog"}, "c" ~~ __]  (* threads over a list of subjects *)
```

```mathematica
In[1]:= StringMatchQ["apple.txt", "*.txt"]  (* the * wildcard inside a literal string *)
```

```mathematica
In[1]:= StringMatchQ["Hello", "H@"]  (* the @ wildcard: one or more characters that are not uppercase *)
```

```mathematica
In[1]:= StringMatchQ["Hello", "hello", IgnoreCase -> True]  (* case-insensitive *)
```

```mathematica
In[1]:= {StringMatchQ["aXa", x_ ~~ "X" ~~ x_], StringMatchQ["aXb", x_ ~~ "X" ~~ x_]}  (* a repeated name must match the same text *)
```

```mathematica
In[1]:= {StringMatchQ["12", x__ /; ToExpression[x] > 10], StringMatchQ["5", x__ /; ToExpression[x] > 10]}  (* a Condition is evaluated during the match *)
```

```mathematica
In[1]:= StringMatchQ["\:00fcber", LetterCharacter..]  (* character classes are Unicode-aware *)
```

```mathematica
In[1]:= Select[{"2024", "v2", "17"}, StringMatchQ[DigitCharacter..]]  (* operator form as a predicate *)
```

### Notes

- `StringMatchQ` asks whether the *whole* string matches. The pattern is wrapped
  `\A(?:...)\z`. Use `StringContainsQ` for a substring test.
- In a literal-string pattern `*` stands for any run of characters and `@` for one
  or more characters that are not uppercase letters. Write `"\\*"` / `"\\@"` for
  the literal characters.
- The full string-pattern language is accepted: `~~`, `|`, `..`, the character
  classes, `Shortest`/`Longest`, `Except`, named patterns as backreferences,
  `Condition`, and an arbitrary `PatternTest`, which runs as a callout into the
  evaluator.
- A list of patterns is a set of alternatives. A non-string subject leaves the
  call unevaluated, while a non-string element of a subject list is passed
  through unchanged.
- Compiled patterns are cached. Mapping one pattern over 40k words dropped from
  556 ms to about 40 ms in v0.456.
- See also: `StringContainsQ`, `StringStartsQ`, `StringFreeQ`, `MatchQ`.
