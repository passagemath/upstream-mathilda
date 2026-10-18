### Worked examples

```mathematica
In[1]:= StringTrim["   hello   "]  (* default: whitespace from both ends *)
```

```mathematica
In[1]:= StringTrim["\n\t text \n"]  (* tabs and newlines are whitespace too *)
```

```mathematica
In[1]:= StringTrim["xxabxx", "x"]  (* one match is removed at each end, as in Mathematica *)
```

```mathematica
In[1]:= StringTrim["xxabxx", "x"..]  (* use a repeated pattern to strip a whole run *)
```

```mathematica
In[1]:= StringTrim["123abc456", DigitCharacter..]  (* any string pattern works *)
```

```mathematica
In[1]:= StringTrim["xxx", "xx"]  (* the two cuts never overlap *)
```

```mathematica
In[1]:= StringTrim["XXhelloxx", "x", IgnoreCase -> True]  (* case-insensitive *)
```

```mathematica
In[1]:= StringTrim[{"  a ", " b  "}]  (* threads over a list of strings *)
```

### Notes

- `StringTrim[s, p]` is `StringReplace[s, {StartOfString ~~ p -> "", p ~~ EndOfString -> ""}]`.
  It removes **one** match at each end and never touches the middle. Before
  v0.455 it stripped each end to a fixed point, so `StringTrim["xxabxx", "x"]`
  wrongly gave `"ab"`.
- The default pattern `Whitespace` is itself a run (`\s+`), so the
  one-argument form still removes all leading and trailing whitespace.
- The front match is compiled anchored and the back match end-anchored, so
  neither end needs a scan of the subject.
- A list holding a non-string emits `StringTrim::strse` and stays unevaluated.
- See also: `StringReplace`, `StringDelete`, `StringPadLeft`.
