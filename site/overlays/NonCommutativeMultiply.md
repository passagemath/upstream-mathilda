### Worked examples

```mathematica
In[1]:= a ** (b ** c)   (* associative: nested products are flattened *)
```

```mathematica
In[1]:= FullForm[(a ** b) ** c]   (* into one flat call *)
```

```mathematica
In[1]:= a ** b === b ** a   (* not commutative *)
```

```mathematica
In[1]:= {0 ** a, 1 ** a, 2 ** 3, NonCommutativeMultiply[a]}   (* no other simplification: zeros, ones, numbers and a single argument are all kept *)
```

```mathematica
In[1]:= FullForm[a b ** c^2]   (* binds tighter than Times and looser than Power *)
```

```mathematica
In[1]:= u ** u ** u ** v /. u ** u -> 1   (* a rule for two operands rewrites a leading pair inside a longer product *)
```

```mathematica
In[1]:= op /: op[i_] ** op[j_] := op[i + j]; op[1] ** op[2] ** op[3]   (* the product is Protected, so new rules go on the operands as UpValues; this one folds the whole product *)
```

```mathematica
In[1]:= op[1] ** op[2] ** x   (* the leading pair is merged and x is kept *)
```

```mathematica
In[1]:= x ** op[1] ** op[2]   (* limitation: a matching pair that is not at the start is not rewritten *)
```

```mathematica
In[1]:= x ** y ** z /. NonCommutativeMultiply -> Times   (* replace the head to get an ordinary product *)
```

### Notes

- A general associative, non-commutative product with no built-in rules. Attributes `{Flat, OneIdentity, Protected}`: only flattening happens automatically, and `NonCommutativeMultiply[a]` is not reduced to `a` (unlike `Times[a]`).
- `Expand`, `Simplify` and `D` do not distribute over it. `Expand[(a + b) ** c]` is unchanged. Write the algebra's rules (linearity, commutation relations) yourself.
- Differences from Mathematica: `Unprotect[NonCommutativeMultiply]` has no effect, because the Protected bit comes from the static attribute table, so DownValues on `**` are refused with `SetDelayed::wrsym`. Attach rules to the operands with `/:` or `^:=` instead. Also, because the head is not `Orderless`, a rule matches only a leading run of operands, so an interior pair is not rewritten.
- Precedence: `^` > `**` > `.` > `*`. So `a ** b^c` is `a ** (b^c)` and `a b ** c` is `a (b ** c)`.
- It is symbolic only, with no packed-array, NDArray or `Compile[]` path.
- See also: `Times`, `Dot`, `Flat`, `OneIdentity`, `TagSetDelayed`.
