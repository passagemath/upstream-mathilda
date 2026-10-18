# NonCommutativeMultiply

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

a \*\* b \*\* c or NonCommutativeMultiply\[a, b, c\] is a general associative, but non-commutative, form of multiplication. Instances are automatically flattened (it has attribute Flat), but no other simplification is performed: a\*\*b is not reordered to b\*\*a, 0\*\*a and 1\*\*a stay as given, and Expand / Simplify do not operate on it. Use it as a generalization of ordinary multiplication for special mathematical objects (e.g. composition in an algebra of differential operators).

## Examples (14)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= a ** b ** c
Out[1]= a ** b ** c

In[2]:= a ** (b ** c) == (a ** b) ** c
Out[2]= True

In[3]:= {a*b == b*a, a**b == b**a}
Out[3]= {True, a ** b == b ** a}

In[4]:= {0 ** a, 1 ** a}
Out[4]= {0 ** a, 1 ** a}
```

### Applications (10)

Associative: nested products are flattened

```mathematica
In[5]:= a ** (b ** c)
Out[5]= a ** b ** c
```

Into one flat call

```mathematica
In[6]:= FullForm[(a ** b) ** c]
Out[6]= NonCommutativeMultiply[a, b, c]
```

Not commutative

```mathematica
In[7]:= a ** b === b ** a
Out[7]= False
```

No other simplification: zeros, ones, numbers and a single argument are all kept

```mathematica
In[8]:= {0 ** a, 1 ** a, 2 ** 3, NonCommutativeMultiply[a]}
Out[8]= {0 ** a, 1 ** a, 2 ** 3, NonCommutativeMultiply[a]}
```

Binds tighter than Times and looser than Power

```mathematica
In[9]:= FullForm[a b ** c^2]
Out[9]= Times[a, NonCommutativeMultiply[b, Power[c, 2]]]
```

A rule for two operands rewrites a leading pair inside a longer product

```mathematica
In[10]:= u ** u ** u ** v /. u ** u -> 1
Out[10]= 1 ** u ** v
```

The product is Protected, so new rules go on the operands as UpValues; this one folds the whole product

```mathematica
In[11]:= op /: op[i_] ** op[j_] := op[i + j]; op[1] ** op[2] ** op[3]
Out[11]= op[6]
```

The leading pair is merged and x is kept

```mathematica
In[12]:= op[1] ** op[2] ** x
Out[12]= op[3] ** x
```

Limitation: a matching pair that is not at the start is not rewritten

```mathematica
In[13]:= x ** op[1] ** op[2]
Out[13]= x ** op[1] ** op[2]
```

Replace the head to get an ordinary product

```mathematica
In[14]:= x ** y ** z /. NonCommutativeMultiply -> Times
Out[14]= x y z
```

## Implementation notes

**Semantics.** `NonCommutativeMultiply` is a purely structural head with no C builtin. Its behaviour comes from its attributes, `Flat`, `OneIdentity` and `Protected`, set in the static table in `src/attr.c`. `Flat` makes the evaluator splice nested `**` calls into one flat call. No other rewriting is done: arguments are not reordered (it is not `Orderless`), numbers are not combined, and identities such as `0 ** a` or `1 ** a` are not applied. A one-argument `NonCommutativeMultiply[a]` stays as it is, because `OneIdentity` only affects pattern matching.

**Parsing and printing.** `src/parse.c` lexes `**` as `OP_NONCOMMULT` at precedence 5900. That is tighter than `.` (`Dot`) and `*` (`Times`) and looser than `^`, so `a b ** c^2` is `Times[a, NonCommutativeMultiply[b, Power[c, 2]]]`. Runs of `**` are flattened when the expression is parsed. `print_standard` prints two or more arguments as infix `a ** b ** c`, and other arities in `FullForm` style.

**Rules and limits.** Since the head is `Flat`, `match_flat_leftover` (`src/match.c`) lets a rule written for two operands rewrite part of a longer product and keep the rest. The head is not `Orderless`, so the matched operands must form a leading run (`op[1] ** op[2] ** x` is rewritten, `x ** op[1] ** op[2]` is not). Its `Protected` bit comes from the static attribute table, which `Unprotect` does not override, so new behaviour is attached through UpValues on the operands (`op /: op[i_] ** op[j_] := ...`) rather than DownValues on `**`. It has no packed-array, NDArray or `Compile[]` path; `**` is symbolic only.

- Attributes `Flat`, `OneIdentity`, `Protected`.
- Instances are automatically flattened (`a ** (b ** c)` and `(a ** b) ** c`
  both become `NonCommutativeMultiply[a, b, c]`), **but no other simplification
  is performed**: arguments are not reordered (`a ** b` stays distinct from
  `b ** a`), `0 ** a` and `1 ** a` are left as given, and `Expand` / `Simplify`
  / `FullSimplify` do not operate on expressions containing it.
- `NonCommutativeMultiply[a]` (one argument) stays unevaluated — unlike
  `Times[a]`, which reduces to `a`.
- Precedence binds tighter than `*` (`Times`) and `.` (`Dot`), and looser than
  `^` (`Power`): `a ** b ^ c` is `a ** (b^c)`, `a b ** c` is `a (b ** c)`.
- No new numeric fast path: `**` is a purely symbolic/structural head, so it has
  no packed/NDArray or `Compile[]` lowering.

**Attributes:** `Flat`, `OneIdentity`, `Protected`.

## References

**See also:** [Flat](../../expression-information/Flat/), [OneIdentity](../../expression-information/OneIdentity/), [Expand](../../algebra/Expand/), [Simplify](../../simplification/Simplify/), [Times](../../arithmetic/Times/), [Dot](../../linear-algebra/Dot/), [Power](../../arithmetic/Power/)

- Source: [`src/parse.c`](https://github.com/stblake/mathilda/blob/main/src/parse.c)
- Specification: [`docs/spec/builtins/arithmetic.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/arithmetic.md)
- Tests: [`tests/test_noncommutativemultiply.c`](https://github.com/stblake/mathilda/blob/main/tests/test_noncommutativemultiply.c)

## Notes & additional examples

### Notes

- A general associative, non-commutative product with no built-in rules. Attributes `{Flat, OneIdentity, Protected}`: only flattening happens automatically, and `NonCommutativeMultiply[a]` is not reduced to `a` (unlike `Times[a]`).
- `Expand`, `Simplify` and `D` do not distribute over it. `Expand[(a + b) ** c]` is unchanged. Write the algebra's rules (linearity, commutation relations) yourself.
- Differences from Mathematica: `Unprotect[NonCommutativeMultiply]` has no effect, because the Protected bit comes from the static attribute table, so DownValues on `**` are refused with `SetDelayed::wrsym`. Attach rules to the operands with `/:` or `^:=` instead. Also, because the head is not `Orderless`, a rule matches only a leading run of operands, so an interior pair is not rewritten.
- Precedence: `^` > `**` > `.` > `*`. So `a ** b^c` is `a ** (b^c)` and `a b ** c` is `a (b ** c)`.
- It is symbolic only, with no packed-array, NDArray or `Compile[]` path.
- See also: `Times`, `Dot`, `Flat`, `OneIdentity`, `TagSetDelayed`.
