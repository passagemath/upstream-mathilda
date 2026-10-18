# StringQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`StringQ[expr]`**

gives True if expr is a string, and False otherwise. The empty

<details>
<summary>Notes</summary>

string "" gives True. Called with any number of arguments other than one it leaves the expression unevaluated (StringQ::argx).

</details>

## Examples (17)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= StringQ["AbC"]
Out[1]= True

In[2]:= StringQ[""]
Out[2]= True

In[3]:= StringQ[123]
Out[3]= False

In[4]:= StringQ[] StringQ::argx: StringQ called with 0 arguments; 1 argument is expected.
```

### Scope (6)

```mathematica
In[5]:= MachineNumberQ[Sin[1000.]]
Out[5]= True
```

Overflows to +inf

```mathematica
In[6]:= MachineNumberQ[Exp[1000.]]
Out[6]= False
```

```mathematica
In[7]:= MachineNumberQ[-29037945.290347]
Out[7]= True
```

MPFR, not machine

```mathematica
In[8]:= MachineNumberQ[N[Pi, 30]]
Out[8]= False
```

```mathematica
In[9]:= MachineNumberQ[1.0 + 2.0 I]
Out[9]= True
```

Exact Gaussian integer

```mathematica
In[10]:= MachineNumberQ[1 + 2 I]
Out[10]= False
```

### Applications (7)

A string

```mathematica
In[11]:= StringQ["AbC"]
Out[11]= True
```

The empty string still counts

```mathematica
In[12]:= StringQ[""]
Out[12]= True
```

A number is not a string

```mathematica
In[13]:= StringQ[123]
Out[13]= False
```

Not Listable: a list is not a string

```mathematica
In[14]:= StringQ[{"a", "b"}]
Out[14]= False
```

Map it for an element-wise test

```mathematica
In[15]:= StringQ /@ {"a", 1, x, "b"}
Out[15]= {True, False, False, True}
```

As a predicate

```mathematica
In[16]:= Select[{1, "a", x, "bc"}, StringQ]
Out[16]= {"a", "bc"}
```

ToString always returns a String

```mathematica
In[17]:= StringQ[ToString[42]]
Out[17]= True
```

## Implementation notes

**Algorithm.** `builtin_stringq` (`src/core.c`) is a one-argument type test: it
returns `True` exactly when `res->data.function.args[0]->type == EXPR_STRING`, and
`False` for every other leaf or compound. The empty string `""` is still an
`EXPR_STRING`, so it gives `True`. `StringQ` is not `Listable`, so a list of
strings is tested as a single object (a `List` is not a string) and gives `False`
rather than threading.

**Data structures.** A single tag check on the argument `Expr`; nothing is
allocated on the common path beyond the `True`/`False` symbol returned.

**Complexity / limits.** `O(1)`. Any arity other than one is a malformed call
shape: the builtin routes a `StringQ::argx` diagnostic through
`builtin_arg_error` (so `Quiet[]`/`Check[]` see it) and leaves the call
unevaluated. Attributes `Protected`.

**Attributes:** `Protected`.

## References

**See also:** [AtomQ](../../expression-information/AtomQ/), [NumberQ](../../expression-information/NumberQ/), [IntegerQ](../../expression-information/IntegerQ/), [MachineNumberQ](../../expression-information/MachineNumberQ/), [Complex](../../arithmetic/Complex/), [ExactNumberQ](../../other-advanced/ExactNumberQ/)

- Source: [`src/core.c`](https://github.com/stblake/mathilda/blob/main/src/core.c)
- Specification: [`docs/spec/builtins/expression-information.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/expression-information.md)
- Tests: [`tests/test_core.c`](https://github.com/stblake/mathilda/blob/main/tests/test_core.c)
- Tests: [`tests/test_graphics.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graphics.c)
- Tests: [`tests/test_ml_classify.c`](https://github.com/stblake/mathilda/blob/main/tests/test_ml_classify.c)
- Tests: [`tests/test_pred_compile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_pred_compile.c)

## Notes & additional examples

### Notes

- `StringQ[expr]` is `True` exactly when `expr` is a string, including `""`, and
  `False` for everything else. The test is a single type-tag check.
- It is deliberately **not** `Listable`, so a list is tested as a single object.
  Map it for an element-wise answer.
- Any arity other than one emits `StringQ::argx` through the message funnel and
  leaves the call unevaluated. `Protected`.
- See also: `StringMatchQ`, `NumberQ`, `MatchQ` with `_String`.
