# StringRiffle

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`StringRiffle[{s1, s2, ...}]`**

Joins the si into a string with spaces between them; nested lists use spaces at the lowest level and increasing numbers of newlines at higher levels. Non-string elements are converted with ToString.

**`StringRiffle[list, sep]`**

Inserts the string sep between the top-level elements.

**`StringRiffle[list, {"left", "sep", "right"}]`**

Joins with sep and wraps the result in the left/right delimiters.

**`StringRiffle[list, sep1, sep2, ...]`**

Inserts separator sep\_i (a string or {left, sep, right}) between elements at level i.

## Examples (12)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= StringRiffle[{"a", "b", "c", "d", "e"}]
Out[1]= "a b c d e"

In[2]:= StringRiffle[{"a", "b", "c", "d", "e"}, ", "]
Out[2]= "a, b, c, d, e"

In[3]:= StringRiffle[{"a", "b", "c", "d", "e"}, {"(", " ", ")"}]
Out[3]= "(a b c d e)"

In[4]:= StringRiffle[{{"a", 27}, {"b", 28}, {"c", 29}}, {"{", ", ", "}"}, ": "]
Out[4]= "{a: 27, b: 28, c: 29}"
```

### Applications (8)

Default scheme: a single space

```mathematica
In[5]:= StringRiffle[{"a", "b", "c"}]
Out[5]= "a b c"
```

A custom separator

```mathematica
In[6]:= StringRiffle[{"2024", "01", "02"}, "-"]
Out[6]= "2024-01-02"
```

A {left, sep, right} triple

```mathematica
In[7]:= StringRiffle[{"x", "y", "z"}, {"(", ", ", ")"}]
Out[7]= "(x, y, z)"
```

Non-string elements are converted as ToString would

```mathematica
In[8]:= StringRiffle[{1, 2, 3}, ", "]
Out[8]= "1, 2, 3"
```

Nested: rows by newline, cells by space

```mathematica
In[9]:= StringRiffle[{{"a", "b"}, {"c", "d"}}]
Out[9]= "a b\nc d"
```

One separator per level, outermost first

```mathematica
In[10]:= StringRiffle[{{1, 2}, {3, 4}}, "; ", ","]
Out[10]= "1,2; 3,4"
```

The inverse of StringSplit

```mathematica
In[11]:= StringRiffle[StringSplit["a  b   c"], "_"]
Out[11]= "a_b_c"
```

The first argument must be a list: StringRiffle::list

```mathematica
In[12]:= StringRiffle[5]
Out[12]= StringRiffle[5]
```

## Algorithm

stringriffle.c - StringRiffle builtin for Mathilda

StringRiffle assembles a string from a (possibly nested) list of elements by inserting separators between them - the inverse of StringSplit. Non-string leaves are converted with expr_to_string (so the integer 27 -> "27"); string leaves are used verbatim.

Forms:

```text
  StringRiffle[list]                       - default separator scheme: a
      single space at the innermost level, one extra newline per level going
      up (2-D: rows by "\n", cells by " "; 3-D: blocks by "\n\n", ...).
  StringRiffle[list, sep]                  - string sep between the top-level
      elements; deeper levels fall back to the default scheme.
  StringRiffle[list, {"l","sep","r"}]      - a 3-string list is a delimiter
      triple: wrap the join with "l"..."r", joining with "sep".
  StringRiffle[list, sep1, sep2, ...]      - sep_i (a string or a delimiter
      triple) between elements at level i (1 = outermost); deeper levels use
      the default scheme.
```

Strings are treated as raw byte arrays (consistent with the rest of the string subsystem); no UTF-8 decoding is performed.

## Implementation notes

**Algorithm.** `builtin_stringriffle` treats the first argument as the data and the rest as separators, one per level, with level 1 outermost. `parse_sep` reads each separator into a `{left, sep, right}` triple. A plain string becomes `{"", sep, ""}`, and a 3-string list is the triple itself. `riffle_build` then recurses over the data. A string leaf is copied verbatim, and any other leaf is rendered in OutputForm (`expr_to_string_outputform`, which is what `ToString` gives). Each level joins its children with its explicit separator. When a level has none, the default scheme applies, chosen by `depth_from_bottom`: one space at the innermost level, and one more newline for each level above it. It is the inverse of `StringSplit`.

**Data structures.** A `SepSpec` array that borrows the argument strings. The whole result is built in one growable `RegexBuf`, so each leaf is copied exactly once whatever the nesting depth, and the buffer is adopted into the result string. This replaced a build-each-level-then-concatenate design in v0.456 (3.6 ms → 1.0 ms on the string benchmark).

**Complexity / limits.** `O(total output)`. No arguments emits `StringRiffle::argm`. A first argument that is neither a list nor a string emits `StringRiffle::list`, and a malformed separator emits `StringRiffle::string`; both leave the call unevaluated. Not `Listable`, so it can see the whole nested structure.

**Attributes:** `Protected`.

## References

**See also:** [StringSplit](../../string-operations/StringSplit/), [ToString](../../expression-information/ToString/)

- Source: [`src/strings/stringriffle.c`](https://github.com/stblake/mathilda/blob/main/src/strings/stringriffle.c)
- Specification: [`docs/spec/builtins/string-operations.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/string-operations.md)
- Tests: [`tests/test_string_review.c`](https://github.com/stblake/mathilda/blob/main/tests/test_string_review.c)
- Tests: [`tests/test_stringriffle.c`](https://github.com/stblake/mathilda/blob/main/tests/test_stringriffle.c)

## Notes & additional examples

### Notes

- `StringRiffle[list, sep1, sep2, ...]` joins a possibly nested list. `sep_i`
  goes between elements at level `i`, counting from the outside. A separator is
  a string or a `{left, sep, right}` triple.
- Levels without an explicit separator use the default scheme: one space at the
  innermost level and one more newline per level above it.
- String leaves are used verbatim, and any other leaf appears in its OutputForm.
- The result is built in one buffer, so each leaf is copied once whatever the
  depth (v0.456: 3.6 ms → 1.0 ms on the string benchmark).
- See also: `StringJoin`, `StringSplit`, `Riffle`.
