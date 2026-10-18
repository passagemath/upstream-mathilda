# Decapitalize

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Decapitalize["string"]`**

Gives "string" with its first character in lowercase.

**`Decapitalize[{s1, s2, ...}]`**

Decapitalizes each string.

## Examples (10)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= Capitalize["hello world"]
Out[1]= "Hello world"

In[2]:= Capitalize["élan"]
Out[2]= "Élan"

In[3]:= Decapitalize[{"Alpha", "Beta"}]
Out[3]= {"alpha", "beta"}
```

### Applications (7)

```mathematica
In[4]:= Decapitalize["Hello World"]
Out[4]= "hello World"
```

Only the first character is lowered

```mathematica
In[5]:= Decapitalize["ABC"]
Out[5]= "aBC"
```

```mathematica
In[6]:= Decapitalize["ÉCOLE"]
Out[6]= "éCOLE"

In[7]:= Decapitalize["Ωmega"]
Out[7]= "ωmega"
```

A list of strings maps

```mathematica
In[8]:= Decapitalize[{"Alpha", "Beta"}]
Out[8]= {"alpha", "beta"}
```

A non-letter first character is unchanged

```mathematica
In[9]:= Decapitalize["123abc"]
Out[9]= "123abc"
```

A list with a non-string: Decapitalize::strse

```mathematica
In[10]:= Decapitalize[{"Alpha", 1}]
Out[10]= Decapitalize[{"Alpha", 1}]
```

## Implementation notes

**Algorithm.** `builtin_decapitalize` calls `capitalize_builtin(res, "Decapitalize", 0)`. A `String` argument goes to the shared case kernel `map_case(s, len, 0, 0, 1)` — the same routine as `ToLowerCase`, restricted by its `from`/`count` window to character 0, so only the first code point is looked up in `k_lower_map` and every later byte is copied verbatim. A flat `List` of strings (checked by `is_string_list`) is mapped element by element into a new `List`; anything else emits `Decapitalize::strse` through `strings_arg_message` (the `mth_message` funnel) and stays unevaluated.

**Data structures.** The generated full-case-mapping table `k_lower_map` (`U8CaseMap` records, Unicode 16.0.0, from `tools/gen_unicode_tables.py`), searched by binary search; ASCII is mapped inline. The first character may expand to up to three code points (`Decapitalize["İ"]` gives `i` plus U+0307`).

**Complexity / limits.** One `O(L)` copy of the string plus one `O(log M)` table lookup. It applies the lowercase mapping, not Unicode titlecasing. Not `Listable`: nested or mixed lists are rejected rather than threaded. No numeric (packed / `Compile[]`) surface.

**Attributes:** `Protected`.

## References

**See also:** [Capitalize](../../string-operations/Capitalize/)

- The Unicode Consortium, *The Unicode Standard, Version 16.0* (2024), §3.13 Default Case Algorithms, and the Unicode Character Database (UnicodeData.txt, SpecialCasing.txt).
- Source: [`src/strings/stringcase.c`](https://github.com/stblake/mathilda/blob/main/src/strings/stringcase.c)
- Specification: [`docs/spec/builtins/string-operations.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/string-operations.md)
- Tests: [`tests/test_string_review.c`](https://github.com/stblake/mathilda/blob/main/tests/test_string_review.c)

## Notes & additional examples

### Notes

- Lowercases the first character with the full Unicode table used by `ToLowerCase`;
  everything after it is copied unchanged.
- Accepts a string or a flat list of strings (not `Listable`); anything else
  emits `Decapitalize::strse` and stays unevaluated.
- Attributes `{Protected}`.
- See also: `Capitalize`, `ToLowerCase`.
