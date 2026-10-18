# Capitalize

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Capitalize["string"]`**

Gives "string" with its first character in uppercase.

**`Capitalize[{s1, s2, ...}]`**

Capitalizes each string.

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

Only the first character changes

```mathematica
In[4]:= Capitalize["hello world"]
Out[4]= "Hello world"
```

```mathematica
In[5]:= Capitalize["élan vital"]
Out[5]= "Élan vital"
```

The first character uses the full mapping, so it can expand

```mathematica
In[6]:= Capitalize["ßx"]
Out[6]= "SSx"
```

A list of strings maps

```mathematica
In[7]:= Capitalize[{"alpha", "beta", "gamma"}]
Out[7]= {"Alpha", "Beta", "Gamma"}
```

```mathematica
In[8]:= Capitalize[""]
Out[8]= ""
```

Title case, word by word

```mathematica
In[9]:= StringRiffle[Capitalize /@ StringSplit["the quick brown fox"]]
Out[9]= "The Quick Brown Fox"
```

Not a string: Capitalize::strse, unevaluated

```mathematica
In[10]:= Capitalize[5]
Out[10]= Capitalize[5]
```

## Implementation notes

**Algorithm.** `builtin_capitalize` calls `capitalize_builtin(res, "Capitalize", 1)`. A `String` argument goes to the shared case kernel `map_case(s, len, 1, 0, 1)` — the same routine as `ToUpperCase`, restricted by its `from`/`count` window to character 0, so only the first code point is looked up in `k_upper_map` and every later byte is copied verbatim. A flat `List` of strings (checked by `is_string_list`) is mapped element by element into a new `List`; anything else emits `Capitalize::strse` through `strings_arg_message` (the `mth_message` funnel) and stays unevaluated.

**Data structures.** The generated full-case-mapping table `k_upper_map` (`U8CaseMap` records, Unicode 16.0.0, from `tools/gen_unicode_tables.py`), searched by binary search; ASCII is mapped inline. The first character may expand to up to three code points (`Capitalize["ßx"]` is `"SSx"`).

**Complexity / limits.** One `O(L)` copy of the string plus one `O(log M)` table lookup. It applies the uppercase mapping, not Unicode titlecasing. Not `Listable`: nested or mixed lists are rejected rather than threaded. No numeric (packed / `Compile[]`) surface.

**Attributes:** `Protected`.

## References

**See also:** [Decapitalize](../../string-operations/Decapitalize/)

- The Unicode Consortium, *The Unicode Standard, Version 16.0* (2024), §3.13 Default Case Algorithms, and the Unicode Character Database (UnicodeData.txt, SpecialCasing.txt).
- Source: [`src/strings/stringcase.c`](https://github.com/stblake/mathilda/blob/main/src/strings/stringcase.c)
- Specification: [`docs/spec/builtins/string-operations.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/string-operations.md)
- Tests: [`tests/test_string_review.c`](https://github.com/stblake/mathilda/blob/main/tests/test_string_review.c)

## Notes & additional examples

### Notes

- Maps only the first character, with the same full Unicode uppercase table as
  `ToUpperCase`; the rest of the string is copied byte for byte.
- It uses the *uppercase*, not the Unicode titlecase, mapping.
- Accepts a string or a flat list of strings; it is not `Listable`, so a mixed
  or nested list emits `Capitalize::strse` and stays unevaluated.
- Attributes `{Protected}`.
- See also: `Decapitalize`, `ToUpperCase`.
