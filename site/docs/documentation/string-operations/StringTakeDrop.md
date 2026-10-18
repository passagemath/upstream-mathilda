# StringTakeDrop

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`StringTakeDrop["string", spec]`**

Gives {StringTake\["string", spec\], StringDrop\["string", spec\]}.

**`StringTakeDrop[{s1, s2, ...}, spec]`**

Gives one pair for each of the si.

## Examples (11)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= StringTakeDrop["abcdef", {2, 4}]
Out[1]= {"bcd", "aef"}

In[2]:= StringTakeDrop["abcdef", 2]
Out[2]= {"ab", "cdef"}

In[3]:= StringTakeDrop[{"abcd", "wxyz"}, -1]
Out[3]= {{"d", "abc"}, {"z", "wxy"}}
```

### Applications (8)

```mathematica
In[4]:= StringTakeDrop["abcdef", 2]
Out[4]= {"ab", "cdef"}
```

The last two characters

```mathematica
In[5]:= StringTakeDrop["abcdef", -2]
Out[5]= {"ef", "abcd"}
```

```mathematica
In[6]:= StringTakeDrop["abcdef", {2, 4}]
Out[6]= {"bcd", "aef"}
```

A stepped span splits odd and even positions

```mathematica
In[7]:= StringTakeDrop["abcdefgh", {1, -1, 2}]
Out[7]= {"aceg", "bdfh"}
```

Clamping with UpTo instead of failing

```mathematica
In[8]:= StringTakeDrop["abcdef", UpTo[10]]
Out[8]= {"abcdef", ""}
```

Positions are characters

```mathematica
In[9]:= StringTakeDrop["héllo", 2]
Out[9]= {"hé", "llo"}
```

One pair per string

```mathematica
In[10]:= StringTakeDrop[{"abcd", "wxyz"}, 2]
Out[10]= {{"ab", "cd"}, {"wx", "yz"}}
```

Out of range: StringTakeDrop::take

```mathematica
In[11]:= StringTakeDrop["abc", 5]
Out[11]= StringTakeDrop["abc", 5]
```

## Algorithm

stringtakedrop.c - StringTakeDrop["string", spec].

Gives {StringTake["string", spec], StringDrop["string", spec]} for one standard sequence specification (n, -n, {n}, {m, n}, {m, n, s}, UpTo[n], All, None), sharing StringTake's / StringDrop's scalar kernels. A list of strings gives one pair per string. Positions count characters.

Diagnostics are StringTakeDrop's own (Mathematica): ::strse for a subject that is not a string or list of strings, ::seqs for an invalid spec, ::take for out-of-range positions -- all leaving the call unevaluated.

## Implementation notes

**Algorithm.** `builtin_stringtakedrop` requires exactly two arguments (else `::argr` via `strings_argr_error`). The subject must be a `String` or a `List` of strings (else `::strse`), and the spec is parsed by `strings_seq_parse` into one `StrSeq` — `n`, `-n`, `{n}`, `{m, n}`, `{m, n, s}`, `UpTo[n]`, `All`, `None` (else `::seqs`). Before producing anything, the spec is resolved against *every* subject with `strings_chars_init` and `strings_seq_resolve`, so an out-of-range position reports `StringTakeDrop::take` rather than the underlying kernels' own messages. Each pair is then `{strings_take_seq(s, q), strings_drop_seq(s, q)}` — exactly the scalar kernels of `StringTake` and `StringDrop` — so the two halves always partition the string's characters.

**Data structures.** `StrSeq` (a kind tag plus `m`, `n`, step `s`) and a `StrChars` view that indexes character positions over the UTF-8 bytes, built only for the prefix the spec needs (`strings_seq_prefix`).

**Complexity / limits.** `O(B)` per subject for the take and the drop. Only a single sequence spec is accepted (no list of specs), matching Mathematica.

**Attributes:** `Protected`.

## References

**See also:** [StringTake](../../string-operations/StringTake/)

- Source: [`src/strings/stringtakedrop.c`](https://github.com/stblake/mathilda/blob/main/src/strings/stringtakedrop.c)
- Specification: [`docs/spec/builtins/string-operations.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/string-operations.md)
- Tests: [`tests/test_string_review.c`](https://github.com/stblake/mathilda/blob/main/tests/test_string_review.c)

## Notes & additional examples

### Notes

- Equivalent to `{StringTake[s, spec], StringDrop[s, spec]}` and built on the same
  scalar kernels (`strings_take_seq` / `strings_drop_seq`), so the two halves always
  partition the characters of `s`.
- Accepts a single sequence spec: `n`, `-n`, `{n}`, `{m, n}`, `{m, n, s}`, `UpTo[n]`,
  `All`, `None`. Positions are validated against every string first, so an
  out-of-range spec reports `StringTakeDrop::take` (not `StringTake::take`).
- Messages: `::strse` (subject), `::seqs` (spec), `::take` (range), `::argr`
  (argument count). Attributes `{Protected}`.
- See also: `StringTake`, `StringDrop`, `TakeDrop`.
