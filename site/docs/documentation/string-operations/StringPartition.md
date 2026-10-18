# StringPartition

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`StringPartition["string", n]`**

Partitions string into non-overlapping substrings of length n.

**`StringPartition["string", n, d]`**

Generates length-n substrings with offset d (all of length n; some trailing or middle characters may be omitted).

**`StringPartition["string", UpTo[n]]`**

Partitions into substrings of length up to n, allowing a shorter final substring so every character appears.

**`StringPartition[{s1, s2, ...}, spec]`**

Threads over a list of strings.

## Examples (11)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= StringPartition["123456789123456789", 9]
Out[1]= {"123456789", "123456789"}

In[2]:= StringPartition["123456789", 2, 1]
Out[2]= {"12", "23", "34", "45", "56", "67", "78", "89"}

In[3]:= StringPartition["123456789", UpTo[2]]
Out[3]= {"12", "34", "56", "78", "9"}

In[4]:= StringPartition["ababababab", 3]
Out[4]= {"aba", "bab", "aba"}
```

### Applications (7)

Non-overlapping blocks; a short tail is dropped

```mathematica
In[5]:= StringPartition["abcdefgh", 3]
Out[5]= {"abc", "def"}
```

Length 2, offset 1: overlapping blocks

```mathematica
In[6]:= StringPartition["abcdef", 2, 1]
Out[6]= {"ab", "bc", "cd", "de", "ef"}
```

An offset larger than the block skips characters

```mathematica
In[7]:= StringPartition["abcdef", 2, 3]
Out[7]= {"ab", "de"}
```

Allow a shorter final block

```mathematica
In[8]:= StringPartition["abcdefg", UpTo[3]]
Out[8]= {"abc", "def", "g"}
```

With an offset, UpTo keeps the clipped tail

```mathematica
In[9]:= StringPartition["abcdefgh", UpTo[3], 2]
Out[9]= {"abc", "cde", "efg", "gh"}
```

Blocks count characters, not bytes

```mathematica
In[10]:= StringPartition["h\:00e9llo", 2]
Out[10]= {"hé", "ll"}
```

As in Mathematica, a list of strings is not threaded

```mathematica
In[11]:= StringPartition[{"ab", "cd"}, 1]
Out[11]= StringPartition[{"ab", "cd"}, 1]
```

## Algorithm

stringpartition.c - StringPartition builtin for Mathilda

```text
StringPartition["string", n]        - non-overlapping length-n blocks
StringPartition["string", n, d]     - length-n blocks starting every d chars
StringPartition["string", UpTo[n]]  - length-<=n blocks; final may be shorter
```

StringPartition["string", UpTo[n], d] - the same, with offset d

Lengths and offsets count characters (code points), not UTF-8 bytes. As in Mathematica, the first argument must be a single string (no threading).

## Implementation notes

**Algorithm.** `builtin_stringpartition` takes 2–3 arguments and gives `::argtu`/`::argt` otherwise. The block length `n` is an integer, in which case every block is exactly `n` characters, or `UpTo[n]`, which clips each block to what remains. The offset `d` defaults to `n`. Block starts run over `0, d, 2d, …` in characters. A fixed `n` stops when a full block no longer fits. `UpTo[n]` takes every start before the end, so with `d < n` the tail repeats shorter and shorter blocks, as in Mathematica. With an explicit `d`, `n = 0` is allowed, and every block is then `""`.

**Data structures.** A `U8View` over the subject, so lengths and offsets count Unicode code points. Each block is created by `expr_new_string_n` straight from the byte span. The block array is bounded by `len/d + 1`.

**Complexity / limits.** `O(output)`. As in Mathematica, the first argument must be a single string. A list of strings is not threaded (a v0.456 change) and gives `StringPartition::string`. A bad block length gives `::intpm`, or `::intnm` once an offset is present. A non-positive or non-integer offset gives `::intp`. Each of these leaves the call unevaluated, and `UpTo[negative]` stays unevaluated without a message.

**Attributes:** `Protected`.

## References

- Source: [`src/strings/stringpartition.c`](https://github.com/stblake/mathilda/blob/main/src/strings/stringpartition.c)
- Specification: [`docs/spec/builtins/string-operations.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/string-operations.md)
- Tests: [`tests/test_string_review.c`](https://github.com/stblake/mathilda/blob/main/tests/test_string_review.c)
- Tests: [`tests/test_stringpartition.c`](https://github.com/stblake/mathilda/blob/main/tests/test_stringpartition.c)

## Notes & additional examples

### Notes

- `StringPartition[s, n, d]` gives the length-`n` substrings that start every `d`
  characters, with `d` defaulting to `n`. A tail too short for a whole block is
  dropped, unless the length is given as `UpTo[n]`.
- Since v0.456 lengths count Unicode code points, and a list first argument is
  rejected with `StringPartition::string` instead of being threaded, matching
  Mathematica. Map it over the list if needed.
- Other messages are `::intpm` (block length not a positive integer), `::intnm`
  (once an offset is given, `n = 0` is allowed) and `::intp` (bad offset). All
  leave the call unevaluated.
- See also: `Partition`, `StringTake`, `Characters`.
