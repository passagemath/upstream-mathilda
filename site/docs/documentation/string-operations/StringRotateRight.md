# StringRotateRight

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`StringRotateRight["string", n]`**

Cycles the characters of "string" n places to the right.

**`StringRotateRight["string"]`**

Cycles one place to the right.

## Examples (9)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= StringRotateLeft["abcdef", 2]
Out[1]= "cdefab"

In[2]:= StringRotateRight["abcdef", 2]
Out[2]= "efabcd"

In[3]:= StringRotateLeft["abcdef", -8]
Out[3]= "efabcd"
```

### Applications (6)

```mathematica
In[4]:= StringRotateRight["abcdef", 2]
Out[4]= "efabcd"
```

One place by default

```mathematica
In[5]:= StringRotateRight["abcdef"]
Out[5]= "fabcde"
```

The two are mirror images

```mathematica
In[6]:= StringRotateRight["abcdef", -2] === StringRotateLeft["abcdef", 2]
Out[6]= True
```

100 mod 3 = 1

```mathematica
In[7]:= StringRotateRight["abc", 100]
Out[7]= "cab"
```

```mathematica
In[8]:= StringRotateRight["αβγδ", 1]
Out[8]= "δαβγ"
```

Not threaded over lists: StringRotateRight::string

```mathematica
In[9]:= StringRotateRight[{"abc"}, 1]
Out[9]= StringRotateRight[{"abc"}, 1]
```

## Algorithm

stringrotate.c - StringRotateLeft / StringRotateRight.

```text
  StringRotateLeft[s, n]   cycles the characters of s n places to the left
  StringRotateLeft[s]      one place
  StringRotateRight[s, n]  n places to the right (StringRotateLeft[s, -n])
```

n may be any integer (negative rotates the other way; |n| > length wraps). Positions count characters (Unicode code points). The result is two memcpys of the subject's bytes.

## Implementation notes

**Algorithm.** `builtin_stringrotateright` calls the shared kernel `rotate(res, "StringRotateRight", 0)`. It accepts one or two arguments: the subject must be a `String` (otherwise `StringRotateRight::string`) and the optional shift a machine `Integer` (otherwise `StringRotateRight::int`), defaulting to 1; both messages route through `strings_arg_message` / `mth_message` and leave the call unevaluated. The character count `m` comes from `u8_length`; the shift is converted to an equivalent left rotation `k = (-n) mod m`, reduced into `[0, m)`, and its byte offset found with `u8_offset`. The result is the tail from that offset followed by the head — two `memcpy`s into one buffer adopted as the result `String`.

**Data structures.** Only the UTF-8 byte buffer of the subject; positions count Unicode code points, so multi-byte characters are never split.

**Complexity / limits.** `O(B)` in the byte length (one length scan, one offset scan, one copy). Any machine integer shift is accepted (it wraps modulo the length); a bignum shift is rejected with `::int`. A list of strings is not threaded.

**Attributes:** `Protected`.

## References

**See also:** [StringRotateLeft](../../string-operations/StringRotateLeft/)

- Source: [`src/strings/stringrotate.c`](https://github.com/stblake/mathilda/blob/main/src/strings/stringrotate.c)
- Specification: [`docs/spec/builtins/string-operations.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/string-operations.md)
- Tests: [`tests/test_string_review.c`](https://github.com/stblake/mathilda/blob/main/tests/test_string_review.c)

## Notes & additional examples

### Notes

- `StringRotateRight[s, n]` is `StringRotateLeft[s, -n]`; both share one kernel
  (`rotate` in `src/strings/stringrotate.c`), reducing `n` modulo the character
  count so any machine integer is accepted.
- A non-string first argument emits `::string`, a non-integer shift `::int`.
- Attributes `{Protected}`.
- See also: `StringRotateLeft`, `RotateRight`.
