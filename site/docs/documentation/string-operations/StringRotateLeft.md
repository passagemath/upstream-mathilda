# StringRotateLeft

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`StringRotateLeft["string", n]`**

Cycles the characters of "string" n places to the left.

**`StringRotateLeft["string"]`**

Cycles one place to the left.

## Examples (11)

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

### Applications (8)

```mathematica
In[4]:= StringRotateLeft["abcdef", 2]
Out[4]= "cdefab"
```

One place by default

```mathematica
In[5]:= StringRotateLeft["abcdef"]
Out[5]= "bcdefa"
```

A negative shift rotates right

```mathematica
In[6]:= StringRotateLeft["abcdef", -1]
Out[6]= "fabcde"
```

Shifts wrap modulo the length

```mathematica
In[7]:= StringRotateLeft["abcdef", 13]
Out[7]= "bcdefa"
```

Positions count characters, not bytes

```mathematica
In[8]:= StringRotateLeft["αβγδ", 1]
Out[8]= "βγδα"
```

```mathematica
In[9]:= Table[StringRotateLeft["abcd", k], {k, 0, 4}]
Out[9]= {"abcd", "bcda", "cdab", "dabc", "abcd"}

In[10]:= StringRotateLeft["", 3]
Out[10]= ""
```

A symbolic shift: StringRotateLeft::int

```mathematica
In[11]:= StringRotateLeft["abc", x]
Out[11]= StringRotateLeft["abc", x]
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

**Algorithm.** `builtin_stringrotateleft` calls the shared kernel `rotate(res, "StringRotateLeft", 1)`. It accepts one or two arguments: the subject must be a `String` (otherwise `StringRotateLeft::string`) and the optional shift a machine `Integer` (otherwise `StringRotateLeft::int`), defaulting to 1; both messages route through `strings_arg_message` / `mth_message` and leave the call unevaluated. The character count `m` comes from `u8_length`; the shift is converted to an equivalent left rotation `k = n mod m`, reduced into `[0, m)`, and its byte offset found with `u8_offset`. The result is the tail from that offset followed by the head — two `memcpy`s into one buffer adopted as the result `String`.

**Data structures.** Only the UTF-8 byte buffer of the subject; positions count Unicode code points, so multi-byte characters are never split.

**Complexity / limits.** `O(B)` in the byte length (one length scan, one offset scan, one copy). Any machine integer shift is accepted (it wraps modulo the length); a bignum shift is rejected with `::int`. A list of strings is not threaded.

**Attributes:** `Protected`.

## References

**See also:** [StringRotateRight](../../string-operations/StringRotateRight/)

- Source: [`src/strings/stringrotate.c`](https://github.com/stblake/mathilda/blob/main/src/strings/stringrotate.c)
- Specification: [`docs/spec/builtins/string-operations.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/string-operations.md)
- Tests: [`tests/test_string_review.c`](https://github.com/stblake/mathilda/blob/main/tests/test_string_review.c)

## Notes & additional examples

### Notes

- The shift is reduced modulo the character count, the cut point located with
  `u8_offset`, and the result assembled with two `memcpy`s of the UTF-8 bytes.
- Takes a single string; a list of strings is not threaded (emits `::string`),
  so use `Map`. A non-integer shift emits `::int`.
- Attributes `{Protected}`.
- See also: `StringRotateRight`, `RotateLeft`, `StringTake`.
