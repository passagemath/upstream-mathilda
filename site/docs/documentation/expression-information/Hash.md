# Hash

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Hash[expr]`**

gives an integer hash code for expr (the "Expression" hash).

**`Hash[expr, type]`**

gives the hash of the given type; only "Expression" is supported.

**`Hash[expr, type, "format"]`**

gives the code in the given format.

<details>
<summary>Notes</summary>

The "Expression" hash is computed from the internal representation of expr and may vary between computer systems and between versions of Mathilda. Supported formats: "Integer" (default), "DecimalString", "HexString", "HexStringLittleEndian", "Base36String". Unsupported hash types or formats stay unevaluated.

</details>

## Examples (15)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (5)

```mathematica
In[1]:= Hash[1 + 1] == Hash[2]
Out[1]= True

In[2]:= Head[Hash[{1, 2, 3}]]
Out[2]= Integer

In[3]:= StringLength[Hash["abc", "Expression", "HexString"]]
Out[3]= 16

In[4]:= FromDigits[Hash[5, "Expression", "HexString"], 16] == Hash[5]
Out[4]= True

In[5]:= Hash[x, "SHA256"]
Out[5]= Hash[x, "SHA256"]
```

### Applications (10)

The argument is evaluated first

```mathematica
In[6]:= Hash[1 + 1] == Hash[2]
Out[6]= True
```

A non-negative integer below 2^64

```mathematica
In[7]:= With[{h = Hash[x^2 + y]}, IntegerQ[h] && h >= 0]
Out[7]= True
```

Equal canonical forms give equal hashes

```mathematica
In[8]:= Hash[a + b] == Hash[b + a]
Out[8]= True
```

The hash depends on type: exact 1 and machine 1., or a string and a symbol, differ

```mathematica
In[9]:= {Hash[1] == Hash[1.], Hash["x"] == Hash[x]}
Out[9]= {False, False}
```

HexString gives 16 lowercase hex digits

```mathematica
In[10]:= StringLength[Hash["abc", "Expression", "HexString"]]
Out[10]= 16
```

Every string format encodes the same integer

```mathematica
In[11]:= FromDigits[Hash["abc", "Expression", "HexString"], 16] == Hash["abc"]
Out[11]= True
```

And so does Base36String

```mathematica
In[12]:= FromDigits[Hash["abc", "Expression", "Base36String"], 36] == Hash["abc"]
Out[12]= True
```

A packed array hashes the same as the unpacked list

```mathematica
In[13]:= Hash[Range[10^5]] == Hash[FromPackedArray[Range[10^5]]]
Out[13]= True
```

No collisions among the first ten thousand integers

```mathematica
In[14]:= Length[Union[Hash /@ Range[10^4]]]
Out[14]= 10000
```

Cryptographic hash types are not implemented, so the call stays unevaluated

```mathematica
In[15]:= Hash[x, "SHA256"]
Out[15]= Hash[x, "SHA256"]
```

## Implementation notes

**Algorithm.** `builtin_hash` (attribute `Protected`) takes 1–3 arguments. It calls `expr_hash` on the evaluated first argument, a 64-bit FNV-1a structural hash: start from the FNV offset basis, fold in the node type, then the payload — the 8 raw bytes of an integer or double, the bytes of a symbol name or string, the sign and GMP limbs of a bigint, and for a function the hashes of its head and each argument in order. With one argument, or with type `"Expression"`, the `uint64` is returned by `hash_make_integer` (a machine integer when it fits in `int64`, otherwise a bigint built with `mpz_import`). A third argument picks the format in `hash_format`: `"Integer"`, `"DecimalString"`, `"HexString"` (16 lowercase digits, big-endian), `"HexStringLittleEndian"` and `"Base36String"`. Any other type (`"SHA256"`, `"MD5"`, `"CRC32"`, …) or format (`"ByteArray"`, `"Base64Encoding"`) returns `NULL`, and the call stays unevaluated with no message.

**Data structures.** `expr_hash` stores its result in the node's `hash_cache` field and clears it on in-place mutation (`expr_invalidate_hash`), so hashing a large tree a second time is O(1). Packed arrays (`EXPR_NDARRAY`) are hashed by `hash_packed_level`, which walks the flat buffer level by level and mixes in the same `List`, `Integer`, `Real` and `Complex` encodings a boxed tree would produce. The packed and unpacked forms of a list therefore hash equally, and the buffer is never unpacked.

**Complexity / limits.** O(n) in the size of the expression on first call. The value depends only on the internal structure: it is the same from run to run of one build, but, as in Mathematica, it is implementation-defined and may change between versions. Exact `1` and machine `1.` hash differently, and so do a string and the symbol of the same name. It is not a cryptographic digest, and none of the cryptographic types are implemented.

- `Protected` (not `Listable`: `Hash[{1, 2, 3}]` hashes the whole list, not element-wise).
- The argument is evaluated first, so `Hash[1 + 1]` equals `Hash[2]`.
- Only the default `"Expression"` hash `type` is implemented: the 64-bit FNV-1a structural
  hash computed from `expr`'s internal representation. Like Mathematica's, this value is
  **implementation-defined** — it may differ between computer systems and between versions of
  Mathilda — so treat it as a stable key within one build, not a portable digest. Other hash
  types (`"SHA256"`, `{"SHA", b}`, `"MD5"`, `"CRC32"`, …) are left unevaluated for now.
- The result is a non-negative integer, returned as a bigint when it exceeds `2^63 - 1`.
- `Hash[expr, "Expression", "format"]` renders the same hash in a string `"format"`:
  - `"Integer"` (default) — the integer hash.
  - `"DecimalString"` — the unsigned value in base 10.
  - `"HexString"` — 16 lowercase hex digits, zero-padded, big-endian.
  - `"HexStringLittleEndian"` — the same eight bytes in little-endian order.
  - `"Base36String"` — the unsigned value in base 36 (digits `0-9`, `a-z`).
  Unsupported formats (`"ByteArray"`, `"Base64Encoding"`, …) are left unevaluated.

**Attributes:** `Protected`.

## References

- G. Fowler, L. C. Noll, K.-P. Vo, D. Eastlake and T. Hansen, *The FNV Non-Cryptographic Hash Algorithm*, IETF Internet-Draft draft-eastlake-fnv.
- Source: [`src/core.c`](https://github.com/stblake/mathilda/blob/main/src/core.c)
- Specification: [`docs/spec/builtins/expression-information.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/expression-information.md)
- Tests: [`tests/test_hash.c`](https://github.com/stblake/mathilda/blob/main/tests/test_hash.c)

## Notes & additional examples

### Notes

- Only the `"Expression"` type is implemented: a 64-bit FNV-1a hash of the internal structure of the expression. It is the same from run to run of one build, but, as in Mathematica, it may change between versions and platforms. Use it as a key within a session, not as a stored digest; the examples above compare hashes rather than print them for that reason.
- Formats: `"Integer"` (default), `"DecimalString"`, `"HexString"`, `"HexStringLittleEndian"`, `"Base36String"`. `"ByteArray"` and `"Base64Encoding"` stay unevaluated, as do `"SHA256"`, `"MD5"`, `"CRC32"` and the other cryptographic types.
- Not `Listable`: `Hash[{1, 2, 3}]` hashes the whole list. A result above `2^63 - 1` is returned as a bigint.
- The hash of a node is cached, so hashing the same large expression again costs O(1). Packed arrays are hashed from their buffer without being unpacked.
- Attributes `{Protected}`.
- See also: `ByteCount`, `LeafCount`, `SameQ`, `Order`.
