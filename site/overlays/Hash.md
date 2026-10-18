### Worked examples

```mathematica
In[1]:= Hash[1 + 1] == Hash[2]   (* the argument is evaluated first *)
```

```mathematica
In[1]:= With[{h = Hash[x^2 + y]}, IntegerQ[h] && h >= 0]   (* a non-negative integer below 2^64 *)
```

```mathematica
In[1]:= Hash[a + b] == Hash[b + a]   (* equal canonical forms give equal hashes *)
```

```mathematica
In[1]:= {Hash[1] == Hash[1.], Hash["x"] == Hash[x]}   (* the hash depends on type: exact 1 and machine 1., or a string and a symbol, differ *)
```

```mathematica
In[1]:= StringLength[Hash["abc", "Expression", "HexString"]]   (* HexString gives 16 lowercase hex digits *)
```

```mathematica
In[1]:= FromDigits[Hash["abc", "Expression", "HexString"], 16] == Hash["abc"]   (* every string format encodes the same integer *)
```

```mathematica
In[1]:= FromDigits[Hash["abc", "Expression", "Base36String"], 36] == Hash["abc"]   (* and so does Base36String *)
```

```mathematica
In[1]:= Hash[Range[10^5]] == Hash[FromPackedArray[Range[10^5]]]   (* a packed array hashes the same as the unpacked list *)
```

```mathematica
In[1]:= Length[Union[Hash /@ Range[10^4]]]   (* no collisions among the first ten thousand integers *)
```

```mathematica
In[1]:= Hash[x, "SHA256"]   (* cryptographic hash types are not implemented, so the call stays unevaluated *)
```

### Notes

- Only the `"Expression"` type is implemented: a 64-bit FNV-1a hash of the internal structure of the expression. It is the same from run to run of one build, but, as in Mathematica, it may change between versions and platforms. Use it as a key within a session, not as a stored digest; the examples above compare hashes rather than print them for that reason.
- Formats: `"Integer"` (default), `"DecimalString"`, `"HexString"`, `"HexStringLittleEndian"`, `"Base36String"`. `"ByteArray"` and `"Base64Encoding"` stay unevaluated, as do `"SHA256"`, `"MD5"`, `"CRC32"` and the other cryptographic types.
- Not `Listable`: `Hash[{1, 2, 3}]` hashes the whole list. A result above `2^63 - 1` is returned as a bigint.
- The hash of a node is cached, so hashing the same large expression again costs O(1). Packed arrays are hashed from their buffer without being unpacked.
- Attributes `{Protected}`.
- See also: `ByteCount`, `LeafCount`, `SameQ`, `Order`.
