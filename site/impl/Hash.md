---
references:
  - "G. Fowler, L. C. Noll, K.-P. Vo, D. Eastlake and T. Hansen, *The FNV Non-Cryptographic Hash Algorithm*, IETF Internet-Draft draft-eastlake-fnv."
source:
  - src/core.c
  - src/expr.c
---
**Algorithm.** `builtin_hash` (attribute `Protected`) takes 1–3 arguments. It calls `expr_hash` on the evaluated first argument, a 64-bit FNV-1a structural hash: start from the FNV offset basis, fold in the node type, then the payload — the 8 raw bytes of an integer or double, the bytes of a symbol name or string, the sign and GMP limbs of a bigint, and for a function the hashes of its head and each argument in order. With one argument, or with type `"Expression"`, the `uint64` is returned by `hash_make_integer` (a machine integer when it fits in `int64`, otherwise a bigint built with `mpz_import`). A third argument picks the format in `hash_format`: `"Integer"`, `"DecimalString"`, `"HexString"` (16 lowercase digits, big-endian), `"HexStringLittleEndian"` and `"Base36String"`. Any other type (`"SHA256"`, `"MD5"`, `"CRC32"`, …) or format (`"ByteArray"`, `"Base64Encoding"`) returns `NULL`, and the call stays unevaluated with no message.

**Data structures.** `expr_hash` stores its result in the node's `hash_cache` field and clears it on in-place mutation (`expr_invalidate_hash`), so hashing a large tree a second time is O(1). Packed arrays (`EXPR_NDARRAY`) are hashed by `hash_packed_level`, which walks the flat buffer level by level and mixes in the same `List`, `Integer`, `Real` and `Complex` encodings a boxed tree would produce. The packed and unpacked forms of a list therefore hash equally, and the buffer is never unpacked.

**Complexity / limits.** O(n) in the size of the expression on first call. The value depends only on the internal structure: it is the same from run to run of one build, but, as in Mathematica, it is implementation-defined and may change between versions. Exact `1` and machine `1.` hash differently, and so do a string and the symbol of the same name. It is not a cryptographic digest, and none of the cryptographic types are implemented.
