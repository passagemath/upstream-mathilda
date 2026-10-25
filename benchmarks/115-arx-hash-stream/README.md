# 115 — ARX hashing and stream ciphers (group F)

The add-rotate-xor primitives on 32-bit words, written in the language:
SHA-256 (FIPS 180-4 padding, message schedule, 64-round compression),
HMAC-SHA256 (RFC 2104), and the ChaCha20 block function (RFC 7539). The native
hash primitive is the `prim:` row. Every app case is written twice, as the
**same algorithm**:

- **idiomatic** uses `BitXor`/`BitAnd`/`BitOr`/`BitShiftLeft`/`BitShiftRight`
  (Python `^ & | << >>`, GP `bitxor`/`bitand`/`bitor`/`shift`). Mathilda has
  no bitwise heads, so these rows are `benchIf`'d and reported ABSENT.
- **(no-bitwise)** does 32-bit XOR/AND as four 8-bit lookups in 256×256 tables
  built once (untimed) from `IntegerDigits[., 2, 8]`, splitting bytes off with
  `Quotient`/`Mod`. Rotations and shifts use `Quotient`/`Mod` by powers of two,
  and NOT is `2³²−1−x`. Python and GP run the identical variant.

| case | what |
|---|---|
| `app: SHA-256 abc x40` / `(no-bitwise)` | FIPS 180-4 one-block vector, 40 times (40 compressions) |
| `app: SHA-256 8x256-byte messages` / `(no-bitwise)` | 8 printable-ASCII LCG messages, 5 blocks each |
| `app: HMAC-SHA256 RFC4231 #2 x12` / `(no-bitwise)` | key "Jefe", message "what do ya want for nothing?" (4 compressions per tag) |
| `app: ChaCha20 keystream 48 blocks` / `(no-bitwise)` | RFC 7539 key 00..1f, nonce 000000090000004a00000000, counters 1..48 |
| `prim: Hash SHA256 8x256-byte messages x100` | WL `Hash[string, "SHA256"]` vs `hashlib.sha256`. Probe-SKIP in Mathilda; PARI has no hash |

Correctness anchors: SHA-256("abc") = `ba7816bf…f20015ad` (CHECK 997922219),
HMAC RFC 4231 #2 = `5bdcc146…64ec3843` (CHECK 626203177), and ChaCha20 block 1
starts `10f1e7e4 d13b5915 500fdd1f …` (RFC 7539 §2.3.2). The Python file
asserts all three against hashlib/hmac/the RFC, plus both variants against each
other. Digests fold as `FromDigits[words, 2³²] mod 10⁹+7`, and keystream blocks
as `FromDigits[bytes, 256]` summed. Messages are printable ASCII, so the same
bytes are also a WL string for `Hash`.

Baselines: Python uses plain ints, with `hashlib` (OpenSSL) for `prim:`.
PARI/GP uses plain GP integers and vectors.

## Findings

Times in ms, one run (Mathematica = wolframscript 14):

| case | Mathilda | Python | PARI | Mathematica |
|---|---:|---:|---:|---:|
| SHA-256 abc x40 | ABSENT | 9.1 | 22 | 82 |
| SHA-256 abc x40 (no-bitwise) | 587 | 66 | 148 | 567 |
| SHA-256 8x256 B | ABSENT | 9.6 | 23 | 81 |
| SHA-256 8x256 B (no-bitwise) | 583 | 65 | 147 | 576 |
| HMAC-SHA256 x12 | ABSENT | 10.9 | 28 | 98 |
| HMAC-SHA256 x12 (no-bitwise) | 699 | 75 | 200 | 679 |
| ChaCha20 48 blocks | ABSENT | 8.1 | 38 | 107 |
| ChaCha20 48 blocks (no-bitwise) | 350 | 32 | 94 | 312 |
| Hash SHA256 x100 | SKIP | 1.3 | SKIP | 337 |

- All checks agree in every system that runs the case. There are no wrong
  answers.
- **ABSENT:** `BitXor`, `BitAnd`, `BitOr`, `BitShiftLeft`, `BitShiftRight`,
  and `BitNot` (4 idiomatic rows). These are the heads that matter here. The
  no-bitwise rows cost about 7× the bitwise ones in every other system (Python
  66 vs 9 ms, Mathematica 567 vs 82 ms). So with native bitwise heads at
  Mathematica's per-op cost, Mathilda's SHA-256 would drop from about 590 ms to
  about 80 ms.
- **Probe-SKIP: `Hash["abc", "SHA256"]` returns unevaluated.** `Hash` exists,
  but it has no SHA-256 (or other cryptographic) method. Reproducer:
  `IntegerQ[Hash["abc", "SHA256"]]` gives `False` (Mathematica: `True`).
- On the no-bitwise rows Mathilda matches Mathematica (within 3–12%) and runs
  about 9–11× slower than pure Python. This is the interpreter cost of
  `Extract` + `Dot` per 32-bit op and `Module` per call.
