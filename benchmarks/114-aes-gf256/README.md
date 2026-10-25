# 114 — AES-128 and GF(2⁸) (group F)

AES-128 (FIPS-197) written in the language, from the field arithmetic up:
GF(2⁸) multiply by repeated `xtime`, the inverse as a²⁵⁴, the S-box as inverse
plus affine map, key expansion, and full block encryption. The experiment also
has a GF(2)[x] polynomial route to the field inverse, and each system's native
AES library where one exists. Every app case is written twice, as the **same
algorithm**:

- **idiomatic** uses `BitXor`/`BitAnd`/`BitOr`/`BitShiftLeft`/`BitShiftRight`
  (Python `^ & | << >>`, GP `bitxor`/`bitand`/`bitor`/`shift`). Mathilda has
  none of the bitwise heads, so these rows are `benchIf`'d and reported ABSENT.
- **(no-bitwise)** XORs bytes through a 256×256 table built once (untimed) from
  `IntegerDigits[., 2, 8]`. `xtime` becomes `If[b >= 128, xor(2b-256, 27), 2b]`
  and rotations become `Mod`/`Quotient`. Python and GP run the identical variant.

| case | what |
|---|---|
| `app: AES S-box generation x4` / `(no-bitwise)` | 256 inverses (a²⁵⁴ by 7 square-and-multiply steps) + affine map, 4 times |
| `app: AES-128 FIPS-197 C.1 x100` / `(no-bitwise)` | key expansion + one block of the FIPS-197 Appendix C.1 vector, 100 times |
| `app: AES-128 ECB encrypt 128 blocks` / `(no-bitwise)` | one LCG key, 128 LCG blocks |
| `prim: GF(2^8) inverse via PolynomialExtendedGCD Modulus->2 x7650` | all 255 inverses in F₂[x]/(x⁸+x⁴+x³+x+1), 30 times |
| `prim: AES-128 ECB library 65536 blocks` | WL `Encrypt` with an AES128/ECB `SymmetricKey` (ABSENT in Mathilda) vs pycryptodome `AES.MODE_ECB`. PARI skips it (no AES) |

Checks fold each byte list as `FromDigits[bytes, 256] mod 10⁹+7`. The FIPS-197
vector must give `69c4e0d86a7b0430d8cdb78070b4c55a`, which is CHECK **413762112**.
Python also asserts the hex. The GF(2)[x] check is Σ a·a⁻¹ mod 10⁹+7. WL's
`Encrypt` pads with PKCS#7, so only the first 16·N bytes are compared.

Baselines: Python uses plain lists and ints, the GF(2)[x] route is a
pure-Python extended Euclid on bit-mask polynomials, and the library is
pycryptodome (C, AES-NI). PARI/GP uses GP vectors, and its GF(2)[x] route is
`Mod(Mod(1,2)*poly, aespoly)^-1`.

## Findings

Times in ms, one run (Mathematica = wolframscript 14):

| case | Mathilda | Python | PARI | Mathematica |
|---|---:|---:|---:|---:|
| S-box x4 | ABSENT | 17.2 | 130 | 380 |
| S-box x4 (no-bitwise) | 534 | 19.6 | 158 | 439 |
| FIPS C.1 x100 | ABSENT | 10.0 | 41 | 94 |
| FIPS C.1 x100 (no-bitwise) | 109 **CHECK-FAIL** | 12.3 | 61 | 180 |
| ECB 128 blocks | ABSENT | 9.2 | 38 | 104 |
| ECB 128 blocks (no-bitwise) | 209 | 11.9 | 62 | 178 |
| GF(2⁸) inverse PolynomialExtendedGCD x7650 | 169 | 28.6 | 32 | 474 |
| ECB library 65536 blocks | ABSENT | 0.22 | SKIP | 10.7 |

- **WRONG ANSWER: `AppendTo`/`PrependTo` on a packed list silently does
  nothing.** The call returns unevaluated and the variable is unchanged, with
  no message. The FIPS key is `Range[0, 15]`, which is packed, so
  `expandKeyN` (`w = Partition[key, 4]; ... AppendTo[w, ...]`) never grows `w`.
  The FIPS row then CHECK-FAILs with a half-symbolic expression. The ECB row
  passes only because its LCG key list happens not to be packed. With
  `w = Append[w, ...]` substituted, the FIPS row checks correctly (413762112)
  in 101 ms. Minimal reproducer:
  ```
  w = Range[8]; AppendTo[w, 0]; w        (* Mathilda: {1,...,8}   WL: {1,...,8,0} *)
  w = Range[8]; PrependTo[w, 0]; w       (* Mathilda: {1,...,8}   WL: {0,1,...,8} *)
  w = Range[3]; AppendTo[w, 0]; w        (* {1,2,3,0}: unpacked short list works *)
  w = 1. Range[8]; AppendTo[w, 0.5]; w   (* real packed list: also a silent no-op *)
  ```
- **ABSENT:** `BitXor`, `BitAnd`, `BitOr`, `BitShiftLeft`, `BitShiftRight`
  (3 idiomatic rows), and `Encrypt`/`SymmetricKey`/`ByteArray` (the library
  row).
- The no-bitwise rows run 10–27× slower than pure Python: interpreter overhead
  on byte-sized scalar code (`MapThread`, `Part` on a 256×256 table, `Module`
  per call). That is comparable to Mathematica, which is about the same speed
  or faster here.
- `PolynomialExtendedGCD[..., Modulus -> 2]` is 2.8× faster than Mathematica
  but 5–6× slower than PARI's native F₂[x] inverse and pure-Python bit-mask
  Euclid.
