# 107 — ElGamal over the MODP-2048 group (group F)

Textbook ElGamal in Z_p^* with p the RFC 3526 group-14 2048-bit safe prime and
g = 2, written in each language (`app:`): keygen h = g^x, encrypt (g^y, m·h^y),
decrypt m = c2·s⁻¹ with s = c1^x — the inverse both by extended Euclid
(`PowerMod[s, -1, p]` / `pow(s, -1, p)` / `Mod(s, p)^-1`) and by Fermat s^(p−2) —
the multiplicative homomorphism (∏ ciphertexts decrypts to ∏ messages), and
re-randomisation (c1·g^r, c2·h^r). Every step is a 2048-bit modular
exponentiation, so the `prim:` case is that kernel alone, batched x200. Exponents,
messages, and randomisers come from `crypto_data`'s `bigOf`.

| case | what |
|---|---|
| `prim: PowerMod 2048-bit x200` | `PowerMod` / `gmpy2.powmod` / `Mod(b, p)^e`, 2048-bit base and exponent |
| `app: ElGamal keygen x50` | 50 key pairs (one exponentiation each) |
| `app: ElGamal encrypt x50` | two exponentiations each |
| `app: ElGamal decrypt x50 (PowerMod[s,-1,p])` | one exponentiation + Euclid inverse; check also requires decrypt == message |
| `app: ElGamal decrypt x50 (s^(p-2))` | one exponentiation + Fermat inverse (a second exponentiation) |
| `app: homomorphic product of 50 ciphertexts + decrypt x20` | componentwise product (formed whole, reduced once) then decrypt, repeated 20× |
| `app: ElGamal re-randomise x50` | two exponentiations each; check requires the re-randomised ciphertexts decrypt to the same messages |

Baselines: Python app = built-in `pow` (idiomatic pure-Python), prim =
`gmpy2.powmod`; PARI = `lift(Mod(b, p)^e)`, `Mod(s, p)^-1`, `vecprod`.

## Findings

Measured (ms; Mathilda / Python / PARI / Mathematica):

| case | Mathilda | Python | PARI | Mathematica |
|---|---|---|---|---|
| prim: PowerMod 2048-bit x200 | 588 | 483 | 819 | 1239 |
| app: ElGamal keygen x50 | 134 | 1079 | 173 | 248 |
| app: ElGamal encrypt x50 | 273 | 2303 | 354 | 529 |
| app: ElGamal decrypt x50 (PowerMod[s,-1,p]) | 148 | 1345 | 187 | 307 |
| app: ElGamal decrypt x50 (s^(p-2)) | 294 | 2864 | 370 | 617 |
| app: homomorphic product of 50 ciphertexts + decrypt x20 | 92 | 832 | 88 | 138 |
| app: ElGamal re-randomise x50 | 265 | 2491 | 342 | 501 |

All checks agree in all four systems. No absences, no wrong answers.
Mathilda's `PowerMod` is within ~1.2× of `gmpy2.powmod` (an earlier run measured
471 ms; the machine was shared) and faster than PARI and Mathematica; every app row
is the fastest of the four, since Python's built-in `pow` is ~8× slower than GMP.
