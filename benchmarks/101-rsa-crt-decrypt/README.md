# 101 — RSA-2048 encrypt / decrypt, plain vs CRT (group F)

The modular-exponentiation workload of textbook RSA once a key exists:
encryption with e = 65537, plain decryption `c^d mod n`, and the CRT (Garner)
decryption real implementations use — two half-size exponentiations mod p and
mod q recombined as `m = m2 + q·(qinv·(m1 − m2) mod p)`. The key is built as in
100-rsa-keygen (incremental search from shared LCG starts, seeds 1101/1102),
outside every timed region; messages are 2000-bit LCG integers (seeds
2001–6000). Checks fold the ciphertexts / recovered plaintexts mod 10⁹+7, so a
wrong decryption cannot pass.

| case | what |
|---|---|
| `prim: PowerMod encrypt e=65537 2048-bit x4000` | `PowerMod[m, e, n]` mapped / `gmpy2.powmod` / `Mod(m, n)^e` |
| `prim: PowerMod[list] encrypt e=65537 2048-bit x4000` | the Listable form `PowerMod[list, e, n]` / `gmpy2.powmod_base_list` / `apply` |
| `prim: PowerMod decrypt c^d mod n 2048-bit x50` | full 2047-bit private exponent, no CRT |
| `app: RSA-2048 CRT decrypt x200` | Garner CRT decryption written in the language |
| `prim: ChineseRemainder recombine 2048-bit x20000` | `ChineseRemainder[{m1, m2}, {p, q}]` / `sympy.ntheory.modular.crt` / `chinese` (4000 pairs x5) |
| `app: Garner recombine 2048-bit x20000` | the same recombination as `m2 + Mod[qinv (m1 - m2), p] q` |
| `app: RSA-2048 encrypt+CRT decrypt round trip x200` | count of messages that survive encrypt→CRT-decrypt (must be 200) |

Baselines: Python "prim" = gmpy2 (`powmod`, `powmod_base_list`) and sympy `crt`;
Python "app" = plain ints and built-in `pow()` (idiomatic pure-Python RSA);
PARI/GP = `Mod(·,·)^e`, `chinese`, plain GP integers for the app cases.

## Findings

No wrong answers; all CHECKs agree in Mathilda, Python, PARI and Mathematica.
`PowerMod` is Listable in Mathilda (as in WL); the list form costs the same as
the mapped form (no batched kernel). Measured (ms, Mathilda / Python / PARI /
Mathematica, 2026-10-21):

| case | Mathilda | Python | PARI | Mathematica |
|---|---|---|---|---|
| encrypt x4000 | 81.9 | 73.5 | 94 | 100.2 |
| encrypt [list] x4000 | 83.3 | 71.9 | 93 | 99.8 |
| decrypt plain x50 | 109.5 | 97.2 | 144 | 223.8 |
| CRT decrypt x200 | 126.5 | 1320.8 | 200 | 403.8 |
| ChineseRemainder x20000 | 194.1 | 410.2 | 144 | 1523.4 |
| Garner recombine x20000 | 96.2 | 101.5 | 21 | 64.8 |
| round trip x200 | 132.3 | 1376.2 | 202 | 459.7 |

Mathilda's modexp is at GMP speed. The gap is per-operation interpreter
overhead on cheap bignum steps: the Garner recombination (one subtract, one
multiply, one `Mod`, one multiply-add on 1024/2048-bit integers per element) is
4.6x slower than PARI, and `ChineseRemainder` is 1.35x slower than `chinese`.
(Python's pure-`pow` CRT is slow because CPython's built-in `pow` is not GMP.)
