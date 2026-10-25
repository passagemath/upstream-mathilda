# 103 — finite-field Diffie–Hellman (group F)

The number theory behind classic mod-p Diffie–Hellman: safe-prime generation
(first odd q ≥ an LCG start with q and 2q+1 both prime, written around each
system's own primality test), a primitive root of a safe prime (builtin, and the
textbook test "least g ≥ 2 with g² ≢ 1 and g^q ≢ 1 mod p" since p−1 = 2q), the
RFC 3526 MODP-2048 exchange with g = 2 and 256-bit private exponents, and
peer-key validation (1 < y < p−1 and y in the order-q subgroup, i.e. a
quadratic residue — Euler's criterion via `PowerMod` vs `JacobiSymbol`). Half
the validation candidates are honest 2^a (2 is a QR as p ≡ 7 mod 8), half
arbitrary LCG integers; 75 of 100 are valid. All inputs from `crypto_data`.

| case | what |
|---|---|
| `app: safe-prime search 256-bit x4` | 4 safe primes from LCG starts (seeds 1–4) |
| `app: safe-prime search 512-bit x2` | 2 safe primes (seeds 9–10), single timed run |
| `prim: PrimitiveRoot 256-bit safe prime x800` | `PrimitiveRoot` / `sympy.primitive_root` / `znprimroot`, 8 safe primes x100 |
| `app: primitive root of safe prime 256-bit x800` | textbook test for p = 2q+1 |
| `prim: PowerMod 2^a mod MODP-2048 (256-bit a) x100` | public-key generation |
| `app: MODP-2048 DH exchange x50` | A=g^a, B=g^b, B^a == A^b; check = sum of secrets |
| `app: subgroup validation y^q mod p MODP-2048 x100` | Euler criterion, 2047-bit exponent |
| `prim: JacobiSymbol validation MODP-2048 x10000` | `JacobiSymbol` / `gmpy2.jacobi` / `kronecker` (100 candidates x100) |

Baselines: Python prim = gmpy2 (`powmod`, `jacobi`), `sympy.primitive_root`;
Python app = plain ints + built-in `pow()`, safe-prime search driven by
`gmpy2.is_prime(n, 25)`; PARI = `ispseudoprime`, `znprimroot`, `Mod()^`,
`kronecker`.

## Findings

No wrong answers or absences; every CHECK agrees in Mathilda, Python, PARI and
Mathematica. Measured ms (Mathilda / Python / PARI / Mathematica):

| case | Mathilda | Python | PARI | Mathematica |
|---|---|---|---|---|
| safe-prime 256-bit x4 | 218.2 | 148.7 | 349 | 1341.2 |
| safe-prime 512-bit x2 | 856.6 | 699.8 | 2111 | 4114.1 |
| PrimitiveRoot x800 | 694.7 | 631.3 | 247 | 488.6 |
| textbook primitive root x800 | 50.9 | 445.7 | 99 | 481.6 |
| PowerMod 2^a x100 | 29.6 | 29.0 | 40 | 48.0 |
| DH exchange x50 | 69.1 | 708.6 | 91 | 118.6 |
| validation y^q x100 | 285.9 | 2788.8 | 377 | 563.4 |
| JacobiSymbol x10000 | 230.1 | 195.8 | 1009 | 3152.1 |

**Slow builtin: `PrimitiveRoot` on a safe prime is 2.8x slower than PARI's
`znprimroot` and 13.6x slower than the four-line textbook test in Mathilda
itself** (0.87 ms vs 0.064 ms per 256-bit safe prime). Breakdown measured
in-process: `FactorInteger[p-1]` ≈ 0.24 ms (2.6x a `PrimeQ` of the 255-bit
cofactor, 0.09 ms), leaving ≈ 0.5 ms in `PrimitiveRoot`'s own candidate
testing. Reproducer: `p = 2 q + 1` with `q = NextPrime`-style 255-bit safe
cofactor; `Do[PrimitiveRoot[p], {100}] // AbsoluteTiming` vs
`Do[FactorInteger[p - 1], {100}]`.

Everything modexp-bound (DH, validation, public keys) is at GMP speed and
ahead of PARI.
