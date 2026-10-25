# 108 — Paillier homomorphic encryption (group F)

Paillier with 1024-bit primes (n ≈ 2048 bits, n² ≈ 4096 bits) and g = n + 1,
written in each language (`app:`): keygen λ = lcm(p−1, q−1), μ = L(g^λ mod n²)⁻¹
mod n with L(u) = (u−1)/n; encrypt c = g^m·r^n mod n² (r^n dominates), and the
g = n + 1 shortcut g^m ≡ 1 + m·n (mod n²) as its own case; decrypt
m = L(c^λ mod n²)·μ mod n; the additive homomorphism (product of 1000 ciphertexts
decrypts to the sum of the messages); and scalar multiplication E(m)^k = E(k·m).
Nearly all of it is exponentiation modulo the 4096-bit n², so the `prim:` case is
that kernel alone. Primes are `primeOf[seed, 1024]` (setup, untimed); the 1000
summed ciphertexts are built with the fast form from 20 precomputed r^n values
(setup, untimed).

| case | what |
|---|---|
| `prim: PowerMod 4096-bit modulus x50` | b^n mod n², `PowerMod` / `gmpy2.powmod` / `Mod(b, n²)^n` |
| `app: Paillier keygen x8 (given p, q)` | λ, μ from 8 precomputed prime pairs |
| `app: Paillier encrypt x20` | generic g^m·r^n (two exponentiations) |
| `app: Paillier encrypt x20 (g=n+1: 1+mn)` | (1 + m·n)·r^n — same ciphertexts, same check |
| `app: Paillier decrypt x20` | check requires decrypt == message |
| `app: homomorphic sum of 1000 ciphertexts + decrypt x10` | running product mod n² then decrypt, repeated 10× (PARI 20 ms floor); check requires result == Σm mod n |
| `app: Paillier scalar mult c^k x20 (k < n)` | c^k mod n², k a full 2048-bit plaintext scalar; check requires decrypt(c₁^k₁) == k₁·m₁ mod n |

Baselines: Python app = built-in `pow`, `math.lcm`, `//`; prim = `gmpy2.powmod`.
PARI = `lift(Mod(b, n²)^e)`, `Mod(x, n)^-1`, `lcm`, `\`.

## Findings

Measured (ms; Mathilda / Python / PARI / Mathematica):

| case | Mathilda | Python | PARI | Mathematica |
|---|---|---|---|---|
| prim: PowerMod 4096-bit modulus x50 | 406 | 386 | 549 | 734 |
| app: Paillier keygen x8 (given p, q) | 67 | 338 | 45 | 64 |
| app: Paillier encrypt x20 | 334 | 2692 | 341 | 387 |
| app: Paillier encrypt x20 (g=n+1: 1+mn) | 165 | 1746 | 221 | 227 |
| app: Paillier decrypt x20 | 166 | 1799 | 221 | 264 |
| app: homomorphic sum of 1000 ciphertexts + decrypt x10 | 159 | 1319 | 170 | 190 |
| app: Paillier scalar mult c^k x20 (k < n) | 167 | 1777 | 215 | 276 |

All checks agree in all four systems. No absences, no wrong answers. Mathilda's
`PowerMod` at a 4096-bit modulus is within 5% of `gmpy2.powmod` and ahead of PARI
and Mathematica; the one row where it trails is keygen (1.5× PARI), whose single
exponentiation per key has a ~2048-bit exponent λ — worth a look but not a gap
in the kernel itself.
