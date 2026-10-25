# 100 — RSA key generation (group F)

The integer kernels an RSA key generator leans on, timed as builtins (`prim:`)
and as a whole keygen written in each language (`app:`).

| case | what |
|---|---|
| `prim: NextPrime 1024-bit x4` | `NextPrime` / `gmpy2.next_prime` / `nextprime` from shared LCG starts |
| `prim: PrimeQ 2048-bit prime x8` | primality of eight 2048-bit primes (the expensive, all-rounds path) |
| `prim: PrimeQ 2048-bit odd x400` | 400 random odd 2048-bit candidates (the early-reject path a search lives on) |
| `prim: PowerMod[e,-1,lambda] 2048-bit x20000` | private-exponent inverse, portable spelling |
| `prim: ModularInverse 2048-bit x20000` | same, WL's `ModularInverse` — **ABSENT** in Mathilda today |
| `app: RSA-1024 keygen x8` / `app: RSA-2048 keygen x2` | incremental odd-candidate search (each system's own primality test) + λ, d, CRT parameters |

Baselines: Python = gmpy2 (`is_prime(n, 25)`, `next_prime`, `invert`), built-in
`pow(e, -1, m)`; PARI/GP = `nextprime`, `ispseudoprime` (BPSW), `Mod(e, l)^-1`.
Every system starts from the same `crypto_data` LCG integers, so all must find
the same primes; checks are sums folded mod 10⁹+7.
