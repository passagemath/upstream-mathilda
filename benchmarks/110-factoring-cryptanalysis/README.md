# 110 — Factoring cryptanalysis (group F)

Three classical integer-factoring algorithms written in each language (`app:`),
each on the family of semiprimes it is designed to break, then the native
factoring primitive (`prim:`) on balanced semiprimes (both factors ~half the
modulus — the honest hard case) at 80/100/120/160 bits.

- **Pollard rho (Brent).** Semiprimes with one small (~28-bit) factor;
  deterministic cycle (start 2, `c=1`) so all systems find the same factor.
- **Pollard p−1 (stage 1).** Semiprimes whose `p−1` is B-smooth: `a=2` raised
  to `pr^e` for each prime `pr ≤ B1` (exponents from a precomputed prime list),
  then `gcd(a−1, n)` reveals `p`.
- **Fermat.** Semiprimes with `|p−q|` small: `a` from `⌈√n⌉` upward until
  `a²−n` is a perfect square.
- **prim.** `FactorInteger` / `sympy.factorint` / PARI `factor` on balanced
  semiprimes.

Each check is the sum of the recovered smaller factors mod `$CheckMod`, plus a
`recovered` count, folded into the check as value + 10⁹+7 × count, (4 / 4 / 4) so a silent miss cannot hide behind a
matching sum.

Baselines: Python = plain-int Pollard rho/p−1/Fermat for `app:`, `sympy.factorint`
for `prim:` (gmpy2 is not a factoriser). PARI/GP = plain-GP rho/p−1/Fermat,
`factor` for `prim:`.

## Findings

- **Mathilda is the fastest of all four systems at `FactorInteger`** on balanced
  (hard) semiprimes — a strong positive result:

  | balanced size | Mathilda | sympy | PARI | Mathematica |
  |---|---|---|---|---|
  | 80-bit ×25 | 59 ms | 5418 ms | 22 ms | 93 ms |
  | 100-bit ×10 | 45 ms | 4766 ms | 30 ms | 316 ms |
  | 120-bit ×3 | 38 ms | 4143 ms | 38 ms | 260 ms |
  | 160-bit ×1 | 83 ms | 51218 ms | 232 ms | 2761 ms |

  At 160-bit Mathilda is **~3× faster than PARI, ~33× faster than Mathematica,
  and ~600× faster than sympy** (whose pure-Python factoring collapses past
  120-bit — 51 s for one 160-bit semiprime). All four agree on every smallest-
  factor sum.

- **WRONG direction, notable slowness: Fermat is ~30× slower in Mathilda than
  in Mathematica, and effectively unmeasurable elsewhere.** The in-language
  Fermat loop (identical algorithm in all systems) over 4 close-prime
  semiprimes (~750 iterations total):

  | case | Mathilda | Python | PARI | Mathematica |
  |---|---|---|---|---|
  | `app: Fermat (close primes) x4` | **3577 ms** | 0.21 ms | ~0 ms | 118 ms |

  Root cause (re-measured by the lead, 2026-10-21): **`Sqrt` of a
  machine-size non-square integer costs ~2.8 ms**. That is the radical
  simplification path that extracts the square part. Here `b2 = a² − n` is about
  46 bits. For comparison:
  - `FactorInteger` of the same `b2` takes 4 µs;
  - `Sqrt` of a 100-digit non-square takes 3 µs;
  - `Ceiling[Sqrt[n]]` of the 80-bit `n` takes 25 µs.

  So the cost is not bignum square roots in general. It is one slow path for
  integers that fit in a machine word, roughly 700× what the factorisation it
  needs costs. (An earlier note here blamed a ~4 ms bignum sqrt on 160-bit
  inputs, which was a misreading.) Reproducer:

  ```
  p = NextPrime[2^39 + 12345]; q = NextPrime[p + 2^25]; n = p q;
  a = Ceiling[Sqrt[n]]; b2s = Table[(a + k)^2 - n, {k, 100}];
  AbsoluteTiming[Do[Sqrt[b], {b, b2s}]]   (* ~0.28 s: 2.8 ms per call *)
  AbsoluteTiming[Do[FactorInteger[b], {b, b2s}]]   (* ~0.4 ms *)
  ```

  PARI `issquare` and Python `isqrt` each take under a microsecond, so no shared
  loop count can put both Mathilda (budget ≤ 5 s) and PARI (≥ 20 ms) in range.
  **The Fermat PARI/Python columns read ~0 ms by necessity, and that gap is the
  finding.** Two improvement targets:
  1. Fix the machine-integer `Sqrt` path.
  2. Add a perfect-square test. Mathilda has no `PerfectSquareQ`, so
     `IntegerQ[Sqrt[·]]` is the only idiom available.

- **Pollard rho and p−1 are correct and competitive.** rho ×4: Mathilda 349 ms
  vs PARI 78 ms / Python 49 ms / Mathematica 545 ms — the in-language modular
  loop is a constant-factor interpreter cost. p−1 stage-1 ×120: Mathilda 104 ms,
  PARI 20 ms, Python 5 ms, Mathematica 29 ms. Every attack recovers 100% of its
  factors in all four systems.

- **No wrong answers, no crashes, no hangs.** All computed checks agree across
  Mathilda, Python, PARI and Mathematica.

Baselines measured on this host (Apple M-series).
