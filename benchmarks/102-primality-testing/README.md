# 102 — primality testing (group F)

Probable-prime testing at cryptographic sizes, as each system's builtin
(`prim:`) and as the two textbook tests written in the language (`app:`):
Miller–Rabin with the first 12 prime bases, and Baillie–PSW (strong base-2
round, perfect-square screen, strong Lucas test with Selfridge method-A
parameters via `JacobiSymbol`). Inputs are deterministic: 30 primes (10 each of
512/1024/2048 bits, `primeOf`), 2000 odd 1024-bit LCG integers, the Mersenne
prime 2⁴⁴²³−1, and 47 composites built to fool weak tests — 8 Chernick
Carmichael numbers (6k+1)(12k+1)(18k+1) searched upward from LCG k, 7 literal
strong pseudoprimes to the first 4…13 prime bases (verified composite and
base-by-base), the 27 composite Mersenne numbers 2ᵖ−1 with 1000 < p < 1200 and
the Fermat numbers F₈…F₁₁ (every composite Mersenne/Fermat number is a base-2
strong pseudoprime, so these exercise the Lucas half of BPSW). Every check is a
count of "prime" verdicts. Expected: 12-base MR is fooled by exactly 2 inputs
(the spsp's to bases 2..37 and 2..41) — the algorithm being honest — while
BPSW and every builtin are fooled by none.

| case | what |
|---|---|
| `prim: PrimeQ 512/1024/2048-bit primes x30` | `PrimeQ` / `gmpy2.is_prime(n, 25)` / `ispseudoprime` |
| `prim: PrimeQ 1024-bit odd x2000` | early-reject path |
| `prim: PrimeQ pseudoprimes (Carmichael/spsp/Mersenne/Fermat) x47` | adversarial composites (count must be 0) |
| `prim: PrimeQ Mersenne 2^4423-1` | one 4423-bit prime |
| `prim: proven primality 256-bit prime x4` | `ProvablePrimeQ` (PrimalityProving\`) — **ABSENT** in Mathilda; PARI `isprime` (APRCL); Python `sympy.isprime` (BPSW — *not* a proof, nothing better in the stack) |
| `app: Miller-Rabin 12 bases primes x30` / `pseudoprimes x47` | natural WL spelling (early exit via `Return[]` inside `While`) |
| `app: Baillie-PSW primes x30` / `pseudoprimes x47` / `1024-bit odd x500` | natural WL spelling |
| `app: ... (no-Return) ...` (4 cases) | the same algorithms with early exits as loop conditions — a Mathilda-runnable variant, identical in Python/GP |

Baselines: Python prim = gmpy2 / sympy; Python app = plain ints, built-in
`pow()`, `sympy.jacobi_symbol`; PARI = `ispseudoprime`, `isprime`, plain GP
integers, `Mod()^`, `kronecker`.

## Findings

**1. WRONG ANSWER — `Return[]` inside `While` inside a function body returns
only from the `While`, not from the function.** Mathematica (and every other
WL reading) returns from the enclosing `SetDelayed` definition. Minimal
reproducer:

```wl
f[s_] := Module[{r = 1}, While[r < s, If[r == 1, Return[True]]; r++]; False];
f[3]   (* Mathilda: False   Mathematica: True *)
f6[s_] := (While[True, If[True, Return[True]]]; False);  f6[3]  (* same *)
```

(`Return` inside `Do` exits only the `Do` in both systems, so that is not the
bug.) It makes the natural Miller–Rabin and Baillie–PSW wrong: 4 CHECK-FAILs —
MR primes 15 vs 30, MR pseudoprimes 0 vs 2, BPSW primes 6 vs 30, BPSW odd 0 vs 3
(BPSW pseudoprimes passes only because "composite" is also the wrong-path
answer). Python, PARI and Mathematica all agree with the expected counts. The
`(no-Return)` variants agree everywhere.

**2. ABSENT — `ProvablePrimeQ`.** Mathematica takes 5.9 s for four 256-bit
proofs; PARI `isprime` 152 ms.

No other disagreement: every `PrimeQ` verdict (incl. all 47 pseudoprimes and
2⁴⁴²³−1) matches. Measured ms (Mathilda / Python / PARI / Mathematica):

| case | Mathilda | Python | PARI | Mathematica |
|---|---|---|---|---|
| PrimeQ primes x30 | 114.8 | 98.6 | 106 | 205.4 |
| PrimeQ odd x2000 | 91.7 | 82.2 | 221 | 275.7 |
| PrimeQ pseudoprimes x47 | 36.3 | 31.8 | 33 | 92.9 |
| PrimeQ 2^4423-1 | 54.4 | 50.5 | 31 | 89.4 |
| proven 256-bit x4 | ABSENT | 0.5 (BPSW) | 152 | 5873.7 |
| MR (no-Return) primes x30 | 318.9 | 3218.5 | 447 | 607.5 |
| MR (no-Return) pseudoprimes x47 | 37.9 | 175.6 | 35 | 72.6 |
| BPSW (no-Return) primes x30 | 517.4 | 1107.7 | 255 | 455.8 |
| BPSW (no-Return) odd x500 | 180.3 | 1708.8 | 280 | 475.6 |
| BPSW pseudoprimes x47 (natural) | 199.3 | 212.6 | 72 | 315.7 |

`PrimeQ` is at gmpy2 speed and 1.7x slower than PARI on the 4423-bit Mersenne
prime. In-language BPSW is 2.0x slower than PARI on primes and 2.8x on the
pseudoprime set: the strong-Lucas ladder is ~5 `Mod` / multiply steps per bit on
1–2 kbit integers, i.e. per-operation interpreter overhead, while the
modexp-dominated MR rows are at parity or better.
