# 113 — lattice attacks (group F)

LLL reduction as a cryptanalyst uses it: breaking a Merkle–Hellman knapsack
public key with the CJLOSS lattice (the improved Lagarias–Odlyzko attack),
recovering a secret from approximate multiples (approximate-GCD / hidden-number
lattice), and plain LLL on square knapsack-type (q-ary) bases. Each `prim:`
case times the whole attack: build basis, reduce, scan the reduced rows. The
legitimate Merkle–Hellman decryption is the `app:` control, so you can compare
the key holder's cost with the attacker's. LLL output is not unique across
implementations, so every check is an invariant: the recovered message or
secret integer, or `|det|` of the reduced basis, which must equal q because
the transform is unimodular.

| case | what |
|---|---|
| `prim: LatticeReduce CJLOSS knapsack n=24 x20` | 20 ciphertexts under one 24-item key (density 0.30), each broken by a 25-dim CJLOSS reduction |
| `prim: LatticeReduce CJLOSS knapsack n=32 x12` | 12 ciphertexts under one 32-item key (density 0.44), 33-dim |
| `app: Merkle-Hellman decrypt n=32 x2000` | key holder: `c·m⁻¹ mod q`, then a greedy superincreasing solve |
| `app: approximate-GCD secret recovery dim=9 x40` | 40 instances of 9 samples `x_i = p q_i + r_i` (128-bit prime p, 16-bit noise); p recovered from the first reduced row |
| `prim: LatticeReduce q-ary basis n=30 300-bit` | plain LLL, basis `(q,0..0),(a_i,e_i)` |
| `prim: LatticeReduce q-ary basis n=50 200-bit` | same, n = 50 |

Baselines: Python uses sympy `DomainMatrix(..., ZZ).lll()`, a pure-Python exact
LLL with δ = 3/4 (no python-flint or fpylll is installed). PARI/GP uses
`qflll(B~, 3)~`: qflll reduces the **columns** of its argument, so the basis is
transposed in and back out. PARI's δ is 0.99. Mathilda's `LatticeReduce` is an
exact-rational (GMP `mpq`) LLL with δ = 3/4. The cases above 1 s use
`benchOnce`.

**Inputs are not raw LCG output.** The `crypto_data` LCG is linear, and LLL is
exactly the tool that detects linearity. With `bigOf[seed + i, bits]` entries
(every word in exact arithmetic progression), the 30-dim q-ary basis reduced to
a vector of squared norm 4. With chunks of one `lcgStates` stream it reduced to
squared norm 120. Uniformly random entries give about 3.8·10⁶ (the Gaussian
heuristic). Each file therefore defines `mixBigs`, which maps every LCG word
through `3^e mod P` before truncating. After that change, reduction behaves like
the random case (squared norms of about 2–9·10⁶). **Other group-F lattice
experiments that feed `bigOf[seed+i]` vectors to LLL are reducing a planted
trivial lattice.**

## Findings

All checks agree in Mathilda, Python, PARI and Mathematica 14. Times are in ms,
from one run (the Mathematica column is wolframscript):

| case | Mathilda | Python (sympy) | PARI | Mathematica |
|---|---:|---:|---:|---:|
| CJLOSS n=24 x20 | 2951 | 1319 | 502 | 588 |
| CJLOSS n=32 x12 | 2939 | 1045 | 80 | 491 |
| MH decrypt x2000 | 128 | 14.5 | 48 | 118 |
| approx-GCD x40 | 3789 | 2831 | 24 | 409 |
| q-ary n=30 300-bit | 2109 | 1500 | 26 | 266 |
| q-ary n=50 200-bit | 4159 | 2244 | 87 | 494 |

- **`LatticeReduce` is about 2× slower than sympy's pure-Python LLL** (which
  uses the same exact δ = 3/4 algorithm and returns the identical reduced basis),
  5–10× slower than Mathematica, and 25–160× slower than PARI. The cause is the
  exact `mpq` Gram–Schmidt in `src/linalg/latticereduce.c`, with no
  floating-point (L²/fplll-style) or integral (de Weger / Cohen 2.6.7) variant.
  For an attack loop of many medium lattices this is the dominant cost. PARI's
  oddly slow n=24 row (502 ms vs 80 ms at n=32) reproduces, so it is a property
  of qflll on this instance.
- **δ = 3/4 vs 0.99 is visible in attack success.** At (n, k) = (32, 24),
  density 0.57, Mathilda and sympy (δ = 3/4) both miss one ciphertext in 20:
  no ±1 vector appears in the reduced basis. PARI's δ = 0.99 recovers it.
  Reproducer: `key = mhKey[11, 32, 24]; m = mhMessage[2015, 32];
  mhBreak[key[[4]], m . key[[4]]]` returns `-1`. Mathilda's `LatticeReduce` has no
  δ option (WL's has none either), so a user cannot strengthen the reduction.
  The sizes above are chosen where every reducer succeeds. Squared norms of
  the shortest reduced q-ary vector also differ: Mathilda/sympy 3.24·10⁶ vs
  PARI 2.21·10⁶ at n=30.
- `app: Merkle-Hellman decrypt` (plain bignum arithmetic and a `Do` loop) runs
  8.8× slower than pure Python. This is interpreter overhead, comparable to
  Mathematica.
