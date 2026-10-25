# 105 — elliptic-curve point arithmetic on secp256k1 (group F)

secp256k1 (y² = x³ + 7 over the 256-bit prime p, from `crypto_data`) point
arithmetic written in each language (`app:`): affine add / double with the slope
inverse via `PowerMod[x, -1, p]` / `pow(x, -1, p)` / `Mod(x, p)^-1`, the
inversion-free Jacobian add (add-1998-cmo-2) and double (dbl-2009-l, a = 0), and
scalar multiplication k·G by left-to-right double-and-add (affine and Jacobian) and
by the Montgomery ladder (Jacobian), for 64 scalars `bigOf[s, 256] mod n`. The
`prim:` cases are a native curve group: PARI's `elladd` / `ellmul` on
`ellinit([0,0,0,0,7], p)`. Every scalar-mult case carries the same check (Σ x(k·G)
mod 10⁹+7), so the three in-language algorithms are cross-verified against PARI's
`ellmul` (and against OpenSSL in Python).

| case | what |
|---|---|
| `app: affine point add x10000` | Q ← Q + G from Q = 2G (one inversion each) |
| `app: Jacobian point add x10000` | same chain, Jacobian (no inversions; one at the end) |
| `app: affine point double x10000` / `app: Jacobian point double x10000` | Q ← 2Q, from G |
| `app: scalar mult affine double-and-add x64` | k·G, affine |
| `app: scalar mult Jacobian double-and-add x64` | k·G, Jacobian |
| `app: scalar mult Jacobian Montgomery ladder x64` | k·G, ladder (one add + one double per bit) |
| `prim: point add x10000` | PARI `elladd`; **ABSENT** in Mathilda (no EC head in WL); Python SKIP (no installed library exposes a bare point add — coincurve not installed) |
| `prim: scalar mult k*G x64` | PARI `ellmul`; Python `cryptography` (OpenSSL) `ec.derive_private_key(k, SECP256K1())` — includes key-object construction; **ABSENT** in Mathilda |

The brief's "x1000 point additions" micro-case is x10000 here (and scalar mult
x64) so every PARI case clears its 20 ms floor. Baselines: Python = plain ints +
built-in `pow`; PARI = plain integers with `%`, `Mod(x,p)^-1`, `binary(k)`.

## Findings

Measured (ms; Mathilda / Python / PARI / Mathematica):

| case | Mathilda | Python | PARI | Mathematica |
|---|---|---|---|---|
| app: affine point add x10000 | 209 | 216 | 43 | 151 |
| app: Jacobian point add x10000 | 662 | 53 | 55 | 369 |
| app: affine point double x10000 | 211 | 216 | 45 | 140 |
| app: Jacobian point double x10000 | 380 | 36 | 39 | 208 |
| app: scalar mult affine double-and-add x64 | 621 | 532 | 112 | 386 |
| app: scalar mult Jacobian double-and-add x64 | 1927 | 104 | 110 | 683 |
| app: scalar mult Jacobian Montgomery ladder x64 | 4087 | 179 | 160 | 1024 |
| prim: point add x10000 | ABSENT | SKIP | 38 | ABSENT |
| prim: scalar mult k*G x64 | ABSENT | 98 | 30 | ABSENT |

All checks agree in all four systems — no wrong answers.

1. **ABSENT: no elliptic-curve group head** (WL has none either). PARI's `ellmul`
   is 20× faster than the best in-language Mathilda scalar mult.
2. **Module temporaries are never reclaimed (memory growth + progressive
   slowdown).** Every `Module` call leaves its `x$NNN` symbols in the symbol table
   for good; Mathematica removes them. The Jacobian code (≈12 locals per add) adds
   hundreds of thousands of symbols per case, and a repeated scalar-mult batch got
   ~10% slower on each re-run in one process (1.12 s → 1.84 s over six runs), which
   is why the Jacobian rows are 4–9× Mathematica while the affine rows (few locals,
   cost dominated by `PowerMod` inverse) are at Python parity. Reproducer:
   ```
   g[x_] := Module[{a = x}, a];
   Do[g[3], {50000}];
   Print[Length[Names["a$*"]], " ", MemoryInUse[]]
   (* Mathilda: 50000 and +19 MB, still growing ~250 B/call; Mathematica: 1, flat *)
   ```
3. **Per-operation interpreter overhead on 256-bit `Mod` arithmetic.** Even
   discounting finding 2, a Jacobian add (≈15 `Mod` products) costs ~66 µs in
   Mathilda against ~5 µs in Python and ~37 µs in Mathematica; the inversion-free
   representation, a 4× win in Python, is a 3× *loss* in Mathilda.
