# 117 — Module-LWE / Kyber-lite (group F)

The integer-array kernels that lattice cryptography leans on, at toy-but-real
sizes (q = 3329, n = 256, R_q = Z_q[x]/(x^256+1), no NTT — 116 covers that):
the negacyclic ring product, plain LWE on a 512×512 matrix, and the
modular-matrix primitives, timed as builtins (`prim:`). Then a whole
Kyber-lite (keygen t = A s + e; encrypt u = Aᵀr + e1, v = t·r + e2 + 1665 m;
decrypt by rounding v − s·u) and a Regev bit-encryption LWE scheme, written in
each language (`app:`). Matrices are uniform mod q and noise is centered
binomial (η = 2), all from `crypto_data`'s `smallList`, built outside the
timed region. App checks are `correct·10^10 + hash` (correct must equal the
total bit count: 2048 for Kyber x8, 256 for Regev).

| case | what |
|---|---|
| `prim: negacyclic polymul n=256 x200` | `ListConvolve[a, b, {1,-1}, 0]` folded with a sign flip / `np.convolve` / PARI `Mod(FpX, x^256+1)` product |
| `prim: LWE Mod[A.s+e,q] 512x512 x200` | matrix–vector mod q on an integer matrix |
| `prim: LWE batch Mod[R.A,q] 256x512x512` | the batched encryption product |
| `prim: Inverse[M,Modulus->q] 128x128 x4` | modular matrix inverse — **probe-SKIP** in Mathilda |
| `prim: LinearSolve[M,b,Modulus->q] 128x128 x16` | modular linear solve — **probe-SKIP** in Mathilda |
| `app: Kyber-lite k=2 / k=3 keygen+enc+dec x8` | rank-2 / rank-3 module LWE, 8 × 256-bit messages; ring sums spelled with `Sum` |
| `app: Kyber-lite k=2 MapThread spelling x8` | the same scheme with the ring sums spelled `Total[MapThread[pmul, …]]` — **CHECK-FAIL** in Mathilda |
| `app: Regev LWE 512x512 enc+dec 256 bits` | one binary r per bit: u = r·A, v = r·b + 1665 μ |

Baselines: Python = numpy int64 (`np.convolve`, `@`, `%`) — how lattice crypto
is prototyped in Python — and sympy `DomainMatrix` over `GF(q)` for the
inverse / solve (pure Python: python-flint is not installed, so those two run
`bench_once`). PARI/GP = POLMOD-over-INTMOD ring elements (`Mod(Pol(Mod(c,q)),
x^256+1)`), integer `matrix * col % q`, `(M*Mod(1,q))^-1`, `matsolve`.
Several PARI rows (polymul, Kyber) are below the 20 ms floor: PARI is ~500×
faster than Mathilda there, and the Mathilda side is already at its 5 s cap.

## Findings

Measured (ms): Mathilda / Python / PARI / Mathematica.

| case | Mathilda | Python | PARI | Mathematica |
|---|---|---|---|---|
| prim: negacyclic polymul n=256 x200 | 5780 | 8.1 | 7 | 10.4 |
| prim: LWE Mod[A.s+e,q] 512x512 x200 | 1078 | 33.5 | 1367 | 9.0 |
| prim: LWE batch Mod[R.A,q] 256x512x512 | 131 | 172 | 334 | 0.65 |
| prim: Inverse[M,Modulus->q] 128x128 x4 | SKIP | 21872 | 18 | 107 |
| prim: LinearSolve[M,b,Modulus->q] 128x128 x16 | SKIP | 19544 | 28 | 185 |
| app: Kyber-lite k=2 x8 | 2844 | 2.7 | 10 | 3.9 |
| app: Kyber-lite k=3 x8 | 4758 | 4.9 | 15 | 5.1 |
| app: Kyber-lite k=2 MapThread x8 | CHECK-FAIL | 2.7 | 9 | 3.0 |
| app: Regev LWE 512x512 256 bits | 1346 | 97.6 | 1613 | 40.2 |

1. **Wrong answer — `MapThread` over a substituted packed array returns
   unevaluated.** When a packed rank-2 array reaches `MapThread` through a
   rule/pure-function/`With` substitution, `MapThread` stays unevaluated (the
   same call on the global symbol works). This breaks the natural
   `Total[MapThread[pmul, {a[[i]], s}]]` Kyber spelling (CHECK-FAIL, check
   emits −1). Minimal reproducer:
   ```
   r = Table[i + j, {i, 2}, {j, 3}];
   (MapThread[g, {#, #}] &)[r]
   (* Mathilda: MapThread[g, {{{2,3,4},{3,4,5}}, {{2,3,4},{3,4,5}}}]
      Mathematica: {g[{2,3,4},{2,3,4}], g[{3,4,5},{3,4,5}]} *)
   ```
   Also `k[x_] := MapThread[g, {x, x}]; k[r]` and
   `With[{x = r}, MapThread[g, {x, x}]]`. Literal `{{2,3,4},{3,4,5}}` works.
2. **Integer `ListConvolve` has no fast path:** 21–29 ms per 256×256 full
   convolution, vs 0.5 ms for the same call on `N[a]`, `N[b]` (40×) and 24 µs
   in Mathematica (~1000×). This single kernel is the whole Kyber gap
   (~900× vs Mathematica). Reproducer:
   `a = Mod[Range[256]^2, 3329]; b = Mod[Range[256]^3, 3329]; AbsoluteTiming[ListConvolve[a, b, {1, -1}, 0];]`.
3. **Absent form — `Inverse[M, Modulus -> p]` and `LinearSolve[M, b, Modulus -> p]`**
   both reject the option (`Inverse::method: Method option value is not one of
   "Automatic", ...`) and stay unevaluated; Mathematica evaluates both. Probed
   in the `.m`, emitted as SKIP.
4. **Packing lost through `crypto_data`'s `smallList`:**
   `Mod[Quotient[lcgStates[...], 2^32], m]` (a list of bignums mapped down to
   machine integers) comes back unpacked in Mathilda (Mathematica: packed), and
   `Partition`, `A + 0`, `Map[Identity, A, {2}]`, `Round[N[A]]` do not repack
   it (only `Table[A[[...]]]` does). `Mod[A.s + e, q]` on that 512×512 matrix
   costs 4.4 ms vs 0.14 ms on a packed one (30×), which is why the LWE rows
   trail Mathematica 100×+. `Mod[R.A, q]` (256×512×512) is 131 ms vs 0.65 ms in
   Mathematica.
5. **`Quiet` leaks messages.** `Quiet[Inverse[{{1,2},{3,4}}, Modulus -> 7]]`
   still prints `Inverse::method` (stderr), and `Quiet[{1,2} + {1,2,3}]` prints
   `Thread::tdlen` **to stdout, twice**. The MapThread case is wrapped in
   `Quiet` and still floods stdout with `Thread::tdlen` lines.
