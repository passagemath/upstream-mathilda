# 119 — Batch GCD over RSA moduli (group F)

Heninger et al.'s weak-RSA-key scan: product tree P = ∏ N_i, remainder tree
z_i = P mod N_i², then g_i = gcd(N_i, z_i / N_i); N_i is weak iff g_i > 1.
Written in each language (`app:`) for M = 256 and 1024 1024-bit moduli, its two
tree phases timed separately at M = 1024, against the naive all-pairs GCD at
M = 256; the big-integer kernels underneath are timed as builtins (`prim:`).
N_i = P_i Q_i with P_i = `primeOf[5000+2i, 512]`, Q_i = `primeOf[5001+2i, 512]`,
except Q_i = P_(i−17) when 37 | i — planting 12 weak moduli in the first 256 and
54 in 1024. Checks are `weak·10^10 + (Σ recovered g_i mod 10⁹+7)`.

| case | what |
|---|---|
| `prim: Times@@ 1024 moduli (1M-bit product) x4` | `Times @@` / `reduce(gmpy2.mul)` / `vecprod` |
| `prim: GCD 1024-bit pairs x1024x8` | the final `GCD[N_i, z_i/N_i]` vector / `gmpy2.gcd` / `gcd` |
| `app: product tree M=1024 x4` | pairwise products bottom-up |
| `app: remainder tree mod N^2 M=1024` | `Mod[parent, child^2]` top-down |
| `app: batch GCD M=256 x4` / `app: batch GCD M=1024` | the whole scan |
| `app: naive pairwise GCD M=256 (32640 GCDs)` | all-pairs GCD for comparison |

Baselines: Python app = plain ints (`math.gcd`, `*`, `%`, `//`); Python prim =
gmpy2. PARI = plain GP integers for the app, `vecprod` / `gcd` for the prims.
The `.m` pairs levels with `Partition[#, 2]` (M is a power of two); Python / GP
use the general odd-carry pairing, which coincides here.

## Findings

Measured (ms): Mathilda / Python / PARI / Mathematica.

| case | Mathilda | Python | PARI | Mathematica |
|---|---|---|---|---|
| prim: Times@@ 1024 moduli (1M-bit product) x4 | 333.5 | 291.4 | 26 | 26.9 |
| prim: GCD 1024-bit pairs x1024x8 | 62.4 | 53.4 | 58 | 62.9 |
| app: product tree M=1024 x4 | 263.4 | 441.4 | 30 | 27.1 |
| app: remainder tree mod N^2 M=1024 | 71.2 | 2716 | 68 | 73.8 |
| app: batch GCD M=256 x4 | 54.0 | 775.9 | 48 | 55.2 |
| app: batch GCD M=1024 | 88.6 | 2833 | 81 | 87.0 |
| app: naive pairwise GCD M=256 | 277.7 | 346.7 | 245 | 279.8 |

1. **Hang / memory blow-up — `Partition[list, 2, 2, 1, {}]` is unevaluated**,
   so the natural ragged product-tree step
   `NestWhileList[Times @@@ Partition[#, 2, 2, 1, {}] &, ns, Length[#] > 1 &]`
   never terminates (`Times @@@` multiplies the arguments of the unevaluated
   `Partition`, and the loop grows numbers without bound; the first draft was
   OOM-killed at 120 s). Reproducer: `Partition[{6, 35, 143, 323}, 2, 2, 1, {}]`
   → Mathilda leaves it unevaluated; Mathematica gives `{{6, 35}, {143, 323}}`
   (and `{{1, 2}, {3}}` for `Partition[{1, 2, 3}, 2, 2, 1, {}]`).
   The `.m` uses `Partition[#, 2]`, valid because M is a power of two.
2. **Slow — `Times @@ list` of big integers is sequential:** 83 ms per
   1024-term product vs 6.5 ms for PARI `vecprod` / Mathematica (13×); a
   balanced product tree in the language is only ~8 ms, so `Times` with many
   big-integer arguments should multiply pairwise.
3. **Slow — reassigning a big ragged structure inside `Do`/`For`.** With `t`
   holding the 1024-leaf product tree, `Do[t = 1, {1}]` costs ~0.26 s, while
   `t = 1` at top level or inside `Table` costs ~0.1 ms. This is why
   `app: product tree M=1024 x4` reads 263 ms (≈66 ms per pass) although one
   tree build is ~8 ms. Reproducer:
   ```
   t = NestWhileList[Times @@@ Partition[#, 2] &, Table[2^1023 + i, {i, 1024}], Length[#] > 1 &];
   AbsoluteTiming[Do[t = 1, {1}]]   (* ~0.27 s; same with For; Table: 0.0001 s *)
   ```
   Flat lists of bignums, a single huge integer, or a single tree level do
   not trigger it.
4. Otherwise correct, and at parity with PARI / Mathematica on the remainder
   tree, the GCD vector and the full scan.
