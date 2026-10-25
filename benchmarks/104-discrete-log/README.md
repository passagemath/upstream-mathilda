# 104 — discrete logarithms (group F)

Generic discrete-log algorithms in (Z/pZ)*, written in the language (`app:`)
and as each system's builtin (`prim:` — WL's `MultiplicativeOrder[g, p, {h}]`,
`sympy.ntheory.discrete_log`, PARI `znlog`). Instances: prime-order-r
subgroups with r = `primeOf(s, bits)` (24 and 32 bits), p = 2kr+1 the first
prime of that shape, g = 2^((p−1)/r), and h = g^x with x = `bigOf(t, 64) mod r`
planted; plus a 131-bit smooth prime p with p−1 = 2³·3·q₁…q₈ (q_i the first
primes above 40000, 45000, …, 75000), g its least primitive root, x =
`bigOf(t, 128) mod (p−1)`. Each check is the sum of recovered exponents mod
10⁹+7 and must equal the planted sum (24-bit 39973122, 32-bit 595305666,
smooth 250044235).

| case | what |
|---|---|
| `app: BSGS 24-bit subgroup (incremental table) x4` | baby steps inserted one at a time: `tab[y] = i` / dict / `mapput` |
| `app: BSGS 24-bit subgroup (bulk table) x4` | table built in one call: `AssociationThread` / `dict(zip)` / `Map(matrix)` |
| `app: BSGS 32-bit subgroup (bulk table) x2` | m = 2¹⁶ baby steps |
| `app: Pollard rho 32-bit subgroup x2` | 3-partition walk (x mod 3), Floyd, restart from g^t on degenerate collision |
| `app: Pohlig-Hellman 131-bit smooth p x10` | per prime power digit-by-digit with bulk BSGS, `ChineseRemainder` / sympy `crt` / `chinese` |
| `prim: discrete log 24-bit subgroup x40` | 4 instances x10; `.m` uses `benchOnce` (>5 s in both CAS) |
| `prim: discrete log 32-bit subgroup x10` | 2 instances x5; `.m` capped by `TimeConstrained[…, 20]`, run once |
| `prim: discrete log 131-bit smooth p x40` | 10 targets x4; `.m` capped likewise |

Baselines: Python app = plain ints, built-in `pow()`, `dict`; Python prim =
`sympy.ntheory.discrete_log`; PARI = `znlog(h, Mod(g, p), r)` (order passed for
the subgroups), `Map`, `chinese`.

## Findings

**1. FIXED in v0.524.** `MultiplicativeOrder[g, p, {h}]` now uses Pohlig–Hellman
with BSGS and Pollard ρ (`src/numbertheory/nt_dlog.c`). All three prim rows
pass their checks and beat PARI on the same machine: 24-bit ×40 3.8 ms (PARI
27), 32-bit ×10 13 ms (PARI 79), 131-bit smooth ×40 32 ms (PARI 54). The
original finding follows.

**1 (original). `MultiplicativeOrder[g, p, {h}]` is an exhaustive search with a cap — it
returns UNEVALUATED for any 32-bit-or-larger order (2 CHECK-FAILs).**
`src/numbertheory/multiplicativeorder.c` scans m = 1, 2, … up to
`MO_SEARCH_CAP` = 10⁸ and declines above. No BSGS, no rho, no Pohlig–Hellman,
so even the 131-bit p whose p−1 is 75011-smooth (PARI: 1.3 ms, sympy: 2.2 ms
per log; the in-language Pohlig–Hellman in Mathilda itself: 6 ms) is left
unevaluated. Reproducer:

```wl
p = 2435453805216900026945200295720533381609;   (* p-1 = 2^3 3 40009 ... 75011 *)
MultiplicativeOrder[17, p, {PowerMod[17, 10^30, p]}]   (* unevaluated; PARI znlog: 10^30 *)
MultiplicativeOrder[64, 18113842987, {7600772646}]     (* unevaluated; 32-bit order subgroup *)
```

Below the cap it is linear: 383 ms per 24-bit log vs 0.6 ms PARI (≈640x).
(Mathematica's `MultiplicativeOrder` is also an exhaustive search — 1.3 s per
24-bit log, 353 s per 32-bit log, no answer within 20 s for the smooth p — so
the `.m` caps those two rows with `TimeConstrained`; Mathematica CHECK-FAILs
them with `$Aborted`.)

**2. Incremental hash-table insertion is super-linear.** `tab[y] = i` into an
`Association` in a loop costs 3170 ms for 4x4097 inserts, vs 33 ms for the same
BSGS with `AssociationThread` (≈96x), 31 ms PARI, 54 ms Mathematica. Measured
scaling: 2 500 inserts 0.34 s, 5 000 1.39 s, 10 000 7.2 s (≈O(n²·³));
`AssociateTo` the same; DownValue memo `dv[i] = i` is also quadratic (10 000:
1.34 s; Mathematica 100 000: 0.10 s). Reproducer:
`tab = <||>; Do[tab[i] = i, {i, 0, 9999}] // AbsoluteTiming` (Mathilda 7.2 s;
Mathematica 100 000 inserts in 0.21 s). The incremental 32-bit BSGS (65 537
inserts) is not run at all because it would take minutes.

**3. Pollard rho in-language is 3.5x PARI, 17x Python** — per-iteration
interpreter cost of `Switch`/`Mod`/list-tuple updates on 34-bit integers.

Measured ms (Mathilda / Python / PARI / Mathematica):

| case | Mathilda | Python | PARI | Mathematica |
|---|---|---|---|---|
| BSGS 24-bit incremental x4 | 3170.5 | 2.8 | 31 | 54.2 |
| BSGS 24-bit bulk x4 | 32.9 | 3.4 | 29 | 21.1 |
| BSGS 32-bit bulk x2 | 304.9 | 40.8 | 336 | 646.3 |
| Pollard rho 32-bit x2 | 1085.2 | 62.4 | 310 | 1597.8 |
| Pohlig-Hellman 131-bit x10 | 62.7 | 17.3 | 32 | 72.3 |
| prim 24-bit x40 | 15316.0 | 39.7 | 24 | 52231.2 |
| prim 32-bit x10 | unevaluated (28.3) | 171.1 | 66 | `$Aborted` at 20 s |
| prim smooth 131-bit x40 | unevaluated (114.4) | 87.3 | 52 | `$Aborted` at 20 s |
