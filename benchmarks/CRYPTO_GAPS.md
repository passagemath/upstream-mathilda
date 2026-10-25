# Group F: cryptographic number theory. Results and improvement roadmap

Twenty experiments (`benchmarks/100-…` to `119-…`) build cryptographic
applications in four systems: Mathilda, Mathematica 13 (the same `.m` file),
Python 3.11 (gmpy2 / sympy / hashlib / pycryptodome) and PARI/GP 2.17. Each
experiment's README describes its cases and records its own findings with
reproducers. This file turns those findings into a ranked work queue.

Run: `HPC_PYTHON=/usr/local/bin/python3.11 python3 benchmarks/run_all.py --only 100-119`
(2026-10-21, ~18 min, Apple M-series). Raw data:
`results/2026-10-21-groupF-partial.json`.

## Headline

149 cases. Coverage of the WL heads the experiments declared: 72.5% (37 of 51).

| class | cases | meaning here |
|---|---|---|
| AHEAD | 52 | Mathilda fastest or within 1.5× of the best of Mathematica / Python / PARI |
| SLOWER | 67 | at least 1.5× behind the best other system |
| ABSENT | 21 | the WL head (or form) a protocol needs does not exist in Mathilda |
| CHECK-FAIL | 9 | **Mathilda computed a different answer.** Every one is a Mathilda bug (C1, C2, C3, C5, C8 below) |
| INCOMPLETE | 0 | |

Median Mathilda time divided by each other system's time, over cases both
timed:

| baseline | median ratio |
|---|---|
| Mathematica | 0.81× |
| Python | 1.12× |
| PARI/GP | 1.50× |

The typical group-F case is faster than Mathematica, level with
gmpy2-backed Python, and 1.5× behind PARI. The tail is long.

**Where Mathilda already stands well.** The GMP-backed kernels match GMP:
- `PowerMod`, `PrimeQ`, `NextPrime`, `JacobiSymbol` and `GCD` are within
  1.0–1.3× of gmpy2.
- `FactorInteger` on balanced semiprimes is the fastest of the four systems
  (160-bit: Mathilda 81 ms, PARI 232 ms, Mathematica 2.7 s, sympy 51 s).
- Whole protocols are bignum-bound rather than interpreter-bound: RSA-CRT,
  DH / MODP-2048, ElGamal, Paillier, BBS, Goldwasser–Micali, batch GCD.
  Mathilda is level with PARI on these or ahead of it, and 5–10× ahead of pure
  Python.

**Where it does not.** There are three clusters:
1. **Correctness:** eight bugs, five of them silent wrong answers.
2. **Absent primitives:** bitwise operations, cryptographic hashes and modular
   linear algebra block whole protocols.
3. **Missing fast paths:**
   - integer `ListConvolve`;
   - incremental `Association` insertion;
   - discrete log;
   - packed integer arrays through `Quotient` / `Map` / `Dot`;
   - LLL;
   - `Sqrt` of machine-size integers.

## Tier 0: correctness. Fix first

All eight reproduce on the current binary (re-checked by the lead on
2026-10-21). Each one silently turns a correct program into a wrong answer or a
hang, so each needs a regression test in `tests/`.

| # | bug | reproducer | found by | where to start |
|---|---|---|---|---|
| C1 | `Return[]` inside `While` / `For` in a function body exits only the loop | `f[s_] := Module[{r = 1}, While[r < s, If[r == 1, Return[True]]; r++]; False]; f[3]` gives `False` (WL gives `True`) | 102: four CHECK-FAILs (Miller–Rabin, BPSW) | `src/iter.c` (`While`), the Return sentinel in `src/eval.c` |
| C2 | `AppendTo` / `PrependTo` on a **packed** list silently does nothing | `w = Range[8]; AppendTo[w, 0]; w` gives `{1,…,8}` | 114: AES key expansion CHECK-FAIL | `src/internal.c` (AppendTo), packed gate in `src/pack.c` |
| C3 | `MapThread` stays unevaluated when a packed array arrives by substitution | `r = Table[i+j,{i,2},{j,3}]; (MapThread[g, {#, #}] &)[r]` | 117: Kyber MapThread CHECK-FAIL | `src/funcprog.c`, packed transparency |
| C4 | `Partition[list, n, d, kL, {}]` (pad-with-nothing overhang) stays unevaluated, so the natural product-tree loop never terminates | `Partition[{6,35,143}, 2, 2, 1, {}]` | 119 | `src/list/partition.c` |
| C5 | `CoefficientList[p, x, n]` (padded length) stays unevaluated | `CoefficientList[1 + 2 x, x, 4]` should be `{1,2,0,0}` | 116 CHECK-FAIL | `src/poly/poly.c` |
| C6 | `Quiet` does not suppress `Inverse::method`, and `Thread::tdlen` goes to **stdout** | `Quiet[Inverse[{{1,2},{3,4}}, Modulus -> 7]]` | 117 | message funnel (`src/message.h`); these sites bypass it. Extend `make check-messages` |
| C7 | `Module` temporaries (`a$NNN`) are never removed: unbounded memory, and a process gets slower the longer it runs | `g[x_] := Module[{a = x}, a]; Do[g[3], {20000}]; Length[Names["a$*"]]` gives 20000 (WL gives ~1) | 105, 106: the identical "tampered verify" batch costs 2× because it runs later | `Module` in `src/internal.c`, symbol-table removal of `Temporary` symbols |
| C8 | **FIXED v0.524.** `MultiplicativeOrder[g, p, {h}]` (discrete log) returns unevaluated once the order exceeds a linear-scan cap | `MultiplicativeOrder[64, 18113842987, {7600772646}]` | 104: two CHECK-FAILs | `src/numbertheory/multiplicativeorder.c` (`MO_SEARCH_CAP`). See P3 |

C7 should go before any interpreter-overhead work (P7). It inflates every
`app:` timing that calls a `Module` function in a loop, so P7 cannot be
measured honestly until C7 is fixed.

## Tier 1: absent heads and forms (feature work, never a ratio)

Ordered by how many cases each one unblocks and how central it is to
cryptographic code.

| # | head / form | cases blocked | best other system | notes |
|---|---|---|---|---|
| A1 | **`BitXor`, `BitAnd`, `BitOr`, `BitNot`, `BitShiftLeft`, `BitShiftRight`** | 7 (`114` AES ×3, `115` SHA-256 / HMAC / ChaCha20 ×4) | Python 7–10 ms; PARI 20–40 ms | In every other system the no-bitwise variants cost about 7× the bitwise ones, so Mathilda's SHA-256 should drop from ~570 ms to under 100 ms. They are integer element-wise heads, so CLAUDE.md requires an exact-int kernel in `src/ndinteger.c` (+ `AWARE`/`INT64_OK`) and a `Compile[]` lowering. Bignum operands via GMP `mpz_xor`/`mpz_and`/… |
| A2 | **`Hash[expr, "SHA256" \| "SHA" \| "MD5" \| "SHA512" …, "HexString" \| "Integer" \| "ByteArray"]`** | 2 (`106`, `115`) | hashlib 1 ms | Only bare `Hash[s]` evaluates today. The same primitive would also feed ECDSA and HMAC. |
| A3 | **`Modulus -> p` for `Inverse`, `LinearSolve`, `RowReduce`, `Det`, `NullSpace`** | 2 (`117`) | PARI 18–28 ms (`matsolvemod`) | Today they fail with `::method` (see C6). Gaussian elimination over F_p on int64 when p < 2³¹, bignum otherwise. Lattice and LWE work needs these. |
| A4 | `ModularInverse[a, m]` | 1 (`100`) | PARI 14 ms | Trivial wrapper (`mpz_invert`). `PowerMod[a, -1, m]` already works. |
| A5 | `KroneckerSymbol[a, n]` | 1 (`118`) | gmpy2 87 ms | `mpz_kronecker`. `JacobiSymbol` is already at GMP speed. |
| A6 | `Convergents[x \| cf]` | 1 (`109`) | PARI 27 ms | Short: `ContinuedFraction` is already fast. |
| A7 | `ProvablePrimeQ` (`PrimalityProving` package) | 1 (`102`) | PARI 141 ms (APR-CL) | ECPP or APR-CL certificate. Larger job. |
| A8 | `16^^ff` base-literal syntax | (parser) | — | Every hex constant in the group F `.m` files goes through `FromDigits[…, 16]` because of this. `src/parse.c`. |
| A9 | `Encrypt` / `GenerateSymmetricKey`, elliptic-curve group heads | 4 (`105`, `106`, `114`) | PARI `ellmul`, pycryptodome | Low priority. Mathematica has no elliptic-curve arithmetic heads either, and `Encrypt` is a large surface. |

`SqrtMod` was probed in `112` and **is not a gap**: it is not a WL function,
and Mathematica leaves it unevaluated too. Mathilda supports the WL spelling,
`PowerMod[a, 1/2, p]` (see P10 for its speed).

## Tier 2: performance, ranked by ratio to the best other system

| # | gap | worst cases (Mathilda vs best) | root cause (as measured) | fix | where |
|---|---|---|---|---|---|
| P1 | **Integer `ListConvolve` has no fast path** | 116 n=1024: 1410 ms vs Mathematica 0.2 ms (**6900×**). 117 Kyber keygen+enc+dec: 2.3 s vs 2.5 ms (**890×**) | Integer inputs fall back to the boxed path, which takes 0.36 s per n=1024 convolution against 3 ms for the same call on reals. The result also comes back unpacked | int64 direct / Karatsuba convolution with an overflow guard (the exact bound \|a\|·\|b\|·n < 2⁶³ is cheap to check up front), and a bignum fallback | `src/convolutions.c`; `AWARE`/`INT64_OK` in `src/pack.c`; Compile lowering |
| P2 | **`Association` insertion `tab[k] = v` is super-linear** | 104 BSGS with an incrementally built table: 2846 ms vs Python 2.9 ms (**970×**) | 2500 inserts take 0.27 s, 5000 take 1.35 s, 10000 take 7.2 s. Each `tab[k] = v` appears to rebuild or copy the association. `AssociateTo` and DownValue memo `dv[i] = i` behave the same | In-place update when the symbol's OwnValue is the only reference (refcount 1), as the O(1) lookup work already does for reads | assignment to `Part`/`Key` of an OwnValue in `src/eval.c` / `src/assoc*.c`. Gate it with the `tests/bench_assoc.c` scaling check |
| P3 | **FIXED v0.524**: Pohlig–Hellman + BSGS/ρ in `nt_dlog.c`; prim rows now 3.8 / 13 / 32 ms vs PARI 27 / 79 / 54. Was: **Discrete log is a capped linear scan** | 104 prim 24-bit: 14.5 s vs PARI 23 ms (**630×**). 32-bit and 131-bit smooth p: unevaluated (C8) | `MultiplicativeOrder[g, p, {h}]` steps g^k one at a time up to `MO_SEARCH_CAP` | Pohlig–Hellman over the factored group order, with BSGS for small prime factors and Pollard ρ for large ones (PARI's `znlog` scheme). Remove the cap | `src/numbertheory/multiplicativeorder.c` |
| P4 | **Packed integer arrays are lost or never formed** | 117 `Mod[A.s+e, q]` 512²: 774 ms vs Mathematica 8.6 ms (**90×**). `R.A` 256×512×512: 95 ms vs 0.6 ms (**150×**). Regev LWE: 32× | `Quotient` / `Mod` over a list of bignums gives an unpacked list even when every result is machine-size. `Map` over a packed matrix hands each row to the function unpacked (`Map[PackedArrayQ, Partition[Range[12], 4]]` gives `{False,False,False}`). That unpacks the 116 NTT loop (937 ms via `/@` against 98 ms via `Table`/`Part`) | Repack a result list that came out all machine integers; pass packed sub-arrays to the function in `Map`; an int64 `Dot` kernel with an overflow check | `src/pack.c`, `src/funcprog.c`, `src/ndarray.c` |
| P5 | **`LatticeReduce` is exact-rational LLL with a fixed δ = 3/4** | 113: 6.5–79× behind PARI `qflll` (5–10× behind Mathematica) | Exact rational Gram–Schmidt. δ cannot be set, so at density 0.57 the attack misses 1 ciphertext in 20 that PARI (δ = 0.99) recovers | Floating-point L² LLL (Nguyen–Stehlé) with exact fallback; accept a δ option (WL has none, so a `Method` suboption) | `src/linalg/latticereduce.c` |
| P6 | **`Sqrt` of a machine-size non-square integer costs ~2.8 ms** | 110 Fermat: 3.6 s vs Mathematica 119 ms (30×; PARI and Python under 1 ms) | The radical-simplification path for a ~46-bit integer costs 700× its own `FactorInteger` (4 µs). A 100-digit non-square takes 3 µs | Find the slow square-part extraction for word-size inputs (probably a trial-division loop); add an `mpz_perfect_square_p` early exit | the integer `Power[n, 1/2]` canonicaliser (`src/power.c`) |
| P7 | **Interpreter overhead in bignum app loops** | 105 Jacobian Montgomery ladder: 3.9 s vs PARI 159 ms (25×). ECDSA sign/verify: 10–20×. Pollard p−1: 21×. NTT n=1024: 47×. Shamir split: 10× | Each modular operation on a 256-bit number costs ~30–60 µs in Mathilda, against ~5 µs in Python and PARI. Avoiding inversions makes Python 4× faster but Mathilda 3× slower, so per-call evaluation overhead dominates. C7 adds to it | Profile `Module` entry/exit, `Set` on locals and `Mod`/`Times` dispatch for small bignums; rule out C7 first. Bignum ops in `Compile[]` would remove the cliff for these loops entirely | `src/eval.c`, `src/internal.c`; see the evaluator-allocation memory notes |
| P8 | **`Times @@ {many bigints}` multiplies left to right** | 119: 329 ms vs PARI 26 ms (12.7×). Product tree: 10× | O(n²) in total limb work. A pairwise product tree written in Mathilda is ~8 ms. Assigning a large product tree to a variable inside `Do` is also slow (~0.27 s vs 0.1 ms at top level) | Balanced binary product (`mpz` tree) when `Times` has more than ~8 bignum arguments. Investigate the `Do`/`Set` copy | `src/times.c`; `Set` in `Do` |
| P9 | **`InterpolatingPolynomial` works over Q only** | 111 n=200: 982 ms vs PARI 45 ms (22×) | Rational denominators grow; there is no `Modulus` option | `Modulus -> p` (Newton divided differences over F_p) | `src/interp.c` |
| P10 | Small-operand bignum builtins | `PowerMod[a, 1/2, p]` 4.5× PARI. `PrimitiveRoot` 2.9× PARI and **13× a four-line Mathilda loop**. Garner recombination 4.3×. `ContinuedFraction` 3.3× (Mathematica) | `PrimitiveRoot` re-factors p−1 and tests candidates slowly. The others are per-call overhead | `PrimitiveRoot`: factor once, test g^((p−1)/q) with early exit. Then profile the call overhead | `src/numbertheory/` |
| P11 | `PolynomialRemainder[p, x^n + 1, x]` is ~n³ | 116: n=128 takes 1.6 s. n=1024 did not finish in 9 minutes | A generic multivariate remainder applied to a sparse binomial modulus | A fast path when the divisor is xⁿ ± c: fold the coefficients directly | `src/poly/` |

## Suggested sequencing

1. **C1–C8**, each with a regression test. C1, C2, C3 and C7 are the most
   likely to break user programs outside cryptography.
2. **A1 bitwise family** (with ND and Compile surfaces) and **A2 `Hash`**.
   Together they unblock 9 cases and the entire symmetric-crypto half of the
   group.
3. **P1 integer `ListConvolve`** and **P4 packing**: the largest ratios, and
   both sit on the packed-array substrate the rest of Mathilda already
   invests in.
4. **P2 `Association` insertion** and **P3 discrete log**. P3 also fixes C8.
5. **A3 modular linear algebra** and **P5 floating-point LLL**: the lattice
   and LWE stack.
6. **P6**, **P8**, **P9**, **P10**, **P11**, then **P7** once C7 is fixed and
   the measurement can be trusted.

After each step, re-run `run_all.py --only 100-119` and diff
`REPORT.partial.md` against this file's numbers.

## Harness notes

- `crypto_data.{m,py,gp}` is a linear LCG, so `bigOf[seed + i]` words form
  exact arithmetic progressions. A **lattice** built from them has a planted
  short vector. `113` passes its inputs through a `3^e mod P` mixer for that
  reason; any future lattice experiment must do the same.
- PARI's `getabstime()` has 1 ms resolution. A few cases in 116 and 117 sit
  under the 20 ms floor because a shared loop count large enough for PARI
  would push Mathilda past its 5 s budget. Their ratios are lower bounds.
- `gp` abandons the rest of a file being read when `parisizemax` is raised, so
  `harness.gp` raises it last.
