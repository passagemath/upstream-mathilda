# 112 — Rabin cryptosystem & modular square roots (group F)

Square roots modulo a prime — the kernel of the Rabin cryptosystem and of
elliptic-curve point decompression. Two classical algorithms written in each
language (`app:` Tonelli–Shanks and Cipolla), both run on the hard `p ≡ 1 mod 8`
case (so `p-1 = Q·2^S` has `S ≥ 3`, past the easy `(p+1)/4` shortcut). Then the
native primitive (`prim:`): `PowerMod[a, 1/2, p]`, which is how a square root
mod p is spelled in WL. Every root is canonicalised to `min(r, p-r)` before the
check, so the four systems — each free to return either root — agree on one
integer. The Rabin case decrypts a ciphertext to its four roots mod `n = p·q`
(`p ≡ q ≡ 3 mod 4`) via CRT and sums their canonical forms.

| case | what |
|---|---|
| `app: Tonelli-Shanks 256/1024-bit p=1mod8` | in-language Tonelli–Shanks |
| `app: Cipolla 256/1024-bit p=1mod8` | in-language Cipolla (arithmetic in `F_p(√w)`) |
| `prim: PowerMod[a,1/2,p] 256/1024-bit` | WL's native modular sqrt spelling |
| `prim: SqrtMod 256-bit x500` | `SqrtMod` head — **ABSENT** in Mathilda (benchIf → SKIP) |
| `app: Rabin decrypt 512-bit n x500` | four roots of `c` mod `n=pq` via CRT, `p≡q≡3 mod 4` |

Baselines: Python = Tonelli/Cipolla as plain ints for `app:`; `sympy.ntheory.sqrt_mod`
for every `prim:` (gmpy2 has no modular square root). PARI/GP = Tonelli/Cipolla as
plain GP; `sqrt(Mod(a,p))` for every `prim:`.

## Findings

- **`SqrtMod` is absent in Mathilda and is not WL either.** The WL spelling of a
  modular square root is `PowerMod[a, 1/2, p]`, which **Mathilda supports and
  returns correctly** — verified root-for-root (canonical `min(r,p-r)`) against
  Tonelli–Shanks, Cipolla, sympy and PARI: all five agree on `214137979`
  (256-bit) and `187074372` (1024-bit). The `prim: SqrtMod` label is benchIf'd
  and reports SKIP in Mathilda. Note: recent Mathematica *lists* `SqrtMod` in
  `Names["SqrtMod"]` but leaves `SqrtMod[a,p]` unevaluated — so the case
  produces a garbage (unevaluated) line in the Mathematica column too; it is a
  non-function in both CASes, and the real spelling is `PowerMod[a,1/2,p]`.

- **`PowerMod[a,1/2,p]` is correct and competitive but trails PARI.**

  | case | Mathilda | Python (sympy) | PARI | Mathematica |
  |---|---|---|---|---|
  | `prim: PowerMod 256-bit x500` | 94 ms | 203 ms | 24 ms | (matches) |
  | `prim: PowerMod 1024-bit x40` | 144 ms | 391 ms | 38 ms | (matches) |

  Mathilda is ~2× faster than sympy's `sqrt_mod` and ~4× slower than PARI's
  `sqrt(Mod)`. No correctness issue.

- **In-language Tonelli–Shanks is fast in Mathilda; Cipolla is ~10–15× slower.**
  Tonelli 256-bit ×100 = 21 ms (fastest of the three; PARI 33 ms, Python 144 ms).
  Cipolla 256-bit ×100 = 306 ms in Mathilda vs 66 ms PARI / 96 ms Python — the
  `F_p(√w)` multiply is written as a 2-vector with four `Mod` bignum multiplies
  per step, and Mathilda's per-`Mod` interpreter overhead dominates the extra
  field arithmetic. A constant-factor interpreter cost on a loop-heavy pure-WL
  inner loop, not a missing capability; the result is correct.

- **Rabin decrypt (CRT, four roots) is correct and within budget** in all four
  systems (Mathilda 27 ms vs PARI 20 ms / Python 202 ms for ×500) — Mathilda is
  competitive here because the per-message work is a handful of `PowerMod`s.

- **No wrong answers, no crashes, no hangs.** Every computed check agrees across
  Mathilda, Python, PARI and Mathematica.
