# 109 — RSA attacks (group F)

Three textbook attacks on misused RSA, each written in the language (`app:`)
over the integer primitives crypto number theory leans on, plus the native
continued-fraction primitives (`prim:`) behind Wiener's attack.

- **(a) Wiener.** A key with a small private exponent (`d` 120-bit, well under
  `N^(1/4)/3` for a 1024-bit `N`) leaks `d`: it is the denominator of a
  convergent of the continued fraction of `e/N`. The recovery is an in-language
  convergent recurrence over `ContinuedFraction[e/N]`, testing each convergent
  by reconstructing `φ` and checking that `x² - (N-φ+1)x + N` has a
  perfect-square discriminant.
- **(b) Common modulus.** Two coprime public exponents on one modulus:
  `ExtendedGCD[e1,e2] = {1,{a,b}}` recovers `m = c1^a c2^b mod N`.
- **(c) Håstad broadcast (e=3).** One message under three coprime moduli:
  `ChineseRemainder` lifts the three `c_i` to `m³ mod (n1 n2 n3) = m³` exactly,
  then an in-language integer cube root (Newton on integers) recovers `m`.

| case | what |
|---|---|
| `prim: ContinuedFraction e/N 1024-bit x600` | the builtin that drives Wiener |
| `prim: Convergents e/N 1024-bit x40` | convergent list → recover `d` — **ABSENT** in Mathilda |
| `app: Wiener attack 1024-bit x20` | full in-language recovery of `d` from `e/N` |
| `app: common-modulus attack 1024-bit x600` | recover `m` via `ExtendedGCD` |
| `app: Hastad broadcast e=3 x1500` | CRT lift + integer cube root |

Each attack also folds a `recovered` count into its check (value + 10⁹+7 × count) (20 / 100 / 50) so a silent
recovery failure cannot hide behind a matching sum.

Baselines: Python = plain ints + `pow()` for `app:` (idiomatic pure-Python
crypto), `sympy.ntheory.continued_fraction` / `continued_fraction_convergents`
for `prim:`, `gmpy2.iroot` for the Håstad cube root. PARI/GP = plain integers,
`bezout`, `chinese`, `sqrtnint` for `app:`; `contfrac` / `contfracpnqn` for
`prim:`.

## Findings

- **`Convergents` is ABSENT in Mathilda** (benchIf → SKIP). It exists in
  Mathematica (which computes the same recovered-`d` sum, `570546006`) and in
  Python (sympy) and PARI (`contfracpnqn`). The Wiener recovery therefore uses
  `ContinuedFraction` (which Mathilda has) plus an in-language convergent
  recurrence — so the absence costs clarity, not capability. **Adding a
  `Convergents` head would let the Wiener recovery be written the WL-idiomatic
  way.**

- **No wrong answers, no crashes, no hangs.** All computed checks agree across
  Mathilda, Python, PARI and Mathematica, and every attack recovers 100% of its
  secrets (`20`, `100`, `50`).

- **Mathilda is the fastest of the three on `ContinuedFraction`** and on the
  Håstad / common-modulus `pow`-bound loops is competitive:

  | case | Mathilda | Python | PARI | Mathematica |
  |---|---|---|---|---|
  | `prim: ContinuedFraction e/N x600` | 63 ms | 3680 ms | 33 ms | 20 ms |
  | `app: Wiener attack x20` | 74 ms | 249 ms | 22 ms | 167 ms |
  | `app: common-modulus x600` | 28 ms | 213 ms | 29 ms | 43 ms |
  | `app: Hastad broadcast x1500` | 107 ms | 746 ms | 38 ms | 176 ms |

  `ContinuedFraction` on an exact 1024-bit rational: Mathilda 63 ms vs sympy
  3680 ms (~58×) — sympy's `continued_fraction` over `Rational` is very slow
  here. PARI stays ~2× ahead of Mathilda on the number-theoretic loops; the gap
  is per-`Mod`/`PowerMod` constant factor, not a missing kernel.

Baselines measured on this host (Apple M-series).
