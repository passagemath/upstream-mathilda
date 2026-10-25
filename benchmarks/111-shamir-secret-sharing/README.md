# 111 — Shamir secret sharing (group F)

A `(t, n)` threshold scheme over the prime field `F_p`, `p` a 256-bit prime
(`primeOf[1111, 256]`). A secret `s` is the constant term of a degree `t-1`
polynomial `f`; share `i` is `(i, f(i) mod p)`; any `t` shares recover
`s = f(0)` by Lagrange interpolation at 0. The `app:` cases are the protocol
written in each language (Horner evaluation mod p, Lagrange-at-0 with one
`PowerMod[d, -1, p]` inverse per term, and the batch form with precomputed
weights). The `prim:` cases are each system's native interpolation of the
shares, reduced mod p.

| case | what |
|---|---|
| `app: split 256-bit n=200 t=50 x20` | 20 secrets, each split into 200 shares by Horner mod p |
| `app: reconstruct Lagrange t=50 x20` | recover `f(0)` from 50 shares, one inverse per term |
| `app: batch reconstruct t=50 x2000` | 2000 secrets on fixed abscissae — weights once, then a dot product mod p each |
| `prim: InterpolatingPolynomial t=50 mod p x20` | each system's interpolation of 50 shares, `f(0) mod p` |
| `prim: InterpolatingPolynomial n=200 mod p x4` | the over-determined degree-199 interpolant, `f(0) mod p` |

Baselines: Python = plain ints + `pow(d, -1, p)` for `app:`; for `prim:`,
sympy's GF(p) Vandermonde solve (`DomainMatrix.lu_solve` over `GF(p)`) —
sympy has no GF(p) interpolation routine, and its Q-route
(`polys.polyfuncs.interpolate`) is unusable at this size (1.7 s at 30 points).
PARI/GP = plain integers + `Mod(d,p)^-1` for `app:`; `polinterpolate` on
`Mod(y, p)` entries for `prim:`, which stays in `F_p`.

## Findings

- **Mathilda's `prim: InterpolatingPolynomial` beats sympy decisively but
  trails PARI by ~10–60×.** `InterpolatingPolynomial` over the shares is the
  WL route a user reaches for (there is no `Interpolation[..., Modulus->p]`).
  It builds the Newton form over Q (exact rational coefficients, growing
  denominators) and we reduce `f(0)` mod p afterward. PARI's `polinterpolate`
  on `Mod(y,p)` entries stays in `F_p` and is a flat finite-field divided
  difference, so it is far cheaper:

  | case | Mathilda | Python (sympy GF(p)) | PARI | Mathematica |
  |---|---|---|---|---|
  | `prim: InterpolatingPolynomial t=50 mod p x20` | 175 ms | 1684 ms | 22 ms | 76 ms |
  | `prim: InterpolatingPolynomial n=200 mod p x4` | 1064 ms | 20808 ms | 51 ms | 355 ms |

  Mathilda's Q-route is ~2.5× faster than Mathematica's `InterpolatingPolynomial`
  at t=50 but ~3× slower at n=200 (the rational denominators blow up with the
  point count). sympy's own Q-route (`interpolate`) is far worse than either
  and was dropped for its GF(p) linear-algebra path, which is still slower than
  both CASes because it forms the full Vandermonde and does an LU solve. **The
  gap to close for Mathilda is a finite-field interpolation primitive** (an
  `InterpolatingPolynomial[..., Modulus->p]` option, or divided differences
  kept in `F_p`); the Q-route is the only option today and it is quadratic in
  the denominators.

- **No wrong answers, no crashes.** All four systems agree on every check.

- **`app:` (hand-written protocol) is correct and well within budget in
  Mathilda** but ~5–10× slower than Python/PARI on the Horner and Lagrange
  loops — the per-`Mod` interpreter overhead on bignum scalars, as expected for
  a loop-heavy pure-WL program (not an absence, a constant-factor interpreter
  cost).

Baselines measured on this host (Apple M-series): see the table above.
