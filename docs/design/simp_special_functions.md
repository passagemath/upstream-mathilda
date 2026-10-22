# Simplify for special functions: the family engine

`src/simp/simp_sf_engine.c` is a table-driven `Simplify` pass. It runs in
`simp_dispatch` before the shape pipelines, and it treats special functions the
way `Together` treats rational functions: it puts them in a canonical form, after
which a vanishing expression visibly vanishes. Each family of functions lives in
its own descriptor file (`simp_bessel.c`, `simp_sf_gamma.c`, …). This document
records, for every identity a descriptor encodes, the **derivation from first
principles** and the **domain of validity**.

No identity here was taken from another system's output. Each one is
re-derived below, and each test in `tests/test_simp_sf_*.c` also checks
**numerically** (30 digits, at two real points and one complex point) that the
identity holds. A mis-derived identity therefore fails its test instead of
being enshrined by it. The numeric checks already caught one evaluator bug this
way: exact `Csc[k Pi/5]` was wrong (fixed in v0.535).

## 1. The model

### Families and lattices

A special function `C[a1, …, an]` has *shift arguments*: the arguments along
which an exact recurrence moves it, such as Bessel's order, the Gamma argument,
and the three parameters of 2F1. Two atoms belong to one **family** when they
share the head, the non-shift arguments, and, for each shift argument, the same
*symbolic part* `s` and *rational part mod 1*. An argument is read as
`s + num/den` by `sf_split`. The members of a family are therefore points of an
integer lattice, at offsets `k ∈ Z^d` from a base point.

### Affine shift matrices

Take the basis vector at lattice point p to be

    v(p) = (C[p], 1)              rank 1   (Gamma, LogGamma, PolyGamma, ...)
    v(p) = (C[p], C[p+e0], 1)     rank 2   (Bessel, Legendre, 2F1, ...)

A descriptor's `up` builder supplies the matrix with `v(p + e_d) = M(p, d) v(p)`.
The trailing `1` carries inhomogeneous terms, as in
`LogGamma[z+1] = LogGamma[z] + Log[z]`. Down steps use `M(p − e_d, d)^(-1)`,
which the engine inverts. A descriptor therefore only ever states the
**raising** relation, which is the one derived below. The engine walks from a
candidate basis point to every member and caches the transfer matrices. Each
member then becomes a linear combination of the basis atoms, and the expression
becomes a rational function in those atoms.

That rewrite is a **canonical form** wherever the basis atoms are algebraically
independent over the coefficient field. For Bessel functions at orders outside
1/2 + Z this is Siegel's theorem. For Gamma it is Hölder's theorem: Γ satisfies
no algebraic differential equation, so Γ(z) is transcendental over C(z), and in
particular so is any Γ(s + c) over Q(s). In those cases the `Together`
numerator vanishes exactly when the expression does.

### Pipeline (`simp_special_functions`)

0. **Prefilter.** One tree walk with pointer-hash probes. It returns `NULL` when
   no registered head is present, which is the common case.
1. **Canonical generators (`canon` hooks).** These apply ratio heads → Gamma,
   reflections, multiplication, and BesselK evenness. A replacement is walked
   again, to a bounded depth.
2. **Families.** Every triggering family is rewritten onto the basis that scores
   best for the *whole* expression. First every family is rewritten
   unconditionally. That is the canonical form, and it is needed when only the
   combination cancels: Pascal's rule spans three Gamma families, and each
   rewrite alone raises the score. If that result does not pass the gate, the
   greedy variant runs instead, adopting a family's rewrite only when it
   improves the current form. A rewrite can introduce atoms of a *lower* family
   (`BarnesG → Gamma`), so the stage repeats while a family still triggers.
   Partner pairs (§2.3) are rewritten jointly and reduced modulo their
   Casoratian.
3. **Passes.** Whole-expression candidates, such as the refolds of Gamma ratios
   into the input's vocabulary.
4. **Gate.** The result must strictly improve on the input (`sf_improves`). It is
   then dispatched once more.

### Normal form of a candidate (`normalize_form`)

`Together` maps into function arguments, so a coefficient kernel such as
`Gamma[x + 1/3]` comes back as `Gamma[1/3 (1 + 3 x)]`. Every candidate therefore
has the **shift arguments** of its atoms expanded. Powers are put in the
spelling that `simp_dispatch` and the evaluator leave on the input:
PrimeRebase, then PowerOneify, then `b^(q + rest) → b^t b^(q − t + rest)`, where
`t` is the integer part of `q` truncated toward zero and the exponent is
expanded. So `2^(1 + 2x)` becomes `2 2^(2x)`. That way a prefactor the engine
builds meets the input's own copy of the same power.

### Acceptance (`sf_improves`)

A candidate is accepted when its score strictly drops. It is also accepted when
its number of special-function **kernels** strictly drops while its score stays
within `s0 + s0/4 + 2`. The leaf score cannot see that eliminating a
transcendental kernel is a simplification: `Gamma[1/3] Gamma[2/3]` and
`2 Pi/Sqrt[3]` tie, and `(n − k)/(1 + k)` outscores
`Binomial[n, k+1]/Binomial[n, k]`. Both measures are bounded below, and every
accepted step lowers one of them, so the rule cannot cycle. Passes additionally
accept ties, because they restore vocabulary at equal score.

`Factor` is among the candidate forms only for rank-1 families and passes. On the
large two-kernel polynomials of a rank-2 rewrite it costs about 6× the whole
pass and never wins.

## 2. Bessel J, Y, I, K (`simp_bessel.c`)

### 2.1 Recurrences

The recurrences come from the series
`J_ν(z) = Σ_k (−1)^k (z/2)^(2k+ν) / (k! Γ(k+ν+1))`. Comparing coefficients
gives `J_{ν−1} + J_{ν+1} = (2ν/z) J_ν` (DLMF 10.6.1). Y obeys the same
recurrence, because it is a combination of `J_{±ν}` with ν-periodic
coefficients. For `I_ν(z) = i^(−ν) J_ν(iz)` the recurrence becomes
`I_{ν−1} − I_{ν+1} = (2ν/z) I_ν`, and K, a combination of `I_{±ν}`, gives
`K_{ν−1} − K_{ν+1} = −(2ν/z) K_ν`. All four have the form
`C[m+1] = p(m) C[m] + q C[m−1]`, and the up matrix is
`[[0, 1, 0], [q, p(m+1), 0]]`.

### 2.2 Reflection

- **K is even in its order, for every order.** This follows from
  `K_ν(z) = ∫_0^∞ e^(−z cosh t) cosh(νt) dt`. The canonical sign must not break
  the family class: a symbolic order flips when its symbolic part reads
  negative. A numeric order flips when `−μ` has the smaller class, that is
  `frac(μ) > 1/2`, or `frac(μ) ∈ {0, 1/2}` with `μ < 0`. When a `BesselI` with
  the same argument is present, K instead follows the class of that I, so that
  the pair relation applies.
- **Integer orders.** The generating function
  `exp(z/2 (t − 1/t)) = Σ J_n t^n` is invariant under `t → −1/t`. That gives
  `J_{−n} = (−1)^n J_n`, and likewise `Y_{−n} = (−1)^n Y_n` and `I_{−n} = I_n`.
  These are applied only when the order is **provably an integer** under the
  assumptions. The sign is written as `(−1)^c (−1)^t`, with `t` the
  positive-reading symbolic part, so that it combines with a `(−1)^n` already
  in the input.

### 2.3 Pairs (Casoratians)

- **J and Y.** Abel's identity for `z² w'' + z w' + (z² − ν²) w = 0` gives
  `W(z) = c/z`. The large-z asymptotics fix `c = 2/π`, so
  `W{J_ν, Y_ν} = 2/(πz)`. Substituting `J' = J_{ν−1} − (ν/z) J` and using the
  recurrence gives the Casoratian `J_{ν+1} Y_ν − J_ν Y_{ν+1} = 2/(πz)`.
- **I and K.** In the same way, `W{I_ν, K_ν} = −1/z` gives
  `I_{ν+1} K_ν + I_ν K_{ν+1} = 1/z`.

Once both families sit on one basis `{b, b+1}`, the expression is reduced modulo
the relation with `PolynomialReduce` in the lex order
`C[b+1] > D[b] > C[b] > D[b+1]`, which eliminates the monomial `C[b+1] D[b]`.
The k-step Casoratian follows:
`J_{n+k} Y_n − J_n Y_{n+k} = B_k · 2/(πz)`, where `B_k` is the k-th Lommel
coefficient, with `B_2 = 2(n+1)/z` and `B_3 = 4(n+1)(n+2)/z² − 1`.

## 3. Gamma and the ratio heads (`simp_sf_gamma.c`)

### 3.1 Recurrences

| head | up matrix | derivation | validity |
|---|---|---|---|
| `Gamma[z]` | `[[z, 0]]` | `Γ(z+1) = ∫_0^∞ t^z e^(−t) dt`, integrated by parts once; continued meromorphically | all z off the poles |
| `Gamma[a, z]` | `[[a, z^a e^(−z)]]` | `Γ(a+1, z) = ∫_z^∞ t^a e^(−t) dt = z^a e^(−z) + a Γ(a, z)`, by parts | all a; z on the principal branch |
| `LogGamma[z]` | `[[1, Log[z]]]` | the logarithm of the Gamma recurrence | C∖(−∞, 0], see below |
| `BarnesG[z]` | `[[Gamma[z], 0]]` | the defining functional equation `G(z+1) = Γ(z) G(z)` | all z |

**LogGamma.** Both sides of `LogGamma(z+1) = LogGamma(z) + log z` are analytic
on the slit plane C∖(−∞, 0], and they agree for real z > 0. By the identity
theorem they agree on the whole slit plane, which is connected. So the
principal-branch identity holds wherever LogGamma is analytic. This was checked
numerically at the complex sample point.

**Bulk shifts.** `Gamma[z + k]/Gamma[z] = Pochhammer[z, k]` for k > 0, and
`1/Pochhammer[z + k, −k]` for k < 0. So `Gamma[x + 500]/Gamma[x + 498]` costs one
step, not 500.

### 3.2 Canonical generators

The ratio heads are rewritten onto Gamma using Euler's Beta integral and the
definitions:

    Pochhammer[a, n] = Γ(a+n)/Γ(a)
    Beta[a, b]       = Γ(a) Γ(b)/Γ(a+b)
    Binomial[n, k]   = Γ(n+1)/(Γ(k+1) Γ(n−k+1))
    Factorial[n]     = Γ(n+1)

These hold as identities of meromorphic functions. Where Binomial is defined by
a limit (negative integer n), the rewrite, like every Simplify rewrite, holds
generically. A rewrite happens only when the expression holds at least two
implied Gamma arguments: a lone `Beta[a, 1 − a]` counts three.
**Factorial-only** expressions are left to `simp_factorial` (the FactorialRules
seed). It groups factorials by offset with or without assumptions and knows
`Factorial2`. The engine converts `Factorial` only to meet another Gamma-type
head.

**Refold passes** turn Gamma ratios back into the heads the input used:
`Γ(v+w)/Γ(v) → Pochhammer[v, w]`, `Γ(a)Γ(b)/Γ(a+b) → Beta`,
`Γ(N)/(Γ(K)Γ(M))` with `N−1 = (K−1) + (M−1)` → `Binomial`, and
`Γ(c + s) → Factorial[c − 1 + s]` for c ≥ 1.

### 3.3 Reflection

Euler's reflection formula `Γ(w) Γ(1 − w) = π / sin(πw)` follows from
`B(w, 1 − w) = ∫_0^∞ t^(w−1)/(1 + t) dt`, evaluated by a keyhole contour, for
0 < Re w < 1, and extends to all non-integer w by analytic continuation.

`Gamma[μ]` is reflected onto `Gamma[1 − μ]` only when two conditions hold:

- an atom of the class of `1 − μ` is present; and
- μ is the less canonical side of the pair, meaning its symbolic part reads
  negative, or, for a numeric μ, `frac(μ) > 1/2`.

The sine is written as `Sin[π(1 − μ)]`, which equals `Sin[πμ]`, so it uses the
partner's own argument and meets the input's trig verbatim. In the
half-integer class, `1 − μ = k + 1/2 + t`, the sine is written `(−1)^k Cos[πt]`.
Shifted pairs follow: `Γ(x+k) Γ(1−x−k) = (−1)^k π / sin(πx)`.

### 3.4 Multiplication

Gauss's multiplication formula is

    Γ(m z) = (2π)^((1−m)/2) m^(m z − 1/2) ∏_{j=0}^{m−1} Γ(z + j/m)

It follows from the Beta-integral product (or from Bohr–Mollerup applied to the
right-hand side), and `m = 2` is Legendre's duplication formula. It is applied to
`Gamma[m s + c]` when a Gamma of symbolic part `s` is present, where `m` is the
**integer content** of the expanded symbolic part (`2a + 2b → m = 2`). For
`m = 2^e m'`, the prefactor is written one power per base:
`2^((1−m)/2 + e(m z − 1/2)) m'^(m z − 1/2) π^((1−m)/2)`. The rewrite is a
heuristic expansion (one Gamma becomes m), not a decision procedure, and is
attempted for `m ≤ 16`.

## 4. PolyGamma, HarmonicNumber, the zeta family (`simp_sf_zeta.c`)

### 4.1 Recurrences (rank 1, affine)

| head | up row | derivation | validity |
|---|---|---|---|
| `PolyGamma[n, z]` | `[[1, (−1)^n n! z^(−n−1)]]` | the n-th derivative of `ψ(z+1) = ψ(z) + 1/z`, which is the log-derivative of `Γ(z+1) = zΓ(z)` | integer `n ≥ 0` only. For a non-integer order, the fractional-derivative PolyGamma does **not** satisfy it (checked numerically), so the descriptor declines. |
| `HarmonicNumber[n]` | `[[1, 1/(n+1)]]` | `H(n) = ψ(n+1) + γ`, then 4.1 at `n = 0` | all n |
| `HarmonicNumber[n, r]` | `[[1, (n+1)^(−r)]]` | the partial sums `Σ_{k≤n} k^(−r)`, continued analytically | all n, r |
| `HurwitzZeta[s, a]` | `[[1, −a^(−s)]]` | split the first term off `Σ_{k≥0} (k+a)^(−s)`; analytic continuation | all a off the poles |
| `Zeta[s, a]` | `[[1, −(a²)^(−s/2)]]` | `Zeta[s, a]` sums `((k+a)²)^(−s/2)`, which is the convention for Re a < 0 and equals `(k+a)^(−s)` when Re a > 0. The term split off is therefore `(a²)^(−s/2)`, checked at a − 1 < 0. | all a |
| `LerchPhi[z, s, a]` | `[[1/z, −a^(−s)/z]]` | `Φ = Σ z^k (k+a)^(−s) = a^(−s) + zΦ(z, s, a+1)` | all a (checked at a < 0) |

### 4.2 Reflection and multiplication

- **PolyGamma reflection.** `ψ(1−w) − ψ(w) = π cot(πw)` is the log-derivative
  of Euler's reflection. Differentiating n times gives
  `ψ^(n)(1−w) = (−1)^n (ψ^(n)(w) + π D^n cot(πw))`. The derivative of the
  cotangent is built with `D` on a dummy symbol.
- **PolyGamma multiplication.** Taking the log-derivative of Gauss's formula,
  and then n more derivatives, gives
  `ψ(mz) = log m + (1/m) Σ_j ψ(z + j/m)` and
  `ψ^(n)(mz) = m^(−n−1) Σ_j ψ^(n)(z + j/m)` for n ≥ 1.
- **`HarmonicNumber[n] → PolyGamma[0, n+1] + EulerGamma`**, applied only when a
  `PolyGamma[0, _]` is present.
- **Riemann functional equation.** It is applied when both `Zeta[w]` and
  `Zeta[1 − w]` are present, rewriting the less canonical side, and it is
  written in the partner's own variable:
  `ζ(1−w) = 2^(1−w) π^(−w) cos(πw/2) Γ(w) ζ(w)`. This is the Riemann form at
  s = 1 − w, using `sin(π(1−w)/2) = cos(πw/2)`. An integer offset in w is folded
  by quarter turns (`cos(π(1+s)/2) = −sin(πs/2)`).

## 5. PolyLog (`simp_sf_zeta.c`)

A `PolyLog[n, v]` atom is rewritten onto a present partner `PolyLog[n, w]`. The
partner is matched by **value** (`Together[v − w] == 0`) and the rewrite is
written with the partner's own argument. Arguments are expanded, so
`PolyLog[2, -(1/10 + z)]` and `PolyLog[2, -1/10 - z]` are one kernel.

| relation | derivation | validity (probed in all four quadrants and on the cuts) |
|---|---|---|
| reflection `Li₂(1−u) = π²/6 − log u log(1−u) − Li₂(u)` | Euler: differentiate both sides, which agree at u = 1/2 | every u, cuts included |
| duplication `Li_n(u²) = 2^(n−1)(Li_n(u) + Li_n(−u))` | `Li_s(u) + Li_s(−u)` keeps the even terms of the series, so it equals `2^(1−s) Li_s(u²)`, for any order | every u |
| Landen `Li₂(u/(u−1)) = −log²(1−u)/2 − Li₂(u)` | differentiate both sides, which agree at u = 0 | **u < 1 proved**: the sides are 5.575 i apart at u = 17/10, so the identity fails on (1, ∞) |
| inversion `Li₂(1/u) = −π²/6 − log²(−u)/2 − Li₂(u)` | the same method | **u < 0 proved**: the identity fails on (0, ∞) |

Simplify identities must hold for **every** value of their variables, the same
reason `Log[1/x]` is not rewritten as `−Log[x]`. An identity that fails on a real
half-line is therefore gated on an assumption that excludes that half-line. It
is not applied on the grounds that the failure set has measure zero in C.

## 6. Error functions, exponential integrals, Fresnel (`simp_sf_errexp.c`)

**Canonical generators.** These hold unconditionally:

- `Erfc = 1 − Erf` and `Erf[z₀, z₁] = Erf[z₁] − Erf[z₀]`, by definition.
- `Erfi[z] = −I Erf[I z]`, by the substitution t → it in the integral. Applied
  when an Erf or Erfc is present.
- `InverseErfc[x] = InverseErf[1 − x]`, applied when an InverseErf is present.
- `Erf[InverseErf[x]] = x`, `Erfc[InverseErfc[x]] = x` and
  `Erf[InverseErfc[x]] = 1 − x`.
- `LogIntegral[x] = ExpIntegralEi[Log[x]]`, which holds on the whole plane, cuts
  included. Applied only to meet that exact atom.
- `Si[I z] = I Shi[z]` and `Shi[I z] = I Si[z]`, applied when the partner head
  is present.
- `FresnelS[I z] = −I FresnelS[z]` and `FresnelC[I z] = I FresnelC[z]`. Under
  t → it the integrands `sin(πt²/2)` and `cos(πt²/2)` become `sin(−πt²/2)` and
  `cos(−πt²/2)`.

A refold pass writes `Erf[z]` back as `1 − Erfc[z]` for every z that the input
spelled with Erfc.

**Odd-function normal form.** For Erf, Erfi, InverseErf, FresnelS, FresnelC, Si
and Shi, the argument is expanded and the sign is pulled out when it reads
negative, including a negative imaginary coefficient: `Erf[−I x] → −Erf[I x]`
and `FresnelS[−1/3 − x] → −FresnelS[1/3 + x]`. The evaluator folds only
`Erf[−x]`.

**Gated identities.** Each of these fails on an open set or on a real
half-line, as checked numerically:

- `InverseErf[Erf[x]] = x` requires x real, because Erf is not injective on C.
- `CosIntegral[−x] = CosIntegral[x] + I π` and
  `CoshIntegral[−x] = CoshIntegral[x] + I π` require x > 0. At a generic complex
  point they are off by 2πi.
- `ExpIntegralE[n, z] = z^(n−1) Γ(1−n, z)` requires z > 0. It follows from the
  substitution t → t/z in the integral, and is applied when an incomplete Gamma
  of z is present.
- `Ci[I z] = Chi[z] + Iπ/2` is not applied: it fails in the second quadrant.

**Recurrence.** `n E_{n+1}(z) = e^(−z) − z E_n(z)`. With
`E_{n+1}(z) = ∫_1^∞ e^(−zt) t^(−n−1) dt` and `dv = t^(−n−1) dt`, integration by
parts gives `e^(−z)/n − (z/n) E_n(z)`. It holds for every n ≠ 0 and every z, as
checked at real, negative and complex z. The up row is `[[−z/n, e^(−z)/n]]`.

## 7. Notes on the candidate normal form

`normalize_form` is the single place where spellings are reconciled. Every rule
in it was added because a correct rewrite failed to cancel against the input's
own spelling of the same value:

- **Shift arguments are expanded.** `Together` writes `Gamma[1/3 (1 + 3 x)]`.
- **Negations are distributed.** `Log[−(−2 + 3x)]` versus `Log[2 − 3x]`.
- **Odd-function arguments.** As in §6.
- **Trig arguments.** A trig argument that is a multiple of π takes the natural
  spelling `π (a + b)`, never `π a + π b`, and quarter turns are folded
  (`Cos[π/2 + π s] → −Sin[π s]`). Simplify's bottom-up pass hands the engine
  *subexpressions*, so the engine's output has to meet the spelling of the
  outer input.
- **Powers.** As in §1.

After the gate, a result that still contains a trig head (the residue of a
reflection) is finished by Simplify's bottom-up pass instead of a single
dispatch, and is then normalised once more, so that `Csc[u] Sin[u]` cancels.

## 8. Out of scope (no builtin)

`GegenbauerC`, `JacobiP`, `HankelH1`/`HankelH2` and `StruveL` have no evaluating
builtin, so they have no descriptors. `MeijerG` has no shift relations beyond
those of the functions it reduces to; `FunctionExpand` handles it.
