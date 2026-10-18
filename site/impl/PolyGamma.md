---
source: src/special_functions/polygamma.c
references:
  - "DLMF §5.15 — the polygamma functions."
  - "M. Abramowitz & I. A. Stegun, Handbook of Mathematical Functions, §6.4.11 — the asymptotic expansion."
  - "W. J. Cody, A. J. Strecok and H. C. Thacher, *Chebyshev approximations for the psi function*, Math. Comp. 27 (1973) 123–127."
  - "DLMF §25.11 — the Hurwitz zeta tail used for psi^(n), n >= 1."
---
**Algorithm.** `builtin_polygamma` evaluates `PolyGamma[n, z] = ψ⁽ⁿ⁾(z)`. It first tries `pg_machine`, the machine fast path. That path applies when `sf_machine_scalar_args` accepts the call and the order is an integer `0 <= n <= 40` (`SFG_POLYGAMMA_MAX_N`); the one-argument form counts as `n = 0`. It dispatches to `sf_machine_gamma.c`:

- **Real digamma** (`sfg_digamma`). On `[1, 2]`, a degree-22 Chebyshev interpolant of `ψ(x)/(x - x₀)` is multiplied back by `(x - x₀)`, with the zero `x₀ = 1.4616…` held as a double-double, so the answer stays relatively accurate right through the zero (the Cody–Strecok–Thacher technique). From `(0, 10)` the kernel reaches `[1, 2]` by integer shifts, which are exact in double. For `x >= 10` it uses the asymptotic series. For `x < 0` it uses reflection with an exactly reduced `π cot(π x)`.
- **Real ψ⁽ⁿ⁾, n >= 1, x > 0** (`sfg_polygamma`). Upward recurrence to `x >= 10 + n/2`, then a Hurwitz-zeta asymptotic tail. Every term has one sign, each gets a first-order correction for the rounding of `x + j`, and the sum is compensated.
- **Complex** (`sfg_polygamma_c`). For ψ, the kernel uses the Taylor tiles shared with `Gamma`/`LogGamma` (`-5 <= Re z < 10`, `|Im z| < 19/4`), a Taylor series about `x₀` within `0.35`, the asymptotic series after shifting to `Re z >= 10`, and reflection when `Re z < -5`. For `n >= 1`, it uses the recurrence plus the tail on `Re z >= 0`.

A machine `PolyGamma[x]` is answered directly and is identical to `PolyGamma[0, x]`. If the kernel declines (a pole, or `n >= 1` with `Re z < 0`), the layered path runs:

1. `PolyGamma[z]` rewrites to `PolyGamma[0, z]`.
2. `n = -1` gives the inert `LogGamma[z]`; `n <= -2` stays unevaluated.
3. Special values: `Infinity` gives `Infinity` (n = 0) or `0` (n >= 1). `ComplexInfinity` and `Indeterminate` give `Indeterminate`. A non-positive integer `z` is a pole, `ComplexInfinity`.
4. Positive-integer `z`: `n = 0` gives `H_{z-1} - EulerGamma`. Odd `n` gives a rational plus a rational multiple of `π^(n+1)` (through ζ(n+1) and Bernoulli numbers). Even `n` stays symbolic.
5. `n = 0` at a positive half-integer: the Gauss digamma closed form, e.g. `-EulerGamma - Log[4]` at `1/2`. Other rationals stay unevaluated.
6. Numeric: an inexact real uses `mpfr_digamma` for `n = 0`, otherwise a recurrence shift plus the A&S 6.4.11 Bernoulli asymptotic in MPFR. An inexact complex `z` uses the same asymptotic in the `pcx` toolkit.

An `Interval` argument with an integer order goes to `interval_polygamma`.

**Data structures.** `Expr`; a process-lifetime exact `mpq` Bernoulli cache; `pcx` (a pair of `mpfr_t`) for the complex asymptotic; `(re, im)` doubles and generated coefficient tables (`tools/sf_accuracy/gen_gamma_tables.py`) in the machine kernel. There are two ND kernels, both packed-aware. The unary `NDKU_PolyGamma = { ndk_PolyGamma_u, ndk_PolyGamma_r }` runs the real digamma or `sfg_polygamma_c(0, ·)`. The binary `NDKB_PolyGamma` runs `sfg_polygamma_c(n, ·)` and declines a non-integer, negative, or `> 40` order. `Compile[]` lowers `PolyGamma[x]` and `PolyGamma[n, x]` through these kernels, for real and complex scalars and for rank-1 arrays. The scalar, packed, `NDArray` and `Compile[]` answers agree bit for bit.

**Complexity / limits.** The machine kernels are `O(1)`, about 0.14 µs real and 0.35 µs complex. Real ψ has a median error of 0.5 ulp and a maximum of 3.5 ulp; complex has a median of 0.7 ulp, with maxima beside ψ's real zero. The MPFR asymptotic needs the argument shifted up to `Re ~ 0.13 wp + n`, so its cost grows with the precision. Limits: exact argument cap `100000`, exact order cap `256`, numeric order cap `1024`. Complex or fractional orders stay unevaluated. Attributes: `Listable`, `NumericFunction`, `Protected`.
