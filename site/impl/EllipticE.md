---
references:
  - "M. Abramowitz and I. A. Stegun, *Handbook of Mathematical Functions* (Dover, 1964), §17.3.36 — the `P(t) − t log(t) Q(t)` form of `E`."
  - "S. L. Moshier, *Cephes Math Library* (1984–2000), `ellpe` — the polynomial-with-logarithm structure (coefficients re-fitted for Mathilda)."
  - "DLMF §19.7.5 — `E(m) = √(1−m) E(m/(m−1))`; §19.8.6 — `E` from the AGM."
  - "D. A. Cox, *The arithmetic-geometric mean of Gauss*, L'Enseignement Math. **30** (1984) 275-330."
  - "B. C. Carlson, *Computing elliptic integrals by duplication*, Numer. Math. **33** (1979) 1-16."
  - "B. C. Carlson, *Numerical computation of real or complex elliptic integrals*, Numer. Algorithms **10** (1995) 13-26."
  - "DLMF §19.25.7 — `E(φ|m) = s R_F(...) − (m/3) s³ R_D(...)`."
source: src/special_functions/elliptic.c
---
**Algorithm.** The head is arity-overloaded and `builtin_elliptice` dispatches on argument
count; a wrong count emits `EllipticE::argt`.

- **Complete, `E[m]`.** The machine fast path `ell_complete_machine` runs first, gated by
  the shared `sf_machine_scalar_args` (machine `Real` or machine `Complex`; exact and
  `EXPR_MPFR` arguments do not pass). A real `m ≤ 1` goes to `sf_machine_ellk_e` in
  `src/special_functions/sf_machine_ellk.c`: Cephes' structure `E = P(t) − t log(t) Q(t)`,
  `t = 1 − m`, with degree-10 `P`, `Q` re-fitted by `tools/sf_accuracy/fit_ellk.py`
  (relative minimax, fit error 5.5·10⁻¹⁹), `t = 0` answered as `E(1) = 1`, and `m < 0`
  folded onto `t' = 1/t` by DLMF 19.7.5. A complex `m`, or a real `m > 1`, goes to
  `sf_machine_ellk_e_c`: the AGM `K = π/(2 M(1, √(1−m)))` with Cox's "right choice" of
  square root, plus the running sum `E = K (1 − Σ 2ⁿ⁻¹ cₙ²)` (DLMF 19.8.6), all in raw
  `double` re/im. The branch cut is `[1, ∞)` and a real `m` on it takes the limit from
  below, so `EllipticE[2.]` is `0.59907 + 0.59907 I`. Only if that kernel declines does
  `elliptice_complete` run: Arb's `acb_elliptic_e` for an inexact argument, then the exact
  reductions `E[0] = π/2`, `E[1] = 1`, `E[∞] = ComplexInfinity`.
- **Incomplete, `E[φ, m]`.** `elliptice_incomplete` runs the numeric path first —
  `elliptic_machine_e_inc` for a real pair, else Arb's `acb_elliptic_e_inc` — then the
  reductions `E[0, m] = 0`, `E[φ, 0] = φ`, `E[π/2, m] = E[m]` and the oddness fold. The
  machine kernel is Carlson's: `elliptic_inc_real(φ, m, want_E = true)` reduces the
  amplitude to the principal strip, evaluates `s R_F(c², 1−ms², 1) − (m/3) s³ R_D(…)`, and
  closes the quasi-period with the complete `E(m)` above.

**`R_F` and `R_D` are one loop.** Every caller that wants `R_D` wants `R_F` at the same
arguments, and the two recurrences walk an identical `xₘ, yₘ, zₘ` sequence, so
`carlson_rf_rd` computes both from one duplication — three square roots per step instead
of six, which pays for the tighter `EC_ERRTOL_RD = 0.0015`.

**`E[φ, 1]` is `Sin[φ]` only on the principal strip.** `E(φ|1) = ∫₀^φ |cos t| dt`, which is
`Sin[φ]` for `|φ| ≤ π/2` and not beyond it. Applied unconditionally the rule made
`EllipticE[2, 1]` answer `Sin[2] = 0.909297` where the value is `2 − Sin[2] = 1.090703`.
The reduction is gated on `ell_in_principal_strip`, so an amplitude off the strip stays
symbolic rather than wrong, and `N[]` routes it to Arb.

**Data structures.** `static const double` coefficient tables and `double` scalars for the
machine kernels, `acb_t` with the bridge's accuracy ladder above machine precision. Both
arities are registered ND kernels: the complete form is `NDKU_EllipticE` (`REG_U`), whose
real arm is `elliptic_machine_e_complete` (delegating to `sf_machine_ellk_e`) and whose
complex arm is `sf_machine_ellk_e_c`, so a complex `NDArray` and `Compile[{{z, _Complex}}, …]`
run the same code as the scalar; the incomplete form is `NDKB_EllipticE` (`REG_B`), real
only. Scalar, `List`, packed, `NDArray` and `Compile[]` agree bit for bit
(`tests/test_elliptic_machine.c`). A real buffer containing an `m > 1` element is
abandoned to the `List` path, whose scalar call then takes the complex kernel.

**Complexity / limits.** `O(1)`, no allocation: about 0.12 µs per real complete call and
0.36 µs per complex one; 1.3 ulp on `[0, 0.95]` and 2.4 ulp on `[−50, 1)` real, 3.9 ulp
complex on `|Re m|, |Im m| ≤ 0.9`. The complete kernels decline only at overflow
(`|Re m|` or `|Im m| > 10³⁰⁰`) or a non-finite result. The incomplete machine kernel is
real-principal-domain only (`1 − m sin²φ ≥ 0`, real `φ`); a complex amplitude or parameter
goes to Arb even at machine precision, and comes back as a 53-bit arbitrary-precision
number rather than a machine `Complex`. `Interval[]` threads the complete form by
certified monotonicity (`E` decreasing in `m` below 1); `Series` at `m = 0` uses the
dedicated kernel `(π/2) Σ aₖ mᵏ/(1−2k)`.
