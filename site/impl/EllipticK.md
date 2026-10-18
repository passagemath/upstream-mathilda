---
references:
  - "M. Abramowitz and I. A. Stegun, *Handbook of Mathematical Functions* (Dover, 1964), §17.3.34 — the `P(t) − log(t) Q(t)` form of `K`."
  - "S. L. Moshier, *Cephes Math Library* (1984–2000), `ellpk` — the polynomial-with-logarithm structure (coefficients re-fitted for Mathilda)."
  - "DLMF §19.7.5 — the imaginary-modulus transformation `K(m) = K(m/(m−1))/√(1−m)`."
  - "DLMF §19.8.5 — `K(m) = π/(2 M(1, √(1−m)))`, the arithmetic-geometric mean."
  - "D. A. Cox, *The arithmetic-geometric mean of Gauss*, L'Enseignement Math. **30** (1984) 275-330 — the \"right choice\" of square root for the complex AGM."
  - "R. L. Smith, *Algorithm 116: Complex division*, Comm. ACM **5** (1962) 435."
  - "B. C. Carlson, *Numerical computation of real or complex elliptic integrals*, Numer. Algorithms **10** (1995) 13-26."
source: src/special_functions/elliptic.c
---
**Algorithm.** `builtin_elliptick` tries, in this order: the pole (`ell_is_one`, which
recognises `1`, `1.` and an `EXPR_MPFR` one alike → `ComplexInfinity`); then the
**machine fast path** `ell_complete_machine`, gated by the shared
`sf_machine_scalar_args` (a machine `Real` or machine `Complex` argument; exact and
`EXPR_MPFR` arguments do not pass); then the exact reductions — `K[0] = π/2`,
`K[∞] = 0` (since `K(m) ~ (π − i log 16m)/(2√m)`), and the two lemniscatic singular values
`m = 1/2` and `m = −1` as closed forms in `Γ` (`ell_k_closed_half`,
`ell_k_closed_minus_one`). Any remaining inexact argument — in practice an arbitrary-precision
one, or a machine one the kernels declined — goes to `flint_num_elliptic_k`, i.e. Arb's
`acb_elliptic_k`, under the bridge's accuracy ladder (`nb_eval1` doubles the working
precision while `acb_rel_accuracy_bits` falls short of the digits requested). An exact
non-special argument stays symbolic.

The machine kernels live in `src/special_functions/sf_machine_ellk.c` and are chosen by the
argument's type and value:

- *real `m < 1`* — `sf_machine_ellk_k`: Cephes' structure `K = P(t) − log(t) Q(t)`,
  `t = 1 − m`, two degree-10 Horner chains and one `log`. The coefficients are not
  Cephes': `tools/sf_accuracy/fit_ellk.py` re-fitted them as a relative-error minimax
  (Lawson IRLS in mpmath, points down to `t = 10⁻³⁰⁰` so that `Q(0)` is pinned against
  `|log t| ≈ 690`), fit error 1.2·10⁻¹⁸. For `m < 0` the DLMF 19.7.5 transformation maps
  the argument to `t' = 1/t`, so the same polynomials serve `(−∞, 1)`;
- *complex `m`, or a real `m > 1`* — `sf_machine_ellk_k_c`: the AGM `π/(2 M(1, √(1−m)))` in
  raw `double` re/im (no `<complex.h>`), taking at every step the square root with
  `Re(aₙ₊₁ conj bₙ₊₁) ≥ 0` (Cox's "right choice"), which yields the principal branch with
  the cut on `[1, ∞)`. A zero imaginary part of either sign is read as `+0`, so on the cut
  the value is the limit from below: `EllipticK[2.]` is `1.31103 − 1.31103 I`. The final
  `(π/2)/M` uses Smith's division. A complex slot holding a real `m < 1` is sent back to the
  real kernel, so `0.9 + 0. I` in a complex buffer and the scalar `0.9` agree bit for bit.

**Data structures.** `static const double` coefficient tables and plain `double`
scalars; nothing is allocated. A complex result is built by `sf_machine_complex_result`,
which hands out one shared `Complex` head node by refcount. The same two kernels are the
`NDKU_EllipticK` unary ND kernel in `src/ndkernels.c` (real arm `elliptic_machine_k`, which
now delegates to `sf_machine_ellk_k`; complex arm `sf_machine_ellk_k_c`), registered with
`REG_U`, which also makes the head `packed_aware`, and `Compile[]` calls those registered
kernels at `_Real`, `_Complex` and rank-1 shapes. So scalar, `List`, packed `List`,
`NDArray` and `Compile[]` agree bit for bit (`tests/test_elliptic_machine.c`). A *real*
buffer stays real-closed: an element with `m > 1` makes the real arm decline, the buffer
is abandoned to the `List` path, and that element's scalar call takes the complex kernel.
The incomplete integrals keep Carlson's `R_F`/`R_D` duplication in `elliptic.c`, and the
quasi-period `2k K(m)` they add off the principal strip goes through the same
`elliptic_machine_k`.

**Complexity / limits.** `O(1)`: about 0.12 µs per real call and 0.37 µs per complex one
(the AGM converges quadratically, 4–6 square roots for `|m| ≤ 1`, at most 40 iterations).
Measured against mpmath the real kernel is within 1.7 ulp on `[0, 0.95]` and 2.4 ulp on
`[−50, 1)`; the complex one within 3.4 ulp on `|Re m|, |Im m| ≤ 0.9` and 4.7 ulp out to 20.
The kernels decline the pole, a non-finite result, and `|Re m|` or `|Im m| > 10³⁰⁰` (where
the AGM products would overflow); those cases, and every request above machine precision,
are Arb's, which is rigorous but slower and single-threaded. `Series` at `m = 0` comes from
a dedicated kernel rather than Taylor-via-`D`, because `d/dm K` at `m = 0` is
`(π/2 − π/2)/0`, and `Interval[]` threads by certified monotonicity below `m = 1`.
