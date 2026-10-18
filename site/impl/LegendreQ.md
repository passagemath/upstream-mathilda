---
source: src/special_functions/legendre.c
references:
  - "DLMF §14 — Legendre functions of the second kind (§14.3, §14.5.22, §14.7, §14.10.3)."
  - "W. Kahan, Branch cuts for complex elementary functions, in The State of the Art in Numerical Analysis (1987) — the cancellation-free Q_0."
---
**Algorithm.** `builtin_legendre_q` handles `LegendreQ[n, x]`,
`LegendreQ[n, m, x]` and `LegendreQ[n, m, a, x]`. For the two-argument form it
first tries the **machine gate** `legendre_machine` (machine-scalar arguments
with at least one machine `Real`, real integer-valued order, `x` not an exact
`1`), which calls `sf_machine_legendre_q_c` in
`src/special_functions/sf_machine_plog.c`. That kernel takes integer `n >= 0`
up to degree 1000 and uses Christoffel's form `Q_n = P_n Q_0 - W_{n-1}`, running
`P_n` and the polynomial `W_{n-1}` up the three-term recurrence. Seeding a
recurrence on `Q` itself with `Q_1 = x Q_0 - 1` would cancel near `|x| ~ 0.8`.
`Q_0` is `atanh(x)` on `(-1, 1)`. Off that interval it uses Kahan's
cancellation-free formulas, with the closed form's branch: a real `x > 1`
reads `Log[1 - x] = log(x - 1) + iπ`, so the value is complex. It declines at
`x = ±1`. Otherwise the exact path runs: an integer `n >= 0` (cap
`LEG_POLY_CAP = 2000`; `n < 0` is singular and stays symbolic) builds
`Q_n(x) = P_n(x) L(x) + v_n(x)` with `L(x) = (1/2)(Log[1+x] - Log[1-x])` and
`v_n` from the same recurrence seeded `v_0 = 0`, `v_1 = -1` (`legendre_q_closed`).
A non-integer order gives the exact special value
`Q_v(0) = -(Sqrt[Pi]/2) Sin[v Pi/2] Gamma[(v+1)/2] / Gamma[v/2+1]` at an exact
zero. For an inexact `x` on the cut `|x| < 1` it evaluates the two Frobenius
series `Q_v(0) 2F1(-v/2, (v+1)/2; 1/2; x^2) + Q_v'(0) x 2F1((1-v)/2, (v+2)/2; 3/2; x^2)`,
built as an `Expr` over `Hypergeometric2F1`/`Gamma` (`legendre_q_numeric`), so it
inherits their machine, MPFR and complex numerics. Associated `Q_n^m` (all three
types, integer `n, m >= 0`) differentiate `Q_n` in a fresh dummy variable, then
substitute. Over a machine buffer the polynomial parts are evaluated by Horner
(`leg_poly_from_coeffs`, shared with `LegendreP`).

**Data structures.** `Expr`; GMP `mpq_t` coefficient arrays; plain re/im doubles
in the machine kernel; the non-integer numeric path is an `Expr` over `2F1`/`Gamma`.
ND: binary kernel `NDKB_LegendreQ = { ndk_LegendreQ_c, ... }`, which calls the same
`sf_machine_legendre_q_c`, registered `REG_B` (so `packed_aware`). It is not
real-closed: a real element with `|x| > 1` returns complex. Attributes: `Listable`,
`NumericFunction`, `Protected`.

**Complexity / limits.** Machine kernel `O(n)` per point up to degree 1000;
exact recurrence `O(n^2)` big rationals, cap `n = 2000`. The non-integer series
converges only inside `|x| < 1` (the result is kept only if it collapsed to a
number), so `LegendreQ[1.5, 2.]` stays unevaluated. `LegendreQ[n, 1.]` is
`Indeterminate`. `D` for symbolic order is in `calculus/deriv.c`. `Compile[]`
lowers `LegendreQ[n, x]` at scalar and rank-1 array shapes (`Compiled -> True`)
through the same kernel.
