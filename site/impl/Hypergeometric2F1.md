---
source: src/special_functions/hypergeopfq.c
references:
  - "DLMF §15 — the Gauss hypergeometric function 2F1."
  - "M. Abramowitz and I. A. Stegun, Handbook of Mathematical Functions (1964), 15.3.3-15.3.12 — the linear transformations and the degenerate logarithmic cases."
  - "R. C. Forrey, Computing the hypergeometric function, J. Comput. Phys. 137 (1997) 79-100."
  - "N. Michel and M. V. Stoitsov, Fast computation of the Gauss hypergeometric function with all its parameters complex, Comput. Phys. Commun. 178 (2008) 535-551."
  - "W. Bühring, An analytic continuation of the hypergeometric series, SIAM J. Math. Anal. 18 (1987) 884-889."
---
**Algorithm.** `builtin_hypergeometric_2f1` (`hypergeopfq.c`) first tries the
machine gate `hg_machine_head` — every argument a machine-scalar leaf
(`sf_machine_scalar_args`), at least one a machine `Real` — and otherwise rewrites
to `HypergeometricPFQ[{a, b}, {c}, z]`, inheriting its exact steps (cancellation,
`z == 0`, termination for non-positive integer `a`/`b`, the closed forms
`-Log[1-z]/z`, `(1-z)^-a` and the `ArcSin[Sqrt[z]]` family). The real machine
kernel `sf_machine_hyper2f1` builds a `Hy2Plan` from the parameters (polynomial,
`(1-z)^-b` / `(1-z)^-a`, or generic) and evaluates: the direct series on
`(0, 0.9]`, Pfaff's `(1-z)^-a F(a, c-b; c; z/(z-1))` on `[-9, 0)`, Gauss's sum
`G(c)G(c-a-b)/(G(c-a)G(c-b))` at `z = 1`, and the shared complex core
`hy2f1_core` elsewhere. That core picks the cheapest of direct / Pfaff series
(`|w| <= 0.75`), the `1-z` and `1/(1-z)` connection formulas (`|w| <= 0.6`,
accepted only while the combination's condition number is at most 2, with the
degenerate logarithmic cases when `c-a-b` is an integer), and — around
`z = e^(±i π/3)` or for ill-conditioned parameters — continuation of the
hypergeometric ODE by Taylor steps from `z0 = 1/2 ± 0.45 i` (`hy2f1_continue`).
A machine-real `z > 1` returns the complex value from below the cut (`zi = -0`).
Complex parameters use the series only, for `|z| <= 0.9`. At arbitrary precision
`|z| < 1` goes to arb's `acb_hypgeom_2f1` (`flint_num_hyp2f1` in
`src/flint_num_bridge.c`), with the MPFR series as fallback.

**Data structures.** Plain C99 `double`, complex values as re/im pairs, with
double-double (`hdd`) re-sums for a series that loses more than one bit or peaks
late; `Hy2Plan` / `Hy2Tab` hold the coefficient recurrence
`c_{k+1} = c_k (p0+k)(p1+k)/((q0+k)(k+1))`, which a buffer tabulates once. When
the terms are one-signed and never grow ("lean", decided from the parameters)
`sf_machine_pfq_vec` sums four buffer elements in lock-step, bit-identical to the
scalar loop. ND: n-ary kernel `NDKN_Hypergeometric2F1` (arity 4, `REG_N`); on
`pack.c`'s `AWARE` list with the family; a buffer containing a cut element comes
back complex. Attributes: `Listable`, `NumericFunction`, `Protected`.

**Complexity / limits.** Scalar about 0.3 µs real and 1.7 µs complex; 10^6 packed
elements about 0.14 s. Term caps 20000 (double) and 6000 (double-double); a
case where nothing converges declines to the exact / MPFR route. The arbitrary
precision path still requires `|z| < 1`, so an exact or MPFR `z` outside the unit
disc stays unevaluated, as does the exact Gauss value at `z = 1`.
`Compile[]` lowers at real and complex scalar shapes (through
`HypergeometricPFQ`'s `OP_KERNN` / `OP_KERNN_C`); the rank-1 array shape reports
`Compiled -> False`, while the packed / `NDArray` path still serves the REPL.
