---
source: src/special_functions/hypergeopfq.c
references:
  - "DLMF §13 — confluent hypergeometric functions (Kummer's 1F1); 13.7.2 for the large-|z| expansion."
  - "M. Abramowitz and I. A. Stegun, Handbook of Mathematical Functions (1964), 13.1.27 — Kummer's transformation."
  - "T. J. Dekker, A floating-point technique for extending the available precision, Numer. Math. 18 (1971) 224-242."
---
**Algorithm.** `builtin_hypergeometric_1f1` first tries the machine gate
`hg_machine_head` (every argument a machine-scalar leaf, at least one a machine
`Real`), which calls the raw-`double` kernel `sf_machine_hyper1f1` /
`sf_machine_hyper1f1_c` in `src/special_functions/sf_machine_hyper.c`; otherwise
it rewrites to `HypergeometricPFQ[{a}, {b}, z]` and inherits that engine's exact
steps (cancellation, so `a == b` gives `Exp[z]`; `z == 0 -> 1`; termination to a
polynomial for non-positive integer `a`; `1F1(1;2;z) = (E^z - 1)/z`) and its
arbitrary-precision MPFR series. The machine kernel tries, in order: the closed
forms for `z = 0` or `a = 0`, `a = b` (`exp`) and `(a, b) = (1, 2)` (`expm1(z)/z`,
with a complex form free of cancellation at small `|z|`); the terminating
polynomial for non-positive integer `a` (also past a lower-parameter pole it never
reaches); for `|z| > 20`, DLMF 13.7.2 with both exponential sectors kept (the
subdominant one with the Stokes-line `cos(π a)` factor), accepted only when its
least term is below half an ulp and the two sectors do not cancel; for
`Re z < 0`, Kummer's transformation `e^z 1F1(b-a; b; -z)` with `b - a` carried in
double-double, so the terms are one-signed; otherwise the plain series. A machine
overflow (e.g. `Hypergeometric1F1[1., 2., 900.]`) is promoted to MPFR by the
`PFQ` builtin.

**Data structures.** Plain C99 `double` with complex values as re/im pairs (no
`<complex.h>`); every series measures its own loss (peak term against the sum)
and re-sums in double-double (`hdd`) when it loses more than a bit. ND: n-ary
kernel `NDKN_Hypergeometric1F1` (arity 3, `REG_N`); a packed or `NDArray` `z`
also runs through `PFQ`'s buffer driver `sf_machine_pfq_vec`, the same kernel as
the scalar call, so the surfaces agree bit for bit. On `pack.c`'s `AWARE` list
with the family. Attributes: `Listable`, `NumericFunction`, `Protected`.

**Complexity / limits.** Entire in `z`. Scalar about 0.15 µs real and 0.35 µs
complex; 10^6 packed elements of `1F1(1;2;x)` about 11 ms. The asymptotic branch
needs `|a|, |b|, |b-a| <= 160`; a declined machine point falls to the exact /
MPFR route. A symbolic call is returned in its `HypergeometricPFQ` form.
`Compile[]` lowers at real and complex scalar shapes (via `HypergeometricPFQ`'s
`OP_KERNN` / `OP_KERNN_C`); the rank-1 array shape reports `Compiled -> False`,
while the buffer path still serves the REPL.
