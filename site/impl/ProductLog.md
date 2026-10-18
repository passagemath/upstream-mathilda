---
source: src/special_functions/productlog.c
references:
  - "DLMF §4.13 — the Lambert W-function."
  - "R. M. Corless, G. H. Gonnet, D. E. G. Hare, D. J. Jeffrey and D. E. Knuth, On the Lambert W function, Adv. Comput. Math. 5 (1996) 329-359."
---
**Algorithm.** `builtin_productlog` evaluates the Lambert W function:
`ProductLog[z]` is the principal branch `W_0`, `ProductLog[k, z]` the `k`-th
branch (`k` an explicit machine `Integer`). It first tries the **machine gate**
`productlog_machine`: when every argument is a machine-scalar leaf with at least
one machine `Real` (`sf_machine_scalar_args`) and the branch is `0`, the raw-double
kernels in `src/special_functions/sf_machine_plog.c` answer. `sf_machine_productlog`
(real `x >= -1/e`) seeds with `x(1 - x + 3/2 x^2)`, `log(1 + x)` or the
asymptotic `L1 - L2 + L2/L1`, and runs Halley's iteration. Near the branch point it
forms `d = x + 1/e` with `1/e` split into a double head and tail (exact by
Sterbenz's lemma), so the real/complex boundary is decided correctly, and it
iterates in `u = w + 1` on `g(u) = (u-1)e^u + 1 - e d`, which cancels the two
`-1/e` terms analytically. `sf_machine_productlog_c` handles `x < -1/e` and complex
`z`: branch-point series, a Padé seed near the origin, or the asymptotic seed,
then complex Halley. Both kernels **verify** the result (`w + log w = log z`
modulo `2πi`, and that `w` lies in the range of `W_0`) and decline otherwise.
Anything the gate does not take goes to `productlog_eval`: exact special values
(`0`, `E -> 1`, `-1/E -> -1` on branches 0 and -1, `-Pi/2 -> I Pi/2`,
`Infinity`, and `ProductLog[k, 0] = -Infinity` for `k != 0`), then the
complex-MPFR Halley core `pl_core` (`ncpx`) for any inexact argument: other
branches at machine precision, and every branch at arbitrary precision, seeded by
the branch-point series in `p = sqrt(2(e z + 1))`, the Maclaurin seed, or
`L1 = log z + 2πik`. A real seed keeps the MPFR iteration exactly real, so a
real-valued branch returns a real leaf. `D` lives in `calculus/deriv.c`; Series at
0, `-1/E` and `Infinity` in `calculus/series.c`.

**Data structures.** `Expr`; raw `double` with complex values as re/im pairs (no
`<complex.h>`) on the machine path; the shared `ncpx` (`mpfr_t` re/im) toolkit with
guard bits on the MPFR path. ND: unary kernel
`NDKU_ProductLog = { sf_machine_productlog_c, ndk_ProductLog_r, ... }`, registered
`REG_U` (so `packed_aware`); a real buffer whose elements fall below `-1/e`
returns complex values. Attributes: `Listable`, `NumericFunction`, `Protected`.

**Complexity / limits.** Halley converges cubically (60 iterations at most on the
machine path, 100 plus a polishing step on the MPFR path). A real leaf is returned
only on the real-valued domains (`W_0` for `x >= -1/e`, `W_{-1}` for
`-1/e <= x < 0`). Exact non-special arguments stay symbolic. The branch index must
be an `Integer`: `N[ProductLog[-1, -1/10]]` numericalizes `k` to `-1.` and stays
unevaluated. `Compile[]` lowers the principal branch at real and complex scalar
shapes and at rank-1 array shape (`Compiled -> True`), calling the same
kernels.
