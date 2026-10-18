---
source: src/special_functions/loggamma.c
references:
  - "DLMF §5.11 — Stirling's series for ln Gamma."
  - "D. E. G. Hare, *Computing the principal branch of log-Gamma*, J. Algorithms 25 (1997) 221–236."
  - "C. Lanczos, *A precision approximation of the gamma function*, SIAM J. Numer. Anal. Ser. B 1 (1964) 86–96."
  - "T. J. Dekker, *A floating-point technique for extending the available precision*, Numer. Math. 18 (1971) 224–242."
---
**Algorithm.** `builtin_loggamma` evaluates the analytic continuation of `log Γ(z)`. This is the continuous branch, with one cut on the negative real axis, and it is not `Log[Gamma[z]]`. `loggamma_one_arg` checks the cases in this order:

1. **Machine real.** For `x > 0` it uses libm `lgamma`; an overflow leaves the call unevaluated. At a non-positive integer it returns `ComplexInfinity`. For `x < 0` it returns `lgamma(x) - i π Ceiling[-x]`, the value approached from above.
2. **Complex with an MPFR part** → `lg_mpfr_complex`: Stirling's series in the file-local `lcx` toolkit. The argument is shifted up to `Re ~ K/3` with `K ~ wp/2` exact `mpq` Bernoulli terms, and the reflection uses a continued `log sin(π z)`.
3. **Machine complex** → `loggamma_machine_complex`, which calls `sfg_loggamma_c` in `sf_machine_gamma.c`. Taylor tiles about 30 centres in the strip `1 <= Re < 5/2` cover `-5 <= Re z < 10`, `|Im z| < 19/4`. On the real-axis row the tiled series is for `lnΓ/((z-1)(z-2))`, so both zeros are carried exactly. Elsewhere the kernel uses Stirling's series after shifting to `Re >= 7`. The shift's `N` logarithms are summed as one log of their product, with the `2πi` windings counted as the product crosses the negative axis (Hare's continuous branch). The kernel declines `Re z < -30`; those inputs fall back to a `g = 7` Lanczos log-gamma with a winding-correct reflection.
4. **Exact** (`lg_exact`). A positive integer gives `Log[(n-1)!]`. A half-integer gives the `Log` of the exact `Sqrt[Pi]` form, plus `-Ceiling[-z] π I` when `z < 0`. A non-positive integer gives `Infinity`.
5. **Symbolic infinities.** `Infinity → Infinity`, `-Infinity → Indeterminate`, and `ComplexInfinity` or a directed complex infinity → `ComplexInfinity`.
6. **Arbitrary-precision real** → `mpfr_lgamma` (which gives `log|Γ|` and its sign), plus the same branch term.

An `Interval` argument routes through `interval_apply_function`.

**Data structures.** `Expr`; `lcx` (a pair of `mpfr_t`) with an exact `mpq` Bernoulli generator for the arbitrary-precision complex path; `(re, im)` doubles in the machine kernel, with the tile constants stored as double-doubles. The ND kernel is `NDKU_LogGamma = { loggamma_machine_complex, ndk_LogGamma_r }`, registered with `REG_U`, so it is packed-aware. The real arm uses `lgamma` and declines `x <= 0`, where the value is complex; a declined buffer drops to the `List` path. The complex arm is the same continued log-gamma the scalar call uses. `Compile[]` lowers `LogGamma` through that kernel for real and complex scalars and for rank-1 arrays, and a compiled complex call agrees bit for bit with the interpreter.

**Complexity / limits.** The machine paths are `O(1)`. The complex kernel has a median error of 0.5 ulp and a maximum of 3.3 ulp on `0.5 <= Re z <= 4`, `|Im z| <= 3`; the Lanczos series it replaced was off by up to 689 ulp near the zeros at 1 and 2. The arbitrary-precision Stirling path costs more as the precision rises. The imaginary part grows without bound as `|Im z|` grows, which is what separates this function from `Log[Gamma[z]]`. One known inconsistency at the poles: the exact path returns `Infinity`, but a machine pole (`LogGamma[0.]`) returns `ComplexInfinity`. Attributes: `Listable`, `NumericFunction`, `Protected`.
