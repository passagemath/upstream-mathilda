---
source:
  - src/special_functions/chebyshev_u.c
  - src/special_functions/chebyshev_common.c
  - src/special_functions/sf_machine.c
references:
  - "DLMF §18 — classical orthogonal polynomials (§18.9 recurrences; U_n(cos θ) = sin((n+1)θ)/sin θ and the 2F1 representation)."
---
**Algorithm.** `builtin_chebyshev_u` dispatches through `chebyshev_u_two_arg`, which mirrors `ChebyshevT` and shares its substrate in `chebyshev_common.c`. (1) A non-integer order at `x = 0` or `x = 1` returns `Cos[n Pi/2]` or `n + 1`. (2) An integer order with an inexact `x` is reflected with `U_{-1} = 0` and `U_{-m} = -U_{m-2}` (`cheb_reflect` returns a degree and a sign). At machine precision with degree ≤ 1000 it goes to `sf_machine_chebyshev_u` / `sf_machine_chebyshev_u_c`, the raw-double real and complex recurrences `U_0 = 1`, `U_1 = 2x`, `U_k = 2x U_{k-1} - U_{k-2}`. Otherwise it runs the same recurrence in MPFR. (3) An exact integer order with an exact or symbolic `x` builds the explicit polynomial from FLINT's `fmpz_poly_chebyshev_u` (the integer recurrence when FLINT is absent). (4) A non-integer or complex order with an inexact argument is built as `(v+1) Hypergeometric2F1[-v, v+2, 3/2, (1-x)/2]` when `|(1-x)/2| < 1`, and as `Sin[(v+1) ArcCos[x]]/Sin[ArcCos[x]]` elsewhere.

**Data structures.** These match `ChebyshevT`: `mpz_t` coefficient arrays for the exact polynomial, MPFR or `ncpx` values for the recurrence, and plain `double`s on the machine path. The binary ND kernel `ndk_ChebyshevU_c` (`src/ndkernels.c`) calls the same helpers and has a complex arm, so a packed complex buffer is not materialised. `Compile[]` lowers the head at scalar and array shapes through that registered kernel.

**Complexity / limits.** The numeric recurrence costs `O(n)` operations. The caps are the same as for `ChebyshevT`: degree 2000 for the explicit polynomial and |n| ≤ 2·10⁶ for the recurrence. An exact fractional order with a symbolic argument stays unevaluated, and symbolic series and half-integer closed forms are not implemented. The derivative `((-1-n) U_{n-1} + n x U_n)/(x² - 1)` is produced by `src/calculus/deriv.c`.
