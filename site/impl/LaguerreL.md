---
source:
  - src/special_functions/laguerre_l.c
  - src/special_functions/sf_machine.c
  - src/ndkernels.c
references:
  - "DLMF §18 — classical orthogonal polynomials (§18.9.13 recurrence, §18.5.12 explicit sum); §13.6.19 for L_n^(α) = Binomial[n+α, n] M(-n, α+1, x)."
---
**Algorithm.** `builtin_laguerre_l` accepts `LaguerreL[n, x]` (with `a = 0`) and `LaguerreL[n, a, x]`. An `NDArray` or packed `x` with scalar `n` and `a` goes to `laguerre_l_ndarray`, which threads the machine helper over the buffer and delists to the scalar path if any element declines. Otherwise `laguerre_l_eval` picks a route by shape. (1) An integer order with an inexact argument, at machine precision with `n ≤ 1000`, calls `sf_machine_laguerre_l_assoc` (or its complex twin), the forward recurrence `(k+1) L_{k+1} = (2k+1+a-x) L_k - (k+a) L_{k-1}` in `double`s. Otherwise `lag_numeric_recurrence` runs the same recurrence in MPFR. (2) An exact order with an exact or symbolic `x` and a numeric `a` gives the explicit polynomial with rational coefficients from that recurrence on `mpq_t` coefficient vectors (`lag_poly_numeric_a`). A symbolic `a` gives the expanded sum `Σ (-1)^i Binomial[n+a, n-i] x^i/i!` (`lag_poly_symbolic_a`). (3) A non-integer or complex order with an inexact argument is built as `Binomial[n+a, n] Hypergeometric1F1[-n, a+1, x]` (`lag_numeric_hyp`), which is valid on the whole plane because `L_n` is entire in `x`. (4) The special values `L_0 = 1` and `L_n^{(a)}(0) = Binomial[n+a, n]` hold for a symbolic order.

**Data structures.** `mpq_t` arrays hold the coefficients for numeric `a`, and `Expr` sums hold them for symbolic `a`. The recurrence uses MPFR or `ncpx` values, and plain `double`s on the machine path. The two-argument form has the binary ND kernel `ndk_LaguerreL_c`, real and complex. The three-argument form has the n-ary kernel `ndk_LaguerreL_n` (`REG_N`), which `Compile[]` uses for the scalar shape.

**Complexity / limits.** The recurrence costs `O(n)` operations. The explicit polynomial is capped at `n ≤ 2000` (`LAG_POLY_CAP`), or `n ≤ 200` when `a` is symbolic (`LAG_SYM_POLY_CAP`), and the numeric recurrence at `2·10⁶`. `Compile[]` does not lower the three-argument form at array shape (`CompileDiagnostics` reports `Compiled -> False`), although the interpreter keeps it on the buffer. An exact fractional order with a symbolic argument stays unevaluated, and symbolic series and the expansion at infinity are not implemented. The derivative `-L_{n-1}^{(a+1)}` is in `src/calculus/deriv.c`.
