---
source:
  - src/special_functions/hypergeometricu.c
  - src/special_functions/sf_machine.c
  - src/special_functions/sf_machine_hyper.c
  - src/flint_num_bridge.c
  - src/calculus/deriv.c
references:
  - "DLMF §13.2 — the Kummer and Tricomi functions; §13.2.9 — U at integer b (the logarithmic case); §13.2.42 — the connection between U and 1F1."
  - "DLMF §13.4.4 — the integral representation U = (1/Γ(a)) ∫₀^∞ e^(−zt) t^(a−1) (1+t)^(b−a−1) dt; §13.7 — the asymptotic expansion; §13.3.22 — the z-derivative."
---
**Algorithm.** `builtin_hypergeometric_u` checks argument shapes in the following order.

1. A buffer z (`NDArray` or packed list) goes to `hypergeometricu_ndarray`. That function never declines a buffer: it either returns a packed result or hands off to `ndarray_delist_and_reeval`.
2. A machine-scalar call goes to `hu_machine`. The call qualifies when a and b are real, z is real or complex, and at least one argument is a machine Real. This gate comes before the exact reductions.
3. `hypergeometricu_eval` applies the exact rules `U(0, b, z) = 1` and `U(a, a+1, z) = z^(−a)`. The second rule tests whether `b − a − 1` evaluates to the literal 0. Any other inexact input goes to FLINT's `acb_hypgeom_u` through `flint_num_hypergeometric_u`. When every inexact argument was machine precision, the result is demoted back to machine doubles.

The real kernel `sf_machine_hyperu` requires z > 0 and tries four regimes in turn:

- `hyperu_asymptotic`, the series z^(−a) 2F0(a, a−b+1; ; −1/z) with TwoSum compensation. It returns exactly when the series terminates. A truncated sum is accepted only if its least term is below half an ulp of the sum.
- `hyperu_integer_b`, the DLMF 13.2.9 log series for integer b ≥ 1. Its three digammas are each seeded once and then advanced with ψ(x+1) = ψ(x) + 1/x. For b ≤ 0 the kernel first applies the fold `U = z^(1−b) U(a−b+1, 2−b, z)`.
- `hyperu_connection`, the Kummer formula with two 1F1 terms, used for non-integer b.
- `sf_machine_hyperu_quad`, an exp-sinh double-exponential quadrature of the integral representation (step 0.125), used for a > 0. For a ≤ 0 it is applied to the Kummer image instead.

The integer-b and connection regimes decline when cancellation amplifies rounding by more than a factor of 4. A real z < 0 and any complex z go to `sf_machine_hyperu_c`, which returns the value on the principal branch.

**Data structures.** Both kernels work on scalar doubles only. `hu_point` dispatches one point between the real and complex kernels. In the buffer loop, a real `NDT_FLOAT64` buffer gets a fast pass. When any point lands on the branch cut, the loop switches to an interleaved re/im scratch buffer and promotes the output to complex64 (or complex32). `REG_N(HypergeometricU)` in `src/ndkernels.c` exposes the real kernel to the Compile VM.

**Complexity / limits.** The asymptotic regime runs up to about z terms (cap 400). The log series is linear in its term count. The quadrature uses 50–90 nodes and gives up when a > 160. `Compile[]` lowers the scalar shape only: `CompileDiagnostics[{{v,_Real,1}}, HypergeometricU[1,2,v]]` reports `Compiled -> False`. The terminating Laguerre case (a a non-positive integer) is not expanded symbolically. There is no Series at the origin and no `FunctionExpand` rule. One defect: when a is a non-positive integer, b is an integer and z is small (for example `HypergeometricU[-1, 2, 1.]`), the asymptotic regime stops at its growing first term. `hyperu_integer_b` then returns 0 because both of its 1/Γ prefactors vanish. That result is wrong.
