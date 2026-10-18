---
source:
  - src/special_functions/expintegrale.c
  - src/special_functions/sf_machine_expint.c
references:
  - "DLMF §8.19 — the generalized exponential integral (8.19.8 power series, 8.19.12 recurrence, 8.19.17 continued fraction)."
  - "W. H. Press et al., *Numerical Recipes*, §6.3 (expint) — the series / continued-fraction split for E_n."
  - "S. L. Moshier, *Methods and Programs for Mathematical Functions* (1989), the Cephes `expn` routine."
---
**Algorithm.** `builtin_expintegrale` works in layers. (1) The **machine-scalar gate** (`sf_machine_scalar_args`) sends a machine-precision call with a real integer order to the `double` kernels in `sf_machine_expint.c`. For real `x > 0`, `sf_machine_expn` handles each case separately. `n = 0` is `e^{-x}/x`, and `n = 1` is the Cody–Thacher `E1` kernel. For `x ≤ 1` and `n ≤ 24` it runs the upward recurrence `E_{k+1} = (e^{-x} - x E_k)/k` from `E1`, which measured about 2 ulp. For `x ≤ 1` and larger `n` it sums the power series (DLMF 8.19.8) with the `ψ(n) - log x` term. For `x > 1` it evaluates the continued fraction 8.19.17 backward from the fixed depth `105/x + 12`, which avoids the error that the Lentz and Wallis forms accumulate. A complex `z`, or a real `x < 0` (where the value lies on the upper side of the cut), goes to `sf_machine_expn_c`. That routine uses the complex continued fraction when `|z| + Re z > 1`, with depth `240/s + 8` capped at 400, and the series elsewhere. (2) Exact rules: `E_0(z) = E^(-z)/z`, and `E_n(0) = 1/(n-1)` for integer `n ≥ 2`. (3) Everything else stays symbolic.

**Data structures.** The kernels use plain `double`s, with complex values spelled out as (re, im) pairs (`ExC`) and no `<complex.h>`. The binary ND kernel `ndk_ExpIntegralE_c` (`src/ndkernels.c`, `REG_B`) calls the same two functions. Packed and `NDArray` buffers therefore stay on the buffer, and `Compile[]` lowers the head at scalar and array shapes. Scalar, packed, `NDArray` and compiled results agree bit for bit. Measured accuracy against mpmath has a median of 0.4 ulp and a maximum of 4.3 ulp on the real line.

**Complexity / limits.** Each evaluation takes `O(1)` to `O(100)` floating-point steps. There is **no arbitrary-precision path**, so `N[ExpIntegralE[2, 3/2], 30]` stays unevaluated. A non-integer order, or an order above 10⁶, is left unevaluated even at machine precision, and so is `ExpIntegralE[1, 0]`. A machine-real underflow returns `0.`. The derivative `-E_{n-1}` lives in `src/calculus/deriv.c`, and the Mellin transform is handled by `Integrate`'s operational-calculus layer.
