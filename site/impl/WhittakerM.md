---
source:
  - src/special_functions/whittaker.c
  - src/special_functions/sf_machine.c
  - src/special_functions/sf_machine_hyper.c
  - src/calculus/deriv.c
  - src/calculus/series.c
references:
  - "DLMF §13.14 — Whittaker functions, M_{k,m}(z) = e^(−z/2) z^(m+1/2) M(m−k+1/2, 1+2m, z) (§13.14.2)."
  - "DLMF §13.4.27 (the derivative recurrence used for D)."
---
**Algorithm.** `builtin_whittaker_m` runs three stages in order.

1. `whit_ndarray` threads the machine kernels over a buffer z, holding k and m constant.
2. `whit_machine` handles a scalar machine call. It applies when k and m are exact or machine reals and z is a machine real or complex number, and it calls `sf_machine_whittaker_m`. That kernel computes `e^(−z/2) z^(m+1/2) · sf_machine_hyper1f1(m−k+1/2, 2m+1, z)` and declines inside the 1F1 when cancellation is too large. For a real z < 0, or a complex z, the call goes to the complex twin `sf_machine_whittaker_m_c`, which takes the principal branch.
3. `whittaker_m_eval` returns the exact `0` at z = 0 when Re m > −1/2. For any other inexact input, `whittaker_m_numeric_hyp` builds `Exp[−z/2] z^(m+1/2) Hypergeometric1F1[m−k+1/2, 2m+1, z]` and evaluates it. The result inherits 1F1's MPFR and complex numerics, and it is kept only if it collapses to a number.

`D` in z comes from DLMF 13.4.27, implemented in `deriv.c`. `series.c` emits the Puiseux origin series `z^(1/2+m)(...)`, and at `Infinity` it emits only the terminating closed form, which arises when `m−k+1/2` is a non-positive integer.

**Data structures.** The machine path keeps no state; it takes doubles in and returns doubles or complex doubles. The function pointer typedefs `WhitReal` and `WhitCplx` let WhittakerM and WhittakerW share the scalar and buffer drivers. A real buffer gets a fast pass that stays real. A point on the cut promotes the whole output buffer to complex. `REG_N(WhittakerM)` registers the n-ary kernel for the Compile VM.

**Complexity / limits.** Machine cost is one 1F1 evaluation plus `exp` and `pow`. The documented timings are 0.2 µs for a real argument and 0.4 µs for a complex one. Precision above machine costs whatever `Hypergeometric1F1` costs. `Compile[]` lowers scalar arguments (real and complex). The rank-1 array shape is not lowered: `CompileDiagnostics[{{v,_Real,1}}, WhittakerM[1,1/2,v]]` reports `Compiled -> False`. The general Stokes asymptotic at infinity is not implemented, and neither is `FunctionExpand`. An `Interval` argument stays unevaluated, because the derivative refers back to WhittakerM.
