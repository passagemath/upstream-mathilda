---
source:
  - src/special_functions/whittaker.c
  - src/special_functions/sf_machine.c
  - src/special_functions/sf_machine_hyper.c
  - src/special_functions/hypergeometricu.c
  - src/calculus/deriv.c
  - src/calculus/series.c
references:
  - "DLMF §13.14 — Whittaker functions, W_{k,m}(z) = e^(−z/2) z^(m+1/2) U(m−k+1/2, 1+2m, z) (§13.14.3)."
  - "DLMF §13.14.31 (the derivative used for D); §13.2.9 (the integer-b log series behind the machine kernel)."
---
**Algorithm.** `builtin_whittaker_w` uses the same three stages as WhittakerM, with the U kernels in place of the 1F1 kernels.

1. `whit_ndarray` threads the kernels over a buffer z.
2. `whit_machine` handles a scalar machine call with `sf_machine_whittaker_w`. That kernel folds m to |m|, which is valid because W is even in m, and so b = 2|m|+1 ≥ 1. It then returns `e^(−z/2) z^(|m|+1/2) · sf_machine_hyperu(|m|−k+1/2, 2|m|+1, z)`, with `sf_machine_whittaker_w_c` for z off the positive axis.
3. `whittaker_w_eval` settles z = 0 from the leading exponent `1/2 − |Re m|`: the value is `0` below 1/2, `ComplexInfinity` above it, and `1/Gamma[1−k]` for a real m = ±1/2. Any other inexact input is evaluated through the template `Exp[−z/2] z^(m+1/2) HypergeometricU[m−k+1/2, 2m+1, z]`, so it inherits the HypergeometricU arb path.

In the common case of integer or half-integer m, the U kernel's log series advances its three digammas by recurrence instead of calling digamma at every term. `D` comes from DLMF 13.14.31. The origin and infinity series in `series.c` cover only the terminating case `U(−n, b, z) = (−1)^n n! LaguerreL[n, b−1, z]`.

**Data structures.** The machine scalar and buffer drivers are shared with WhittakerM through the `WhitReal` and `WhitCplx` function pointers. Buffers stay real when every point is real and are promoted to complex otherwise. `REG_N(WhittakerW)` registers the kernel for the Compile VM.

**Complexity / limits.** Machine cost is one `sf_machine_hyperu` call plus the prefactor. The docs report the packed path running about 1.8× faster than SciPy's `exp·z^(m+½)·hyperu`. `Compile[]` lowers the scalar shape but not rank-1 arrays: `CompileDiagnostics[{{v,_Real,1}}, WhittakerW[1,1/2,v]]` reports `Compiled -> False`. The general two-branch origin series and the general asymptotic are not produced. The HypergeometricU terminating-case defect carries over: `WhittakerW[2, 1/2, 1.]` returns `0.0`, while the correct value, −e^(−1/2) = −0.606531, is obtained with `N[WhittakerW[2, 1/2, 1]]`.
