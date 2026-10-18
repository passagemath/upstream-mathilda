---
source:
  - src/complex.c
  - src/ndkernels.c
  - src/compile/compile_infer.c
  - src/compile/compile_emit_mathfn.c
---
**Algorithm.** `builtin_arg` returns the principal argument in `(-Pi, Pi]`, trying the following cases in order. (1) A pure MPFR real gives the exact `0` or `Pi` by sign. (2) An exact real leaf (an Integer, bigint or Rational, bignum parts included) gives `0` or `Pi` from its exact sign, so even `Arg[-1/10^400]` is `Pi`. (3) **Root-of-unity constants.** `arg_turns_of` reads a product of real-rational powers and pure-imaginary Gaussians, such as `(-1)^(1/6)`, `-(-1)^(5/6)` or `I (-1)^(1/4)`. It adds their half-turns mod 2 (`[c<0] + Σ eᵢ[bᵢ<0]`), with denominators capped at 2^40, and `arg_from_turns` reduces the total into `(-Pi, Pi]` as an exact rational multiple of `Pi`. (4) A `Complex[re, im]`: with an MPFR part it uses `mpfr_atan2` at the combined precision. With exact parts it returns exact multiples of `Pi` for the axes and diagonals (`0`, `Pi`, `±Pi/2`, `±Pi/4`, `±3 Pi/4`), and otherwise the symbolic two-argument `ArcTan[re, im]`. With a machine part it uses libm `atan2(im, re)`. (5) A machine Real gives `atan2(0, x)`. Anything else, including symbols, returns `NULL`.

**Data structures.** The routine allocates no structures. A rational turn count, an int64 pair, carries the root-of-unity case. On the packed and NDArray surfaces `Arg` is a projection kernel (`UK_PROJ(Arg, atan2(ai, ar))`) that writes a real buffer. In `Compile[]` it lowers at scalar and array shapes: an integer argument stays integral, and the `Pi` case is handed back to the interpreter.

**Complexity / limits.** Each case is O(1), and the root-of-unity reader is linear in the number of factors. Generic exact Gaussians stay as `ArcTan[re, im]`, which Mathematica would print as `ArcTan[3/2]` for `Arg[2 + 3 I]`; `N` of either gives the same number. A gap: an exact real constant that is not rational, such as `Arg[Pi - 4]`, stays unevaluated (Mathematica gives `Pi`). A product such as `-Sqrt[2]` is covered by the power-product reader.
