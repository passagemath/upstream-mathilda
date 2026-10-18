---
source:
  - src/complex.c
  - src/ndkernels.c
  - src/compile/compile_infer.c
  - src/compile/compile_emit_mathfn.c
---
**Algorithm.** `builtin_conjugate` first checks its arity. Any argument count other than one emits `Conjugate::argx` (through `mth_message`) and leaves the call unevaluated. Then, in order: `Conjugate[Conjugate[z]]` gives `z` (it is an involution). A call to a head that is real by construction (`Re`, `Im`, `Abs`, `Arg`) is a fixed point. A `Complex[re, im]` literal gives `Complex[re, -im]`. An Integer, machine Real or Rational is returned unchanged. So is any expression that `is_numeric_real` accepts, meaning it reduces to a real number under machine `N`; this covers `Sqrt[2] Log[3]` and `Pi - 4`. Finally, an expression that `complex_decompose` splits into numeric real and imaginary parts is rebuilt as `re - I im` and evaluated, which is how `Conjugate[(2 + I)/(1 - 3 I)]` becomes `-1/10 - 7/10 I`. Anything else returns `NULL`.

**Data structures.** The routine allocates no state beyond the decomposition's fresh `re`/`im` trees. On the packed and NDArray surfaces `Conjugate` is a closed kernel (`UK_CLOSED(Conjugate, x, conj(z))`): the identity on real buffers and `conj` on complex ones. In `Compile[]` its result type is its argument's type, at scalar and array shapes.

**Complexity / limits.** The leaf cases are O(1), and the reality test costs one machine `N`. Symbols are not assumed real, so `Conjugate[z]` stays unevaluated, as in Mathematica. Note that the reality test for non-literal expressions is numerical, at machine precision.
