---
source:
  - src/complex.c
  - src/poly/flint_qqbar.c
  - src/ndkernels.c
  - src/compile/compile_infer.c
  - src/compile/compile_emit_mathfn.c
---
**Algorithm.** `builtin_im` mirrors `builtin_re`. It returns `0` for a call to a head that is real by construction (`Re`, `Im`, `Abs`, `Arg`), for a real `AlgebraicNumber` or `Root` (FLINT `qqbar` reality test), and for any real number leaf (an Integer, bigint, machine Real, MPFR, or Rational with int64 or bignum parts). A `Complex[re, im]` literal gives `im`. Otherwise `complex_decompose`, a recursive Plus/Times walk that multiplies Complex factors through, supplies the imaginary part when both decomposed parts are concretely numeric (`is_numeric_real`). This is how `Im[Sqrt[-4]]` is `2` and `Im[Log[-1]]` is `Pi`: those arguments have already evaluated to `2 I` and `I Pi`. Anything else returns `NULL`, including special-function values with no closed form such as `Im[Gamma[1 + I]]`, whose `N` still evaluates.

**Data structures.** The routine works on fresh `re`/`im` trees from the decomposition. On the packed and NDArray surfaces the kernel `NDKU_ImInt` narrows: a real or int64 buffer gives an int64 buffer of exact zeros (matching the List path's exact `0`), and a complex buffer gives the real buffer of imaginary parts. In `Compile[]` `Im` lowers at scalar and array shapes, and it is `0` over integers.

**Complexity / limits.** The leaf cases are O(1), and the decomposition is linear in the size of the expression. Symbols are not assumed real (`Im[x]` stays unevaluated). A gap: an exact real constant with no Complex literal, such as `Im[Pi - 4]`, stays unevaluated at present (Mathematica gives `0`).
