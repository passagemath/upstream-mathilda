---
source:
  - src/complex.c
  - src/poly/flint_qqbar.c
  - src/ndkernels.c
  - src/compile/compile_infer.c
  - src/compile/compile_emit_mathfn.c
---
**Algorithm.** `builtin_re` tries the following cases in order. A call to a head that is real by construction (`Re`, `Im`, `Abs`, `Arg`) is returned unchanged, so `Re[Re[z]]` is `Re[z]`. A real `AlgebraicNumber` or `Root` (tested with FLINT `qqbar`, `flint_qqbar_is_real`) is returned as itself; a non-real one stays wrapped, as in Mathematica. A `Complex[re, im]` literal gives `re`. A real number leaf (`is_real_number_leaf`: an Integer, bigint, machine Real, MPFR, or Rational with int64 or bignum parts) gives itself. Otherwise `complex_decompose` walks the Plus/Times structure, multiplying Complex factors through. When both parts it finds are concretely numeric (`is_numeric_real`, judged by machine `N`), the real part is returned; this is what reduces `Re[1/(2 + 3 I)]` to `2/13`. Anything else returns `NULL` and stays unevaluated.

**Data structures.** The decomposition builds and evaluates fresh `re`/`im` trees. No other state is kept. On the packed and NDArray surfaces `Re` is a projection kernel (`UK_PROJ(Re, ar)`): a complex buffer gives a real buffer, and real data passes through. In `Compile[]` it lowers at scalar and array shapes, and over integers it is the identity.

**Complexity / limits.** The leaf cases are O(1). The decomposition is linear in the size of the expression, plus one machine-`N` reality check per part. Symbols are not assumed real, so `Re[x + I y]` stays unevaluated (use `ComplexExpand` for that convention). A gap: an exact real constant with no Complex literal inside, such as `Re[Pi]`, is not yet recognised and stays unevaluated (Mathematica gives `Pi`).
