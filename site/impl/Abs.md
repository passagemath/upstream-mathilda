---
source:
  - src/complex.c
  - src/interval.c
  - src/ndkernels.c
  - src/compile/compile_infer.c
  - src/compile/compile_emit_mathfn.c
---
**Algorithm.** `builtin_abs` tries the following cases in order. (1) An `Interval` goes to `interval_apply_function` (`iv_abs`), so `Abs[Interval[{-2, 1}]]` is `Interval[{0, 2}]`. (2) A real `AlgebraicNumber` or `Root` object is compared with 0 through FLINT `qqbar` (`flint_qqbar_compare`) and returned as itself or its negation, so it stays an exact algebraic object. (3) A `Complex[re, im]` literal, or an expression that `complex_decompose` splits into concretely numeric real and imaginary parts (a recursive Plus/Times walk that propagates Complex factors), becomes the modulus. When either part is MPFR this goes through `mpfr_hypot` at the combined precision, which is also stable when the two magnitudes are very different. Otherwise it builds `(re^2 + im^2)^(1/2)` exactly, so `Abs[3 + 4 I]` is `5` and `Abs[Sqrt[2] + Sqrt[3] I]` is `Sqrt[5]`. (4) Real leaves: `mpz_abs` for a bigint, a sign flip for an int64 (`INT64_MIN` is promoted to the bignum 2^63 instead of wrapping), `fabs` for a machine Real, `mpfr_abs` at the input precision, and `|n|/d` for a Rational, including one with bignum parts. Anything else returns `NULL`, so the call stays unevaluated.

**Data structures.** The routine works directly on `Expr` leaves. MPFR scratch values are sized by `numeric_combined_bits`. The packed and NDArray surfaces use `NDKU_AbsInt` (`ndkernels.c`): an exact int64 arm (it abandons the whole buffer on `INT64_MIN`, so the List path answers with a bignum), and a `hypot` arm for real and complex buffers that writes a real buffer. `Abs` is on `pack.c`'s `AWARE` list. In `Compile[]`, `Abs` of an Integer stays an Integer, of a Real a Real, and of a Complex a Real, at scalar and array shapes.

**Complexity / limits.** Each case costs O(1) numeric work, except the `complex_decompose` walk, which is linear in the size of the expression. Symbolic arguments stay unevaluated (`Abs[x]`), as in Mathematica. A gap: a *real* exact constant that is neither a rational nor a `Complex` literal, such as `Abs[Sqrt[2]]`, `Abs[-Pi]` or `Abs[Pi - 4]`, also stays unevaluated at present, whereas Mathematica returns `Sqrt[2]`, `Pi` and `4 - Pi`. Apply `N` first for a numeric value.
