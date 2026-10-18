---
source:
  - src/complex.c
  - src/interval.c
  - src/ndkernels.c
  - src/pack.c
  - src/compile/compile_infer.c
  - src/compile/compile_emit_mathfn.c
---
**Algorithm.** `builtin_sign` maps an `Interval` monotonically through `interval_thread_monotone`, so `Sign[Interval[{-1, 2}]]` is `Interval[{-1, 1}]`. For a real number it returns the **exact** integer -1, 0 or 1 whatever the input type: a direct comparison for an int64 or a machine Real (so `Sign[-2.5]` is `-1`, not `-1.`), `mpz_sgn` for a bigint, `mpfr_sgn` for MPFR, and `expr_numeric_sign` for a Rational, including one with bignum parts such as `1/10^25`. For a `Complex[re, im]` with numeric parts it returns the unit vector `z/Abs[z]`, and `0` for `0 + 0 I`. Exact parts build `z · Abs[z]^-1` and let the evaluator fold it, so `Sign[3 + 4 I]` is `3/5 + 4/5 I` and `Sign[2 - 2 I]` is `(1/2 - 1/2 I) Sqrt[2]`. MPFR parts are divided by `mpfr_hypot` at the combined precision. Anything else returns `NULL`.

**Data structures.** The routine allocates no temporaries beyond the MPFR scratch values. On the packed and NDArray surfaces `NDKU_Sign` is a *narrowing* kernel: a real or int64 buffer gives an int64 buffer of -1, 0 and 1 (`-0.` gives 0; a non-finite element abandons the buffer), and a complex buffer gives `z/|z|`. `Sign` is on `pack.c`'s `AWARE` and `INT64_OK` lists. In `Compile[]` the result type is Integer for an Integer or Real argument and Complex for a Complex one, at scalar and array shapes.

**Complexity / limits.** Each case is O(1). `Sign[x]` stays unevaluated for symbolic `x`. A gap: an exact real constant that is not a rational, such as `Sign[Pi]`, `Sign[Sqrt[2]]` or `Sign[Pi - 3]`, also stays unevaluated at present, whereas Mathematica returns `1`. Use `Sign[N[expr]]`, which is still exact.
