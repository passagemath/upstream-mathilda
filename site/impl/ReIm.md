---
source: src/complex.c
---
**Algorithm.** `builtin_reim` returns `{Re[z], Im[z]}` as a two-element `List` and follows the same order of cases as `builtin_re` and `builtin_im`. A call to a head that is real by construction (`Re`, `Im`, `Abs`, `Arg`) gives `{f[z], 0}`. A `Complex[re, im]` literal gives `{re, im}`. A real number leaf (an Integer, bigint, machine Real, MPFR, or Rational with int64 or bignum parts) gives `{x, 0}`. Otherwise the shared `complex_decompose` Plus/Times walk supplies both parts when they are concretely numeric. Anything else returns `NULL`. Unlike `Re` and `Im`, `ReIm` has no special case for real `AlgebraicNumber` or `Root` objects.

**Data structures.** The result is a fresh `List` of the two parts. Because `ReIm` is `Listable`, a list argument gives a list of pairs.

**Complexity / limits.** The cost is the same as one `Re` plus one `Im`, and the decomposition is done only once. `ReIm` has no packed or NDArray kernel and no `Compile[]` lowering. `ReIm[NDArray[...]]` stays unevaluated, a packed complex list is unpacked before the call, and a compiled body containing `ReIm` falls back to the interpreter. The head is recorded in the fast-path sweep's `OFF_BUFFER` backlog (`tools/nd_fastpath_sweep.py`). Use `{Re[v], Im[v]}` inside `Compile[]` or on large arrays. Symbolic `z` stays unevaluated, as do exact real constants such as `ReIm[Pi]`.
