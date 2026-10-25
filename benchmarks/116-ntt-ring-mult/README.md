# 116 — NTT ring multiplication (group F)

Negacyclic polynomial multiplication in R_q = Z_q[x]/(xⁿ+1), the inner loop
of every Ring/Module-LWE scheme. Two rings are used: NewHope's q = 12289,
n = 1024, and a Kyber-sized ring with a **complete** NTT, q = 7681, n = 256.
Kyber's own q = 3329 has 512 ∤ q−1, so it admits only an incomplete NTT. The
`app:` rows are the number-theoretic transform written in the language: an
iterative radix-2 Cooley–Tukey NTT over a bit-reversed input with precomputed
twiddles, and the negacyclic wrap done by ψ-pre-multiplication (ω = ψ²),
pointwise product, inverse NTT and n⁻¹ψ⁻ⁱ post-multiplication. In WL each
butterfly stage is the idiomatic whole-array form (`Partition` / `Transpose` /
`Mod` over the packed list). Python and GP run the same butterflies as loops.
The `prim:` rows are each system's builtins on the same products.

| case | what |
|---|---|
| `app: NTT negacyclic mult q=12289 n=1024 x20` | 20 products, 3 NTTs each |
| `prim: ListConvolve negacyclic q=12289 n=1024 x4` | `ListConvolve[a, b, {1,-1}, 0]` folded by xⁿ = −1, `Mod q`; numpy `convolve` (int64); PARI integer `Polrev` product |
| `app: NTT negacyclic mult q=7681 n=256 x100` | 100 products |
| `prim: ListConvolve negacyclic q=7681 n=256 x100` | same batch as the row above, so the checks coincide |
| `prim: PolynomialRemainder negacyclic q=7681 n=128 x2` | `PolynomialMod[PolynomialRemainder[Expand[A B], x^n+1, x], q]` → `CoefficientList[., x, n]`; sympy integer `Poly(...).rem`; PARI `Mod(1,q)`-coefficient polynomials `% (x^n+1)` |

Inputs come from `smallList` (coefficients in [0, q)). The check is
Σ_batch Σᵢ i·cᵢ mod 10⁹+7, weighted so that a permuted result cannot pass.
The Python file also asserts that convolution equals NTT and that the sympy
route equals convolution.

## Findings

Times in ms, one run (Mathematica = wolframscript 14):

| case | Mathilda | Python | PARI | Mathematica |
|---|---:|---:|---:|---:|
| NTT n=1024 x20 | 880 | 110 (pure Python) | 688 | 19.4 |
| ListConvolve n=1024 x4 | 1436 | 1.9 (numpy) | 3 | 0.21 |
| NTT n=256 x100 | 1009 | 103 | 658 | 51.6 |
| ListConvolve n=256 x100 | 2263 | 3.7 | 21 | 2.9 |
| PolynomialRemainder n=128 x2 | 3559 **CHECK-FAIL** | 40 (sympy) | <1 | 19.2 |

- **WRONG ANSWER: `CoefficientList[p, x, n]` (length argument) returns
  unevaluated.** The PolynomialRemainder row therefore CHECK-FAILs. Reproducer:
  `CoefficientList[1 + 2 x, x, 4]` gives `{1, 2, 0, 0}` in WL and is left
  unevaluated by Mathilda. The truncating form `CoefficientList[1 + 2x + 3x^5,
  x, 4]` (WL: `{1, 2, 0, 0}`) and the multivariate form
  `CoefficientList[1 + 2 x y, {x, y}, {2, 3}]` are also unevaluated.
- **`PolynomialRemainder[p, x^n + 1, x]` is ~O(n³).** For dense p of degree
  2n−2 (q = 7681 coefficients) it takes 0.028 s at n = 32, 0.20 s at 64, and
  1.6 s at 128. At n = 1024 two products did not finish in 9 minutes, so the
  row runs at n = 128. `Expand[A B]` (0.0013 s at n = 128) and `PolynomialMod`
  are cheap; the division is the cost. Dividing by a monic binomial needs only
  O(n) per quotient term (O(n²) overall, or just a fold).
- **Integer `ListConvolve` has no packed fast path.** It takes 337 ms per
  n = 1024 convolution, about 1600× Mathematica and about 180× numpy, and it
  returns an **unpacked** list (`PackedArrayQ` → `False`). The same call on
  reals (`N[a]`) takes 6 ms. The cost is quadratic with a large constant:
  ~330 ns per multiply-add.
- **`Map` over a packed array hands each row to the function UNPACKED.**
  This is what makes the NTT rows slow. `pairsA` (a packed 20×2×1024 tensor) is
  processed by the natural `nttMul[#[[1]], #[[2]], ctxA] & /@ pairsA`, and
  every `#` arrives as a boxed list. All the vectorised arithmetic inside the
  NTT then runs on unpacked lists, which is ~100× slower per op (`Mod[r r, q]`
  takes 0.51 ms vs 0.005 ms). The same 20 products cost **937 ms via `Map`**
  and **98 ms via `Table[nttMul[pairsA[[k, 1]], pairsA[[k, 2]], ctxA], {k, 20}]`**
  (`Part` keeps the packing). One NTT on a packed input is 1.3 ms, against
  ~0.3 ms for Mathematica. Reproducer:
  ```
  Map[PackedArrayQ, Partition[Range[12], 4]]     (* Mathilda {False,False,False}; WL {True,True,True} *)
  m = Partition[Range[12], 4]; Table[PackedArrayQ[m[[k]]], {k, 3}]      (* Mathilda {True,True,True} *)
  ```
  A related, smaller point: `smallList` itself returns an unpacked list,
  because `Quotient[lcgStates[...], 2^32]` of 64-bit bigints yields machine-size
  values that are never auto-packed (`PackedArrayQ` → `False`, and the same
  after `Mod[., q]`). `Table` repacks them when building `pairsA`.
  Inside the vectorised stage, the one slow op is `Join[u, v, 2]` on many short
  rows (0.27 ms at h = 1 vs Mathematica's 0.002 ms).
- **Packed-array status:** on a packed input, `Mod` over the packed list stays
  packed, and so do `Partition`, the `Transpose`·twiddle product,
  `Join[., ., 2]`, `Flatten` and the final `nttMul` result. `ListConvolve`
  unpacks its result, and `Map` unpacks its argument rows (see above).
  Note that `PackedArrayQ` lives in ``System` `` in Mathilda;
  ``Developer`PackedArrayQ`` is absent there, although it is the WL spelling.
- PARI's GP-loop NTT (≈660 ms) is slower than pure Python (≈105 ms). Two PARI
  rows fall below the 20 ms floor (ListConvolve n=1024 3 ms, PolynomialRemainder
  <1 ms). They cannot be looped up, because the same loop would put Mathilda
  far over its 5 s budget.
