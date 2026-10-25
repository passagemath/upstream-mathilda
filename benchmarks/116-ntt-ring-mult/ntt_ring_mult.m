(* ==========================================================================
   Experiment 116 -- NTT ring multiplication (group F, lattice crypto core)
   ==========================================================================
   WHAT IT MEASURES.  Negacyclic polynomial multiplication in
   R_q = Z_q[x]/(x^n + 1), the inner loop of every Ring/Module-LWE scheme
   (NewHope: q = 12289, n = 1024; a Kyber-sized ring with a COMPLETE NTT:
   q = 7681, n = 256 -- Kyber's own q = 3329 only admits an incomplete NTT).

   "app:"  the number-theoretic transform written in the language: an
           iterative radix-2 Cooley-Tukey NTT over a bit-reversed input with
           precomputed twiddles; negacyclic wrap by pre-multiplying with
           powers of a primitive 2n-th root psi (omega = psi^2), pointwise
           product, inverse NTT, post-multiply by n^-1 psi^-i.  The WL form is
           the idiomatic whole-array one: each butterfly stage is a vectorised
           Partition / Transpose / Mod over the packed integer list.
   "prim:" the system's builtins on the same products: the schoolbook
           convolution ListConvolve[a, b, {1, -1}, 0] folded by x^n = -1 and
           reduced mod q; and the symbolic route
           PolynomialMod[PolynomialRemainder[A B, x^n + 1, x], q].

   Inputs: N pairs of coefficient vectors in [0, q) from crypto_data.m's
   smallList.  Check: sum over the batch of sum_i i * c_i (weighted, so a
   permuted result cannot pass), mod $CheckMod.  Every route computes the
   same unique product (the n = 256 NTT and ListConvolve rows share their
   batch, so their checks coincide).  PolynomialRemainder is cubic in n in
   Mathilda, so its row runs in the smaller ring n = 128 -- see README.
   ========================================================================== *)

Get["../harness.m"];
Get["../crypto_data.m"];

require[{"ListConvolve", "PolynomialRemainder", "PolynomialMod", "PrimitiveRoot",
         "PowerMod", "IntegerDigits", "Partition", "Transpose"}];

(* ---- ring setup (untimed): roots, twiddles, bit reversal ---------------- *)

nttSetup[q_, n_] := Module[{g = PrimitiveRoot[q], psi, psiInv, w, wInv, k = Log2[n]},
  psi = PowerMod[g, (q - 1)/(2 n), q]; psiInv = PowerMod[psi, -1, q];
  w = PowerMod[psi, 2, q]; wInv = PowerMod[w, -1, q];
  <|"q" -> q, "n" -> n,
    "perm" -> Table[FromDigits[Reverse[IntegerDigits[i, 2, k]], 2] + 1, {i, 0, n - 1}],
    "tw" -> PowerMod[w, Range[0, n/2 - 1], q],
    "twInv" -> PowerMod[wInv, Range[0, n/2 - 1], q],
    "psiPow" -> PowerMod[psi, Range[0, n - 1], q],
    "postPow" -> Mod[PowerMod[n, -1, q] PowerMod[psiInv, Range[0, n - 1], q], q]|>];

(* Cyclic NTT, iterative Cooley-Tukey on a bit-reversed copy.  Stage with
   half-width h: blocks of 2h, v-half scaled by omega^(j n/(2h)), j < h. *)
ntt[a_, tw_, ctx_] := Module[{q = ctx["q"], n = ctx["n"], x = a[[ctx["perm"]]], h = 1,
                             blk, u, v, t},
  While[h < n,
    blk = Partition[x, 2 h];
    t = tw[[1 ;; n/2 ;; n/(2 h)]];
    u = blk[[All, 1 ;; h]];
    v = Mod[Transpose[Transpose[blk[[All, h + 1 ;; 2 h]]] t], q];
    x = Flatten[Join[Mod[u + v, q], Mod[u - v, q], 2]];
    h *= 2];
  x];

nttMul[a_, b_, ctx_] := Module[{q = ctx["q"], fa, fb},
  fa = ntt[Mod[a ctx["psiPow"], q], ctx["tw"], ctx];
  fb = ntt[Mod[b ctx["psiPow"], q], ctx["tw"], ctx];
  Mod[ntt[Mod[fa fb, q], ctx["twInv"], ctx] ctx["postPow"], q]];

(* Schoolbook: the linear convolution, folded by x^n = -1. *)
convMul[a_, b_, q_, n_] := Module[{full = ListConvolve[a, b, {1, -1}, 0]},
  Mod[full[[1 ;; n]] - Append[full[[n + 1 ;; 2 n - 1]], 0], q]];

(* Symbolic: as polynomials in x. *)
polyMul[a_, b_, q_, n_] := Module[{pa = a . x^Range[0, n - 1], pb = b . x^Range[0, n - 1]},
  CoefficientList[PolynomialMod[PolynomialRemainder[Expand[pa pb], x^n + 1, x], q], x, n]];

digest[cs_] := Mod[Total[(# . Range[Length[#]]) & /@ cs], $CheckMod];

(* ---- inputs ------------------------------------------------------------- *)

ctxA = nttSetup[12289, 1024];
ctxB = nttSetup[7681, 256];
pairsA = Table[{smallList[1600 + 2 k, 1024, 12289], smallList[1601 + 2 k, 1024, 12289]}, {k, 1, 20}];
pairsB = Table[{smallList[1700 + 2 k, 256, 7681], smallList[1701 + 2 k, 256, 7681]}, {k, 1, 100}];

(* ---- q = 12289, n = 1024 ------------------------------------------------ *)

bench["app: NTT negacyclic mult q=12289 n=1024 x20",
      cA1 = nttMul[#[[1]], #[[2]], ctxA] & /@ pairsA;];
check["app: NTT negacyclic mult q=12289 n=1024 x20", digest[cA1]];

bench["prim: ListConvolve negacyclic q=12289 n=1024 x4",
      cA2 = convMul[#[[1]], #[[2]], 12289, 1024] & /@ pairsA[[1 ;; 4]];];
check["prim: ListConvolve negacyclic q=12289 n=1024 x4", digest[cA2]];

(* ---- q = 7681, n = 256 -------------------------------------------------- *)

bench["app: NTT negacyclic mult q=7681 n=256 x100",
      cB1 = nttMul[#[[1]], #[[2]], ctxB] & /@ pairsB;];
check["app: NTT negacyclic mult q=7681 n=256 x100", digest[cB1]];

bench["prim: ListConvolve negacyclic q=7681 n=256 x100",
      cB2 = convMul[#[[1]], #[[2]], 7681, 256] & /@ pairsB;];
check["prim: ListConvolve negacyclic q=7681 n=256 x100", digest[cB2]];

(* ---- the symbolic route, at a size it finishes ------------------------- *)

(* PolynomialRemainder[., x^n + 1, x] scales ~n^3 in Mathilda (0.03 s at
   n = 32, 0.2 s at 64, 1.6 s at 128; n = 1024 x2 did not finish in 9 min),
   so this case runs in the smaller ring Z_7681[x]/(x^128 + 1).  See README. *)
pairsC = Table[{smallList[1900 + 2 k, 128, 7681], smallList[1901 + 2 k, 128, 7681]}, {k, 1, 2}];
benchOnce["prim: PolynomialRemainder negacyclic q=7681 n=128 x2",
      cC3 = polyMul[#[[1]], #[[2]], 7681, 128] & /@ pairsC;];
check["prim: PolynomialRemainder negacyclic q=7681 n=128 x2", digest[cC3]];
