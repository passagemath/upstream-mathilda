(* ==========================================================================
   Experiment 119 -- batch GCD over RSA moduli (group F)
   ==========================================================================
   WHAT IT MEASURES.  Heninger et al.'s weak-key scan ("Mining your Ps and
   Qs", 2012): given M RSA moduli, find every modulus that shares a prime
   with another, in quasi-linear time:

       product tree     P = N_1 N_2 ... N_M          (pairwise, bottom-up)
       remainder tree   z_i = P mod N_i^2            (top-down)
       z_i / N_i        exact, then  g_i = GCD[N_i, z_i / N_i]

   N_i is weak iff g_i > 1.  Written in the language ("app:") for M = 256
   and M = 1024 1024-bit moduli, its two tree phases timed separately at
   M = 1024, and compared with the naive all-pairs GCD at M = 256.  The
   big-integer kernels underneath are timed as builtins ("prim:"): the
   M-fold product (Times @@, the big-multiplication kernel) and the final
   vector of GCDs.

   Moduli: N_i = P_i Q_i with P_i = primeOf[5000 + 2 i, 512] and
   Q_i = primeOf[5001 + 2 i, 512], except that for every i divisible by 37
   the second factor is replaced by P_(i - 17): that plants a shared prime
   between moduli i and i - 17 (6 pairs = 12 weak moduli in the first 256,
   27 pairs = 54 weak in 1024).

   CHECK.  weak * 10^10 + (sum of the recovered g_i over weak moduli
   mod $CheckMod).  Each shared prime is counted once per modulus it
   breaks, in every system.
   ========================================================================== *)

Get["../harness.m"];
Get["../crypto_data.m"];

require[{"Times", "GCD", "Mod", "Quotient", "NextPrime"}];

$M = 1024;
bgP = Table[primeOf[5000 + 2 i, 512], {i, $M}];
bgQ = Table[If[Mod[i, 37] == 0, bgP[[i - 17]], primeOf[5001 + 2 i, 512]], {i, $M}];
moduli = bgP bgQ;
moduli256 = Take[moduli, 256];

(* ---- the algorithm, written in the language ---------------------------- *)

(* levels bottom-up: {ns, pair products, ..., {P}}.  M is a power of two, so
   plain Partition[#, 2] pairs every level; the general ragged spelling
   Partition[#, 2, 2, 1, {}] is unevaluated in Mathilda and sends this
   NestWhileList into an unbounded loop (see README, Findings). *)
productTree[ns_] := NestWhileList[Times @@@ Partition[#, 2] &, ns, Length[#] > 1 &];

(* z_i = P mod N_i^2, walking the tree top-down *)
remainderTree[tree_] := Fold[
  Function[{r, level}, Table[Mod[r[[Ceiling[j/2]]], level[[j]]^2], {j, Length[level]}]],
  Last[tree], Reverse[Most[tree]]];

batchGCD[ns_] := Module[{z = remainderTree[productTree[ns]]},
  MapThread[GCD[#1, Quotient[#2, #1]] &, {ns, z}]];

weakDigest[g_] := Module[{w = Select[g, # > 1 &]}, Length[w] 10^10 + Mod[Total[w], $CheckMod]];

(* ---- prim: the big-integer kernels -------------------------------------- *)

bench["prim: Times@@ 1024 moduli (1M-bit product) x4",
      Do[bigProd = Times @@ moduli, {4}];];
check["prim: Times@@ 1024 moduli (1M-bit product) x4", Mod[bigProd, $CheckMod]];

tree1024 = productTree[moduli];
z1024 = remainderTree[tree1024];
zq1024 = Quotient[z1024, moduli];
bench["prim: GCD 1024-bit pairs x1024x8",
      Do[gFinal = MapThread[GCD, {moduli, zq1024}], {8}];];
check["prim: GCD 1024-bit pairs x1024x8", weakDigest[gFinal]];

(* ---- app: the two tree phases, then the whole scan ---------------------- *)

bench["app: product tree M=1024 x4", Do[tr = productTree[moduli], {4}];];
check["app: product tree M=1024 x4", Mod[Last[tr][[1]], $CheckMod] + 10^10 Length[tr]];

bench["app: remainder tree mod N^2 M=1024", zz = remainderTree[tree1024];];
check["app: remainder tree mod N^2 M=1024", Mod[Total[zz], $CheckMod]];

bench["app: batch GCD M=256 x4", Do[g256 = batchGCD[moduli256], {4}];];
check["app: batch GCD M=256 x4", weakDigest[g256]];

bench["app: batch GCD M=1024", g1024 = batchGCD[moduli];];
check["app: batch GCD M=1024", weakDigest[g1024]];

(* all pairs: g_i is the product of every nontrivial GCD found for modulus i *)
naiveGCD[ns_] := Module[{m = Length[ns], g = ConstantArray[1, Length[ns]], d},
  Do[d = GCD[ns[[i]], ns[[j]]];
     If[d > 1, g[[i]] *= d; g[[j]] *= d],
     {i, m - 1}, {j, i + 1, m}];
  g];

benchOnce["app: naive pairwise GCD M=256 (32640 GCDs)", gNaive = naiveGCD[moduli256];];
check["app: naive pairwise GCD M=256 (32640 GCDs)", weakDigest[gNaive]];
