(* ==========================================================================
   Experiment 120 -- SparseArray construction, arithmetic and element-wise
   functions against scipy.sparse
   ==========================================================================
   WHAT IT MEASURES.  The sparse-array substrate a sparse linear-algebra or
   graph workload sits on, on a 10^5 x 10^5 matrix with 10^6 explicit
   entries (density 10^-4):

       build      coordinate list -> compressed rows   (SparseArray[pos -> v])
       arith      A + B, A B (element-wise), 2.5 A, A^2, A + A
       functions  Sin[A], Sqrt[A]                      (sparsity-preserving)
       access     one row, one element, the explicit entries
       dense      a 2000 x 2000 dense <-> sparse round trip

   Positions are a Lehmer LCG, i_k = (48271 k mod (2^31 - 1)) mod n + 1 and
   j_k = (69621 k mod (2^31 - 1)) mod n + 1 (B shifts j by 7 before the
   reduction for odd k only, so A and B share half their positions), so
   every system builds the same matrices; no position repeats within a
   matrix, so first-rule-wins and sum-duplicates agree.  Values are
   (k mod 97)/10 + 0.1 > 0.

   CHECK.  Round[1000 x (sum of the explicit values)] of each result -- an
   Integer, so the systems compare it exactly.
   ========================================================================== *)

Get["../harness.m"];

require[{"SparseArray", "ArrayRules", "Normal"}];

$n = 100000; $nnz = 1000000;
kk = Range[$nnz];
ii = Mod[Mod[48271 kk, 2147483647], $n] + 1;
jj = Mod[Mod[69621 kk, 2147483647], $n] + 1;
jj2 = Mod[Mod[69621 kk + 7 Mod[kk, 2], 2147483647], $n] + 1;
vv = N[Mod[kk, 97]]/10. + 0.1;
posA = Transpose[{ii, jj}];
posB = Transpose[{ii, jj2}];

digest[s_] := Round[1000 Total[s["NonzeroValues"]]];

(* ---- build -------------------------------------------------------------- *)

bench["build: 10^6 entries, 10^5 x 10^5 (coordinates -> rows)",
      A = SparseArray[posA -> vv, {$n, $n}]];
check["build: 10^6 entries, 10^5 x 10^5 (coordinates -> rows)", digest[A]];
B = SparseArray[posB -> vv, {$n, $n}];

(* ---- arithmetic --------------------------------------------------------- *)

bench["arith: A + B", r = A + B];
check["arith: A + B", digest[r]];
bench["arith: A B (element-wise)", r = A B];
check["arith: A B (element-wise)", digest[r]];
bench["arith: 2.5 A", r = 2.5 A];
check["arith: 2.5 A", digest[r]];
bench["arith: A^2 (element-wise)", r = A^2];
check["arith: A^2 (element-wise)", digest[r]];
bench["arith: A + A", r = A + A];
check["arith: A + A", digest[r]];

(* ---- element-wise functions --------------------------------------------- *)

bench["fn: Sin[A]", r = Sin[A]];
check["fn: Sin[A]", digest[r]];
bench["fn: Sqrt[A]", r = Sqrt[A]];
check["fn: Sqrt[A]", digest[r]];

(* ---- access ------------------------------------------------------------- *)

bench["access: row 500 x100", Do[r = A[[500]], {100}]];
check["access: row 500 x100", Round[1000 Total[Normal[A[[500]]]]]];
bench["access: element (500, 600) x1000", Do[r = A[[500, 600]], {1000}]];
bench["access: explicit positions and values", (p = A["NonzeroPositions"]; v = A["NonzeroValues"])];
check["access: explicit positions and values", Total[p, 2]];

(* ---- dense <-> sparse --------------------------------------------------- *)

dense = Table[If[Mod[i + 3 j, 97] == 0, 1.5 + Mod[i j, 7], 0.], {i, 2000}, {j, 2000}];
bench["dense: 2000 x 2000 dense -> sparse", Dm = SparseArray[dense]];
check["dense: 2000 x 2000 dense -> sparse", digest[Dm]];
bench["dense: 2000 x 2000 sparse -> dense", r = Normal[Dm]];
check["dense: 2000 x 2000 sparse -> dense", Round[1000 Total[r, 2]]];
