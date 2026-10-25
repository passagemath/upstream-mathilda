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

(* ---- products ------------------------------------------------------------ *)

x = N[Mod[Range[$n], 13]]/7.;
bench["dot: A . x (SpMV) x100", Do[r = A . x, {100}]];
check["dot: A . x (SpMV) x100", Round[1000 Total[A . x]]];
bench["dot: x . A x10", Do[r = x . A, {10}]];
check["dot: x . A x10", Round[1000 Total[x . A]]];
bench["dot: A . A (SpGEMM)", r = A . A];
check["dot: A . A (SpGEMM)", digest[r]];
bench["struct: Transpose[A]", At = Transpose[A]];
check["struct: Transpose[A]", digest[At]];
bench["dot: A . Transpose[A] (SpGEMM)", r = A . At];
check["dot: A . Transpose[A] (SpGEMM)", digest[r]];

(* ---- reductions ----------------------------------------------------------- *)

bench["reduce: Total[A] (column sums)", r = Total[A]];
check["reduce: Total[A] (column sums)", Round[1000 Total[r]]];
bench["reduce: Total[A, {2}] (row sums)", r = Total[A, {2}]];
check["reduce: Total[A, {2}] (row sums)", Round[1000 Total[r]]];
bench["reduce: Total[A, 2]", r = Total[A, 2]];
check["reduce: Total[A, 2]", Round[1000 r]];
bench["reduce: Max[A]", r = Max[A]];
check["reduce: Max[A]", Round[1000 r]];
bench["reduce: Norm[A, Frobenius]", r = Norm[A, "Frobenius"]];
check["reduce: Norm[A, Frobenius]", Round[1000 r]];

(* ---- assembly ------------------------------------------------------------- *)

bench["struct: Join[A, B] (vstack)", r = Join[A, B]];
check["struct: Join[A, B] (vstack)", digest[r]];
bench["struct: ArrayFlatten 2x2 blocks", r = ArrayFlatten[{{A, B}, {B, A}}]];
check["struct: ArrayFlatten 2x2 blocks", digest[r]];
bench["struct: IdentityMatrix[10^6, SparseArray]", r = IdentityMatrix[10^6, SparseArray]];
check["struct: IdentityMatrix[10^6, SparseArray]", Round[Total[r["NonzeroValues"]]]];

(* ---- 2-D Poisson: 5-point Laplacian on a 300 x 300 grid ------------------- *)

$m = 300;
bench["poisson: build L (Band, KroneckerProduct)",
      (T1 = SparseArray[{Band[{1, 1}] -> 2., Band[{1, 2}] -> -1., Band[{2, 1}] -> -1.}, {$m, $m}];
       I1 = IdentityMatrix[$m, SparseArray];
       L = KroneckerProduct[T1, I1] + KroneckerProduct[I1, T1])];
check["poisson: build L (Band, KroneckerProduct)", digest[L]];
bb = N[Mod[Range[$m^2], 7]] - 3.;
bench["poisson: LinearSolve[L, b] direct", sol = LinearSolve[L, bb]];
check["poisson: LinearSolve[L, b] direct", Round[1000 Total[sol]]];
bench["poisson: LinearSolve[L, b] CG 1e-10",
      solcg = LinearSolve[L, bb, Method -> {"Krylov", Method -> "ConjugateGradient",
                                            "Preconditioner" -> None, Tolerance -> 10^-10}]];
check["poisson: LinearSolve[L, b] CG 1e-10", Round[Total[solcg]]];
bench["poisson: 6 largest eigenvalues", ev = Eigenvalues[L, 6]];
check["poisson: 6 largest eigenvalues", Round[10^6 Total[ev]]];
bench["poisson: 6 smallest eigenvalues", ev = Eigenvalues[L, -6]];
check["poisson: 6 smallest eigenvalues", Round[10^6 Total[ev]]];
