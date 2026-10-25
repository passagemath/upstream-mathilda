(* ==========================================================================
   Experiment 117 -- Module-LWE / Kyber-like encryption (group F, lattices)
   ==========================================================================
   WHAT IT MEASURES.  The integer-array kernels lattice cryptography leans
   on, at toy-but-real sizes (q = 3329, the Kyber modulus):

     * the negacyclic product in R_q = Z_q[x]/(x^256 + 1), spelled the way an
       array programmer writes it without an NTT (experiment 116 covers the
       NTT): a full linear ListConvolve, then fold the top half back with a
       sign flip;
     * plain (Regev) LWE on a 512 x 512 packed integer matrix:
       Mod[A . s + e, q], and the batched R . A of an encryption;
     * the modular-matrix primitives Inverse[M, Modulus -> q] and
       LinearSolve[M, b, Modulus -> q];

   then whole protocols written in the language ("app:" cases):

     * Kyber-lite (no NTT, no compression), k = 2 and k = 3:
         keygen    t = A s + e
         encrypt   u = A^T r + e1,  v = t . r + e2 + Round[q/2] m
         decrypt   m' = round-to-{0, q/2} of v - s . u
       for N 256-bit messages;
     * Regev LWE bit encryption, one 512-vector r per bit.

   All matrices, secrets and noise come from crypto_data.m's LCG
   (smallList): uniform entries mod q, and centered-binomial noise with
   eta = 2 (four LCG bits per coefficient, b1 + b2 - b3 - b4).  The data is
   built OUTSIDE the timed region; the timed region is the ring arithmetic.

   CHECKS.  An app case checks  correct * 10^10 + hash,  where correct is the
   number of decrypted bits that match (must equal the total) and hash is the
   position-weighted sum of every ciphertext coefficient mod $CheckMod.

   Inverse / LinearSolve with Modulus -> q currently error in Mathilda
   (Inverse::method), so those two cases are guarded by an explicit probe and
   emit SKIP (brief rule 6).
   ========================================================================== *)

Get["../harness.m"];
Get["../crypto_data.m"];

require[{"ListConvolve", "Dot", "Mod", "Quotient", "Inverse", "LinearSolve"}];

q = 3329;
$n = 256;
$half = 1665;                              (* Round[q/2], spelled exactly *)

(* centered binomial, eta = 2, from 4 LCG bits per coefficient *)
cbd[seed_, len_] := Partition[smallList[seed, 4 len, 2], 4] . {1, 1, -1, -1};
uniform[seed_, len_] := smallList[seed, len, q];

(* negacyclic product in Z_q[x]/(x^256 + 1) *)
pmul[a_, b_] := Module[{c = ListConvolve[a, b, {1, -1}, 0]},
  Mod[c[[;; $n]] - Append[c[[$n + 1 ;;]], 0], q]];

(* bit recovery: 1 iff Mod[w, q] lies in [q/4, 3q/4) *)
toBits[w_] := Mod[Quotient[4 Mod[w, q] + q, 2 q], 2];

(* position-weighted fold of a nested integer array *)
foldHash[x_] := Module[{f = Flatten[x]}, Mod[Range[Length[f]] . f, $CheckMod]];

(* ---- prim: the ring product ------------------------------------------- *)

pa = uniform[11, $n]; pb = cbd[12, $n];
benchOnce["prim: negacyclic polymul n=256 x200",
          Do[pc = pmul[pa, pb], {200}];];
check["prim: negacyclic polymul n=256 x200", foldHash[pc]];

(* ---- prim: plain LWE on a 512 x 512 packed matrix ---------------------- *)

lweA = Partition[uniform[21, 512 512], 512];
lweS = cbd[22, 512]; lweE = cbd[23, 512];
bench["prim: LWE Mod[A.s+e,q] 512x512 x200",
      Do[lweB = Mod[lweA . lweS + lweE, q], {200}];];
check["prim: LWE Mod[A.s+e,q] 512x512 x200", foldHash[lweB]];

lweR = Partition[smallList[24, 256 512, 2], 512];
bench["prim: LWE batch Mod[R.A,q] 256x512x512",
      lweU = Mod[lweR . lweA, q];];
check["prim: LWE batch Mod[R.A,q] 256x512x512", foldHash[lweU]];

(* ---- prim: modular matrix inverse / solve ------------------------------ *)

invM = Partition[uniform[31, 128 128], 128];
invB = uniform[32, 128];
okModInv = MatrixQ[Quiet[Inverse[{{1, 2}, {3, 4}}, Modulus -> 7]], IntegerQ];
okModSolve = VectorQ[Quiet[LinearSolve[{{1, 2}, {3, 4}}, {1, 1}, Modulus -> 7]], IntegerQ];
If[okModInv,
  bench["prim: Inverse[M,Modulus->q] 128x128 x4",
        Do[invMi = Inverse[invM, Modulus -> q], {4}];];
  check["prim: Inverse[M,Modulus->q] 128x128 x4", foldHash[invMi]],
  Print["SKIP\t", "prim: Inverse[M,Modulus->q] 128x128 x4", "\tInverse[Modulus]"]];
If[okModSolve,
  bench["prim: LinearSolve[M,b,Modulus->q] 128x128 x16",
        Do[invX = LinearSolve[invM, invB, Modulus -> q], {16}];];
  check["prim: LinearSolve[M,b,Modulus->q] 128x128 x16", foldHash[invX]],
  Print["SKIP\t", "prim: LinearSolve[M,b,Modulus->q] 128x128 x16", "\tLinearSolve[Modulus]"]];

(* ---- app: Kyber-lite --------------------------------------------------- *)

(* Inputs for rank k and N messages: seeds are offset by 1000 k. *)
kyberData[k_, nMsg_] := Module[{b = 1000 k},
  <|"A" -> Table[uniform[b + 10 i + j, $n], {i, k}, {j, k}],
    "s" -> Table[cbd[b + 100 + i, $n], {i, k}],
    "e" -> Table[cbd[b + 150 + i, $n], {i, k}],
    "r" -> Table[cbd[b + 200 + 10 t + i, $n], {t, nMsg}, {i, k}],
    "e1" -> Table[cbd[b + 500 + 10 t + i, $n], {t, nMsg}, {i, k}],
    "e2" -> Table[cbd[b + 800 + t, $n], {t, nMsg}],
    "m" -> Table[smallList[b + 900 + t, $n, 2], {t, nMsg}]|>];

(* The ring-vector sums are spelled with Sum.  The equally natural MapThread
   spelling is kept as its own case below: it currently returns unevaluated in
   Mathilda (see README, Findings). *)
kyberKeygen[a_, s_, e_] := Module[{k = Length[s]},
  Table[Mod[Sum[pmul[a[[i, j]], s[[j]]], {j, k}] + e[[i]], q], {i, k}]];

kyberEncrypt[a_, t_, r_, e1_, e2_, m_] := Module[{k = Length[r], u, v},
  u = Table[Mod[Sum[pmul[a[[j, i]], r[[j]]], {j, k}] + e1[[i]], q], {i, k}];
  v = Mod[Sum[pmul[t[[j]], r[[j]]], {j, k}] + e2 + $half m, q];
  Append[u, v]];

kyberDecrypt[s_, ct_] := Module[{k = Length[s]},
  toBits[ct[[k + 1]] - Sum[pmul[s[[j]], ct[[j]]], {j, k}]]];

(* the same three steps with MapThread over the ring vectors *)
mtKeygen[a_, s_, e_] :=
  Table[Mod[Total[MapThread[pmul, {a[[i]], s}]] + e[[i]], q], {i, Length[s]}];

mtEncrypt[a_, t_, r_, e1_, e2_, m_] := Module[{at = Transpose[a], u, v},
  u = Table[Mod[Total[MapThread[pmul, {at[[i]], r}]] + e1[[i]], q], {i, Length[r]}];
  v = Mod[Total[MapThread[pmul, {t, r}]] + e2 + $half m, q];
  Append[u, v]];

mtDecrypt[s_, ct_] :=
  toBits[ct[[-1]] - Total[MapThread[pmul, {s, Most[ct]}]]];

kyberRun[d_, kg_, enc_, dec_] := Module[{t, cts, msg},
  t = kg[d["A"], d["s"], d["e"]];
  cts = Table[enc[d["A"], t, d["r"][[j]], d["e1"][[j]], d["e2"][[j]], d["m"][[j]]],
              {j, Length[d["m"]]}];
  msg = dec[d["s"], #] & /@ cts;
  Count[Flatten[msg - d["m"]], 0] 10^10 + foldHash[cts]];

ky2 = kyberData[2, 8];
benchOnce["app: Kyber-lite k=2 keygen+enc+dec x8", kres2 = kyberRun[ky2, kyberKeygen, kyberEncrypt, kyberDecrypt];];
check["app: Kyber-lite k=2 keygen+enc+dec x8", kres2];

ky3 = kyberData[3, 8];
benchOnce["app: Kyber-lite k=3 keygen+enc+dec x8", kres3 = kyberRun[ky3, kyberKeygen, kyberEncrypt, kyberDecrypt];];
check["app: Kyber-lite k=3 keygen+enc+dec x8", kres3];

benchOnce["app: Kyber-lite k=2 MapThread spelling x8",
          kres2m = Quiet[kyberRun[ky2, mtKeygen, mtEncrypt, mtDecrypt]];];
check["app: Kyber-lite k=2 MapThread spelling x8", If[IntegerQ[kres2m], kres2m, -1]];

(* ---- app: Regev LWE, one bit per 512-vector ---------------------------- *)

lweMsg = smallList[25, 256, 2];
lweRun[] := Module[{b, cts, dec},
  b = Mod[lweA . lweS + lweE, q];
  cts = Table[{Mod[lweR[[t]] . lweA, q], Mod[lweR[[t]] . b + $half lweMsg[[t]], q]},
              {t, 256}];
  dec = Table[toBits[c[[2]] - lweS . c[[1]]], {c, cts}];
  Count[dec - lweMsg, 0] 10^10 + foldHash[cts]];

bench["app: Regev LWE 512x512 enc+dec 256 bits", lres = lweRun[];];
check["app: Regev LWE 512x512 enc+dec 256 bits", lres];
