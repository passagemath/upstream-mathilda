(* ==========================================================================
   Experiment 106 -- ECDSA sign / verify over secp256k1
   (group F, cryptographic number theory)
   ==========================================================================
   WHAT IT MEASURES.  ECDSA written in the language ("app:" cases) on the
   in-language Jacobian point arithmetic of experiment 105 (copied below,
   not imported -- each experiment is self-contained):

       sign:    R = k G,  r = x(R) mod n,  s = k^-1 (z + r d) mod n
       verify:  w = s^-1, u1 = z w, u2 = r w,  X = u1 G + u2 Q  (Shamir's
                trick: one joint double-and-add over the bits of u1, u2 with
                G + Q precomputed),  valid iff x(X) mod n == r

   The nonce k is NOT RFC 6979: it is a deterministic value from
   crypto_data's LCG (bigOf[s, 256] mod n), so all systems sign identically.
   The message hash z is an integer from bigOf (portable, no hashing).  The
   separate "prim: SHA-256 of 1000 messages" case times the hash a real
   signer applies first, via Hash[string, "SHA256"]; Mathilda's Hash returns
   that form unevaluated today, so the case is guarded by an explicit probe
   (it is not caught by benchIf -- the head exists).

   "prim: ECDSA sign / verify" are a native curve group (PARI ellmul /
   elladd; Python's OpenSSL verify via the `cryptography` package); WL has no
   elliptic-curve head, so they are benchIf'd on "EllipticCurve" (ABSENT).

   Checks: sign -> Sum(r + s) mod $CheckMod; verify -> count of valid
   signatures (48); tampered (z + 1) -> count of valid (0).
   ========================================================================== *)

Get["../harness.m"];
Get["../crypto_data.m"];

require[{"PowerMod", "IntegerDigits", "Hash", "EllipticCurve"}];

ecP = $Secp256k1P;
ecN = $Secp256k1N;
ecG = {$Secp256k1Gx, $Secp256k1Gy};

(* ---- Jacobian point arithmetic (as in 105-ec-arithmetic) -------------- *)

toJ[{}] := {1, 1, 0};
toJ[{x_, y_}] := {x, y, 1};
fromJ[{_, _, 0}] := {};
fromJ[{x_, y_, z_}] := Module[{zi = PowerMod[z, -1, ecP]},
  {Mod[x zi^2, ecP], Mod[y zi^3, ecP]}];

ecDblJ[{x_, y_, z_}] := If[z == 0 || y == 0, {1, 1, 0},
  Module[{a = Mod[x^2, ecP], b = Mod[y^2, ecP], c, d, e, x3},
    c = Mod[b^2, ecP];
    d = Mod[2 ((x + b)^2 - a - c), ecP];
    e = 3 a;
    x3 = Mod[e^2 - 2 d, ecP];
    {x3, Mod[e (d - x3) - 8 c, ecP], Mod[2 y z, ecP]}]];

ecAddJ[p1_, p2_] := Which[
  p1[[3]] == 0, p2,
  p2[[3]] == 0, p1,
  True, Module[{x1, y1, z1, x2, y2, z2, z1z1, z2z2, u1, u2, s1, s2, h, r,
                h2, h3, u1h2, x3},
    {x1, y1, z1} = p1; {x2, y2, z2} = p2;
    z1z1 = Mod[z1^2, ecP]; z2z2 = Mod[z2^2, ecP];
    u1 = Mod[x1 z2z2, ecP]; u2 = Mod[x2 z1z1, ecP];
    s1 = Mod[y1 z2 z2z2, ecP]; s2 = Mod[y2 z1 z1z1, ecP];
    h = Mod[u2 - u1, ecP]; r = Mod[s2 - s1, ecP];
    If[h == 0,
      If[r == 0, ecDblJ[p1], {1, 1, 0}],
      h2 = Mod[h^2, ecP]; h3 = Mod[h h2, ecP]; u1h2 = Mod[u1 h2, ecP];
      x3 = Mod[r^2 - h3 - 2 u1h2, ecP];
      {x3, Mod[r (u1h2 - x3) - s1 h3, ecP], Mod[z1 z2 h, ecP]}]]];

ecMulJ[k_, pt_] := Module[{r = {1, 1, 0}, pj = toJ[pt]},
  Do[r = ecDblJ[r]; If[b == 1, r = ecAddJ[r, pj]], {b, IntegerDigits[k, 2]}];
  fromJ[r]];

(* Shamir's trick: u1 G + u2 Q in one pass, G + Q precomputed. *)
shamirJ[u1_, u2_, g_, q_] := Module[{gj = toJ[g], qj = toJ[q], gq, len, b1, b2,
                                     r = {1, 1, 0}},
  gq = ecAddJ[gj, qj];
  len = Max[BitLength[u1], BitLength[u2]];
  b1 = IntegerDigits[u1, 2, len]; b2 = IntegerDigits[u2, 2, len];
  Do[r = ecDblJ[r];
     Which[b1[[i]] == 1 && b2[[i]] == 1, r = ecAddJ[r, gq],
           b1[[i]] == 1, r = ecAddJ[r, gj],
           b2[[i]] == 1, r = ecAddJ[r, qj]],
     {i, len}];
  fromJ[r]];

(* ---- ECDSA ------------------------------------------------------------ *)

ecdsaSign[z_, d_, k_] := Module[{r, s},
  r = Mod[First[ecMulJ[k, ecG]], ecN];
  s = Mod[PowerMod[k, -1, ecN] (z + r d), ecN];
  {r, s}];

ecdsaVerify[z_, {r_, s_}, q_] := If[! (1 <= r < ecN && 1 <= s < ecN), 0,
  Module[{w, x},
    w = PowerMod[s, -1, ecN];
    x = shamirJ[Mod[z w, ecN], Mod[r w, ecN], ecG, q];
    If[x === {}, 0, If[Mod[First[x], ecN] == r, 1, 0]]]];

privD = Mod[bigOf[701, 256], ecN];
pubQ = ecMulJ[privD, ecG];

zs = Table[bigOf[s, 256], {s, 801, 848}];
nonces = Table[Mod[bigOf[s, 256], ecN], {s, 1801, 1848}];

bench["app: ECDSA sign x48",
      sigs = MapThread[ecdsaSign[#1, privD, #2] &, {zs, nonces}];];
check["app: ECDSA sign x48", Mod[Total[Flatten[sigs]], $CheckMod]];

bench["app: ECDSA verify x48 (Shamir)",
      nValid = Total[MapThread[ecdsaVerify[#1, #2, pubQ] &, {zs, sigs}]];];
check["app: ECDSA verify x48 (Shamir)", nValid];

bench["app: ECDSA verify tampered x48 (Shamir)",
      nBad = Total[MapThread[ecdsaVerify[#1 + 1, #2, pubQ] &, {zs, sigs}]];];
check["app: ECDSA verify tampered x48 (Shamir)", nBad];

(* ---- prim: SHA-256 of the message (probe: Hash may not evaluate) ------- *)

msgs = Table["message " <> ToString[i], {i, 1, 1000}];
shaOkQ = IntegerQ[Hash["abc", "SHA256"]];
If[shaOkQ,
  bench["prim: SHA-256 of 1000 messages", hs = Hash[#, "SHA256"] & /@ msgs;];
  check["prim: SHA-256 of 1000 messages", Mod[Total[hs], $CheckMod]],
  Print["SKIP\tprim: SHA-256 of 1000 messages\tHash[SHA256]"]];

(* ---- prim: a native curve group (none in WL) --------------------------- *)

benchIf["prim: ECDSA sign x48", "EllipticCurve", Null];
benchIf["prim: ECDSA verify x48", "EllipticCurve", Null];
