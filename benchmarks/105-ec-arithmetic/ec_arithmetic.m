(* ==========================================================================
   Experiment 105 -- elliptic-curve point arithmetic on secp256k1
   (group F, cryptographic number theory)
   ==========================================================================
   WHAT IT MEASURES.  The 256-bit modular arithmetic an elliptic-curve
   implementation lives on, written in the language ("app:" cases): affine
   add / double with the slope's inverse taken by PowerMod[x, -1, p], the
   inversion-free Jacobian-coordinate add / double, and scalar multiplication
   k*G by left-to-right double-and-add (affine and Jacobian) and by the
   Montgomery ladder (Jacobian).  Each "prim:" case is a system's native
   elliptic-curve builtin: PARI's ellmul / elladd on ellinit([0,0,0,0,7], p).

   secp256k1 (y^2 = x^3 + 7 over F_p) and the scalars come from
   crypto_data.m; scalars are bigOf[s, 256] mod n.  Every scalar-mult check is
   the sum of the affine x-coordinates of k*G folded mod $CheckMod, so the
   three in-language algorithms and PARI's ellmul are cross-verified against
   each other by carrying the SAME check value.

   WL (and so Mathilda) has no elliptic-curve-group head, so the "prim:" cases
   are benchIf'd on "EllipticCurve" and report ABSENT.

   The point at infinity is {} in affine coordinates and {1, 1, 0} (Z = 0) in
   Jacobian coordinates.
   ========================================================================== *)

Get["../harness.m"];
Get["../crypto_data.m"];

require[{"PowerMod", "IntegerDigits", "EllipticCurve"}];

ecP = $Secp256k1P;
ecN = $Secp256k1N;
ecG = {$Secp256k1Gx, $Secp256k1Gy};

(* ---- affine coordinates (one inversion per add / double) ------------- *)

ecDblA[{}] := {};
ecDblA[{x_, y_}] := If[y == 0, {},
  Module[{l = Mod[3 x^2 PowerMod[2 y, -1, ecP], ecP], x3},
    x3 = Mod[l^2 - 2 x, ecP];
    {x3, Mod[l (x - x3) - y, ecP]}]];

ecAddA[{}, q_] := q;
ecAddA[p_, {}] := p;
ecAddA[{x1_, y1_}, {x2_, y2_}] := Which[
  x1 == x2, If[Mod[y1 + y2, ecP] == 0, {}, ecDblA[{x1, y1}]],
  True, Module[{l = Mod[(y2 - y1) PowerMod[x2 - x1, -1, ecP], ecP], x3},
    x3 = Mod[l^2 - x1 - x2, ecP];
    {x3, Mod[l (x1 - x3) - y1, ecP]}]];

(* left-to-right double-and-add *)
ecMulA[k_, pt_] := Module[{r = {}},
  Do[r = ecDblA[r]; If[b == 1, r = ecAddA[r, pt]], {b, IntegerDigits[k, 2]}];
  r];

(* ---- Jacobian coordinates: (X, Y, Z) <-> (X/Z^2, Y/Z^3), no inversions - *)

toJ[{}] := {1, 1, 0};
toJ[{x_, y_}] := {x, y, 1};
fromJ[{_, _, 0}] := {};
fromJ[{x_, y_, z_}] := Module[{zi = PowerMod[z, -1, ecP]},
  {Mod[x zi^2, ecP], Mod[y zi^3, ecP]}];

(* dbl-2009-l (a = 0) *)
ecDblJ[{x_, y_, z_}] := If[z == 0 || y == 0, {1, 1, 0},
  Module[{a = Mod[x^2, ecP], b = Mod[y^2, ecP], c, d, e, x3},
    c = Mod[b^2, ecP];
    d = Mod[2 ((x + b)^2 - a - c), ecP];
    e = 3 a;
    x3 = Mod[e^2 - 2 d, ecP];
    {x3, Mod[e (d - x3) - 8 c, ecP], Mod[2 y z, ecP]}]];

(* add-1998-cmo-2 *)
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

(* Montgomery ladder, Jacobian: R1 - R0 = P is invariant *)
ecLadderJ[k_, pt_] := Module[{r0 = {1, 1, 0}, r1 = toJ[pt]},
  Do[If[b == 0,
       r1 = ecAddJ[r0, r1]; r0 = ecDblJ[r0],
       r0 = ecAddJ[r0, r1]; r1 = ecDblJ[r1]],
     {b, IntegerDigits[k, 2]}];
  fromJ[r0]];

xDigest[pts_] := Mod[Total[First /@ pts], $CheckMod];

ks = Table[Mod[bigOf[s, 256], ecN], {s, 501, 564}];

(* ---- app: single operations, batched ---------------------------------- *)

g2 = ecDblA[ecG];

bench["app: affine point add x10000",
      qa = g2; Do[qa = ecAddA[qa, ecG], {10000}];];
check["app: affine point add x10000", Mod[First[qa], $CheckMod]];

gJ = toJ[ecG];
bench["app: Jacobian point add x10000",
      qj = toJ[g2]; Do[qj = ecAddJ[qj, gJ], {10000}];];
check["app: Jacobian point add x10000", Mod[First[fromJ[qj]], $CheckMod]];

bench["app: affine point double x10000",
      da = ecG; Do[da = ecDblA[da], {10000}];];
check["app: affine point double x10000", Mod[First[da], $CheckMod]];

bench["app: Jacobian point double x10000",
      dj = gJ; Do[dj = ecDblJ[dj], {10000}];];
check["app: Jacobian point double x10000", Mod[First[fromJ[dj]], $CheckMod]];

(* ---- app: scalar multiplication k*G ----------------------------------- *)

bench["app: scalar mult affine double-and-add x64",
      ma = ecMulA[#, ecG] & /@ ks;];
check["app: scalar mult affine double-and-add x64", xDigest[ma]];

bench["app: scalar mult Jacobian double-and-add x64",
      mj = ecMulJ[#, ecG] & /@ ks;];
check["app: scalar mult Jacobian double-and-add x64", xDigest[mj]];

bench["app: scalar mult Jacobian Montgomery ladder x64",
      ml = ecLadderJ[#, ecG] & /@ ks;];
check["app: scalar mult Jacobian Montgomery ladder x64", xDigest[ml]];

(* ---- prim: a native elliptic-curve group (none in WL) ----------------- *)

benchIf["prim: point add x10000", "EllipticCurve", Null];
benchIf["prim: scalar mult k*G x64", "EllipticCurve", Null];
